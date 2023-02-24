###############################################################################
###	 @author: Alexander Sandercock                                          ###
###        September 2022                                                   ###
###                                               							###
### need to use absolute path for files                                     ###
####                                                                        ###
####
#### Refer to https://github.com/alex-sandercock/American_chestnut_landscape_genomics
###############################################################################


######### Packages needed ##############################
#pip install numpy
#pip3 install fuc
import pandas as pd
from fuc import pyvcf #This requires fuc version >= 0.31.0
import gzip
import argparse
import statistics as stat
import numpy as np
import scipy.stats as st

############ 1. Getting the info from VCF file and subsetting population of interest ###############

def get_sample_list(sample_list:str):
    #Importing sample name list and converting to list
    
    # opening the file in read mode
    sample_file = open(sample_list, "r")
  
    # reading the file
    names = sample_file.read()
    
    #Making comma separated list from each new line
    names_list = names.split("\n")
    sample_file.close()
    
    return(names_list)

def get_subset_vcf(vcf_file,pop_list):
    #This function will import the VCF file and subset it based on the user supplied list
    vf = pyvcf.VcfFrame.from_file(str(vcf_file))
    if pop_list is not None: #If no pop list supplied, use all samples in vcf file
        vf_new = vf.subset(pop_list)
    else:
        vf_new = vf
        
    return(vf_new)


########### 2. Get allele frequency values for each allele across all samples ############

def get_pop_allele_freq(vcf_subset):
    #Calculate the allele frequency (AF) for each SNP for the subset samples
    
    vf_new = vcf_subset.compute_info('AF')
    
    return(vf_new)

def pop_AF_values(vcf_subset):
    #Need to get list of AF values that the bootstrap sampling will reference
    #These are the AF values when all individuals from a population of interest are included.
    
    #Check to make sure that the input vcf is a vcfFrame file and not a pandas dataframe, if not, convert.
    try:
        af_values = vcf_subset.extract_info('#AF')
    
    except AttributeError:
        af_values = pyvcf.VcfFrame(['##fileformat=VCFv4.3'], vcf_subset).extract_info('#AF')
        
     
    return(af_values)


############ 3. Perform bootstrap sampling ###########################

def bootstrap_batch_sample_regression(filtered_vcf,batch:int,target_values,r2_threshold):
    #Perform bootstrap sampling until threshold of diversity captured is met based on R2 value
    
    #Initialize new df that randomly sampled genotypes will be added to while retaining info columns
    df_merged = filtered_vcf.df.loc[:, filtered_vcf.df.columns.isin(['CHROM', 'POS','ID','REF','ALT','QUAL','FILTER','INFO','FORMAT'])]
    
    sampling_round = 0 #initialize recording of number of sampling rounds (while loops)
    
    diversity_captured = False #Initialize the diversity captured as False
    
    #Continue bootstrap sampling until the AF value for every allele is >= the threshold*original AF value
    while diversity_captured == False:
        #Add 1 to the number of sampling performed
        sampling_round += 1
    
        #Randomly select X number of individuals for the bootstrap sample size
        df = filtered_vcf.df.loc[:, ~filtered_vcf.df.columns.isin(['CHROM', 'POS','ID','REF','ALT','QUAL','FILTER','INFO','FORMAT'])].sample(n=batch,axis='columns')
    
        #Add the columns from each into a common df and erase the INFO column
        df_merged = pd.concat([df_merged, df], axis=1)
        df_merged['INFO'] = '.'
    
        #Convert to vcfFrame file again so that the AF can be calculated
        df_frame = pyvcf.VcfFrame(['##fileformat=VCFv4.3'], df_merged)
    
        #Calculate AF from only subset samples
        df_frame = df_frame.compute_info('AF')
        df_AF_values = list(df_frame.extract_info('#AF')) #retrieve those values and make list
    
        #Perform linear regression between the two lists of allele frequencies
        slope, intercept, r_value, p_value, std_err = st.linregress(df_AF_values, target_values)
        r2 = (r_value**2)
        
        #End if r2 value is greater than user defined threshold
        if r2 >= r2_threshold:
            break

    return(sampling_round)


############### 4. Calculate statistics ################################

def calculate_statistics(sampling_round_list,batch:int):
    
    #Convert sampling rounds to number of trees
    #Can I just multiply the list by the batch number?
    sampling_round_list = [x * batch for x in sampling_round_list]
    
    #Trees to sample from population of interest
    mean_trees = np.mean(sampling_round_list)
    
    #I am using population standard deviation since the list of numbers that I am evaluating are all of the iterations and not a subset of them: https://stats.stackexchange.com/questions/485326/confused-when-to-use-population-vs-sample-standard-deviation-in-engineering-test#:~:text=The%20two%20forms,calculating%20standard%20deviation%3F
    stand_dev = stat.pstdev(sampling_round_list)
    
    #create 95% confidence interval for population mean weight
    CI = st.t.interval(alpha=0.95, df=len(sampling_round_list)-1, loc=np.mean(sampling_round_list), scale=st.sem(sampling_round_list))
    
    return(print('Number of trees to sample =',mean_trees,'\nPopulation Standard Deviation =',stand_dev,'\n95% Confidence Intervals =',CI,'\nIterations performed =',len(sampling_round_list)))



######## 5. Define main() and user arguments #############################

def main():
    
    ########User arguments
    # Create the parser
	parser = argparse.ArgumentParser()

	# Add arguments
	parser.add_argument('--batch','-b', type=int, nargs='?', default=5, help = 'This is the number of samples for the bootstrap sampling to pull. Default is 5.')
	parser.add_argument('--vcf','-v', type=str, required=True, help = 'This is the full path to the adaptive SNP vcf file.')
	parser.add_argument('--iterations','-i', type=int, nargs='?', default=100, help = 'This is the number iterations to perform. Default is 100')
	parser.add_argument('--sample_file','-f', type=str, help = 'Provide a text file of sample names for a population of interest, with one sample on each line, that will be used to subset the vcf file.')
	parser.add_argument('--diversity_captured','-d', type=float, default=0.9, help = 'This is the threshold value of the percent of diversity to capture through sampling. Default is to capture 90% (0.9) genomic diversity of a population through sampling.')
	#parser.add_argument('--samples','-s', type=list, help = 'Provide a space separated list of sample names for the population of interest that will be used to subset the vcf file.')

	# Parse the arguments
	args = parser.parse_args()

	#assign variable names
	vcf_file = args.vcf
	sample_list = args.sample_file
	batch = args.batch
	iters = args.iterations
	threshold_ratio = args.diversity_captured


	#Get sample names from population of interest
	if sample_list is not None:
		name_list = get_sample_list(sample_list)
        name_list = [x for x in name_list if x != ''] #This removes any empty strings from name list

	else:
		name_list = None

	#Subset VCF file if list of sample provided
	subset_vcf = get_subset_vcf(vcf_file,name_list)

	#Calculate AF for each allele across all samples in VCF file
	subset_vcf = get_pop_allele_freq(subset_vcf)

	#Get the AF values for all samples as a list
	target_AF_values = pop_AF_values(subset_vcf)

	#Perform boostrap sampling over user input or default iterations
	sample_round_list = []

	for i in range(iters):
		sample_value = bootstrap_batch_sample_regression(subset_vcf,batch,target_AF_values,threshold_ratio)
		sample_round_list.append(sample_value)

	#Calculate ouput statistics from analysis
	calculate_statistics(sample_round_list,batch)


if __name__ == "__main__":
    main()
