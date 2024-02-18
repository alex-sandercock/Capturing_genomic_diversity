## Background

This script was developed for the American chestnut landscape genomics manuscript. The aim was to develop a method for estimating the number of trees to sample from a population in order to capture a predefined percentage of diversity within a breeding population. Specifically, we sought to capture the adaptive diversity from Adaptive Units within TACF chestnut breeding program to develop locally adapted, blight-resistent American chestnut trees.

The preprint for the article can be found at bioRxiv: https://doi.org/10.1101/2023.05.30.542850

## Getting Started

A vcf file should be made beforehand that contains only the putatively adaptive alleles that were previously identified.

## Usage

The script will be executed in the command line. This script can be used to estimate the number of trees to sample from a population in order to capture a desired percentage of the genomic diversity. 

#### Required arguments

- ```--vcf```, ```-v```: This specifies the path to the adaptive SNP vcf file. The vcf can be gzipped.

#### Optional arguments

- ```--batch```, ```-b```: This is the number of genotypes to be randomly selected for each bootstrap sampling. Default is 1.

- ```--iteration```, ```-i```: This is the number of iterations to perform for the analysis. Default is 100.

- ```--sample_file```, ```-f```: This specifies a .txt file containing the list of sample names to subset the vcf file by. There should be one sample name per line. By default, all samples will be used in the vcf file.

- ```--diversity_captured```,```-d```: This is the r2 threshold of the desired percentage of diversity to capture. Default is 0.9.

#### Example:

```

python capturing_diversity.py --vcf path_to_file/adaptive_diversity_SNPs.vcf.gz -i 999

```

## Output data

The output will be provided in the terminal following the completion of all iterations:

```

Number of trees to sample = 22.0 

95% Confidence Intervals = (16.378127142052662, 27.621872857947338) 

Iterations performed = 5

```

## Additional help

Refer to the jupyter notebook with a guided example in the /doc folder

