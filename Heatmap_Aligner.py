#!/usr/bin/env python3

#import sys
#import os
#import biopython

'''
Script to take multiple PyLipID outputs and create a sequence alignment coloured by lipid statistics
Can be used with a single PyLipID output just to make a sequence-based heatmap
'''

############################
### Defining input files ###
############################

# Where files are
dir = '/sansom/s156a/bioc1535/MraY/EC_other_protein/Data'

# List of input csv files
csvfile_list = ["%s/Interaction_UDP1_DEDA/Dataset_UDP1/Dataset.csv" % dir, 
		"%s/Interaction_UDP1_UPTA/Dataset_UDP1/Dataset.csv" % dir ]


#########################
### Setting variables ###
#########################



#######################################
### Don't touch anything below here ###
#######################################

for csvfile in csvfile_list:
	print(csvfile)
