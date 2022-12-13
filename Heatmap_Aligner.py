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
csvfile_list = [ "%s/Interaction_UDP1_DEDA/Dataset_UDP1/Dataset.csv" % dir,
		 "%s/Interaction_UDP1_UPTA/Dataset_UDP1/Dataset.csv" % dir ]

# List of input coordinate files - can be CG or AT
coord_files = [	"%s/Interaction_UDP1_UPTA/Coordinate_UDP1/Coordinate_UDP1_Occupancy.pdb" % dir,
		"%s/Interaction_UDP1_UPTA/Coordinate_UDP1/Coordinate_UDP1_Occupancy.pdb" % dir ]

##########################
### Set variables here ###
##########################



#######################################
### Don't touch anything below here ###
#######################################

def get_number_of_systems():
	csv_number = len(csvfile_list)
	coord_number = len(coord_files)
	if csv_number == coord_number:
		
		print('%s input files being read' % csv_number)
	else:
		print('mismatch in input files, %s csv files and %s coordinate files' % (csv_number, coord_number))

def sequence_alignment():
	pass

def get_lipid_heatmap():
	pass

def plot_heatmap_on_alignment():
	pass

for csvfile in csvfile_list:
	#print(csvfile)
	pass

######################
### Run everything ###
######################

get_number_of_systems()
