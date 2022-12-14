#!/usr/bin/env python3

import os
import numpy as np
import csv
import re
from subprocess import Popen, PIPE, STDOUT
# biopython needed: pip install biopython 
# still need to work out how best to align these sequences
#from Bio import AlignIO,SeqIO
#from Bio import pairwise2
#from Bio.pairwise2 import format_alignment

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

# List of names for sequnece alignment (optional)
system_names = [ "DEDA", "UPTA" ]

# List of input coordinate files - can be CG or AT
#coord_files = [	"%s/Interaction_UDP1_UPTA/Coordinate_UDP1/Coordinate_UDP1_Occupancy.pdb" % dir,
#		"%s/Interaction_UDP1_UPTA/Coordinate_UDP1/Coordinate_UDP1_Occupancy.pdb" % dir ]

##########################
### Set variables here ###
##########################



#######################################
### Don't touch anything below here ###
#######################################

def get_number_of_systems():
	# not currently needed
	csv_number = len(csvfile_list)
	coord_number = len(coord_files)
	if csv_number == coord_number:
		print('%s input files being read' % csv_number)
		return csv_number
	else:
		print('mismatch in input files, %s csv files and %s coordinate files' % (csv_number, coord_number))

def get_sequence(csvfile,column):
	#sequence, heatmap = np.loadtxt(delimiter=',',fname=csvfile, usecols=(0, 1), skiprows=1, unpack=True)
	with open (csvfile) as inf:
		reader = csv.reader(inf, delimiter=",")
		sequence = list(zip(*reader))[0] 
	with open (csvfile) as inf:
		reader = csv.reader(inf, delimiter=",")
		heatmap = list(zip(*reader))[column]
	return sequence[1:],heatmap[1:], heatmap[:1]

def write_fasta(sequence,csvfile):
	f = open('HeatmapAlignment/Sequences.txt','a')
	f.write('> sequence from %s\n' % csvfile) 
	for res in sequence:
		# strip residue number
		output = re.sub(r'\d+', '', res)
		f.write(res_dict[output]) 
	f.write('\n\n')

def write_data(heatmap,attribute,csvfile):
	f = open('HeatmapAlignment/%s.txt' % attribute,'a')
	f.write('> %s from %s\n' % (attribute, csvfile))
	for res in heatmap:
		num = float(res)
		f.write('%.2f ' % num)
	f.write('\n\n')

def sequence_alignment():
	# build and perform mafft
	command = 'mafft --retree 2 --inputorder --inputorder "%s/HeatmapAlignment/Sequences.txt" > "%s/HeatmapAlignment/Alignment.txt"' % (dir, dir)
	process = Popen(command, shell=True, stdout=PIPE, stderr=STDOUT)
	# write mafft output to log file, not stdout
	f = open('%s/HeatmapAlignment/Alignment.log' % dir, 'a')
	with process.stdout:
		for line in iter(process.stdout.readline, b''):
			f.write(line.decode("utf-8"))

def plot_heatmap_on_alignment():
	pass

######################
### Run everything ###
######################

# build directories and file structure

try:
	os.mkdir('%s/HeatmapAlignment' % dir)
except OSError as error:
	pass

f = open('%s/HeatmapAlignment/Sequences.txt' % dir,'w')
f = open('%s/HeatmapAlignment/Occupancy.txt' % dir,'w')
f = open('%s/HeatmapAlignment/Alignment.log' % dir,'w')

# define dictionary for AA codes - allows customisation

res_dict = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
	'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N', 
	'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W', 
	'ALA': 'A', 'VAL':'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M',
	'HSD': 'H', 'HSE': 'H', 'GLU0': 'E', 'ASP0': 'D'}

print('Loading Data')

# loop through inputs and get sequences and heatmaps

for i in np.arange(len(csvfile_list)):
	# get AA sequence and occupancy values
	sequence,heatmap,attribute = get_sequence(csvfile_list[i], 4)
	# write to file for later analysis/posterity
	write_fasta(sequence,csvfile_list[i])
	write_data(heatmap,attribute[0],csvfile_list[i])

print('Done. Files written to %s/HeatmapAlignment' % dir)
	
# run sequence alignment on all sequneces

print('Running alignment')

sequence_alignment()

print('Done. Alignment file written to %s/HeatmapAlignment' % dir)

# combine all into a lovely plot

print('Making plot')

plot_heatmap_on_alignment()

print('Done')
