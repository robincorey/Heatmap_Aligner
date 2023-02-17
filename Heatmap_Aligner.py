#!/usr/bin/env python3

import os
import string
import matplotlib
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import csv
import re
import argparse
from subprocess import Popen, PIPE, STDOUT

####
# Script to take multiple PyLipID outputs (csv) and create a sequence alignment coloured by lipid statistics, currently just occupancy
# Can be used with a single PyLipID output (csv) to make a sequence-based heatmap
#
# Note - mafft needed: conda install -c bioconda mafft
#
# Written by Robin Corey
# Current to-do:
#   - add options for other metrics other than occupancy
#   - test on a different env (i.e. other user)
####

#######################################
### Don't touch anything below here ###
#######################################

if __name__ == '__main__':
        parser = argparse.ArgumentParser(description='Script to take multiple PyLipID outputs and create a sequence alignment coloured by lipid statistics', formatter_class=argparse.RawTextHelpFormatter)
        group_input = parser.add_argument_group('INPUT arguments')
        group_input.add_argument("-i", "--input", nargs='+', metavar='filename', type=str, required=True, help='Paths to input files')
        group_input.add_argument("-d", "--dir", type=str, required=False, help='Directory where jobs to be run. Default = current', default='./')
        group_input.add_argument("-n", "--names_list", nargs='+', metavar='filename', type=str, required=False, help='Names of systems', default=[])
	group_input.add_argument("-c", "--cmap", type=str, required=False, help='Which matplotlib color range to use. Default = Reds', default='Reds')

        group_output = parser.add_argument_group('OUTPUT arguments')
        group_output.add_argument("-o", "--output", metavar='filename', required=False, type=str, help="Output path. Default = ./Heatmap_Aligner", default='./Heatmap_Aligner')

        # Parse arguments from command line
        args = parser.parse_args()


########################
### Define functions ###
########################

# this gets the AA sequence from the input csv files
def get_sequence(csvfile,column):
	# read sequence
	with open (csvfile) as inf:
		reader = csv.reader(inf, delimiter=",")
		sequence = list(zip(*reader))[0] 
	# read occupancy > can be any attribute, although not supported later in script just yet
	with open (csvfile) as inf:
		reader = csv.reader(inf, delimiter=",")
		heatmap = list(zip(*reader))[column]
	return sequence[1:],heatmap[1:], heatmap[:1]

# this coverts the AA sequence to fasta file, which it writes to a new file
def write_fasta(sequence,csvfile):
	f = open('%s/HeatmapAlignment/Sequences.txt' % args.dir, 'a')
	f.write('> sequence from %s\n' % csvfile) 
	for res in sequence:
		# strip residue number
		output = re.sub(r'\d+', '', res)
		f.write(res_dict[output]) 
	f.write('\n\n')

# this gets the occupancy data from the csv and writes as text file
# attribute here can be changed to other features - to implement
def write_data(heatmap,attribute,csvfile):
	f = open('%s/HeatmapAlignment/%s.txt' % attribute, 'a')
	f.write('> %s from %s\n' % (attribute, csvfile))
	for res in heatmap:
		num = float(res)
		f.write('%.2f ' % num)
	f.write('\n\n')

# this run sequence alignment between all of the input sequences 
# mafft is needed here conda install -c bioconda mafft
def sequence_alignment():
	# build and perform mafft
	command = 'mafft --retree 2 --inputorder --inputorder "%s/HeatmapAlignment/Sequences.txt" > "%s/HeatmapAlignment/Alignment.txt"' % (args.dir, args.dir)
	process = Popen(command, shell=True, stdout=PIPE, stderr=STDOUT)
	# write mafft output to log file, not stdout
	f = open('%s/HeatmapAlignment/Alignment.log' % args.dir, 'a')
	with process.stdout:
		for line in iter(process.stdout.readline, b''):
			f.write(line.decode("utf-8"))

# functions after here initially worked through in ipynb notebook

# this reads the newly poducted alignment file into NumPy array

table = str.maketrans('', '', string.ascii_lowercase)

def reorder_alignment_for_plot(num):
    new_sequence = []
    with open('%s/HeatmapAlignment/Sequences.txt' % args.dir, "r") as ifile:
        # this count keeps track of total lines
        count = 0
        # define one data series for each sequences used
        for line_number,line in enumerate(ifile):
            if 'sequence' in line:
                # this bit simply removes non-sequence chars from the line
                line = (line.translate(table).strip('\n '))
                # save line with revised line number
                line_array = [ line, count ]
                new_sequence.append(line.translate(table).strip('\n '))
                count=count+1
    # count/num gives the lines per sequence
    return count/num, new_sequence

# this reformats the file with attribute date to match the alignment file
# writes a new file - need to check conversion from "occupancy" works
def attribute_to_new_file ():
    # this list called occupancy due to legacy
    occupancy = []
    with open('%s/HeatmapAlignment/%s.txt' % attribute,'r') as ifile:
            for line in ifile:
                if '>' not in line:
                    if '0' in line:
                        occupancy.append(line.split( ))
    f = open('%s/HeatmapAlignment/%s_reformatted.txt' % attribute ,'w')
    for system in np.arange(0,count):
        residue_count = 0
        for line in np.array(range(int(lines))):
            sequence = new_sequence[int(line)+((int(lines)-1)*system)+system]
            f = open('%s/HeatmapAlignment/%s_reformatted.txt','a')
            for aa,char in enumerate(sequence):
                if char != '-':
                    f.write('%s ' % occupancy[system][residue_count]) 
                    residue_count = residue_count+1
                else:
                    f.write(' -1 ')
            f.write('\n')

# this reorders the attribute file as a prelude to plotting, and stores in a heatmap
def get_occupancy_reordered(sequence, system, alignment_line, occupancy_count_in):
    occupancy = []
    with open('%s/HeatmapAlignment/%s_reformatted.txt','r') as ifile:
        for line in ifile:
                occupancy.append(line.split( ))
    heatmap = []
    sequence_array = []
    occupancy_count = occupancy_count_in
    for aa,char in enumerate(sequence):
        if char != '-1':
            heatmap.append(float(occupancy[alignment_line][aa]))
            sequence_array.append(char)
            occupancy_count = occupancy_count+1
        else:
            heatmap.append(-1)
            sequence_array.append(char)
    return heatmap, sequence_array, occupancy_count

# this plots the attribute data as a heatmap 
def plot_sequence(heatmap_array,alignment_line,line, sys1, sys2):
    #takes an array of single letter AA codes and residue attribute values
    axs[line].set_yticklabels(['',sys])
    axs[line].set_yticks([-1,0])
    #axs[line].set_xticklabels(np.arange(-10,60,step=10))
    axs[line].set_xticklabels('')
    #axs[line].set_xticks([])
    im = axs[line].imshow(heatmap_array, cmap=args.cmap, vmin=0)

def add_labels(sequence_array,alignment_line,line):
    for x,y in enumerate(np.arange(len(sequence_array))):
        axs[line].text(x-0.3,0.3, sequence_array[x], fontsize=8) # was 0.25 and 10 for 2, hm

######################
### Run everything ###
######################

# build args.directories and file structure

try:
	os.mkdir('%s/HeatmapAlignment' % args.dir)
except OSError as error:
	pass

f = open('%s/HeatmapAlignment/Sequences.txt' % args.dir,'w')
f = open('%s/HeatmapAlignment/Occupancy.txt' % args.dir,'w')
f = open('%s/HeatmapAlignment/Alignment.log' % args.dir,'w')

# define dictionary for AA codes - allows customisation

res_dict = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
	'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N', 
	'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W', 
	'ALA': 'A', 'VAL':'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M',
	'HSD': 'H', 'HSE': 'H', 'GLU0': 'E', 'ASP0': 'D'}

print('Loading Data')

# loop through inputs and get sequences and heatmaps

for i in np.arange(len(args.input)):
	# get AA sequence and occupancy values
	sequence,heatmap,attribute = get_sequence(args.input[i], 4)
	# write to file for later analysis/posterity
	write_fasta(sequence,args.input[i])
	write_data(heatmap,attribute[0],args.input[i])

print('Done. Files written to %s/HeatmapAlignment' % args.dir)
	
# run sequence alignment on all sequneces

print('Running alignment')

sequence_alignment()

print('Done. Alignment file written to %s/HeatmapAlignment' % args.dir)

# reformat data for plotting

print('Reformating files for plotting')

lines, new_sequence = reorder_alignment_for_plot(count)
occupancy_to_new_file()

print('Done. Reformatted files written to %s/HeatmapAlignment' % args.dir)

# combine all into a lovely plot

print('Making plot')

fig = matplotlib.pyplot.gcf()
fig.set_size_inches(lines*(count+1), 6)
gs = fig.add_gridspec(int(lines),1, hspace=0)
axs = gs.subplots(sharex=True)

a = np.array(range(int(lines)))
occupancy_count = 0 
a_line = 0
for system in np.arange(0,count):
    for line in np.arange(2):
        heatmap_array = []
        text_array = []
        alignment_line = line*count+system
        heatmap_values, sequence_array, occupancy_count = get_occupancy_reordered(
            new_sequence[alignment_line], system, a_line, occupancy_count)
        heatmap_array.append(heatmap_values)
        text_array.append(sequence_array)
        a_line = a_line + 1
        plot_sequence(heatmap_array, alignment_line, line, sys1, sys2)
    add_labels(text_array[0],text_array[1],alignment_line, line)

plt.savefig('%s/HeatmapAlignment/Heatmap.png' % args.dir, bbox_inches='tight', dpi=600 )
print('Plot finished, exiting program')
