#!/usr/bin/env python3

import os
import string
import matplotlib
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import csv
import re
from subprocess import Popen, PIPE, STDOUT

'''
Script to take multiple PyLipID outputs (csv) and create a sequence alignment coloured by lipid statistics, currently just occupancy

Can be used with a single PyLipID output (csv) to make a sequence-based heatmap

Written by Robin Corey
'''

if __name__ == '__main__':
        parser = argparse.ArgumentParser(description='Script to take multiple PyLipID outputs and create a sequence alignment coloured by lipid statistics', formatter_class=argparse.RawTextHelpFormatter)
        group_input = parser.add_argument_group('INPUT arguments')
        group_input.add_argument("-i", "--input", nargs='+', metavar='filename', type=str, required=True, help='Paths to input files')
        group_input.add_argument("-d", "--dir", metavar='filename', type=str, required=False, help='Directory where jobs to be run. Default = current', default='./')
        group_input.add_argument("-n", "--names_list", nargs='+', metavar='filename', type=str, required=False, help='Names of systems', default=[])

        group_output = parser.add_argument_group('OUTPUT arguments')
        group_output.add_argument("-o", "--output", metavar='filename', required=False, type=str, help="Output path. Default = ./Heatmap_Aligner", default='./Heatmap_Aligner')

        # Parse arguments from command line
        args = parser.parse_args()

#######################################
### Don't touch anything below here ###
#######################################

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

# read alignment file into NumPy array
table = str.maketrans('', '', string.ascii_lowercase)
def reorder_alignment_for_plot(num):
    new_sequence = []
    with open('test.txt', "r") as ifile:
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
                #print(new_sequence)'''
                count=count+1
    # count/num gives the lines per sequence
    return count/num, new_sequence

# reformat occupancy file
def occupancy_to_new_file ():
    occupancy = []
    with open('Occupancy.txt','r') as ifile:
            for line in ifile:
                if '>' not in line:
                    if '0' in line:
                        occupancy.append(line.split( ))
    f = open('Occupancy_reformatted.txt','w')
    for system in np.arange(0,count):
        residue_count = 0
        for line in np.array(range(int(lines))):
            sequence = new_sequence[system+(line*count)]
            f = open('Occupancy_reformatted.txt','a')
            for aa,char in enumerate(sequence):
                if char != '-':
                    f.write('%s ' % occupancy[system][residue_count])
                    residue_count = residue_count+1
                else:
                    f.write(' -1 ')
            f.write('\n')

def get_occupancy_reordered(sequence, system, alignment_line, occupancy_count_in):
    occupancy = []
    with open('Occupancy_reformatted.txt','r') as ifile:
        for line in ifile:
                occupancy.append(line.split( ))
    heatmap = []
    sequence_array = []
    occupancy_count = occupancy_count_in
    residue_count = occupancy_count_in
    for aa,char in enumerate(sequence):
        if char != '-1':
            heatmap.append(float(occupancy[alignment_line][aa]))
            sequence_array.append(char)
            occupancy_count = occupancy_count+1
            residue_count = residue_count+1
        else:
            heatmap.append(-1)
            sequence_array.append(char)
    return heatmap, sequence_array, occupancy_count

def plot_sequence(heatmap_array,alignment_line,line, sys1, sys2):
    #takes an array of single letter AA codes and residue attribute values
    axs[line].set_yticklabels(['',sys1,sys2])
    axs[line].set_xticklabels('')
    im = axs[line].imshow(heatmap_array, cmap='Reds', vmin=0, vmax=100 )

def add_labels(text_array_1,text_array_2,alignment_line,line):
    for x,y in enumerate(np.arange(len(text_array_1))):
        axs[line].text(x-0.25,0.25, text_array_1[x], fontsize=10)
    for x,y in enumerate(np.arange(len(text_array_2))):
        axs[line].text(x-0.25,(1.25), text_array_2[x], fontsize=10)

def add_ticks(offset,line):
    for x,y in enumerate(np.arange(0,60,step=10)):
        axs[line].text(y-0.4,2.5, y+offset, fontsize=10)

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

# reformat data for plotting

print('Reformating files for plotting')

lines, new_sequence = reorder_alignment_for_plot(count)
occupancy_to_new_file()

print('Done. Reformatted files written to %s/HeatmapAlignment' % dir)

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
        print(new_sequence[alignment_line], system, occupancy_count)
        heatmap_values, sequence_array, occupancy_count = get_occupancy_reordered(
            new_sequence[alignment_line], system, a_line, occupancy_count)
        heatmap_array.append(heatmap_values)
        text_array.append(sequence_array)
        a_line = a_line + 1
        print(heatmap_array)
        plot_sequence(heatmap_array, alignment_line, line, sys1, sys2)
    add_labels(text_array[0],text_array[1],alignment_line, line)

print('Plot finished, exiting program')
