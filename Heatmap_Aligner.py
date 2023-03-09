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
from PIL import Image
import warnings
warnings.filterwarnings("ignore")

####
# Script to take multiple PyLipID outputs (csv) and create a sequence alignment coloured by lipid statistics, currently just occupancy
# Can be used with a single PyLipID output (csv) to make a sequence-based heatmap
#
# Note - mafft needed: conda install -c bioconda mafft
#
# Written by Robin Corey
# Current to-do:
#   - add options for other metrics other than occupancy, like duration and residence time
#   - test on a different env (i.e. other user)
#   - residue numbers perhaps - add to end of line
#   - ensure visual consistency between different numbers of files/image scaling
####

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Script to take multiple PyLipID outputs and create a sequence alignment coloured by lipid statistics', formatter_class=argparse.RawTextHelpFormatter)
    group_input = parser.add_argument_group('INPUT arguments')
    group_input.add_argument("-i", "--input", nargs='+', metavar='filename', type=str, required=True, help='Paths to input files, e.g. ./data/A.csv ./data/B.csv ./data/C.csv')
    group_input.add_argument("-d", "--dir", type=str, required=False, help='Directory where jobs to be run. Default = current', default='.')
    group_input.add_argument("-n", "--names_list", nargs='+', metavar='filename', type=str, required=True, help='Names of systems, e.g. ABCD DEFG HIJK')
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
    f = open('%s/HeatmapAlignment/%s.txt' % (args.dir, attribute), 'a')
    f.write('> %s from %s\n' % (attribute, csvfile))
    for res in heatmap:
        num = float(res)
        f.write('%.2f ' % num)
    f.write('\n\n')

# this run sequence alignment between all of the input sequences 
# mafft is needed here conda install -c bioconda mafft
def sequence_alignment():
    # build and perform mafft
    command = 'mafft --retree 2 --inputorder --inputorder "%s/HeatmapAlignment/Sequences.txt" > \
		 "%s/HeatmapAlignment/Alignment.txt"' % (args.dir, args.dir)
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
    with open('%s/HeatmapAlignment/Alignment.txt' % args.dir, "r") as ifile:
        # this count keeps track of total lines
        count_lines = 0
        # define one data series for each sequences used
        for line_number,line in enumerate(ifile):
            if '>' not in line:
                new_sequence.append(line.strip('\n '))
                count_lines=count_lines+1
    return count_lines/num, new_sequence

# this reformats the file with attribute date to match the alignment file
# writes a new file - need to check conversion from "occupancy" works
def attribute_to_new_file(count, attribute):
    # this list called occupancy due to legacy
    occupancy = []
    with open('%s/HeatmapAlignment/%s.txt' % (args.dir, attribute),'r') as ifile:
            for line in ifile:
                if '>' not in line:
                    if '0' in line:
                        occupancy.append(line.split( ))
    f = open('%s/HeatmapAlignment/%s_reformatted.txt' % (args.dir, attribute) ,'w')
    for system in np.arange(0,count):
        residue_count = 0
        for line in np.array(range(int(lines))):
            sequence = new_sequence[int(line)+((int(lines)-1)*system)+system]
            f = open('%s/HeatmapAlignment/%s_reformatted.txt' % (args.dir, attribute) ,'a')
            for aa,char in enumerate(sequence):
                if char != '-':
                    f.write('%s ' % occupancy[system][residue_count]) 
                    residue_count = residue_count+1
                else:
                    f.write(' -1 ')
            f.write('\n')
    max_occ = np.max(occupancy)
    max_occ = np.array(max_occ, dtype=np.float32)
    return np.max(max_occ)

# this reorders the attribute file as a prelude to plotting, and stores in a heatmap
def get_occupancy_reordered(sequence, system, alignment_line, occupancy_count_in, attribute):
    occupancy = []
    with open('%s/HeatmapAlignment/%s_reformatted.txt' % (args.dir, attribute) ,'r') as ifile:
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
def plot_sequence(heatmap_array,alignment_line,line, sys, max_occ ):
    #takes an array of single letter AA codes and residue attribute values
    axs[line].set_yticklabels(['', sys])
    axs[line].set_yticks([-1,0])
    axs[line].set_xticklabels('')
    im = axs[line].imshow(heatmap_array, cmap=args.cmap, vmin=0, vmax=max_occ+15)
    return im
    # here, could add numbers to the end of each line. Interesting idea.

def add_labels(sequence_array,alignment_line,line):
    for x,y in enumerate(np.arange(len(sequence_array))):
        axs[line].text(x-0.3,0.3, sequence_array[x], fontsize=8) # was 0.25 and 10 for 2, hm
    # work out how to optimise the spacing and fontsize.

def plot_each_line(count):
    # for each line, make an array of sequences and plot
    a_line = 0
    for system in np.arange(0,count):
        heatmap_array = []
        occupancy_count = 0 # currentl not doing anything, but might be useful
        alignment_line = int(line)+((int(lines)-1)*system)+system
        heatmap_values, sequence_array, occupancy_count = get_occupancy_reordered(
            new_sequence[alignment_line], system, alignment_line, occupancy_count, attribute [0])
        heatmap_array.append(heatmap_values)
        im = plot_sequence(heatmap_array, alignment_line, a_line, args.names_list[system] ,max_occ+15)
        add_labels(sequence_array,alignment_line, a_line)
        a_line = a_line + 1
    return len(sequence_array), im

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

# define dictionary for AA codes - this allows for customisation if needed

res_dict = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
    'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N', 
    'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W', 
    'ALA': 'A', 'VAL':'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M',
    'HSD': 'H', 'HSE': 'H', 'GLU0': 'E', 'ASP0': 'D'}

print('Loading Data')

# loop through inputs and get sequences and heatmaps

for i in np.arange(len(args.input)):
    # get AA sequence and occupancy values
    # Here 2 = occupancy, can change
    sequence, heatmap, attribute = get_sequence(args.input[i], 2)
    # write to file for later analysis/posterity
    write_fasta(sequence,args.input[i])
    write_data(heatmap,attribute[0],args.input[i])

print('Done. Files written to %s/HeatmapAlignment' % args.dir)
    
# run sequence alignment on all sequences

print('Running alignment')

sequence_alignment()

print('Done. Alignment file written to %s/HeatmapAlignment' % args.dir)

# reformat data for plotting

print('Reformating files for plotting')

lines, new_sequence = reorder_alignment_for_plot(len(args.input))
max_occ = attribute_to_new_file(len(args.input),attribute[0])

print('Done. Reformatted files written to %s/HeatmapAlignment' % args.dir)

# combine all into a lovely plot

print('Making plot')

# this makes a separate plot for each line of the alignment
list_im = []
for line in np.arange(0,int(lines)): 
    fig = matplotlib.pyplot.gcf()
    #fig.set_size_inches(72/(len(args.input)), 1) # check scalability here 
    fig.set_size_inches(12, 1)
    gs = fig.add_gridspec(len(args.input),1,hspace=0) #lines
    axs = gs.subplots(sharex=True)
    fig.patch.set_facecolor('w')
    # call the main plotting function
    line_length, im = plot_each_line(len(args.input))
    pad = (line_length/60)
    if pad == 1.0:
        plt.subplots_adjust(left=0.1, right=0.9, top=0.9, bottom=0.2)
    else:
        plt.subplots_adjust(left=0.1, right=0.9*(pad), top=0.9, bottom=0.2)
    plt.savefig('%s/HeatmapAlignment/alignment_line.%s.png' % (args.dir, line ) )
    list_im.append('%s/HeatmapAlignment/alignment_line.%s.png' % (args.dir, line ))
    plt.close()

# plot colorbar separately
a = np.array([[0,max_occ+15]])
#plt.figure(figsize=(4, 0.9))
img = plt.imshow(a, cmap=args.cmap)
plt.gca().set_visible(False)
#cax = plt.axes([0.1, 0.1, 0.6, 0.2])
plt.colorbar(orientation='horizontal') #, cax=cax)
plt.savefig('%s/HeatmapAlignment/colorbar.png' % (args.dir) )

# finally, combine the individial plots into a final plot
imgs = [ Image.open(i) for i in list_im ]
min_shape = sorted( [(np.sum(i.size), i.size ) for i in imgs])[0][1]
imgs_comb = np.vstack([i.resize(min_shape) for i in imgs])
imgs_comb = np.vstack(imgs)
imgs_comb = Image.fromarray(imgs_comb)
imgs_comb.save( '%s/HeatmapAlignment/full_alignment.png' % args.dir )

print('Plot finished, see %s/HeatmapAlignment/full_alignment.png' % args.dir ) 
print('Color bar written separately, see %s/HeatmapAlignment/colorbar.png' % args.dir )
