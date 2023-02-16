# from ipynb, need to bring here

import string
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import pandas as pd 

# these bits are set earlier in the script
dsys = ['YQED' , 'YNGC']
count = 2 ### THIS NEEDS TO BE SET FROM FIRST PART OF SCRIPT

# reads alignment into a NumPy array

table = str.maketrans('', '', string.ascii_lowercase)
def reorder_alignment_for_plot(num):
    new_sequence = []
    with open('Alignment.txt', "r") as ifile:
        # this count keeps track of total lines
        count = 0
        # define one data series for each sequences used
        for line_number,line in enumerate(ifile):
            if '>' not in line:
                # this bit simply removes non-sequence chars from the line
                #line = (line.translate(table).strip('\n '))
                # save line with revised line number
                #line_array = [ line, count ]
                #new_sequence.append(line.translate(table).strip('\n '))
                new_sequence.append(line.strip('\n '))
                #print(new_sequence)'''
                count=count+1
    # count/num gives the lines per sequence
    return count/num, new_sequence
#print(alignment)

lines, new_sequence = reorder_alignment_for_plot(count)

# reformatting occupancy file
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
            sequence = new_sequence[int(line)+((int(lines)-1)*system)+system]
            f = open('Occupancy_reformatted.txt','a')
            for aa,char in enumerate(sequence):
                if char != '-':
                    f.write('%s ' % occupancy[system][residue_count]) 
                    residue_count = residue_count+1
                else:
                    f.write(' -1 ')
            f.write('\n')

occupancy_to_new_file()

import warnings
warnings.filterwarnings("ignore")

def get_occupancy_reordered(sequence, system, alignment_line, occupancy_count_in):
    #print(sequence, system, alignment_line, occupancy_count_in)
    occupancy = []
    with open('Occupancy_reformatted.txt','r') as ifile:
        x = ifile.readlines()[alignment_line]
        occupancy.append(x.split( ))
    heatmap = []
    sequence_array = []
    occupancy_count = occupancy_count_in
    residue_count = occupancy_count_in
    for aa,char in enumerate(sequence):
        if char != '-1':
            heatmap.append(float(occupancy[0][aa]))
            sequence_array.append(char)
            occupancy_count = occupancy_count+1
            residue_count = residue_count+1
        else:
            heatmap.append(-1)
            sequence_array.append(char)
    return heatmap, sequence_array, occupancy_count

def plot_sequence(heatmap_array,alignment_line,line, sys):
    #takes an array of single letter AA codes and residue attribute values
    axs[line].set_yticklabels(['',sys])
    axs[line].set_xticklabels(np.arange(-10,60,step=10))
    #axs[line].set_xticks([])
    im = axs[line].imshow(heatmap_array, cmap='Reds', vmin=0, vmax=100 )
    
def add_labels(sequence_array,alignment_line,line):
    for x,y in enumerate(np.arange(len(sequence_array))):
        axs[line].text(x-0.25,0.25, sequence_array[x], fontsize=10)
        
#def add_ticks(offset,line):
#    for x,y in enumerate(np.arange(0,60,step=10)):
#        axs[line].text(y-0.4,2.5, y+offset, fontsize=10)

def plot_each_line():
    # for each line, make an array of sequences and plot
    a_line = 0
    for system in np.arange(0,count): #count):
        heatmap_array = []
        occupancy_count = 0 # remove
        #alignment_line = int(line)*count+system #line*count+system
        alignment_line = int(line)+((int(lines)-1)*system)+system
        #print(line, alignment_line, new_sequence[alignment_line], system, occupancy_count)
        heatmap_values, sequence_array, occupancy_count = get_occupancy_reordered(
            new_sequence[alignment_line], system, alignment_line, occupancy_count)
        heatmap_array.append(heatmap_values)
        plot_sequence(heatmap_array, alignment_line, a_line, sys[system])
        add_labels(sequence_array,alignment_line, a_line)
        a_line = a_line + 1

        # for each line, make an array of sequences and plot
#a = np.array(range(int(lines)))
#occupancy_count = 0 

for line in np.arange(0,1): 
    fig = matplotlib.pyplot.gcf()
    fig.set_size_inches(18, 1)
    gs = fig.add_gridspec(count,1,hspace=0)
    axs = gs.subplots(sharex=True)
    plot_each_line()

