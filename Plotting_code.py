import string
import numpy as np
import matplotlib.pyplot as plt

#reorder alignment
table = str.maketrans('', '', string.ascii_lowercase)
def reorder_alignment_for_plot(num):
    new_sequence = []
    with open('test.txt', "r") as ifile:
        # this count keeps track of total lines
        count = 0
        for line in ifile:
            if 'sequence' in line:
                # this bit simply removes non-sequence chars from the line
                new_sequence.append(line.translate(table).strip('\n '))
                count=count+1
    # count/num gives the lines per sequence
    return count/num, new_sequence
#print(alignment)

lines, new_sequence = reorder_alignment_for_plot(2)

def get_occupancy(sequence,lines, system, occupancy_count_in):
    occupancy = []
    with open('Occupancy.txt','r') as ifile:
        for line in ifile:
            if '>' not in line:
                if '0' in line:
                    occupancy.append(line.split( ))
    heatmap = []
    sequence_array = []
    # to keep track of progress through occupancy file
    occupancy_count = occupancy_count_in
    residue_count = occupancy_count_in
    for aa,char in enumerate(sequence):
        if char != '-':
            heatmap.append(float(occupancy[system][residue_count]))
            sequence_array.append(char)
            occupancy_count = occupancy_count+1
            residue_count = residue_count+1
        else:
            heatmap.append(-1)
            sequence_array.append(char)
    return heatmap, sequence_array, occupancy_count

def plot_sequence(text_array_1,text_array_2,heatmap_array,alignment_line,line):
    #takes an array of single letter AA codes and residue attribute values
    #fig.set_size_inches(30, 6)
    axs[line].set(yticks=[0,1],label=['A','B'])
    im = axs[line].imshow(heatmap_array, cmap='Reds', vmin=0, vmax=100 )
    for x,y in enumerate(np.arange(len(text_array_1))):
        axs[line].text(x-0.25,alignment_line-0.75, text_array_1[x], fontsize=15)
    for x,y in enumerate(np.arange(len(text_array_2))):
        axs[line].text(x-0.25,alignment_line+0.25, text_array_2[x], fontsize=15)
    #print(sequence_array,heatmap_values)
    #print(np.mean(heatmap_values))
    #plt.plot(np.arange(len(heatmap_values)),heatmap_values)

count = 2 ### THIS NEEDS TO BE SET FROM FIRST PART OF SCRIPT
lines = 3

fig = matplotlib.pyplot.gcf()
fig.set_size_inches(20, 6)
gs = fig.add_gridspec(lines,1, hspace=0)
axs = gs.subplots(sharex=True)
    
# for each line, make an array of sequences and plot
for line in np.arange(0,int(lines)): 
    # loop through systems, i.e. input files
    heatmap_array = []
    text_array = []
    for system in np.arange(0,count):
        occupancy_count = 0 # this will be wrong
        alignment_line = line*count+system
        heatmap_values, sequence_array, occupancy_count = get_occupancy(
            new_sequence[alignment_line], lines, system, occupancy_count )
        heatmap_array.append(heatmap_values)
        text_array.append(sequence_array)
    plot_sequence(text_array[0],text_array[1],heatmap_array,alignment_line, line)
    
# Need to resolve: 
# alignment of text with image
# space between plots
# more than 2 seqs
# occupancy counting for new order

