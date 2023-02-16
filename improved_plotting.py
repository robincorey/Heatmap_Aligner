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
