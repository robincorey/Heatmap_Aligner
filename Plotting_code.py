import string
table = str.maketrans('', '', string.ascii_lowercase)
def reorder_alignment_for_plot(num):
    new_sequence = []
    with open('test.txt', "r") as ifile:
        count = 0
        for line in ifile:
            if 'sequence' in line:
                new_sequence.append(line.translate(table).strip('\n '))
                count=count+1
    #print(new_sequence)
    return count/num, new_sequence
#print(alignment)

lines, new_sequence = reorder_alignment_for_plot(2)

def get_occupancy(sequence,lines, i):
    residue_count = 0
    occupancy = []
    with open('Occupancy.txt','r') as ifile:
        for line in ifile:
            if '>' not in line:
                if '0' in line:
                    occupancy.append(line.split( ))
    #print(sequence,occupancy[i]) # need to make occ a list 
    heatmap = []
    for aa,char in enumerate(sequence):
        if char != '-':
            #print(char,occupancy[0][residue_count])
            residue_count = residue_count+1
            heatmap.append(occupancy[0][residue_count])
        else:
            heatmap.append('NA') # maybe set to 0
    return heatmap
