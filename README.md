About
====

A Python script to take the output of multiple PyLipID runs (csv), and plot onto a sequence alignment.

Currently most of the way there.

Required dependencies
====

NumPy
Matplotlib

Current structure
====

Part 1 in Heatmap_Aligner.py

1. Read in single file from PyLipID and extract sequence and Occupancy (duration easy to add)
2. Write sequence as a fasta, and write occupancy as separate file.
3. Perform sequence alignment with mafft.

Need to check that this still works as intended.

Part 1 currently in Plotting_code.py. Still in progress, but close. Have tried a lot of things here.
1.  Reorder alignment for plot - i.e. take from fasta to a plotting-friendly format
2.  get the occuancy from the file, and reformat to match the above - including gaps.
3.  plot occupancy vs residue number
4.  add labels which is the sequnce data

Still needs to be finalised and tested for 1 or >2 sequences. Plotting code being tested in test_plotting.ipynb
