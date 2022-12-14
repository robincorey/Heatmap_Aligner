About
====

A Python script to take the output of multiple PyLipID runs (csv), and plot onto a sequence alignment.

Pseudo-code:

1. Get sequence and occupancy/duration info from PyLipID
2. Perform sequence alignment with mafft.
3. Somehow format the occupancy/duration data to match the sequence alignment. Looping over sequence probably.
4. Plot the occupancy/duration data as a grid. Label cells with sequence data.

Tried biopython, biotite, dash. Very hard to see how to use custom colours. I think it has to be manual.
