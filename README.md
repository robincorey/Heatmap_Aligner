About
====

A Python script to take the output of multiple PyLipID runs (csv), and plot onto a sequence alignment. Currently runnable with test data for any number of seqences, still needs testing on different datasets/systems

Produces several txt output files as part of the plotting process, as well as png files for each line of the alignment ('alignment_line.X.png') and the full alignment ('full_alignment.png'). A colorbar is plotted separately ('colorbar.png'). 

Required dependencies
====

MAFFT: 
```
conda install -c bioconda mafft
```

To run using test data
====

Download Heatmap_Aligner.py to local system. To see the options, run:

```
python Heatmap_Aligner.py -h
```

Then download test data. Make sure MAFFT installed before running on data.

Then run Heatmap_Aligner.py using Python 3

```
python Heatmap_Aligner.py -i test_data/YNGC/Dataset.csv test_data/YABI/Dataset.csv test_data/YQED/Dataset.csv -n YNGC YABI YQED
```

Note that -i and -n have to be defined. Directory (-d) and colour scheme (-c) are optional variable.

Next steps
====
* check documention
* consider CI
* consider adding optional residue numbers
