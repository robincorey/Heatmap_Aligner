About
====

A Python script to take the output of multiple PyLipID runs (csv), and plot onto a sequence alignment.

Currently runnable with the test data, still needs testing on different datasets/systems

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
* add suggested new features
