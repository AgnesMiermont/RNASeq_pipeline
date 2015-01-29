Set of scripts for RNA-Seq data processing, in particular differential expression analysis


# Requirements

R packages and software:

- R >= 3.1.2, use
```
/cluster/project8/vyp/vincent/Software/R-3.1.2/bin/R
```
- DESeq  >= 1.18
- DEXSeq >= 1.10.8
- BiocParallel >= 1.0.1 

Notes for installation of DEXSeq:

- Need to add /share/apps/binutils-2.25/bin to the PATH
- type 
```
scl enable devtoolset-1.1 'bash'
```
and then run R

# Input format

The key input is a table in plain text format that contains one row per sample, with the header line:
```
sample f1 f2 condition
```
Assuming that condition is what you want to run the differential expression analysis on.
You also need to specify the input folder, so that the fastq can be found at ${iFolder}/${f1} and ${iFolder}/${f2}.
Note that f1 and f2 can specify subfolders themselves. Also a species parameter and and output folder.