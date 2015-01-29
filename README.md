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