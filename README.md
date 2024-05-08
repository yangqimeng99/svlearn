# SVLearn
A machine learning-based genotyping tool for structural variation of short reads

# Installation
## Requirements
The following dependency software needs to be installed:
1. [RepeatMasker>=4.1.5](https://www.repeatmasker.org/RepeatMasker/)
2. [trf=4.09](https://tandem.bu.edu/trf/trf.html)
3. [GenMap=1.3.0](https://github.com/cpockrandt/genmap)
4. [BISER=1.4](https://github.com/0xTCG/biser)

## Python environment
```
conda create -n svlearn python=3.9 pysam=0.22.0 polars=0.20.15 pandas=2.2.1 scikit-learn=1.3.0 pyfaidx pyarrow bioconda::pybedtools conda-forge::intervaltree
conda activate svlearn
```

## Build
