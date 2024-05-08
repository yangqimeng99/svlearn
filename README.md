# SVLearn
A machine learning-based genotyping tool for structural variation of short reads

# Installation
## Requirements
The following dependency software needs to be installed:
1. [RepeatMasker>=4.1.5](https://www.repeatmasker.org/RepeatMasker/)
2. [trf=4.09](https://tandem.bu.edu/trf/trf.html)
3. [GenMap=1.3.0](https://github.com/cpockrandt/genmap)
4. [BISER=1.4](https://github.com/0xTCG/biser)
5. [bwa-mem2=2.2.1](https://github.com/bwa-mem2/bwa-mem2)
6. [samtools>=1.17](https://github.com/samtools/samtools)
7. [sambamba>=1.0.1](https://github.com/biod/sambamba)

## Python environment
```
conda create -n svlearn python=3.9 pysam=0.22.0 polars=0.20.15 pandas=2.2.1 scikit-learn=1.3.0 pyfaidx pyarrow bioconda::pybedtools conda-forge::intervaltree
conda activate svlearn
```

## Build
Download the [Release](https://github.com/yangqimeng99/svlearn/releases) 
```
tar xvzf svlearn-0.0.1.tar.gz
cd svlearn-0.0.1
bash install.sh
```
or:
```
git clone https://github.com/yangqimeng99/svlearn.git
cd svlearn
bash install.sh
```


