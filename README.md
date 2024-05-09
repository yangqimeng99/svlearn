# SVLearn
A machine learning-based genotyping tool for structural variation of short reads

## Installation
### Requirements
The following dependency software needs to be installed:
1. [RepeatMasker>=4.1.5](https://www.repeatmasker.org/RepeatMasker/)
2. [trf=4.09](https://tandem.bu.edu/trf/trf.html)
3. [GenMap=1.3.0](https://github.com/cpockrandt/genmap)
4. [BISER=1.4](https://github.com/0xTCG/biser)
5. [bwa-mem2=2.2.1](https://github.com/bwa-mem2/bwa-mem2)
6. [samtools>=1.17](https://github.com/samtools/samtools)
7. [sambamba>=1.0.1](https://github.com/biod/sambamba)

### Python environment
```
conda create -n svlearn python=3.9 pysam=0.22.0 polars=0.20.15 pandas=2.2.1 scikit-learn=1.3.0 pyfaidx pyarrow bioconda::pybedtools conda-forge::intervaltree
conda activate svlearn
```

### Build
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

### Download the trained model


## Usage
Before starting the SVLearn workflow, please ensure that all Requirements and SVLearn are configured in your environment.
Please adjust the **Input files** paths according to the actual situation.

### 1. Create Alt Genome
In this step, the input VCF file will be formatted and used to generate an alternative genome relative to the reference genome.

**Input files:**
```
ref.fasta # Reference Genome
sv.vcf    # SV Set
```
Please note that the input file `sv.vcf` should adhere to the following requirements in VCF 4.0+ format:
  * The VCF file should include complete REF and ALT allele sequences in 4th and 5th columns(like Indel).
  * The INFO field must contain SVTYPE.
  * Each SV should have a unique ID in 3th columns.

**Running:**
```
svlearn prepareAlt --ref_fasta ref.fasta --ref_sv_vcf sv.vcf --out 01.prepareAlt_output
```

**Output files:**
```
$tree 01.prepareAlt_output
01.prepareAlt_output/
├── alt.fasta
├── alt.fasta.fai
├── alt_sorted_format_filtered_sv.bed
├── alt_sorted_format_filtered_sv.vcf
├── ref_sorted_format_close_sv.vcf
├── ref_sorted_format_filtered_sv.vcf
├── ref_sorted_format.vcf
└── ref_sorted.vcf
```
 * `ref_sorted.vcf`: The output of the input `sv.vcf` after sorting
 * `ref_sorted_format.vcf`: Further formatting of `ref_sorted.vcf`
 * `ref_sorted_format_filtered_sv.vcf`: Filtered SV set used for subsequent analysis
 * `ref_sorted_format_close_sv.vcf`: SV set that has been filtered out
 * `alt.fasta`: Alternative genome created based on `ref_sorted_format_filtered_sv.vcf` and `ref.fasta`
 * `alt_sorted_format_filtered_sv.vcf`: SV records based on the alt genome, changing DEL to INS, INS to DEL
 * `alt_sorted_format_filtered_sv.bed`: Positions of SV’s ALT sequences in the alt genome

### 2. Extracting each SV feature
In this step, four external software tools are employed separately to extract the features of each SV from both the `ref.fasta` and `alt.fasta`. Afterward, the svFeature module of svlearn is used to integrate these extracted features.

**Input files:**
`ref.fasta` and `alt.fasta`

**Running:**
```
mkdir 02.SV.feature; cd 02.SV.feature

# 01.RepeatMasker, Please adjust the `-pa` and `-species` parameters based on the actual circumstances.
RepeatMasker -pa 64 -engine ncbi -species human -xsmall -s -no_is -cutoff 255 -frag 20000 -dir ./ -gff ref.fasta
RepeatMasker -pa 64 -engine ncbi -species human -xsmall -s -no_is -cutoff 255 -frag 20000 -dir ./ -gff alt.fasta

# 02.trf
trf ref.fasta 2 7 7 80 10 50 500 -f -d -h
trf alt.fasta 2 7 7 80 10 50 500 -f -d -h
python /software/svlearn/script/TRF2GFF.py -d ref.fasta.2.7.7.80.10.50.500.dat -o ref.trf.gff
python /software/svlearn/script/TRF2GFF.py -d alt.fasta.2.7.7.80.10.50.500.dat -o alt.trf.gff

# 03.GenMap
bash /software/svlearn/script/genmap.sh ref.fasta ref.fasta.genmapK50E1
bash /software/svlearn/script/genmap.sh alt.fasta alt.fasta.genmapK50E1

# 04.BISER, ref.fasta.masked and alt.fasta.masked are the softmasked FASTA files obtained from the output of 01.RepeatMasker
biser -o ref.fasta.masked.out -t 2 --gc-heap 32G ref.fasta.masked
biser -o alt.fasta.masked.out -t 2 --gc-heap 32G alt.fasta.masked

# 05.svlearn svFeature
svlearn svFeature 
        --ref_sv_vcf ../01.prepareAlt_output/ref_sorted_format_filtered_sv.vcf \
        --alt_sv_bed ../01.prepareAlt_output/alt_sorted_format_filtered_sv.bed \
        --ref_rm ref.fasta.out \ # 01.RepeatMasker output in ref.fasta
        --alt_rm alt.fasta.out \ # 01.RepeatMasker output in alt.fasta
        --ref_trf ref.trf.gff \
        --alt_trf alt.trf.gff \
        --ref_genmap ref.fasta.genmapK50E1.txt \
        --alt_genmap alt.fasta.genmapK50E1.txt \
        --ref_biser ref.fasta.masked.out \
        --alt_biser alt.fasta.masked.out \
        --out sv_feature.tsv

cd ..
```
**Output files:**
`02.SV.feature/sv_feature.tsv`: The feature matrix of the SV set

### 3. Reads mapping
In this step, each sample’s short reads (fastq files) are separately linearly aligned to ref.fasta and alt.fasta to obtain each sample’s ref_bam and alt_bam.

**Input files:**
```
# alignment reference
ref.fasta
alt.fasta

# samples short reads
sample1_R1.fastq.gz
sample1_R2.fastq.gz

sample2_R1.fastq.gz
sample2_R2.fastq.gz
...
```

**Running:**
```
mkdir 03.mapping;cd 03.mapping
mkdir 01.RefGenome 02.AltGenome


# Index
cd 01.RefGenome
bwa-mem2.avx512bw index ref.fasta
cd ../02.AltGenome
bwa-mem2.avx512bw index alt.fasta

# Reads mapping, only one sample is shown as an example, but multiple samples can be mapped in parallel
cd 01.RefGenome
bash /software/svlearn/script/genmap.sh ref.fasta ref.fasta.genmapK50E1





```


