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
 * `ref.fasta`: Reference Genome
 * `sv.vcf`   : SV Set

**Please note:** The input file `sv.vcf` should adhere to the following requirements in VCF 4.0+ format:
   (1) The VCF file should include complete REF and ALT allele sequences in 4th and 5th columns(like Indel).
   (2) The INFO field must contain SVTYPE.
   (3) Each SV should have a unique ID in 3th columns.

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
 * `ref.fasta` and `alt.fasta`
 * `ref_sorted_format_filtered_sv.vcf` and `alt_sorted_format_filtered_sv.bed`

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
svlearn svFeature \
        --ref_sv_vcf ../01.prepareAlt_output/ref_sorted_format_filtered_sv.vcf \ # output of 1. Create Alt Genome
        --alt_sv_bed ../01.prepareAlt_output/alt_sorted_format_filtered_sv.bed \ # output of 1. Create Alt Genome
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
`02.SV.feature/sv_feature.tsv`: The **SV feature** matrix file

**Note:**
The TRF2GFF.py script in the script/ directory comes from [Adam Taranto’s open-source project](https://github.com/Adamtaranto/TRF2GFF).

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
mkdir 03.mapping;cd 03.mapping;mkdir 01.RefGenome 02.AltGenome

# Index
cd 01.RefGenome
bwa-mem2.avx512bw index ref.fasta
cd ../02.AltGenome
bwa-mem2.avx512bw index alt.fasta

# Reads mapping, only one sample is shown as an example, but multiple samples can be mapped in parallel
cd 01.RefGenome
bash /software/svlearn/script/bwa_dedup.sh sample1_R1.fastq.gz sample1_R2.fastq.gz ref.fasta sample1
cd ../02.AltGenome
bash /software/svlearn/script/bwa_dedup.sh sample1_R1.fastq.gz sample1_R2.fastq.gz alt.fasta sample1_alt

cd ../..
```

**Output files:**
 * `03.mapping/01.RefGenome/sample1.dedup.sort.bam`: ref_bam
 * `03.mapping/02.AltGenome/sample1_alt.dedup.sort.bam`: alt_bam

### 4. Extract the Alignment feature
In this step, combine ref_bam and alt_bam to run `svlearn alignFeature` to obtain the alignment feature of each SV site for the genotyping sample.

**Input files:**
 * `ref.fasta` and `alt.fasta`
 * `ref_sorted_format_filtered_sv.vcf` and `alt_sorted_format_filtered_sv.bed`
 * `sample1.dedup.sort.bam` and `sample1_alt.dedup.sort.bam`

**Running:**
```
mkdir 04.align_feature;cd 04.align_feature

svlearn alignFeature \
        --ref_fasta ref.fasta \
        --alt_fasta alt.fasta \
        --ref_sv_vcf ../01.prepareAlt_output/ref_sorted_format_filtered_sv.vcf \
        --alt_sv_bed ../01.prepareAlt_output/alt_sorted_format_filtered_sv.bed \
        --ref_bam ../03.mapping/01.RefGenome/sample1.dedup.sort.bam \
        --alt_bam ../03.mapping/02.AltGenome/sample1_alt.dedup.sort.bam \
        --threads 1 \
        --out sample1

cd ..
```

**Output files:**
```
$tree 04.align_feature/sample1/
04.align_feature/sample1/
├── BP_2Bam_result.tsv
├── BreakPoint_ReadDepth_2Bam_feature.tsv
└── RD_2Bam_result.tsv
```
 * `BP_2Bam_result.tsv`: Alignment feature based on breakpoint information
 * `RD_2Bam_result.tsv`: Alignment feature based on read depth information
 * `BreakPoint_ReadDepth_2Bam_feature.tsv`: The **alignment feature** matrix of the sample1

**Note:**
In the `svlearn alignFeature`, the [mosdepth](https://github.com/brentp/mosdepth) is called to calculate read depth. The SVLearn package includes a binary distribution of [mosdepth 0.3.7](https://github.com/brentp/mosdepth/releases/tag/v0.3.7), and we haven’t made any changes to its code.

### 5. Extract the Paragraph feature (optional)
In this step, run another genotyping tool [paragraph](https://github.com/Illumina/paragraph) separately on the ref bam and alt bam to obtain the paragraph feature of each SV site for the genotyping sample.

**Input files:**
 * `ref.fasta` and `alt.fasta`
 * `ref_sorted_format_filtered_sv.vcf` and `alt_sorted_format_filtered_sv.vcf`
 * `sample1.dedup.sort.bam` and `sample1_alt.dedup.sort.bam`

**Running:**
```
mkdir 05.para_feature;cd 05.para_feature

svlearn runParagraph \
        --ref_fasta ref.fasta \
        --alt_fasta alt.fasta \
        --ref_sv_vcf ../01.prepareAlt_output/ref_sorted_format_filtered_sv.vcf \
        --alt_sv_vcf ../01.prepareAlt_output/alt_sorted_format_filtered_sv.vcf \
        --ref_bam ../03.mapping/01.RefGenome/sample1.dedup.sort.bam \
        --alt_bam ../03.mapping/02.AltGenome/sample1_alt.dedup.sort.bam \
        --threads 20 \
        --out sample1

cd ..
```

**Output files:**
```
$tree 05.para_feature/sample1/
05.para_feature/sample1/
├── alt.idxdepth.json
├── alt.idxdepth.json.log
├── alt.manifest
├── alt_paragraph_out
│   ├── genotypes.json.gz
│   ├── genotypes.vcf.gz
│   ├── grmpy.log
│   ├── variants.json.gz
│   └── variants.vcf.gz
├── para_feature.tsv
├── ref.idxdepth.json
├── ref.idxdepth.json.idxdepth.log
├── ref.manifest
└── ref_paragraph_out
    ├── genotypes.json.gz
    ├── genotypes.vcf.gz
    ├── grmpy.log
    ├── variants.json.gz
    └── variants.vcf.gz
```
 * `para_feature.tsv`: The **paragraph feature** matrix of the sample1
 * other files: Paragraph’s output

**Note:**
In the `svlearn runParagraph`, the [Paragraph](https://github.com/Illumina/paragraph) is called twice to extract six paragraph features. The SVLearn package includes a binary distribution of [Paragraph v2.4a](https://github.com/Illumina/paragraph/releases/tag/v2.4a), and we haven’t made any changes to its code.

### 6. Genotype
In this step, the feature matrixs obtained from previous steps is integrated, and the [trained model]() is called to obtain the SV genotypes.

**Input files:**
 * `model.joblib`: trained model file
 * `ref_sorted_format_filtered_sv.vcf`: output of **1. Create Alt Genome**
 * `sv_feature.tsv`: output of **2. Extracting each SV feature**
 * `BreakPoint_ReadDepth_2Bam_feature.tsv`: output of **4. Extract the Alignment feature**
 * (optional) `para_feature.tsv`: output of **5. Extract the Paragraph feature**

**Running:**
```
mkdir 06.genotype;cd 06.genotype

# If you haven’t extracted paragraph features, please call the 18-feature model
svlearn genotype \
        --model model.joblib \
        --sv_feature ../02.SV.feature/sv_feature.tsv \
        --align_feature ../04.align_feature/sample1/BreakPoint_ReadDepth_2Bam_feature.tsv \
        --name sample1 \
        --out sample1_svlearn_18feature_genotype.vcf

# If you have extracted paragraph features, please call the 24-feature model
svlearn genotype \
        --model model.joblib \
        --sv_feature ../02.SV.feature/sv_feature.tsv \
        --align_feature ../04.align_feature/sample1/BreakPoint_ReadDepth_2Bam_feature.tsv \
        --paragraph_feature ../05.para_feature/sample1/para_feature.tsv \
        --name sample1 \
        --out sample1_svlearn_24feature_genotype.vcf

cd ..
```

**Output files:**
`06.genotype/sample1_svlearn_18feature_genotype.vcf` and `06.genotype/sample1_svlearn_24feature_genotype.vcf`: SV genotyping result file for sample1 using different models.

## Advanced usage
### Training new model

### Benchmark
