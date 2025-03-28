# SVLearn
A machine learning-based genotyping tool for structural variation of short reads

## Update (0.0.5)
In the latest **0.0.5** release, SVLearn has introduced the `--no-filter-overlaps` option to the `prepareAlt` subcommand. When building the ALT genome, this option disables the filtering of overlapping SVs, ensuring that every input SV is included in the genotyping analysis.


## Installation
SVLearn version 0.0.5 has undergone stability testing on Linux, including Ubuntu 20.04.3 and Red Hat Enterprise Linux 8.3 (Ootpa).

Installing SVLearn and the required Python environment using conda/mamba is convenient and efficient. Typically, the installation should not take more than half an hour.

### Requirements
1. python=3.9
2. pysam=0.22.0
3. polars=0.20.15
4. pandas=2.2.1
5. scikit-learn=1.3.0
6. pyfaidx
7. pyarrow
8. pybedtools
9. intervaltree

The following mapping software may need to be installed separately:

10. [bwa-mem2=2.2.1](https://github.com/bwa-mem2/bwa-mem2) 
11. [samtools>=1.17](https://github.com/samtools/samtools)
12. [sambamba>=1.0.1](https://github.com/biod/sambamba)


### Python environment
```
conda create -n svlearn python=3.9 bioconda::pysam=0.22.0 polars=0.20.15 pandas=2.2.1 scikit-learn=1.3.0 bioconda::pyfaidx pyarrow bioconda::pybedtools conda-forge::intervaltree
conda activate svlearn
```
or:
```
mamba create -n svlearn python=3.9 bioconda::pysam=0.22.0 polars=0.20.15 pandas=2.2.1 scikit-learn=1.3.0 bioconda::pyfaidx pyarrow bioconda::pybedtools conda-forge::intervaltree
mamba activate svlearn
```

### Build
Download the [Release](https://github.com/yangqimeng99/svlearn/releases) .

Ensure that `install.sh` is run within an activated conda environment.
```
tar xvzf svlearn-0.0.5.tar.gz
cd svlearn-0.0.5
wget https://github.com/Illumina/paragraph/releases/download/v2.4a/paragraph-v2.4a-binary.zip
unzip -q paragraph-v2.4a-binary.zip -d bin/paragraph-v2.4a
bash install.sh 
```
or:
```
git clone https://github.com/yangqimeng99/svlearn.git
cd svlearn
wget https://github.com/Illumina/paragraph/releases/download/v2.4a/paragraph-v2.4a-binary.zip
unzip -q paragraph-v2.4a-binary.zip -d bin/paragraph-v2.4a
bash install.sh
```

### Download the [trained model](https://doi.org/10.5281/zenodo.11144997)
Please select the corresponding coverage genotyping model to achieve the best genotyping results.


### Extra Requirements
SVLearn has provided available SV datasets and their corresponding SV features for [three species](https://doi.org/10.5281/zenodo.13309024). If you need to extract SV features from new SV datasets, the following additional tools need to be installed:
1. [RepeatMasker>=4.1.5](https://www.repeatmasker.org/RepeatMasker/)
2. [trf=4.09](https://tandem.bu.edu/trf/trf.html)
3. [GenMap=1.3.0](https://github.com/cpockrandt/genmap)
4. [BISER=1.4](https://github.com/0xTCG/biser)


## Test
Below are the steps to test whether SVLearn has been installed successfully.

If everything is working correctly, you will see the final output: `SVLearn test successful!`
```
cd test
bash test.sh
```

## Demo
We provide one sample dataset for each of the three species: human, cattle, and sheep, which can be used for demonstration and validation. You can download them from [here](https://doi.org/10.5281/zenodo.13309024). Before testing the demo, the [trained model](https://doi.org/10.5281/zenodo.11144997) should also be downloaded to the working directory.


Here is a demonstration using one example from each of the three species. Each demonstration is expected to complete in approximately 5 minutes using a single CPU core and 16GB of memory.

#### Demo 1: Human 24 feature
```
## 24 feature
svlearn genotype \
      -v 02.SV-set_Genomes/human_ref_sorted_format_filtered_sv.vcf \
      -m Human_30x_24feature_RandomForest_model.joblib \
      -s 04.validation-dataset/human/Human_SV_Set-sv_feature.tsv \
      -a 04.validation-dataset/human/HG002_30x/BreakPoint_ReadDepth_2Bam_feature.tsv \
      -p 04.validation-dataset/human/HG002_30x/para_feature.tsv \
      -o HG002_30_svlearn_genotype.vcf \
      -n HG002_30

svlearn benchmark \
      -b 04.validation-dataset/human/Human_SV_Set-HG002_True_GT.vcf \
      -c HG002_30_svlearn_genotype.vcf \
      -o HG002_30_svlearn.bench.tsv \
      -n1 HG002 -n2 HG002_30
```
The above demonstration will generate two files:

`HG002_30_svlearn_genotype.vcf`: The SV genotyping file output by SVLearn.

`HG002_30_svlearn.bench.tsv`: The genotyping performance of `HG002_30_svlearn_genotype.vcf` compared with the truth set, with the specific results as follows:
```
sv_set_number                 38613
genotyped_sv_number           37021
genotype_rate                 0.9588
accuracy_genotyped_sv_number  32956
precision_GT                  0.8292
recall_GT                     0.7585
f1_GT                         0.7922
precision                     0.9209
recall                        0.8424
f1                            0.8799
conc_00                       0.95
conc_01                       0.7534
conc_11                       0.8575
wgc                           0.8537
```
**Note:**
The `precision_GT`, `recall_GT`, and `f1_GT` metrics emphasize the genotyping performance when distinguishing between heterozygous and homozygous genotypes. 
The `precision`, `recall`, and `f1` metrics focus on the presence or absence of variants in the genotyping results, without distinguishing between heterozygous and homozygous variants.

#### Demo 2: Cattle 18 feature
```
## 18 feature
svlearn genotype \
      -v 02.SV-set_Genomes/cattle_ref_sorted_format_filtered_sv.vcf \
      -m Cattle_18feature_RandomForest_model.joblib \
      -s 04.validation-dataset/cattle/Cattle_SV_Set-sv_feature.tsv \
      -a 04.validation-dataset/cattle/Charolais_30x/BreakPoint_ReadDepth_2Bam_feature.tsv \
      -o Charolais_30_svlearn_genotype.vcf \
      -n Charolais_30

svlearn benchmark \
      -b 04.validation-dataset/cattle/Cattle_SV_Set-Charolais_True_GT.vcf \
      -c Charolais_30_svlearn_genotype.vcf \
      -o Charolais_30_svlearn.bench.tsv \
      -n1 Charolais -n2 Charolais_30
```
`Charolais_30_svlearn.bench.tsv`: The genotyping performance of `Charolais_30_svlearn_genotype.vcf` compared with the truth set, with the specific results as follows:
```
sv_set_number                 121435
genotyped_sv_number           120151
genotype_rate                 0.9894
accuracy_genotyped_sv_number  115362
precision_GT                  0.8619
recall_GT                     0.7924
f1_GT                         0.8257
precision                     0.9104
recall                        0.837
f1                            0.8721
conc_00                       0.9854
conc_01                       0.7866
conc_11                       0.8621
wgc                           0.878
```

#### Demo 3: Sheep 18 feature
```
## 18 feature
svlearn genotype \
      -v 02.SV-set_Genomes/sheep_ref_sorted_format_filtered_sv.vcf \
      -m Sheep_18feature_RandomForest_model.joblib \
      -s 04.validation-dataset/sheep/Sheep_SV_Set-sv_feature.tsv \
      -a 04.validation-dataset/sheep/Romanov_30x/BreakPoint_ReadDepth_2Bam_feature.tsv \
      -o Romanov_30_svlearn_genotype.vcf \
      -n Romanov_30

svlearn benchmark \
      -b 04.validation-dataset/sheep/Sheep_SV_Set-Romanov_True_GT.vcf \
      -c Romanov_30_svlearn_genotype.vcf \
      -o Romanov_30_svlearn.bench.tsv \
      -n1 Romanov -n2 Romanov_30
```
`Romanov_30_svlearn.bench.tsv`: The genotyping performance of `Romanov_30_svlearn_genotype.vcf` compared with the truth set, with the specific results as follows:
```
sv_set_number                 113042
genotyped_sv_number           110163
genotype_rate                 0.9745
accuracy_genotyped_sv_number  101289
precison_GT                   0.8773
recall_GT                     0.8189
f1_GT                         0.8471
precison                      0.9451
recall                        0.8822
f1                            0.9126
conc_00                       0.9658
conc_01                       0.8187
conc_11                       0.8928
wgc                           0.8924
```

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
Filter overlapping SV:
svlearn prepareAlt --ref_fasta ref.fasta --ref_sv_vcf sv.vcf --out 01.prepareAlt_output

Retain overlapping SVs:
svlearn prepareAlt --ref_fasta ref.fasta --ref_sv_vcf sv.vcf --out 01.prepareAlt_output --no-filter-overlaps
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

**Note:**
With `--no-filter-overlaps`: Overlapping SVs are retained, and files with `*_filtered_sv*` are not generated. Instead, `ref_sorted_format.vcf` and `alt.bed` serve as the primary inputs for downstream analysis.


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
In the `svlearn runParagraph`, the [Paragraph](https://github.com/Illumina/paragraph) is called twice to extract six paragraph features.

### 6. Genotyping
In this step, the feature matrixs obtained from previous steps is integrated, and the [trained model](https://doi.org/10.5281/zenodo.11144997) is called to obtain the SV genotypes.

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
        --ref_sv_vcf ../01.prepareAlt_output/ref_sorted_format_filtered_sv.vcf \
        --model model.joblib \
        --sv_feature ../02.SV.feature/sv_feature.tsv \
        --align_feature ../04.align_feature/sample1/BreakPoint_ReadDepth_2Bam_feature.tsv \
        --name sample1 \
        --out sample1_svlearn_18feature_genotype.vcf

# If you have extracted paragraph features, please call the 24-feature model
svlearn genotype \
        --ref_sv_vcf ../01.prepareAlt_output/ref_sorted_format_filtered_sv.vcf \
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

## Training new model
Download the [prepared training dataset](https://doi.org/10.5281/zenodo.13309024) for human, cattle, or sheep.

### 1. Prepare the training label file, tab-separated
training_sample.true.gt.tsv
```
sv_id   GT_true
INS.0   0/1
DEL.0   0/0
INS.1   1/1
.....   .....
```

### 2. Generate the training matrix

```
# 18 features
csvtk join -t -f 'sv_id' training_sample1.true.gt.tsv sv_feature.tsv training_sample1.BreakPoint_ReadDepth_2Bam_feature.tsv > training_sample1.18feature.training.tsv
csvtk join -t -f 'sv_id' training_sample2.true.gt.tsv sv_feature.tsv training_sample2.BreakPoint_ReadDepth_2Bam_feature.tsv > training_sample2.18feature.training.tsv

csvtk concat training_sample1.18feature.training.tsv training_sample2.18feature.training.tsv > sample1_2.18feature.training.tsv

# 24 features
csvtk join -t -f 'sv_id' training_sample1.true.gt.tsv sv_feature.tsv training_sample1.BreakPoint_ReadDepth_2Bam_feature.tsv training_sample1.para_feature.tsv > training_sample1.24feature.training.tsv
csvtk join -t -f 'sv_id' training_sample2.true.gt.tsv sv_feature.tsv training_sample2.BreakPoint_ReadDepth_2Bam_feature.tsv training_sample2.para_feature.tsv > training_sample2.24feature.training.tsv

csvtk concat training_sample1.24feature.training.tsv training_sample2.24feature.training.tsv > sample1_2.24feature.training.tsv
```

### 3. Training Model
```
# 18 features
svlearn trainingModel \
      --train_set sample1_2.18feature.training.tsv \
      --train_model RandomForest \
      --out RandomForest.18feature.model
      --threads 64

# 24 features
svlearn trainingModel \
      --train_set sample1_2.24feature.training.tsv \
      --other_feature paragraph \
      --train_model RandomForest \
      --out RandomForest.24feature.model
      --threads 64
```
**Output files:**
`RandomForest.18feature.model` and `RandomForest.24feature.model` folders contain the .joblib files, which are the generated genotyping models.

## Benchmark
SVLearn benchmark compares SV genotypes based on the unique IDs in the 3th columns of the true set and call set VCF files.

Download the [prepared validation dataset](https://doi.org/10.5281/zenodo.13309024) for human, cattle, or sheep.
```
svlearn benchmark
usage: svlearn benchmark [-h] -b file -c file -n1 str -n2 str [-o file]

-h, --help                     show this help message and exit
-b, --base_set file            True sv vcf set
-c, --call_set file            Call sv vcf set
-n1, --base_sample_name str    Sample name in true sv vcf set
-n2, --call_sample_name str    Sample name in call sv vcf set
-o, --out file                 The out of bechmark result, Default: benchmark.result.tsv
```

## Cite 
Yang, Q., Sun, J., Wang, X. *et al.* SVLearn: a dual-reference machine learning approach enables accurate cross-species genotyping of structural variants. *Nat Commun* **16**, 2406 (2025). https://doi.org/10.1038/s41467-025-57756-z
