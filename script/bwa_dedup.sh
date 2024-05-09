#!/bin/sh

R1=$1
R2=$2
Ref=$3
outprefix=$4

mkdir -p ./tmp/${outprefix}
bwa-mem2.avx512bw mem -t 24 -R "@RG\tID:${outprefix}\tLB:${outprefix}\tSM:${outprefix}\tPL:ILLUMINA" \
	${Ref} ${R1} ${R2} | samtools sort --threads 24 -T ./tmp/${outprefix} -o ${outprefix}.bam

samtools index ${outprefix}.bam
rm -rf ./tmp/${outprefix}

sambamba markdup -r -t 48 --tmpdir ./tmp/ ${outprefix}.bam ${outprefix}.dedup.sort.bam
rm -rf ${outprefix}.bam ${outprefix}.bam.bai
