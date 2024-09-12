#!/bin/bash
../bin/svlearn prepareAlt --ref_fasta ref.fasta --ref_sv_vcf raw.vcf --out test_output

file1_a="expected_result/alt.fasta"
file1_b="test_output/alt.fasta"

file2_a="expected_result/ref_sorted_format_filtered_sv.vcf"
file2_b="test_output/ref_sorted_format_filtered_sv.vcf"

file3_a="expected_result/alt_sorted_format_filtered_sv.bed"
file3_b="test_output/alt_sorted_format_filtered_sv.bed"

if diff "$file1_a" "$file1_b" > /dev/null && \
   diff "$file2_a" "$file2_b" > /dev/null && \
   diff "$file3_a" "$file3_b" > /dev/null; then
    echo "SVLearn test successful!"
else
    echo "The test output is inconsistent with the expected result, please check."
fi
