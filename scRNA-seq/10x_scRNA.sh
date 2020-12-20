#!/bin/bash

rawdata1=../../rawdata/scRNA_seq_2020_11_1
rawdata2=../../rawdata/scRNA_seq_2020_11_7
reference=../../rawdata/reference/HeLa_s1_32_and_mm10
result1=../../result/scRNA_seq_2020_11_1
result2=../../result/scRNA_seq_2020_11_7
cellranger=/work/bio-chenr/software/cellranger-3.1.0/cellranger

primary_analysis_1()
{
    echo "process will start at:"
    date
    ${cellranger} count --id=${1} \
                   --transcriptome=${reference} \
                   --fastqs=${rawdata1}/${1} \
                   --expect-cells=13000 \
                   --nosecondary \
                   --r1-length=28 \
                   --jobmode=local --localcores=40 --localmem=150
    echo "process end at:"
    date
}

primary_analysis_2()
{
    echo "process will start at:"
    date
    ${cellranger} count --id=${1} \
                   --transcriptome=${reference} \
                   --fastqs=${rawdata2}/${1} \
                   --expect-cells=13000 \
                   --nosecondary \
                   --jobmode=local --localcores=40 --localmem=150
    echo "process end at:"
    date
}

# process technology repeat 1
for prefix in A1A2 C1C2 B1A3; do 
{
    primary_analysis_1 "${prefix}" > "${result1}/${prefix}.primary_analysis.log" 2>&1
} & done

# process technology repeat 2
for prefix in A1A2 C1C2 B1A3; do 
{
    primary_analysis_2 "${prefix}" > "${result2}/${prefix}.primary_analysis.log" 2>&1
} & done
