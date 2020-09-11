#!/bin/bash

rawdata=../rawdata/AJCHS2200107003
metadata=../metadata/chip-seq-2
reference=../metadata/reference/hg19_insert
result=../result/chip-seq-2

pre_mapping()
{
    fastp -w 4 -c -h ${metadata}/${1}.fastp.html -j ${metadata}/${1}.fastp.json \
    -i ${rawdata}/${1}_R1.fastq.gz -o ${metadata}/${1}_r1.fq.gz \
    -I ${rawdata}/${1}_R2.fastq.gz -O ${metadata}/${1}_r2.fq.gz 
}

mapping()
{
    bowtie2 -p 5 -t --mm -x ${reference} \
    -1 ${metadata}/${1}_r1.fq.gz \
    -2 ${metadata}/${1}_r2.fq.gz | samtools view -bS - > ${metadata}/${1}.bam
    samtools sort -n -o ${metadata}/${1}_nsort.bam ${metadata}/${1}.bam
    samtools fixmate -m ${metadata}/${1}_nsort.bam ${metadata}/${1}_nsort_fix.bam
    samtools sort -o ${metadata}/${1}_fix_sort.bam ${metadata}/${1}_nsort_fix.bam
    samtools index ${metadata}/${1}_fix_sort.bam
}

post_mapping()
{
    picard MarkDuplicates I=${metadata}/${1}_fix_sort.bam \
    O=${metadata}/${1}_sort_mdup.bam M=${metadata}/${1}_mdup_metrics.txt
    samtools view -F 1804 -f 2 -q 1 -u ${metadata}/${1}_sort_mdup.bam | samtools sort  \
    -o ${metadata}/${1}_rmdup_sort.bam -
    samtools index ${metadata}/${1}_rmdup_sort.bam
    samtools view -h ${metadata}/${1}_rmdup_sort.bam | grep -v -P '\tchrM\t' | samtools view \
    -b - > ${metadata}/${1}_sort_rm.bam
    samtools index ${metadata}/${1}_sort_rm.bam
}

peak_calling()
{
    macs2 callpeak -t ${metadata}/${1}_sort_rm.bam -n ${1} --outdir ${metadata} -g hs -f BAMPE -B --SPMR
    bedGraphToBigWig ${metadata}/${1}_treat_pileup.bdg ../metadata/reference/hg19_insert.chromsize ${metadata}/${1}.bw
}

quality_report()
{
    ataqv --threads 5 --peak-file ${metadata}/${1}_peaks.narrowPeak human ${metadata}/${1}_sort_mdup.bam --metrics-file ${metadata}/${1}.ataqv.json.gz
    mkarv ${metadata}/${1} ${metadata}/${1}.ataqv.json.gz
}

cross_correlation()
{
    Rscript phantompeakqualtools/run_spp_nodups.R -rf -c=${metadata}/${1}_sort_rm.bam \
    -savp=${metadata}/${1}_cross_correlation.pdf \
    -savd=${metadata}/${1}_cross_correlation.Rdata \
    -out=${metadata}/${1}_cross_correlation.txt
}

scale_factor()
{
    multiBamSummary bins -b ${metadata}/*_sort_rm.bam --smartLabels -p 20 \
    --scalingFactors ${metadata}/scalingFactors.txt -o ${metadata}/readCounts.npz 
}

get_coverage()
{
    bamCoverage -b ${metadata}/1-H3K27AC_combined_sort_rm.bam \
     -o ${metadata}/1-H3K27AC_combined.SeqDepthNorm.bw -e -p 5 --scaleFactor 1.8695 \
     --binSize 25
    bamCoverage -b ${metadata}/2-H3K27AC_combined_sort_rm.bam \
     -o ${metadata}/2-H3K27AC_combined.SeqDepthNorm.bw -e -p 5 --scaleFactor 0.8715 \
     --binSize 25
    bamCoverage -b ${metadata}/3-H3K27AC_combined_sort_rm.bam \
     -o ${metadata}/3-H3K27AC_combined.SeqDepthNorm.bw -e -p 5 --scaleFactor 1.3995 \
     --binSize 25
}
 
for prefix in 1-H3K27AC_combined 2-H3K27AC_combined 3-H3K27AC_combined; do 
{
    pre_mapping "${prefix}" > "${result}/${prefix}.pre_mapping.log" 2>&1
    mapping "${prefix}" > "${result}/${prefix}.mapping.log" 2>&1
    post_mapping "${prefix}" > "${result}/${prefix}.post_mapping.log" 2>&1
    peak_calling "${prefix}" > "${result}/${prefix}.peak_calling.log" 2>&1
    quality_report "${prefix}" > "${result}/${prefix}.quality_report.log" 2>&1
    cross_correlation "${prefix}" > "${result}/${prefix}.cross_correlation.log" 2>&1
} & done

scale_factor > "${result}/scale_factor.log" 2>&1
get_coverage > "${result}/get_coverage.log" 2>&1
