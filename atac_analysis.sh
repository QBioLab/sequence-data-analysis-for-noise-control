#!/bin/bash

rawdata=../rawdata/atac-seq
metadata=../metadata/atac-seq
reference=../metadata/reference/hg19_insert
result=../result/atac-seq

fastq_merge()
{
    cat ${rawdata}/${1}_L3_R1.fq.gz ${rawdata}/${1}_L4_R1.fq.gz > ${metadata}/${1}_R1.fastq.gz
    cat ${rawdata}/${1}_L3_R2.fq.gz ${rawdata}/${1}_L4_R2.fq.gz > ${metadata}/${1}_R2.fastq.gz
}

pre_mapping()
{
    fastp -w 4 -c -h ${metadata}/${1}.fastp.html -j ${metadata}/${1}.fastp.json \
    -i ${metadata}/${1}_R1.fastq.gz -o ${metadata}/${1}_r1.fq.gz \
    -I ${metadata}/${1}_R2.fastq.gz -O ${metadata}/${1}_r2.fq.gz 
}

mapping()
{
    bowtie2 -p 20 -t --mm --very-sensitive -x ${reference} \
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
    macs2 callpeak -t ${metadata}/${1}_sort_rm.bam -n ${1} --outdir ${metadata} -g hs -f BAM --nomodel --shift -100 --extsize 200 -B --SPMR
    bedGraphToBigWig ${metadata}/${1}_treat_pileup.bdg ../metadata/reference/hg19_insert.chromsize ${metadata}/${1}.bw
}

quality_report()
{
    ataqv --threads 20 --peak-file ${metadata}/${1}_peaks.narrowPeak human ${metadata}/${1}_sort_mdup.bam --metrics-file ${metadata}/${1}.ataqv.json.gz
    mkarv ${metadata}/${1} ${metadata}/${1}.ataqv.json.gz
}

scale_factor()
{
    multiBamSummary bins -b ${metadata}/*_sort_rm.bam --smartLabels -p 20 \
    --scalingFactors ${metadata}/scalingFactors.txt -o ${metadata}/readCounts.npz 
}

get_coverage()
{
    bamCoverage -b ${metadata}/high_sort_rm.bam \
     -o ${metadata}/high.SeqDepthNorm.bw -e -p 5 --scaleFactor 1.0903 \
     --binSize 25
    bamCoverage -b ${metadata}/low_sort_rm.bam \
     -o ${metadata}/low.SeqDepthNorm.bw -e -p 5 --scaleFactor 1.1836 \
     --binSize 25
    bamCoverage -b ${metadata}/A1_sort_rm.bam \
     -o ${metadata}/A1.SeqDepthNorm.bw -e -p 5 --scaleFactor 0.7303 \
     --binSize 25
}

for prefix in high low A1; do 
{
    fastq_merge "${prefix}" > "${result}/${prefix}.fastq_merge.log" 2>&1
    pre_mapping "${prefix}" > "${result}/${prefix}.pre_mapping.log" 2>&1
    mapping "${prefix}" > "${result}/${prefix}.mapping.log" 2>&1
    post_mapping "${prefix}" > "${result}/${prefix}.post_mapping.log" 2>&1
    peak_calling "${prefix}" > "${result}/${prefix}.peak_calling.log" 2>&1
    quality_report "${prefix}" > "${result}/${prefix}.quality_report.log" 2>&1
} & done

scale_factor > "${result}/scale_factor.log" 2>&1
get_coverage > "${result}/get_coverage.log" 2>&1
