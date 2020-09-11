#!/bin/bash
rawdata1=../rawdata/scatac-seq/HGC20191230001-0003_lane7
rawdata2=../rawdata/scatac-seq/HGC20191230001-0003_lane8
name1=${rawdata1}/L7_md5sum.check.out
name2=${rawdata2}/L8_md5sum.check.out
metadata=../metadata/scatac-seq
reference=../metadata/reference/hg19_insert
result=../result/scatac-seq

pre_mapping()
{
    fastp -w 4 -c -h ${metadata}/${1}.fastp.html -j ${metadata}/${1}.fastp.json \
    -i ${2}/${1}_R1_001.fastq.gz -o ${metadata}/${1}_r1.fq.gz \
    -I ${2}/${1}_R2_001.fastq.gz -O ${metadata}/${1}_r2.fq.gz 
}

mapping()
{
    bowtie2 -p 20 -t --mm --very-sensitive -x ${reference} \
    -1 ${metadata}/${1}_r1.fq.gz \
    -2 ${metadata}/${1}_r2.fq.gz | samtools view -bS - > ${metadata}/${1}.bam
    samtools view -f 2 ${metadata}/${1}.bam | samtools sort -n -o ${metadata}/${1}_nsort.bam -
    samtools fixmate -rm ${metadata}/${1}_nsort.bam ${metadata}/${1}_nsort_fix.bam
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

get_coverage()
{
    bamCoverage -b ${metadata}/${1}_sort_rm.bam \
    -o ${metadata}/${1}.SeqDepthNorm.bw -e -p 20 \
    --binSize 25
}

single_cell()
{
    for line in `cat ${1}|grep HW`
    do
        if [[ ${line} =~ R1 ]]
        then
            echo ${line%_R1_001.fastq*}
            prefix=${line%_R1_001.fastq*}
            pre_mapping "${prefix}" "${2}"
            mapping "${prefix}"
            post_mapping "${prefix}"
            get_coverage "${prefix}"
        fi
    done
}

bulk_fastq_merge()
{
    ulimit -Sn 10000
    cat ${rawdata1}/*_R1_*.fastq.gz ${rawdata2}/*_R1_*.fastq.gz > ${metadata}/bulk_R1.fastq.gz
    cat ${rawdata1}/*_R2_*.fastq.gz ${rawdata2}/*_R2_*.fastq.gz > ${metadata}/bulk_R2.fastq.gz
}

bulk_bam_merge()
{
    ulimit -Sn 10000
    samtools merge ${metadata}/*HW*_sort.bam > ${metadata}/bulk_sort.bam
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

single_cell "${name1}" "${rawdata1}" > "${result}/L7_single_cell.log" 2>&1 &
single_cell "${name2}" "${rawdata2}" > "${result}/L8_single_cell.log" 2>&1 &
bulk_bam_merge > "${result}/bulk_bam_merge.log" 2>&1 &
post_mapping "bulk" > "${result}/bulk_post_mapping.log" 2>&1 
peak_calling "bulk" > "${result}/bulk_peak_calling.log" 2>&1
quality_report "bulk" > "${result}/quality_report.log" 2>&1 &
