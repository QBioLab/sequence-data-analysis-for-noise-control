#!/bin/bash

build_reference()
{
    python ./hg19_insert.py
    bowtie2-build ../metadata/reference/hg19_insert.fa ../metadata/reference/hg19_insert
    cut -f1,2 hg19_insert.fa.fai > hg19_insert.chromsize
}

build_reference > "../result/build-reference.log" 2>&1 &