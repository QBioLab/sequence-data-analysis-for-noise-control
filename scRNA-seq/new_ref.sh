#!/bin/bash

cat mRuby-2A-bla.fa >> genome.fa
cat mRuby-2A-bla.gtf >> genes.gtf
cellranger mkref --genome=HeLa_s1_32 --fasta=genome.fa --genes=genes.gtf