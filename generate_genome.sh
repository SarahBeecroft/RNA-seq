#!/bin/bash
##generates genome indices for use with STAR fastq alignment.

STAR --runThreadN 4 --runMode genomeGenerate \
--genomeDir . \
--genomeFastaFiles GRCh37.p13.genome.fa
