#!/bin/bash
##generates genome indices for use with STAR fastq alignment. assumed 150bp unpaired reads 

data_dir=/data/RNAseq/Broad_gtex

##if nor already done- need to do this step first:
#zcat gencode.v19.annotation.gtf.gz | \
#    sed 's/chrM/chrMT/;s/chr//' > gencode.v19.annotation.patched_contigs.gtf

STAR \
--runMode genomeGenerate \
--genomeDir /data/RNAseq/Broad_gtex/genome_index \
--genomeFastaFiles /data/RNAseq/Broad_gtex/Homo_sapiens_assembly19.fasta \
--sjdbGTFfile /data/RNAseq/Broad_gtex/gencode.v19.annotation.patched_contigs.gtf \
--sjdbOverhang 149 \
--runThreadN 4

rsem-prepare-reference \
$data_dir/Homo_sapiens_assembly19.fasta \
$data_dir/rsem_reference/rsem_reference \
--gtf $data_dir/gencode.v19.annotation.patched_contigs.gtf \
--num-threads 4
