#!/bin/bash
base_dir=/data/RNAseq/Broad_gtex
data_dir=$base_dir/data

set -x
 echo "$1 file"
  python3 run_rnaseqc.py \
  $1.Deduped.Aligned.sortedByCoord.out.bam \
  $base_dir/gencode.v19.annotation.patched_contigs.gtf \
  $base_dir/rsem_reference/rsem_reference.transcripts.fa \
  $1 \
  --output_dir $data_dir
