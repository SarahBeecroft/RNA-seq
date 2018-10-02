#!/bin/bash
set -x
echo "$1 file"

samtools sort  -o $1.Aligned.sortedByCoord.out.bam $1.Aligned.out.bam