#!/bin/bash
set -x
echo "$1 file"

samtools index $1.Aligned.sortedByCoord.out.bam