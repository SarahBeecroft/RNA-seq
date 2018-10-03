#!/bin/bash
set -x
echo "$1 file"

samtools sort  -o sortedByCoord.$1 $1
