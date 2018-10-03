#!/bin/bash
set -x
echo "$1 file"

samtools sort  -o Sorted.$1 $1
