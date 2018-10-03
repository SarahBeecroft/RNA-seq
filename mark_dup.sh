#!/bin/bash
set -x
echo "$1 file"
 java '-Djava.io.tmpdir=/data/tmp' -jar ~/bin/picard.jar MarkDuplicates \
        INPUT=$1 \
        OUTPUT=$1.Deduped.sortedByCoord.Aligned.out.bam \
        METRICS_FILE=$1.sorted.deduped.metrics \
        CREATE_INDEX=TRUE \
        SORTING_COLLECTION_SIZE_RATIO=0.1 \
        ASSUME_SORTED=TRUE \
        VALIDATION_STRINGENCY=LENIENT
