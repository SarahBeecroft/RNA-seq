#!/bin/bash
set -x
echo "$1 file"
 java '-Djava.io.tmpdir=/data/tmp' -jar ~/bin/picard.jar MarkDuplicates \
        INPUT=$1 \
        OUTPUT=Deduped.$1 \
        METRICS_FILE=$1.MarkDup.metrics \
        CREATE_INDEX=TRUE \
        SORTING_COLLECTION_SIZE_RATIO=0.1 \
        ASSUME_SORTED=TRUE \
        VALIDATION_STRINGENCY=LENIENT
