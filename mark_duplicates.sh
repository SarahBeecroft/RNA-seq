#!/bin/bash

java -Djava.io.tmpdir=/data/tmp -jar ~/bin/picard.jar MarkDuplicates \
    INPUT=$sample_id.bam \
    OUTPUT=$sample_id.sorted.deduped.bam \
    METRICS_FILE=$sample_id.sorted.deduped.metrics \
    CREATE_INDEX=TRUE \
    SORTING_COLLECTION_SIZE_RATIO=0.1 \
    ASSUME_SORTED=TRUE \
    VALIDATION_STRINGENCY=LENIENT
