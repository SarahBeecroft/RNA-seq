#!/bin/bash
Gencode_genomeDir=/data/RNAseq/data/gencode_genome
thread_N=4
script_dir=/data/RNAseq/MendelianRNA-seq/Analysis
processes_N=12
ID_list=$(ls | grep 'fastq.gz' | cut -d'_' -f1 | uniq)
transcript_model_dir=/data/RNAseq/MendelianRNA-seq

echo 'samples included in this processing batch are' $ID_list
echo 'start STAR processing'
for sample in $(echo $ID_list)
do 
fastq_input=$(ls | grep $sample | paste -s -d, -)
echo 'start' $sample 'Two pass mode processing'
STAR --genomeDir $Gencode_genomeDir \
--readFilesCommand zcat \
--twopassMode Basic \
--readFilesIn $fastq_input \
--sjdbGTFfile /data/RNAseq/data/gencode_genome/gencode.v19.annotation.gtf \
--outFileNamePrefix $sample. \
--sjdbOverhang 149 \
--alignSoftClipAtReferenceEnds No \
--runThreadN $thread_N \
--limitBAMsortRAM 40000000000 \
--outSAMtype BAM SortedByCoordinate \
--outSAMunmapped Within \
--outSAMmapqUnique 60 \
--outSAMattributes NH HI AS NM MD \
--outFilterMultimapNmax 20 \
--outFilterMismatchNmax 999 \
--alignIntronMin 20 \
--alignIntronMax 1000000 \
--alignMatesGapMax 1000000 \
--limitSjdbInsertNsj 3000000 \
-outFilterMismatchNoverReadLmax 0.05 \
--alignSJoverhangMin 8 \
--alignSJDBoverhangMin 1 \
--sjdbScore 1  \
--outFilterScoreMinOverLread 0.33 \
--outFilterMatchNminOverLread 0.33 \
--chimJunctionOverhangMin 15
done

echo 'stage mark_duplicates'

#for sample in $(echo $ID_list)
#do
#	echo 'start' $sample 'mark duplicates processing'
#	java '-Djava.io.tmpdir=/data/tmp' -jar ~/bin/picard.jar MarkDuplicates \
#	INPUT=$sample.Aligned.sortedByCoord.out.bam \
#	OUTPUT=$sample.sorted.deduped.bam \
#	METRICS_FILE=$sample.sorted.deduped.metrics \
#	CREATE_INDEX=TRUE \
#	SORTING_COLLECTION_SIZE_RATIO=0.1 \
#	ASSUME_SORTED=TRUE \
#	VALIDATION_STRINGENCY=LENIENT
#done

echo $ID_list | parallel java '-Djava.io.tmpdir=/data/tmp' -jar ~/bin/picard.jar MarkDuplicates \
	INPUT={.}.Aligned.sortedByCoord.out.bam \
	OUTPUT={.}.sorted.deduped.bam \
	METRICS_FILE={.}.sorted.deduped.metrics \
	CREATE_INDEX=TRUE \
	SORTING_COLLECTION_SIZE_RATIO=0.1 \
	ASSUME_SORTED=TRUE \
	VALIDATION_STRINGENCY=LENIENT

ls *.deduped.bam | grep '' > bamlist.list

python3 $script_dir/SpliceJunctionDiscovery.py -transcriptFile=transcript_file.bed -bamList=bamlist.list -processes=$processes_N

python $script_dir/NormalizeSpliceJunctionValues.py -transcript_model=transcript_model_dir/gencode.comprehensive.splice.junctions.txt -splice_file=All.transcript_file.bed.splicing.list --normalize > $(date +"%Y%d%b%T").normalised.junctions
