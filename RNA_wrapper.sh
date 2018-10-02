#!/bin/bash

base_dir=/data/RNAseq/Broad_gtex
data_dir=$base_dir/data
Gencode_genomeDir=$base_dir/genome_index
transcript_model_dir=/data/RNAseq/MendelianRNA-seq
ID_list=$(ls $data_dir | grep 'fastq.gz' | cut -d'_' -f1 | uniq)
thread_N=3
set -x

echo 'samples included in this processing batch are' $ID_list
echo 'start STAR alignment'

for sample in $(echo $ID_list) ; do
fastq_input=$(ls $data_dir | grep $sample | paste -s -d, -)
echo 'samples included in this processing batch are' $ID_list
echo 'start' $sample 'Two pass mode processing'

STAR --genomeDir $Gencode_genomeDir \
--readFilesCommand zcat \
--twopassMode Basic \
--readFilesIn $fastq_input \
--outFileNamePrefix $sample. \
--sjdbOverhang 149 \
--alignSoftClipAtReferenceEnds No \
--runThreadN $thread_N \
--limitBAMsortRAM 38000000000 \
--outSAMtype BAM Unsorted \
--outSAMunmapped Within \
--outSAMmapqUnique 60 \
--outSAMattributes NH HI AS NM MD \
--outFilterMultimapNmax 20 \
--outFilterMismatchNmax 999 \
--alignIntronMin 20 \
--alignIntronMax 1000000 \
--alignMatesGapMax 1000000 \
--limitSjdbInsertNsj 3000000 \
--outFilterMismatchNoverReadLmax 0.05 \
--alignSJoverhangMin 8 \
--alignSJDBoverhangMin 1 \
--sjdbScore 1  \
--outFilterScoreMinOverLread 0.33 \
--outFilterMatchNminOverLread 0.33 \
--chimJunctionOverhangMin 15
done

# Sort BAM (use samtools to get around the memory gluttony of STAR) in parallel

parallel -j5 --no-notice \
        ./sort_bam.sh \
        ::: *.

# Index BAM in parallel

parallel -j5 --no-notice \
        ./index_bam.sh \
        ::: *.

# Mark Duplicates in parallel

parallel -j5 --no-notice \
        ./mark_dup.sh \
        ::: *.

#Mark duplicates
#for sample in $(echo $ID_list)
#do
#echo 'marking duplicates for' $sample
#java -Djava.io.tmpdir=/data/tmp -jar ~/bin/picard.jar MarkDuplicates \
#    INPUT=$data_dir/$sample.Aligned.sortedByCoord.out.bam \
#    OUTPUT=$data_dir/$sample.Aligned.sortedByCoord.out.patched.md.bam \
#    METRICS_FILE=data_dir/$sample.sorted.deduped.metrics \
#    CREATE_INDEX=TRUE \
#    SORTING_COLLECTION_SIZE_RATIO=0.1 \
#    ASSUME_SORTED=TRUE \
#    VALIDATION_STRINGENCY=LENIENT
#done

# RNA-SeQCin parallel
parallel -j5 --no-notice \
        ./RNASeQC.sh \
        ::: *.Deduped.Aligned.sortedByCoord.out.bam

#for sample in $(echo $ID_list)
#do
#echo 'RNA-SeQC for ' $sample
#python run_rnaseqc.py \
#    $data_dir/$sample.Aligned.sortedByCoord.out.patched.md.bam \
#    $base_dir/gencode.v19.annotation.patched_contigs.gtf \
#    $base_dir/rsem_reference/rsem_reference.transcripts.fa \
#    $sample \
#    --output_dir $data_dir
#done

# RSEM transcript quantification in parallel
parallel -j5 --no-notice \
        ./RSEM_quant.sh \
        ::: *.Deduped.Aligned.sortedByCoord.out.bam


#for sample in $(echo $ID_list)
#do
#echo 'RSEM quantification for ' $sample
#python run_RSEM.py \
#        $base_dir/rsem_reference \
#        $data_dir/star_out/$sample.Aligned.toTranscriptome.out.bam \
#        $data_dir/$sample \
#        --threads 4
#done

ls *.Deduped.Aligned.sortedByCoord.out.bam > bamlist.list

echo 'starting Splice Junction Discovery'
python3 $base_dir/SpliceJunctionDiscovery.py \
	-transcriptFile=$base_dir/transcript_file.bed \
	-bamList=bamlist.list \
	-processes=12

echo 'Starting Normalization of splice junction values'
python $base_dir/NormalizeSpliceJunctionValues.py \
	-transcript_model=$transcript_model_dir/gencode.comprehensive.splice.junctions.txt \
	-splice_file=$base_dir/All.transcript_file.bed.splicing.list \
	--normalize > $(date +"%Y%d%b%T").normalised.junctions
