#!/bin/bash
##wrapper around Broad's run_STAR.py script downloaded from their git repo and used to process the Gtex samples

base_dir=/data/RNAseq/Broad_gtex
data_dir=$base_dir/data
Gencode_genomeDir=$base_dir/genome_index
ID_list=$(ls $data_dir | grep 'fastq.gz' | cut -d'_' -f1 | uniq)


echo 'samples included in this processing batch are' $ID_list
echo 'start STAR alignment'
for sample in $(echo $ID_list)
do

fastq_input=$(ls $data_dir | grep $sample | paste -s -d, -)

python3 $base_dir/run_STAR.py \
        $data_dir/$fastq_input \
        $sample \
        --threads 3 \
        --output_dir $data_dir/star_out \
        && mv /tmp/star_out /data/star_out
done

#Mark duplicates in parallel
ID_list=$(echo $ID_list)
parallel -j5 --no-notice \
        ./mark_dup.sh \
        ::: $ID_list

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
ID_list=$(echo $ID_list)
parallel -j5 --no-notice \
        ./RNASeQC.sh \
        ::: $ID_list

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
ID_list=$(echo $ID_list)
parallel -j5 --no-notice \
        ./RNASeQC.sh \
        ::: $ID_list


#for sample in $(echo $ID_list)
#do
#echo 'RSEM quantification for ' $sample
#python run_RSEM.py \
#        $base_dir/rsem_reference \
#        $data_dir/star_out/$sample.Aligned.toTranscriptome.out.bam \
#        $data_dir/$sample \
#        --threads 4
#done


echo 'starting Splice Junction Discovery'
python3 $base_dir/SpliceJunctionDiscovery.py \
	-transcriptFile=transcript_file.bed \
	-bamList=bamlist.list\ 
	-processes=12

echo 'Starting Normalization of splice junction values'
python $base_dir/NormalizeSpliceJunctionValues.py \
	-transcript_model=$base_dir/gencode.comprehensive.splice.junctions.txt \
	-splice_file=$base_dir/All.transcript_file.bed.splicing.list \
	--normalize > $(date +"%Y%d%b%T").normalised.junctions
