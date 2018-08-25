# RNA-seq #
RNA-seq pipeline based on Beryl Cummings' and Dennis Kao's methods.

This pipeline is designed to be flexible to the user's directory structure, and require minimal user input once configured. It will use STAR to align Illumina RNA-seq fastq files using 2-pass mode, which is more sensitive to novel splice junction discovery. It is written with the assumption of 150bp unpaired reads. It will then deduplicate your samples (picard tools), and use Dennis Kao's update of SpliceJunctionDiscovery.py and SpliceJunctionNormalization.py. the deduplication step is parallelised with GNU parallel. 

# Usage:
Add STAR_to_norm.sh to your PATH. Alternatively, add it to /usr/bin or equivilent, or specify the path when using the script.
Place all fastq files into one directory. This should not contain any text files (ie foo.txt). Run the script from the dir containing the fastq files. 

This assumes you have generated the appropriate genome index files with STAR. if you have not done this already then you will need to download ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_19/gencode.v19.annotation.gtf.gz and ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_19/GRCh37.p13.genome.fa.gz unzip them, and run the genome_generate.sh script. 

## Note ##
STAR is RAM intensive, and requires ~32GB for a human alignment. It will crash if you have less than that available. Also, if you are likely to change the number of threads, keep in mind that >4 threads will use >32GB of RAM. Do not increase thread number unless you have sufficient RAM available. This also means being careful about running other processes concurrently! I have a VM for RNA-seq alone. 
you also need enough space for output files, ~100gb as a minimum
