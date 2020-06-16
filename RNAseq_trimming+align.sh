#!/usr/bin/env bash

TRIMMOMATIC="$HOME/miniconda3/bin/trimmomatic"
BWA="$HOME/miniconda3/bin/bwa"
STAR="$HOME/miniconda3/bin/STAR"

THREAD_NUM=6

BWA_DB="/work/SeqRef/bwaRef/RAT/rat98_longest_bwa"
BWA_DB_NAME=$(basename $BWA_DB)
STAR_DB_DIR="/work/SeqRef/starRef/RAT"
STAR_DB_NAME="RAT_ens98_dna_rm"

for FASTQ in $(ls *.fastq.gz)
do
        FASTQ_NM=$(basename $FASTQ .fastq.gz)
        TRIM_FASTQ=$FASTQ_NM".fastq.trimmed.gz"
        TRIM_LOG=$FASTQ_NM".trimmed.log"
        if [ ! -e $TRIM_FASTQ ]; then
	
		$TRIMMOMATIC SE -threads $THREAD_NUM -trimlog $TRIM_LOG $FASTQ $TRIM_FASTQ ILLUMINACLIP:NEB_adapter.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:50
	
	fi
	
	SAM_NAME=$(basename $TRIM_FASTQ .fastq.gz)
        BWA_SAM=$SAM_NAME"."$BWA_DB_NAME".bwa_mem.sam"

#        if [ -e $TRIM_FASTQ ] && [! -e $BWA_SAM ]; then
#
#                echo "Make $SAM by BWA" 
#                $BWA mem -t $THREAD_NUM $BWA_DB $TRIM_FASTQ > $BWA_SAM 2>log
#
#        fi
	
	STAR_SAM=$(basename $SAM_NAME .gz)"."$STAR_DB_NAME".STAR."
	STAR_SAM_DIR="$FASTQ_NM"
	if [ -e $TRIM_FASTQ ] && [ ! -e $STAR_SAM ]; then
		
		echo "Make $SAM by STAR"
                if [ ! -d $STAR_SAM_DIR ]; then
	        
			mkdir $STAR_SAM_DIR

	        fi		
		$STAR --genomeDir $STAR_DB_DIR --runThreadN $THREAD_NUM --readFilesIn $TRIM_FASTQ --readFilesCommand zcat --quantMode TranscriptomeSAM GeneCounts \
			--outFileNamePrefix ./$STAR_SAM_DIR/$STAR_SAM  
		
	fi
done
