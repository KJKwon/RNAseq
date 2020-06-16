#!/usr/bin/env bash

SAMTOOLS="$HOME/miniconda3/bin/samtools"
for SAM in $(ls *.sam)
do
#Mapping statistics
        STAT=$SAM".flagstat.txt"
        if [ ! -e $STAT ]; then
                echo "Calculating $SAM reads statistic..."
                $SAMTOOLS flagstat $SAM > $STAT
        fi
#SAM filtering for downstream process
        SAM_BEST=${SAM/.sam/.best.sam}
        echo $SAM_BEST  
        if [ ! -e $SAM_BEST ]; then
                echo "Filtering $SAM to $SAM_BEST..."
                python sam-to-sam_best.py $SAM
        fi
#Paired distance check for QC
        echo "Paired distance calculation from $SAM_BEST..."
        python sam_best-to-distance.py $SAM_BEST
#Total count, cpm, rpkm table
        echo "Calculating count data from $SAM_BEST..."
        python sam-to-count+cpm+rpkm.py $SAM_BEST

        echo "Process $SAM done!"
done

#Total mapping statistics
echo "Combine reads statistics to clean table..."
python Seq_stat2table.py
