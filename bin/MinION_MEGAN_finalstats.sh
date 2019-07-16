#!/bin/bash
#SBATCH -t 01:00:00

####################################################
## FINAL STATS FILES
####################################################
## STATS FOR TRIMMING WITH PRINSEQ-LITE:
for i in $(grep -l 'ERROR' $WDIR/1_Rawdata/trim_logs/*.log); do
	grep -vwE 'ERROR:' $i > tmp.log
	mv tmp.log $i
done
rm tmp.log

for i in $(ls $WDIR/1_Rawdata/trim_logs/*.log); do
	echo $i | xargs -n1 basename > Trimcol.tmp
	awk -F': ' 'FNR>=4 && FNR<=12 {print $2}' $i >> Trimcol.tmp
	awk -F': ' 'FNR>=14, FNR==15 {print $2}' $i >> Trimcol.tmp
	perl -00lpe 's/\n/ /g' Trimcol.tmp
done > $WDIR/Minion_trim.tsv
rm Trimcol.tmp

sed -i '/^$/d; s/ /\t/g' $WDIR/Minion_trim.tsv
sed -i $'1 i\\\nOutName\tInput_sequences\tInput_bases\tInput_meanlength\tGood_seqs\tGood_seqs_%\tGood_bases\tGood_meanlength\tBad_seqs\tBad_seqs_%\tBad_bases\tBad_meanlength\tmin_len\tmin_qual_mean' $WDIR/Minion_trim.tsv

## STATS FOR DIAMOND TIME
for i in $(ls $WDIR/6_logs_dmnd/DIAMOND/*.out); do
	grep -E 'Input file:|Total time =|^Reported|queries aligned.$' $i | sed -e 's/^.*:\s*//; s/^.*=\s*//' | perl -00lpe 's/\n/\t/g'
done | sed -e '/^$/d' > $WDIR/Minion_dmnd.tsv
sed -i $'1 i\\ID\tTotal time\tPairwise alignments\tAligned queries' $WDIR/Minion_dmnd.tsv

## STATS FOR MEGAN TIME
for i in $(ls $WDIR/6_logs_dmnd/MEGAN/*.out); do
	grep -E 'Input file:|Total time:|Peak memory:' $i | sed -e 's/^.*:\s*//' | perl -00lpe 's/\n/\t/g'
done | sed -e '/^$/d' > Minion_megan.tsv
sed -i $'1 i\\ID\tTotal time\tUsed memory' Minion_megan.tsv

mv logs_dmnd $WDIR
mv Minion_*.tsv $WDIR
