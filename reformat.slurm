#!/bin/bash
# module load plink2


#cp /home/maihofer/freeze3_gwas/sumstats/*"$study"*.assoc.gz $TMPDIR/.
cp /home/maihofer/freeze3_gwas/sumstats/"$infile" $TMPDIR/.


echo "working on $study"

for chr in {1..23} X
do
 zcat "$TMPDIR"/$infile | awk -v chr=$chr -v colno=$colno -v freqcol=$freqcol '{if (NR == 1 || ($colno==chr && $freqcol >= 0.01 && $freqcol <= 0.99) ) print}'  | gzip > /home/maihofer/freeze3_gwas/sumstats/bychr/"$infile"_"$chr".gz
done

