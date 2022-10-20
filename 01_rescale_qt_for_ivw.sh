
#all PLINK studies
IFS=$'\n'
for studycomb in $(cat conversion_factors_ivw_males.csv  |   tail -n+2 | grep -v ukbb | grep -v vetsa | grep ncmh )
do
 echo converting $studycomb
 study=$(echo $studycomb | awk 'BEGIN{ FS=","}{ print $1}')
 study2=$(echo $study | sed 's/.gz//g')
 
 #subset=$(echo $studycomb | awk 'BEGIN{ FS=","}{ print $2}')
 fraction=$(echo $studycomb | awk 'BEGIN{ FS=","}{ print $4}')

  for chr in {1..22}
  do
   zcat "$study2".chr"$chr".gz | awk -v multiplier=$fraction '{if (NR>1) {$12=$12*multiplier; $13=$13*multiplier}; print}'  | gzip > "$study"_"$chr".rescale.gz
  done

 done
 #note: betr is now c/c
 #check on ftca, wach. 
 #grac, nhs2, wach. are a warning to ignore, those dont exist in males
 #results_cat/ftca_ftca_eur_pcs_males.Lifetime_PTSD_Continuous.assoc.chr1.gz

 
 
 #special handling for UKBB and VETSA (BOLT LMM STUDIES) - done for males
 for studycomb in $(cat /home/maihofer/freeze3_gwas/conversion_factors_ivw_males.csv  |  tail -n 2 | grep ukbb)
do
  echo converting $studycomb
 study=$(echo $studycomb | awk 'BEGIN{ FS=","}{ print $1}')
 #subset=$(echo $studycomb | awk 'BEGIN{ FS=","}{ print $2}')
 fraction=$(echo $studycomb | awk 'BEGIN{ FS=","}{ print $4}')

 for chr in {1..22}
 do
  zcat "$study"_"$chr".gz | awk -v multiplier=$fraction '{if (NR>1) {$11=$11*multiplier; $12=$12*multiplier}; print}'  | gzip > "$study"_"$chr".rescale.gz
 done
done
