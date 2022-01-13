#!/bin/bash
# module load plink2

#Only has minor modifications to be run for f3, namely changing datadir, reslength -le, and looking for chr10_000_020 to account for chunk size
#have also added sex chromosome. Sex chromosomes are imputed separately by sex so I do some stuff to account for this, namely dosage list command
#Don't really worry about PAR, ricopili doesn't really impute it (some overlap, have to worry about this in filtering).
# see this for https://www.cog-genomics.org/plink/1.9/assoc x chromosome handling. may have to reduce to plink 1.9 and produce output separately!!
# https://www.cog-genomics.org/plink/2.0/assoc#glm ill only hae to do this if I can't figure out how to combine male and female x imputation files!!
#may just want to combine them  prior to any of this! if plink2 doesnt take dosage-lsit that sucks!
#I think this may be the best option...

echo $study
echo "--export=study="$study",phenolevel="$phenolevel",pheno="$pheno",cov="$cov",study_2="$study_2",ancgroup="$ancgroup""


datadir=/home/pgcdacpt/dac/pts/wave3/v1/"$study"/qc1/ 

mkdir "$TMPDIR"/"$study"_"$study_2"

#for files in $( ls $datadir | grep .gz$ | sed 's/.gz//g' | grep  chr10_000_020) TEST CASE
for files in $( ls $datadir | grep .gz$ | sed 's/.gz//g' | grep -v "chrX.mal" | grep -v "chrX.fem")
do
 outname=$(echo $files | awk 'BEGIN{FS="ch.fl.chr"}{print $2}'  | sed 's/.out.dosage//g') #get output filename...
 echo $outname
 
if [ $sex == "males" ]
then
 #plink_flag="--filter-males"
 awk '{ if ($5 == 1) print $1,$2}' "$datadir"/"$files".fam > "$TMPDIR"/"$files".keepm
 plink_flag="--keep "$TMPDIR"/"$files".keepm"
 extension="_males"
fi

if [ $sex == "females" ]
then
 #plink_flag="--filter-females"
 awk '{ if ($5 == 2) print $1,$2}' "$datadir"/"$files".fam > "$TMPDIR"/"$files".keepf
 plink_flag="--keep "$TMPDIR"/"$files".keepf"
 extension="_females"
fi

    echo "Assuming $plink_flag"
    
   /home/maihofer/trauma_gwas/plink2 --import-dosage "$datadir"/"$files".gz --fam "$datadir"/"$files".fam --map "$datadir"/"$files".map --pheno pheno/p2_"$study"_"$study_2"_"$ancgroup".pheno   \
    --covar pheno/p2_"$study"_"$ancgroup"_"$study_2"_"$cov".cov --maf 0.01 \
    --glm cols=chrom,pos,ref,alt,ax,a1freq,machr2,a1freqcc,firth,test,nobs,orbeta,se,ci,tz,p   hide-covar firth-fallback  \
    --covar-variance-standardize --threads 11 --ci 0.95 $plink_flag \
    --out "$TMPDIR"/"$study"_"$study_2"/"$study"_"$study_2"_"$ancgroup"_chr"$outname"_"$cov""$extension"
done


  #Auto detect pheno levels evaluated and loop over these... 
  for phenoset in $(ls "$TMPDIR"/"$study"_"$study_2" | grep chr10_000_020 | awk 'BEGIN{FS="."}{print $2}' | grep -v log | grep -v PHENO1  )
  do
   reslength=$(ls "$TMPDIR"/"$study"_"$study_2"/"$study"_"$study_2"_"$ancgroup"_*_"$cov""$extension"."$phenoset".glm.* | wc -l) #Get linecount to make sure thats good!!
   echo "Working on "$study"_"$study_2"_"$ancgroup"_*_"$cov""$extension"."$phenoset""
   if [ $reslength -le 130 ]
   then
    echo "Fewer than 130 results files --  $reslength -- something maybe be wrong "
   fi
  
  if [ $( echo $phenoset | grep -c "Case") -eq 1 ]
   then
    cat "$TMPDIR"/"$study"_"$study_2"/"$study"_"$study_2"_"$ancgroup"_*_"$cov""$extension"."$phenoset".glm.* | awk '{if (NR == 1 || $3 != "ID") print}'  | sed 's/#//g'  | awk '{ if (NR==1 || ($8 >= 0.01 && $8 <= 0.99 && $11 > 0.6)) print }' | gzip > $TMPDIR/"$study"_"$study_2"_"$ancgroup"_"$cov""$extension"."$phenoset".assoc.gz 
   else
    cat "$TMPDIR"/"$study"_"$study_2"/"$study"_"$study_2"_"$ancgroup"_*_"$cov""$extension"."$phenoset".glm.* | awk '{if (NR == 1 || $3 != "ID") print}'  | sed 's/#//g'  | awk '{ if (NR==1 || ($8 >= 0.01 && $8 <= 0.99 && $9 > 0.6)) print }' | gzip > $TMPDIR/"$study"_"$study_2"_"$ancgroup"_"$cov""$extension"."$phenoset".assoc.gz
   fi
   
   #Now put it in the home dir
   cp "$TMPDIR"/"$study"_"$study_2"_"$ancgroup"_"$cov""$extension"."$phenoset".assoc.gz results_cat/"$study"_"$study_2"_"$ancgroup"_"$cov""$extension"."$phenoset".assoc.gz
   #And make a split one to make this easy for us...
   for chr in {1..22} X 
   do
    zcat "$TMPDIR"/"$study"_"$study_2"_"$ancgroup"_"$cov""$extension"."$phenoset".assoc.gz | awk -v chr=$chr '{if (NR == 1 || $1==chr) print}'  | gzip > results_cat/"$study"_"$study_2"_"$ancgroup"_"$cov""$extension"."$phenoset".assoc.gz_"$chr".gz
   done
   
   echo "completed"
   #rm -f results/"$study"_*_"$phenolevel"_"$ancgroup"_"$study_2"_"$cov"."$pheno".glm.linear 
  done




 
	