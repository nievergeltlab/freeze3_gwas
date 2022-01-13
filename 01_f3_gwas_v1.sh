
### 1)Study Level Analysis steps: 

##Group 1 and 2 GWAS only: GWAS

cov=pcs
sex=all # females males
working_dir=$(pwd)

#may need to do gtpc as well
#Make sure time codes are correct
IFS=$'\n'
for line in $(cat dosage_locations_f3_aam_cont.csv | head -n 14  |  tail -n2 )
do
 study_1=$(echo $line | awk 'BEGIN{FS=","}  {print $1}')
 study_2=$(echo $line | awk 'BEGIN{FS=","}  {print $2}')

 ancgroup=$(echo $line | awk 'BEGIN{FS=","} {print $3}')
 timecode=$(echo $line | awk 'BEGIN{FS=","} {print $4}')
 exclude=$(echo $line | awk 'BEGIN{FS=","}  {print $5}')

 timecode=02:45:00
 
#analysis of phenotypes coded into .pheno files

 #if [ $exclude == "0" ] 
 #then
  echo gwas for $study_1 $study_2 $timecode 
  sbatch --time=$timecode --error errandout/"$study_1"_"$study_2"_"$cov"_"$ancgroup"_"$sex".e --output errandout/"$study_1"_"$study_2"_"$cov"_"$ancgroup"_"$sex".o  \
   --export=ALL,study="$study_1",cov="$cov",study_2="$study_2",ancgroup="$ancgroup",sex="$sex" run_trauma_gwas_v2_freeze3.slurm -D $working_dir 
# fi


#Case/Control analysis of built in phenotypes (For Grotzinger) - beware, phenotypes in files may not be right  or at least not optimal - such as for army starrs data!
 # if [ $exclude == "1" ] 
 # then
  # echo gwas for $study_1 $study_2 $timecode 
  # sbatch --time=$timecode --error errandout/"$study_1"_"$study_2"_"$cov"_"$ancgroup"_"$sex".e --output errandout/"$study_1"_"$study_2"_"$cov"_"$ancgroup"_"$sex".o  \
   # --export=ALL,study="$study_1",cov="$cov",study_2="$study_2",ancgroup="$ancgroup",sex="$sex" run_trauma_gwas_v2_freeze3_case_control.slurm -D $working_dir 
 # fi

done


##Group 3 and 4 and 5 only: Summary stat conversion from genome-wide to split by chromosome
 #Caution: these data may need to be reformatted first for meta-analysis. Do this PRIOR to running this script
 
 #To reformat by chromosome, supply a list .gz files. Each row should be a file name followed by the column number for the chromosome
 IFS=$'\n'
 for line in $(cat sumstat_reformat.csv  )
 do
  infile=$(echo $line | awk ' {print $1}')
  colno=$(echo $line | awk '{print $2}')
  freqcol=$(echo $line | awk ' {print $3}')
  echo $infile $colno $freqcol
  sbatch --time=00:25:00 --error errandout/"$infile"_reformat.e --output errandout/"$infile"_reformat.o  \
   --export=ALL,colno="$colno",freqcol="$freqcol",infile="$infile" reformat.slurm -D $working_dir 
 done
  
### 2)Meta analysis

#User: Make a meta-analysis script. Follow format of eur_ptsd_m0.mi
#User: Studies should be weighted by effective sample size. Weights are set in the meta-script
#00_weight_checkr.r contains code to take input phenotypes and estimate a weighting factor


#Adjust meta-analysis script for each chromosome (replaces "XXX" with a chromosome number)
 for chr in {1..22} #X
 do 
  #Basic analysis of all data
 # sed s/XXX/$chr.gz/g f3_gwas_may13_2021.mi  >   metal_scripts/f3_gwas_may13_2021.mi_$chr
 #   sed s/XXX/$chr.gz/g f3_gwas_may13_2021_ehr_only.mi  >   metal_scripts/f3_gwas_may13_2021_ehr_only.mi_$chr
#  sed s/XXX/$chr.gz/g f3_gwas_may13_2021_pgc_only.mi  >   metal_scripts/f3_gwas_may13_2021_pgc_only.mi_$chr
  # sed s/XXX/$chr.gz/g f3_gwas_may13_2021_nomvp.mi  >   metal_scripts/f3_gwas_may13_2021_nomvp.mi_$chr
  #   sed s/XXX/$chr.gz/g f3_gwas_may13_2021_noukbb.mi  >   metal_scripts/f3_gwas_may13_2021_noukbb.mi_$chr
     # sed s/XXX/$chr.gz/g f3_gwas_may13_2021_nomvpukbb.mi  >   metal_scripts/f3_gwas_may13_2021_nomvpukbb.mi_$chr
     
    sed s/XXX/$chr.gz/g f3_gwas_may13_2021_aam.mi  >   metal_scripts/f3_gwas_may13_2021_aam.mi_$chr
   #sed s/XXX/$chr.gz/g f3_gwas_may13_2021_casecontrol.mi  >   metal_scripts/f3_gwas_may13_2021_casecontrol.mi_$chr

 
 done


#List all chromosome meta-analysis files into metafilelist.txt file
 ls  metal_scripts/f3_gwas_may13_2021.mi_* > metafilelist.txt
 ls  metal_scripts/f3_gwas_may13_2021_ehr_only.mi_* > metafilelist.txt #EHR studies only
 ls  metal_scripts/f3_gwas_may13_2021_pgc_only.mi_* > metafilelist.txt #PGC studies only


 ls  metal_scripts/f3_gwas_may13_2021_nomvp.mi_* > metafilelist.txt #EHR studies only
 ls  metal_scripts/f3_gwas_may13_2021_noukbb.mi_* >> metafilelist.txt #PGC studies only
 ls  metal_scripts/f3_gwas_may13_2021_nomvpukbb.mi_* > metafilelist.txt #PGC studies only

 ls  metal_scripts/f3_gwas_may13_2021_casecontrol.mi_* > metafilelist.txt #PGC studies only
 ls  metal_scripts/f3_gwas_may13_2021_casecontrol_aam.mi_* > metafilelist.txt #PGC studies only
   ls  metal_scripts/f3_gwas_may13_2021_aam.mi_* > metafilelist.txt #PGC studies only
 
  #make sure all input files exist
 grep PROCESS  metal_scripts/f3_gwas_may13_2021_casecontrol_aam.mi_1 |  sed 's/PROCESS/head/g' > filecheck.txt
 

#User: Give a name for the error log files
 dataset=PTSD_F3_v5_jan4_2022_aam

 sbatch -t 00:50:00  --error errandout/"$dataset".e --output errandout/"$dataset".o   --export=ALL,metafile=metafilelist.txt -D /home/maihofer/freeze3_gwas run_meta_v2_loo_v2.slurm
 
 
#Concatenate meta-analysis results
 cat metal_results/eur_ptsd_pcs_v4_aug3_2021_*.gz1.tbl | awk '{if (NR == 1 || ($3 != "MarkerName" )) print}' > metal_results/eur_ptsd_pcs_v4_aug3_2021.gz1.tbl
 cat metal_results/aam_ptsd_pcs_v5_jan4_2021_*.gz1.tbl | awk '{if (NR == 1 || ($3 != "MarkerName" )) print}' > metal_results/aam_ptsd_pcs_v5_jan4_2022.gz1.tbl


 #Then run 00_translate_meta_pluses_to_samplesize.sh to get Ns
 
 
#Format meta results for fuma #Note: COLUMNS WILL CHANGE if ANALYZE HET IS ON!
#User: find out what the max N is make sure only decently covered loci are returned (by default, 90% of total N)
#TBD: Don't do this for the x chromosome! N will definitely be less!
 totalN=641553
 percentN=0.8
 cat metal_results/eur_ptsd_pcs_v4_aug3_2021_*.gz1.tbl   |  awk  -v totalN=$totalN -v percentN=$percentN   '{ if (NR ==1 || ($6 >= 0.01 && $6 <= 0.99 && $3 != "MarkerName" && $10 >= percentN*totalN)) print $1,$2,$3,$4,$5,$6,$10,$11,$12}'   | grep -v : | sort -g -k 9 | gzip > results_filtered/eur_ptsd_pcs_v4_aug3_2021.fuma.gz

 totalN=42804
 percentN=0.8
 cat metal_results/aam_ptsd_pcs_v5_jan4_2022.gz1.tbl  |  awk  -v totalN=$totalN -v percentN=$percentN   '{ if (NR ==1 || ($6 >= 0.01 && $6 <= 0.99 && $3 != "MarkerName" && $10 >= percentN*totalN)) print $1,$2,$3,$4,$5,$6,$10,$11,$12}'   | grep -v : | sort -g -k 9 | gzip > results_filtered/aam_ptsd_pcs_v5_jan4_2022.fuma.gz
 

#Format case/control analysis: First combine all files
 cat metal_results/eur_ptsdcasecontrol_pcs_v4_aug3_2021_*.gz1.tbl  |  awk   '{ if (NR ==1 || ($6 >= 0.01 && $6 <= 0.99 && $3 != "MarkerName" )) print $0}'   | grep -v : | sort -g -k 12 | gzip > eur_ptsdcasecontrol_pcs_v4_aug3_2021.temp1.tbl 
 cat metal_results/eur_ptsdcasecontrol_pcs_v4_aug3_2021_*.gz1.tbl  |  awk   '{ if (NR ==1 || ($6 >= 0.01 && $6 <= 0.99 && $3 != "MarkerName" )) print $0}'   | grep -v : | sort -g -k 12 | gzip > eur_ptsdcasecontrol_pcs_v4_aug3_2021.temp1.tbl 
 
#Next, apply  the same size counting script


###LDSC analysis
 conda activate ldsc

 #Filter data down to just LDSC SNPs then munge
 LC_ALL=C join -1 1 -2 3  <(awk '{print $1}' ~/trauma_gwas/ldsc-master/ldsc_data/eur_w_ld_chr/w_hm3.snplist.sorted | LC_ALL=C sort -k1b,1 ) <(cat metal_results/eur_ptsd_pcs_v4_aug3_2021.gz1.tbl | awk '{if (NR == 1) $3="SNP"; print}'   | LC_ALL=C sort -k3b,3 ) | gzip > results_filtered/eur_ptsd_pcs_v4_aug3_2021.gz1.tbl.premunge.gz 
 python2 ~/trauma_gwas/ldsc-master/munge_sumstats.py --sumstats results_filtered/eur_ptsd_pcs_v4_aug3_2021.gz1.tbl.premunge.gz  --N-col Weight --out results_filtered/eur_ptsd_pcs_v4_aug3_2021.gz1.tbl.munge.gz #add --N-col OBS_CT for the sample size
 
 #Get format for LDhub
 zcat results_filtered/eur_ptsd_pcs_v4_aug3_2021.fuma.gz | awk '{if (NR==1) print "CHR","BP","SNP","A1","A2","FRQ","N","Z","P"; else print $0}' > results_filtered/eur_ptsd_pcs_v4_aug3_2021.gz1.tbl.munge.txt
 cd results_filtered
 zip eur_ptsd_pcs_v4_aug3_2021.gz1.tbl.munge.txt.zip eur_ptsd_pcs_v4_aug3_2021.gz1.tbl.munge.txt
 
 #Run LDSC
  python2 /home/maihofer/trauma_gwas/ldsc-master/ldsc.py \
 --h2 results_filtered/eur_ptsd_pcs_v4_aug3_2021.gz1.tbl.munge.gz.sumstats.gz \
 --ref-ld-chr /home/maihofer/trauma_gwas/ldsc-master/ldsc_data/eur_w_ld_chr/ \
 --w-ld-chr  /home/maihofer/trauma_gwas/ldsc-master/ldsc_data/eur_w_ld_chr/ \
 --out results_filtered/eur_ptsd_pcs_v4_aug3_2021.gz1.tbl.munge.gz.tbl.ldsc

#extra analyses

 totalN=469191.00
 percentN=0.8
 cat metal_results/eur_ptsd_pcs_v4_aug3_2021_noukbb_*.gz1.tbl   |  awk  -v totalN=$totalN -v percentN=$percentN   '{ if (NR ==1 || ($6 >= 0.01 && $6 <= 0.99 && $3 != "MarkerName" && $10 >= percentN*totalN)) print $1,$2,$3,$4,$5,$6,$10,$11,$12}'   | grep -v : | sort -g -k 9 | gzip > results_filtered/eur_ptsd_pcs_v4_aug3_2021_noukbb.fuma.gz


 totalN=454844.00
 percentN=0.8
 cat metal_results/eur_ptsd_pcs_v4_aug3_2021_nomvp_*.gz1.tbl   |  awk  -v totalN=$totalN -v percentN=$percentN   '{ if (NR ==1 || ($6 >= 0.01 && $6 <= 0.99 && $3 != "MarkerName" && $10 >= percentN*totalN)) print $1,$2,$3,$4,$5,$6,$10,$11,$12}'   | grep -v : | sort -g -k 9 | gzip > results_filtered/eur_ptsd_pcs_v4_aug3_2021_nomvp.fuma.gz



 totalN=279432.00
 percentN=0.8
 cat metal_results/eur_ptsd_pcs_v4_aug3_2021_nomvpukbb_*.gz1.tbl   |  awk  -v totalN=$totalN -v percentN=$percentN   '{ if (NR ==1 || ($6 >= 0.01 && $6 <= 0.99 && $3 != "MarkerName" && $10 >= percentN*totalN)) print $1,$2,$3,$4,$5,$6,$10,$11,$12}'   | grep -v : | sort -g -k 9 | gzip > results_filtered/eur_ptsd_pcs_v4_aug3_2021_nomvpukbb.fuma.gz



#Format case/control analysis: First combine all files
 cat metal_results/eur_ptsdcasecontrol_pcs_v5_jan4_2021_*.gz1.tbl  |  awk   '{ if (NR ==1 || ($6 >= 0.01 && $6 <= 0.99 && $3 != "MarkerName" )) print }'   | grep -v : | sort -g -k 12 | gzip > eur_ptsdcasecontrol_pcs_v5_jan4_2022.temp1.tbl  
 cat metal_results/aam_ptsdcasecontrol_pcs_v5_jan4_2021_*.gz1.tbl  |  awk   '{ if (NR ==1 || ($3 != "MarkerName" )) print }'   | grep -v : | grep -v '???????'  | gzip > aam_ptsdcasecontrol_pcs_v5_jan4_2022.temp1.tbl 


#Next, apply  the same size counting script

#Then filter on 80% neff
 totalN=463244
 percentN=0.8
 cat eur_ptsdcasecontrol_pcs_v5_jan4_2022.temp1.tbl.neffX |  awk  -v totalN=$totalN -v percentN=$percentN   '{ if (NR ==1 || ($6 >= 0.01 && $6 <= 0.99 && $3 != "MarkerName" && $16 >= percentN*totalN)) print $1,$2,$3,$4,$5,$6,$10,$11,$12,$14,$15,$16,$17}'   | grep -v : | gzip > eur_ptsdcasecontrol_pcs_v5_jan4_2022X.fuma.gz
 #using eur_ptsdcasecontrol_pcs_v5_jan4_2022.fuma.gz and not eur_ptsdcasecontrol_pcs_v5_jan4_2022X.fuma.gz - eur_ptsdcasecontrol_pcs_v5_jan4_2022X.fuma.gz contains alterations to huckins sumstats. dont bother.
 
 totalN=51252
 percentN=0.8
 cat aam_ptsdcasecontrol_pcs_v5_jan4_2022.temp1.tbl.neff |  awk  -v totalN=$totalN -v percentN=$percentN   '{ if (NR ==1 || ($6 >= 0.01 && $6 <= 0.99 && $3 != "MarkerName" && $16 >= percentN*totalN)) print $1,$2,$3,$4,$5,$6,$10,$11,$12,$14,$15,$16,$17}'  | sort -g -k 9 | grep -v :  | gzip > aam_ptsdcasecontrol_pcs_v5_jan4_2022.fuma.gz


 #Filter data down to just LDSC SNPs then munge
 LC_ALL=C join -1 1 -2 3  <(awk '{print $1}' ~/trauma_gwas/ldsc-master/ldsc_data/eur_w_ld_chr/w_hm3.snplist.sorted | LC_ALL=C sort -k1b,1 ) <(zcat eur_ptsdcasecontrol_pcs_v5_jan4_2022X.fuma.gz  | awk '{if (NR == 1) $3="SNP"; print}'   | LC_ALL=C sort -k3b,3 ) | gzip > results_filtered/eur_ptsdcasecontrol_pcs_v5_jan4_2022.fuma.gz.premunge.gz 
 python2 ~/trauma_gwas/ldsc-master/munge_sumstats.py --sumstats results_filtered/eur_ptsdcasecontrol_pcs_v5_jan4_2022.fuma.gz.premunge.gz  --N-col ntot --out results_filtered/eur_ptsdcasecontrol_pcs_v5_jan4_2022.fuma.gz.munge.gz #add --N-col OBS_CT for the sample size
 
  LC_ALL=C join -1 1 -2 3  <(awk '{print $1}' ~/trauma_gwas/ldsc-master/ldsc_data/eur_w_ld_chr/w_hm3.snplist.sorted | LC_ALL=C sort -k1b,1 ) <(zcat results_filtered/aam_ptsdcasecontrol_pcs_v5_jan4_2022.fuma.gz  | awk '{if (NR == 1) $3="SNP"; print}'   | LC_ALL=C sort -k3b,3 ) | gzip > results_filtered/aam_ptsdcasecontrol_pcs_v5_jan4_2022.fuma.gz.premunge.gz 
 python2 ~/trauma_gwas/ldsc-master/munge_sumstats.py --sumstats results_filtered/aam_ptsdcasecontrol_pcs_v5_jan4_2022.fuma.gz.premunge.gz  --N-col ntot --out results_filtered/aam_ptsdcasecontrol_pcs_v5_jan4_2022.fuma.gz.munge.gz #add --N-col OBS_CT for the sample size
 
   LC_ALL=C join -1 1 -2 3  <(awk '{print $1}' ~/trauma_gwas/ldsc-master/ldsc_data/eur_w_ld_chr/w_hm3.snplist.sorted | LC_ALL=C sort -k1b,1 ) <(zcat results_filtered/aam_ptsd_pcs_v5_jan4_2022.fuma.gz  | awk '{if (NR == 1) $3="SNP"; print}'   | LC_ALL=C sort -k3b,3 ) | gzip > results_filtered/aam_ptsd_pcs_v5_jan4_2022.fuma.gz.premunge.gz 
 python2 ~/trauma_gwas/ldsc-master/munge_sumstats.py --sumstats results_filtered/aam_ptsd_pcs_v5_jan4_2022.fuma.gz.premunge.gz    --out results_filtered/aam_ptsd_pcs_v5_jan4_2022.fuma.gz.munge.gz #add --N-col OBS_CT for the sample size
 
 
 
 #Run LDSC
  python2 /home/maihofer/trauma_gwas/ldsc-master/ldsc.py \
 --h2 results_filtered/aam_ptsd_pcs_v5_jan4_2022.fuma.gz.munge.gz.sumstats.gz \
 --samp-prev 0.11 \
--pop-prev 0.20 \
 --ref-ld-chr /home/maihofer/trauma_gwas/ldsc-master/ldsc_data/eur_w_ld_chr/ \
 --w-ld-chr  /home/maihofer/trauma_gwas/ldsc-master/ldsc_data/eur_w_ld_chr/ \
 --out results_filtered/aam_ptsd_pcs_v5_jan4_2022.fuma.gz.tbl.ldsc2


    python2 /home/maihofer/trauma_gwas/ldsc-master/ldsc.py \
 --rg results_filtered/eur_ptsdcasecontrol_pcs_v5_jan4_2022.fuma.gz.munge.gz.sumstats.gz,results_filtered/eur_ptsd_pcs_v4_aug3_2021.gz1.tbl.munge.gz.sumstats.gz \
 --ref-ld-chr /home/maihofer/trauma_gwas/ldsc-master/ldsc_data/eur_w_ld_chr/ \
 --w-ld-chr  /home/maihofer/trauma_gwas/ldsc-master/ldsc_data/eur_w_ld_chr/ \
 --out results_filtered/rg_eur_ptsd_pcs_v4_aug3_2021.gz1.tbl.munge.gz.sumstats.gz_eur_ptsdcasecontrol_pcs_v5_jan4_2022.fuma.gz.tbl.ldsc2


    python2 /home/maihofer/trauma_gwas/ldsc-master/ldsc.py \
 --rg results_filtered/eur_ptsdcasecontrol_pcs_v5_jan4_2022.fuma.gz.munge.gz.sumstats.gz,results_filtered/eur_ptsdcasecontrol_pcs_v4_aug3_2021.fuma.gz.munge.gz.sumstats.gz \
 --ref-ld-chr /home/maihofer/trauma_gwas/ldsc-master/ldsc_data/eur_w_ld_chr/ \
 --w-ld-chr  /home/maihofer/trauma_gwas/ldsc-master/ldsc_data/eur_w_ld_chr/ \
 --out results_filtered/rg_eur_ptsdcasecontrol_pcs_v4_aug3_2021.fuma.gz.munge.gz.sumstats.gz_eur_ptsdcasecontrol_pcs_v5_jan4_2022.fuma.gz.tbl.ldsc2




