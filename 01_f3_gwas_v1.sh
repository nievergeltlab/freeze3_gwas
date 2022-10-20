
### 1)Study Level Analysis steps: 

##Group 1 and 2 GWAS only: GWAS

cov=pcs
sex=all # females males
working_dir=$(pwd)

#may need to do gtpc as well
#Make sure time codes are correct
IFS=$'\n'
for line in $(cat dosage_locations_f3_aam_cont.csv  | grep saf2)
do
 study_1=$(echo $line | awk 'BEGIN{FS=","}  {print $1}')
 study_2=$(echo $line | awk 'BEGIN{FS=","}  {print $2}')

 ancgroup=$(echo $line | awk 'BEGIN{FS=","} {print $3}')
 timecode=$(echo $line | awk 'BEGIN{FS=","} {print $4}')
 exclude=$(echo $line | awk 'BEGIN{FS=","}  {print $5}')

 
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
 
 #Some X chromosome data needs to be converted, from showing 23 to showing X
 for files in ptsd_qt_vetsa_may12_2021_related_filtered.imputed.stats.gz_23.gz pts_ukbb_may13_2021_unrelated_chrX.bgen.stats.gz # daner_iPSYCH2015_PTSDbroad_chrX_HRC_MAF01.gz PTSD_EHR_Broad_2020_12_10_ReGENIE_Step2_With_X_WG.txt_Broad.regenie.gz_23.gz PTSD_broad_EstBB_GWAS_saige_chr23.txt.noinfo.gz
 do
 # zcat sumstats/bychr/$files | sed 's/^23[[:blank:]]/X /g' | gzip > sumstats/bychr/"$files".23.gz
    zcat sumstats/bychr/$files | sed 's/[[:blank:]]23[[:blank:]]/ X /g' | gzip > sumstats/bychr/"$files".23.gz
  
  zcat sumstats/bychr/"$files".23.gz | head
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
 #sed s/XXX/$chr.gz/g f3_gwas_may13_2021_ehr.mi  >   metal_scripts/f3_gwas_may13_2021_ehr.mi_$chr
 #sed s/XXX/$chr.gz/g f3_gwas_may13_2021_pgc.mi  >   metal_scripts/f3_gwas_may13_2021_pgc.mi_$chr
 
  #sed s/XXX/$chr.gz/g f3_gwas_may13_2021_lat.mi  >   metal_scripts/f3_gwas_may13_2021_lat.mi_$chr
 
 #f3_gwas_may13_2021_casecontrol_nomvp.mi
 
 #   sed s/XXX/$chr.gz/g f3_gwas_may13_2021_ehr_only.mi  >   metal_scripts/f3_gwas_may13_2021_ehr_only.mi_$chr
#  sed s/XXX/$chr.gz/g f3_gwas_may13_2021_pgc_only.mi  >   metal_scripts/f3_gwas_may13_2021_pgc_only.mi_$chr
  # sed s/XXX/$chr.gz/g f3_gwas_may13_2021_nomvp.mi  >   metal_scripts/f3_gwas_may13_2021_nomvp.mi_$chr
  #   sed s/XXX/$chr.gz/g f3_gwas_may13_2021_noukbb.mi  >   metal_scripts/f3_gwas_may13_2021_noukbb.mi_$chr
     # sed s/XXX/$chr.gz/g f3_gwas_may13_2021_nomvpukbb.mi  >   metal_scripts/f3_gwas_may13_2021_nomvpukbb.mi_$chr
      sed s/XXX/$chr.gz/g f3_gwas_may13_2021_no25.mi  >   metal_scripts/f3_gwas_may13_2021_no25.mi_$chr 
   # sed s/XXX/$chr.gz/g f3_gwas_may13_2021_aam_nogtp.mi  >   metal_scripts/f3_gwas_may13_2021_aam_nogtp.mi_$chr
  # sed s/XXX/$chr.gz/g f3_gwas_may13_2021_casecontrol_nomvp.mi  >   metal_scripts/f3_gwas_may13_2021_casecontrol_nomvp.mi_$chr
    # sed s/XXX/$chr/g f3_gwas_males_may13_2021_qtonly.mi  >   metal_scripts/f3_gwas_males_may13_2021_qtonly.mi_$chr
  # sed s/XXX/$chr/g f3_gwas_females_may13_2021_qtonly.mi  >   metal_scripts/f3_gwas_females_may13_2021_qtonly.mi_$chr

    # sed s/XXX/$chr.gz/g f3_gwas_males_may13_2021_cconly.mi  >   metal_scripts/f3_gwas_males_may13_2021_cconly.mi_$chr
       # sed s/XXX/$chr.gz/g f3_gwas_females_may13_2021_cconly.mi  >   metal_scripts/f3_gwas_females_may13_2021_cconly.mi_$chr
     
 done


#List all chromosome meta-analysis files into metafilelist.txt file
 ls  metal_scripts/f3_gwas_may13_2021.mi_* > metafilelist.txt
 ls  metal_scripts/f3_gwas_may13_2021_ehr.mi_* > metafilelist.txt #EHR studies only
 ls  metal_scripts/f3_gwas_may13_2021_pgc.mi_* >> metafilelist.txt #PGC studies only

 ls  metal_scripts/f3_gwas_males_may13_2021_qtonly.mi_* > metafilelist.txt #PGC studies only
 ls  metal_scripts/f3_gwas_females_may13_2021_qtonly.mi_* > metafilelist.txt #PGC studies only 
 ls  metal_scripts/f3_gwas_males_may13_2021_cconly.mi_* > metafilelist.txt #PGC studies only
 ls  metal_scripts/f3_gwas_females_may13_2021_cconly.mi_* > metafilelist.txt #PGC studies only


 ls  metal_scripts/f3_gwas_may13_2021_nomvp.mi_* > metafilelist.txt #EHR studies only
 
 ls  metal_scripts/f3_gwas_may13_2021_noukbb.mi_* >> metafilelist.txt #PGC studies only
 ls  metal_scripts/f3_gwas_may13_2021_nomvpukbb.mi_* > metafilelist.txt #PGC studies only

 ls  metal_scripts/f3_gwas_may13_2021_no25.mi_* > metafilelist.txt #EHR studies only
 
 ls  metal_scripts/f3_gwas_may13_2021_casecontrol.mi_* > metafilelist.txt #PGC studies only
  ls  metal_scripts/f3_gwas_may13_2021_casecontrol_nomvp.mi_* > metafilelist.txt #PGC studies only
 
 ls  metal_scripts/f3_gwas_may13_2021_casecontrol_aam.mi_* > metafilelist.txt #PGC studies only
 ls  metal_scripts/f3_gwas_may13_2021_aam.mi_* > metafilelist.txt #PGC studies only
 
  ls  metal_scripts/f3_gwas_may13_2021_lat.mi_* > metafilelist.txt #PGC studies only
 
  ls  metal_scripts/f3_gwas_may13_2021_aam_nogtp.mi_* > metafilelist.txt #PGC studies only
 
 echo f3_gwas_may13_2021_X.mi > metafilelist.txt
 
  echo f3_gwas_may13_2021_aam_x.mi > metafilelist.txt
  echo f3_gwas_may13_2021_lat_X.mi > metafilelist.txt
 
 
 f3_gwas_may13_2021_lat.mi
 
  #make sure all input files exist

 grep PROCESS  metal_scripts/f3_gwas_may13_2021_casecontrol_aam.mi_1 |  sed 's/PROCESS/head/g' > filecheck.txt
 
   #COUNT SAMPLE SIZE using this, i.e. 
 grep PROCESS metal_scripts/f3_gwas_females_may13_2021_qtonly.mi_9 | sed 's/PROCESS/zcat/g' | sed 's/rescale.gz/rescale.gz | head/g'> f25_femalecountqt.txt
 # in the output file, append code to get the Nth column where N is listed.
 #    zcat results_cat/onga_onga_eur_pcs_females.Lifetime_PTSD_Continuous.assoc.gz_9.rescale.gz | head |   tail -n 1 | awk '{print $11}' >> n_qt_females.txt

#User: Give a name for the error log files
 dataset=PTSD_3_no25

 sbatch -t 01:50:00  --error errandout/"$dataset".e --output errandout/"$dataset".o   --export=ALL,metafile=metafilelist.txt -D /home/maihofer/freeze3_gwas run_meta_v2_loo_v2.slurm
 
  sbatch -t 01:50:00  --error TEST.e --output TEST.o   --export=ALL,metafile=metafilelist.txt -D /home/maihofer/freeze3_gwas run_meta_v2_loo_v2.slurm
 
#Concatenate meta-analysis results

 

#Format meta results for fuma #Note: COLUMNS WILL CHANGE if ANALYZE HET IS ON!

#User: find out what the max N is make sure only decently covered loci are returned (by default, 90% of total Neff)
#TBD: X chromosome needs its own cutoffs

 #Main analysis
  #european
   cat metal_results/eur_ptsd_pcs_v4_aug3_2021_*.gz1.tbl | awk '{if (NR == 1 || ($3 != "MarkerName" )) print}' > metal_results/eur_ptsd_pcs_v4_aug3_2021.gz1.tbl
   totalN=641553
   percentN=0.8
   cat metal_results/eur_ptsd_pcs_v4_aug3_2021_*.gz1.tbl   |  awk  -v totalN=$totalN -v percentN=$percentN   '{ if (NR ==1 || ($6 >= 0.01 && $6 <= 0.99 && $3 != "MarkerName" && $10 >= percentN*totalN)) print $1,$2,$3,$4,$5,$6,$10,$11,$12}'   | grep -v : | sort -g -k 9 | gzip > results_filtered/eur_ptsd_pcs_v4_aug3_2021.fuma.gz

#Test: low cutoff
   totalN=641553
   percentN=0.25
  cat metal_results/eur_ptsd_pcs_v4_aug3_2021_[0-9].gz1.tbl metal_results/eur_ptsd_pcs_v4_aug3_2021_[0-9][0-9].gz1.tbl   |  awk  -v totalN=$totalN -v percentN=$percentN   '{ if (NR ==1 || ($6 >= 0.01 && $6 <= 0.99 && $3 != "MarkerName" && $10 >= percentN*totalN)) print $1,$2,$3,$4,$5,$6,$10,$11,$12}'   | grep -v : | sort -g -k 9 | gzip > results_filtered/eur_ptsd_pcs_v4_aug3_2021_LOWCUTOFF.fuma.gz

  #Add the X (if starting from scratch, beware the previous *, which may have already added the X. In that case, you have to code an if statement
   totalN=424500.00
   percentN=0.8
   cat metal_results/eur_ptsd_pcs_v4_aug3_2021_chrX.gz1.tbl   |  awk  -v totalN=$totalN -v percentN=$percentN   '{ if (NR ==1 || ($6 >= 0.01 && $6 <= 0.99 && $3 != "MarkerName" && $10 >= percentN*totalN)) print $1,$2,$3,$4,$5,$6,$10,$11,$12}'   | grep -v : | sort -g -k 9 | gzip > results_filtered/eur_ptsd_pcs_v4_aug3_2021.X.fuma.gz

   zcat results_filtered/eur_ptsd_pcs_v4_aug3_2021.fuma.gz results_filtered/eur_ptsd_pcs_v4_aug3_2021.X.fuma.gz | sort -g -k 9 | awk '{if (NR == 1 || $1 != "Chromosome") print}' | gzip > results_filtered/eur_ptsd_pcs_v4_aug3_2021.allchr.fuma.gz

  #AAM
   cat metal_results/aam_ptsd_pcs_v5_jan4_2021_*.gz1.tbl | awk '{if (NR == 1 || ($3 != "MarkerName" )) print}' > metal_results/aam_ptsd_pcs_v5_jan4_2022.gz1.tbl
   totalN=42804
   percentN=0.8
   cat metal_results/aam_ptsd_pcs_v5_jan4_2022.gz1.tbl  |  awk  -v totalN=$totalN -v percentN=$percentN   '{ if (NR ==1 || ($6 >= 0.01 && $6 <= 0.99 && $3 != "MarkerName" && $10 >= percentN*totalN)) print $1,$2,$3,$4,$5,$6,$10,$11,$12}'   | grep -v : | sort -g -k 9 | gzip > results_filtered/aam_ptsd_pcs_v5_jan4_2022.fuma.gz
  
   totalN=13246.00
   percentN=0.8
   cat metal_results/aam_ptsd_pcs_v5_jan4_2021_X.gz1.tbl  |  awk  -v totalN=$totalN -v percentN=$percentN   '{ if (NR ==1 || ($6 >= 0.01 && $6 <= 0.99 && $3 != "MarkerName" && $10 >= percentN*totalN)) print $1,$2,$3,$4,$5,$6,$10,$11,$12}'   | grep -v : | sort -g -k 9 | gzip > results_filtered/aam_ptsd_pcs_v5_jan4_2021_X.fuma.gz

   zcat results_filtered/aam_ptsd_pcs_v5_jan4_2022.fuma.gz results_filtered/aam_ptsd_pcs_v5_jan4_2021_X.fuma.gz | sort -g -k 9 | awk '{if (NR == 1 || $1 != "Chromosome") print}' | gzip > results_filtered/aam_ptsd_pcs_v5_jan4_2022.allchr.fuma.gz
  
   #LAT
   cat metal_results/hna_ptsd_pcs_v4_aug3_2021_*.gz1.tbl | awk '{if (NR == 1 || ($3 != "MarkerName" )) print}' > metal_results/hna_ptsd_pcs_v4_aug3_2021.gz1.tbl
   totalN=6530
   percentN=0.8
   cat metal_results/hna_ptsd_pcs_v4_aug3_2021.gz1.tbl  |  awk  -v totalN=$totalN -v percentN=$percentN   '{ if (NR ==1 || ($6 >= 0.01 && $6 <= 0.99 && $3 != "MarkerName" && $10 >= percentN*totalN)) print $1,$2,$3,$4,$5,$6,$10,$11,$12}'   | grep -v : | sort -g -k 9 | gzip > results_filtered/hna_ptsd_pcs_v4_aug3_2021.fuma.gz
  
   totalN=6530
   percentN=0.8
   cat metal_results/hna_ptsd_pcs_v4_aug3_2021_X.gz1.tbl  |  awk  -v totalN=$totalN -v percentN=$percentN   '{ if (NR ==1 || ($6 >= 0.01 && $6 <= 0.99 && $3 != "MarkerName" && $10 >= percentN*totalN)) print $1,$2,$3,$4,$5,$6,$10,$11,$12}'   | grep -v : | sort -g -k 9 | gzip > results_filtered/hna_ptsd_pcs_v4_aug3_2021_X.fuma.gz
  
   zcat results_filtered/hna_ptsd_pcs_v4_aug3_2021.fuma.gz results_filtered/hna_ptsd_pcs_v4_aug3_2021_X.fuma.gz | sort -g -k 9 | awk '{if (NR == 1 || $1 != "Chromosome") print}' | gzip > results_filtered/hna_ptsd_pcs_v4_aug3_2021.allchr.fuma.gz
    
   
   
   #TRANS
   totalN=42804
   percentN=0.8
   cat metal_results/trans_ptsd_pcs_v5_jan4_20221.tbl  |  awk  -v totalN=$totalN -v percentN=$percentN   '{ if (NR ==1 || ($6 >= 0.01 && $6 <= 0.99 && $3 != "MarkerName" && $10 >= percentN*totalN)) print $1,$2,$3,$4,$5,$6,$10,$11,$12}'   | grep -v : | sort -g -k 9 | gzip > results_filtered/trans_ptsd_pcs_v5_jan4_2022.fuma.gz
   
    
   
   
 #Main, PGC datasets only
   totalN=259695.00
   percentN=0.8
   cat metal_results/eur_ptsd_ehronly_pcs_v4_aug3_2021_*.gz1.tbl   |  awk  -v totalN=$totalN -v percentN=$percentN   '{ if (NR ==1 || ($6 >= 0.01 && $6 <= 0.99 && $3 != "MarkerName" && $10 >= percentN*totalN)) print $1,$2,$3,$4,$5,$6,$10,$11,$12}'   | grep -v : | sort -g -k 9 | gzip > results_filtered/eur_ptsd_ehronly_pcs_v4_aug3_2021.fuma.gz

   totalN=191915.00
   percentN=0.8
   cat metal_results/eur_ptsd_pgconly_pcs_v4_aug3_2021_*.gz1.tbl   |  awk  -v totalN=$totalN -v percentN=$percentN   '{ if (NR ==1 || ($6 >= 0.01 && $6 <= 0.99 && $3 != "MarkerName" && $10 >= percentN*totalN)) print $1,$2,$3,$4,$5,$6,$10,$11,$12}'   | grep -v : | sort -g -k 9 | gzip > results_filtered/eur_ptsd_pgconly_pcs_v4_aug3_2021.fuma.gz

 #Main, no Freeze 2.5 datasets (EHR + MVP ONLY)
     totalN=409843.00
   percentN=0.8
   cat metal_results/eur_ptsd_pcs_v4_aug3_2021_no25_*.gz1.tbl   |  awk  -v totalN=$totalN -v percentN=$percentN   '{ if (NR ==1 || ($6 >= 0.01 && $6 <= 0.99 && $3 != "MarkerName" && $10 >= percentN*totalN)) print $1,$2,$3,$4,$5,$6,$10,$11,$12}'   | grep -v : | sort -g -k 9 | gzip > results_filtered/eur_ptsd_pcs_v4_aug3_2021_no25.fuma.gz



 #Main, EHR datasets only

 
 #Case control
  #Run 00_translate_meta_pluses_to_samplesize.sh to get exact Ns
  
  #Format case/control analysis: First combine all files

  #European
  cat metal_results/eur_ptsdcasecontrol_pcs_v4_aug3_2021_*.gz1.tbl  |  awk   '{ if (NR ==1 || ($ >= 0.01 && $6 <= 0.99 && $3 != "MarkerName" )) print }'   | grep -v : | sort -g -k 12 | gzip > eur_ptsdcasecontrol_pcs_v4_aug3_2021.temp1.tbl 
  totalN=463244
  percentN=0.8
  cat eur_ptsdcasecontrol_pcs_v5_jan4_2022.temp1.tbl.neffX |  awk  -v totalN=$totalN -v percentN=$percentN   '{ if (NR ==1 || ($6 >= 0.01 && $6 <= 0.99 && $3 != "MarkerName" && $16 >= percentN*totalN)) print $1,$2,$3,$4,$5,$6,$10,$11,$12,$14,$15,$16,$17}'   | grep -v : | gzip > eur_ptsdcasecontrol_pcs_v5_jan4_2022X.fuma.gz
  
    cat metal_results/eur_ptsdcasecontrolnomvp_pcs_v5_jan4_2021*.gz1.tbl  |  awk   '{ if (NR ==1 || ($6 >= 0.01 && $6 <= 0.99 && $3 != "MarkerName" )) print }'   | grep -v : | sort -g -k 12 | gzip > eur_ptsdcasecontrolnomvp_pcs_v4_aug3_2021.temp1.tbl 
  totalN=342624
  percentN=0.8
  cat eur_ptsdcasecontrolnomvp_pcs_v5_jan4_2022.tbl.neff  |  awk  -v totalN=$totalN -v percentN=$percentN   '{ if (NR ==1 || ($6 >= 0.01 && $6 <= 0.99 && $3 != "MarkerName" && $16 >= percentN*totalN)) print $1,$2,$3,$4,$5,$6,$10,$11,$12,$14,$15,$16,$17}'   | grep -v : | gzip > eur_ptsdcasecontrolnomvp_pcs_v5_jan4_2022.tbl.neff.fuma.gz
  
  #using eur_ptsdcasecontrol_pcs_v5_jan4_2022.fuma.gz and not eur_ptsdcasecontrol_pcs_v5_jan4_2022X.fuma.gz - eur_ptsdcasecontrol_pcs_v5_jan4_2022X.fuma.gz contains minor alterations to huckins sumstats. dont bother, result is worse somehow.
  
  
  
  #AAM
  cat metal_results/aam_ptsdcasecontrol_pcs_v5_jan4_2021_*.gz1.tbl  |  awk   '{ if (NR ==1 || ($3 != "MarkerName" )) print }'   | grep -v : | grep -v '???????'  | gzip > aam_ptsdcasecontrol_pcs_v5_jan4_2022.temp1.tbl 

  totalN=51252
  percentN=0.8
  cat aam_ptsdcasecontrol_pcs_v5_jan4_2022.temp1.tbl.neff |  awk  -v totalN=$totalN -v percentN=$percentN   '{ if (NR ==1 || ($6 >= 0.01 && $6 <= 0.99 && $3 != "MarkerName" && $16 >= percentN*totalN)) print $1,$2,$3,$4,$5,$6,$10,$11,$12,$14,$15,$16,$17}'  | sort -g -k 9 | grep -v :  | gzip > aam_ptsdcasecontrol_pcs_v5_jan4_2022.fuma.gz

#Special stuff:

 # AAM, no GTP or WACH EGHS ETC
 
   #AAM
   cat metal_results/aam_ptsd_pcs_v5_jan4_2021_nogtp_*.gz1.tbl | awk '{if (NR == 1 || ($3 != "MarkerName" )) print}' > metal_results/aam_ptsd_pcs_v5_jan4_2022_nogtp.gz1.tbl
   totalN=37201.00
   percentN=0.8
   cat metal_results/aam_ptsd_pcs_v5_jan4_2022_nogtp.gz1.tbl  |  awk  -v totalN=$totalN -v percentN=$percentN   '{ if (NR ==1 || ($6 >= 0.01 && $6 <= 0.99 && $3 != "MarkerName" && $10 >= percentN*totalN)) print $1,$2,$3,$4,$5,$6,$10,$11,$12}'   | grep -v : | sort -g -k 9 | gzip > results_filtered/aam_ptsd_pcs_v5_jan4_2022_nogtp.fuma.gz
  
  #F2.5 males QT and CC
   cat metal_results/eur_ptsd_pcs_males_f25qtrescale_v4_aug3_2021_*.tbl | awk '{if (NR == 1 || ($3 != "MarkerName" )) print}' > metal_results/eur_ptsd_pcs_males_f25qtrescale_v4_aug3_2021.gz1.tbl
   cat metal_results/eur_ptsd_pcs_males_f25qtrescale_v4_aug3_2021.gz1.tbl  | awk '{qcol=$13; gsub("?","",qcol); if (NR== 1|| length(qcol) >= 10) print}' | sort -g -k 12  |  awk  '{ if (NR ==1 || ($6 >= 0.01 && $6 <= 0.99 && $3 != "MarkerName")) print }'   | gzip > results_filtered/eur_ptsd_pcs_males_f25qtrescale_v4_aug3_2021.fuma.gz
  
   cat metal_results/eur_ptsd_pcs_females_f25qtrescale_v4_aug3_2021_*.tbl | awk '{if (NR == 1 || ($3 != "MarkerName" )) print}' > metal_results/eur_ptsd_pcs_females_f25qtrescale_v4_aug3_2021.gz1.tbl
   cat metal_results/eur_ptsd_pcs_females_f25qtrescale_v4_aug3_2021.gz1.tbl  | awk '{qcol=$13; gsub("?","",qcol); if (NR== 1|| length(qcol) >= 10) print}' | sort -g -k 12  |  awk  '{ if (NR ==1 || ($6 >= 0.01 && $6 <= 0.99 && $3 != "MarkerName")) print }'   | gzip > results_filtered/eur_ptsd_pcs_females_f25qtrescale_v4_aug3_2021.fuma.gz
  
   cat metal_results/eur_ptsd_pcs_males_f25cc_v4_aug3_2021_*.tbl | awk '{if (NR == 1 || ($3 != "MarkerName" )) print}' > metal_results/eur_ptsd_pcs_males_f25cc_v4_aug3_2021.gz1.tbl
   cat metal_results/eur_ptsd_pcs_males_f25cc_v4_aug3_2021.gz1.tbl  | awk '{qcol=$13; gsub("?","",qcol); if (NR== 1|| length(qcol) >= 5) print}' | sort -g -k 12  |  awk  '{ if (NR ==1 || ($6 >= 0.01 && $6 <= 0.99 && $3 != "MarkerName")) print }'   | gzip > results_filtered/eur_ptsd_pcs_males_f25cc_v4_aug3_2021.fuma.gz
  
   cat metal_results/eur_ptsd_pcs_females_f25cc_v4_aug3_2021_*.tbl | awk '{if (NR == 1 || ($3 != "MarkerName" )) print}' > metal_results/eur_ptsd_pcs_females_f25cc_v4_aug3_2021.gz1.tbl
   cat metal_results/eur_ptsd_pcs_females_f25cc_v4_aug3_2021.gz1.tbl  | awk '{qcol=$13; gsub("?","",qcol); if (NR== 1|| length(qcol) >= 5) print}' | sort -g -k 12  |  awk  '{ if (NR ==1 || ($6 >= 0.01 && $6 <= 0.99 && $3 != "MarkerName")) print }'   | gzip > results_filtered/eur_ptsd_pcs_females_f25cc_v4_aug3_2021.fuma.gz
   
  
  
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


#Reformat X chromosome data to allow for column tracking...


#Next, apply  the same size counting script

#Then filter on 80% neff



 #Filter data down to just LDSC SNPs then munge
 
 #euro case control
 LC_ALL=C join -1 1 -2 3  <(awk '{print $1}' ~/trauma_gwas/ldsc-master/ldsc_data/eur_w_ld_chr/w_hm3.snplist.sorted | LC_ALL=C sort -k1b,1 ) <(zcat eur_ptsdcasecontrol_pcs_v5_jan4_2022X.fuma.gz  | awk '{if (NR == 1) $3="SNP"; print}'   | LC_ALL=C sort -k3b,3 ) | gzip > results_filtered/eur_ptsdcasecontrol_pcs_v5_jan4_2022.fuma.gz.premunge.gz 
 python2 ~/trauma_gwas/ldsc-master/munge_sumstats.py --sumstats results_filtered/eur_ptsdcasecontrol_pcs_v5_jan4_2022.fuma.gz.premunge.gz  --N-col ntot --out results_filtered/eur_ptsdcasecontrol_pcs_v5_jan4_2022.fuma.gz.munge.gz #add --N-col OBS_CT for the sample size
 
 #aam case control
  LC_ALL=C join -1 1 -2 3  <(awk '{print $1}' ~/trauma_gwas/ldsc-master/ldsc_data/eur_w_ld_chr/w_hm3.snplist.sorted | LC_ALL=C sort -k1b,1 ) <(zcat results_filtered/aam_ptsdcasecontrol_pcs_v5_jan4_2022.fuma.gz  | awk '{if (NR == 1) $3="SNP"; print}'   | LC_ALL=C sort -k3b,3 ) | gzip > results_filtered/aam_ptsdcasecontrol_pcs_v5_jan4_2022.fuma.gz.premunge.gz 
 python2 ~/trauma_gwas/ldsc-master/munge_sumstats.py --sumstats results_filtered/aam_ptsdcasecontrol_pcs_v5_jan4_2022.fuma.gz.premunge.gz  --N-col ntot --out results_filtered/aam_ptsdcasecontrol_pcs_v5_jan4_2022.fuma.gz.munge.gz #add --N-col OBS_CT for the sample size
 
 #aam qt
   LC_ALL=C join -1 1 -2 3  <(awk '{print $1}' ~/trauma_gwas/ldsc-master/ldsc_data/eur_w_ld_chr/w_hm3.snplist.sorted | LC_ALL=C sort -k1b,1 ) <(zcat results_filtered/aam_ptsd_pcs_v5_jan4_2022.fuma.gz  | awk '{if (NR == 1) $3="SNP"; print}'   | LC_ALL=C sort -k3b,3 ) | gzip > results_filtered/aam_ptsd_pcs_v5_jan4_2022.fuma.gz.premunge.gz 
 python2 ~/trauma_gwas/ldsc-master/munge_sumstats.py --sumstats results_filtered/aam_ptsd_pcs_v5_jan4_2022.fuma.gz.premunge.gz    --out results_filtered/aam_ptsd_pcs_v5_jan4_2022.fuma.gz.munge.gz #add --N-col OBS_CT for the sample size
 
   LC_ALL=C join -1 1 -2 3  <(awk '{print $1}' ~/trauma_gwas/ldsc-master/ldsc_data/eur_w_ld_chr/w_hm3.snplist.sorted | LC_ALL=C sort -k1b,1 ) <(cat pgc-mdd2022-no23andMe-eur-v3.49.24.09.pgc.gz2 | awk '{if (NR == 1) $3="SNP"; print}'   | LC_ALL=C sort -k3b,3 ) | gzip > results_filtered/pgc-mdd2022-no23andMe-eur-v3.49.24.09.pgc.gz2.premunge.gz 
 python2 ~/trauma_gwas/ldsc-master/munge_sumstats.py --sumstats results_filtered/pgc-mdd2022-no23andMe-eur-v3.49.24.09.pgc.gz2.premunge.gz  --N-cas-col NCAS --N-con-col NCON  --out results_filtered/pgc-mdd2022-no23andMe-eur-v3.49.24.09.pgc.gz2.munge.gz #add --N-col OBS_CT for the sample size
 
 #eur 2.5
   LC_ALL=C join -1 1 -2 3  <(awk '{print $1}' ~/trauma_gwas/ldsc-master/ldsc_data/eur_w_ld_chr/w_hm3.snplist.sorted | LC_ALL=C sort -k1b,1 ) <(zcat results_filtered/eur_ptsd_pgconly_pcs_v4_aug3_2021.fuma.gz  | awk '{if (NR == 1) $3="SNP"; print}'   | LC_ALL=C sort -k3b,3 ) | gzip > results_filtered/eur_ptsd_pgconly_pcs_v4_aug3_2021.fuma.gz.premunge.gz 
 python2 ~/trauma_gwas/ldsc-master/munge_sumstats.py --sumstats results_filtered/eur_ptsd_pgconly_pcs_v4_aug3_2021.fuma.gz.premunge.gz    --out results_filtered/eur_ptsd_pgconly_pcs_v4_aug3_2021.fuma.gz.munge.gz #add --N-col OBS_CT for the sample size
 
  #eur EHR
   LC_ALL=C join -1 1 -2 3  <(awk '{print $1}' ~/trauma_gwas/ldsc-master/ldsc_data/eur_w_ld_chr/w_hm3.snplist.sorted | LC_ALL=C sort -k1b,1 ) <(zcat results_filtered/eur_ptsd_ehronly_pcs_v4_aug3_2021.fuma.gz  | awk '{if (NR == 1) $3="SNP"; print}'   | LC_ALL=C sort -k3b,3 ) | gzip > results_filtered/eur_ptsd_ehronly_pcs_v4_aug3_2021.fuma.gz.premunge.gz 
 python2 ~/trauma_gwas/ldsc-master/munge_sumstats.py --sumstats results_filtered/eur_ptsd_ehronly_pcs_v4_aug3_2021.fuma.gz.premunge.gz    --out results_filtered/eur_ptsd_ehronly_pcs_v4_aug3_2021.fuma.gz.munge.gz #add --N-col OBS_CT for the sample size
 
 #MVP PCL
 LC_ALL=C join -1 1 -2 1  <(awk '{print $1}' ~/trauma_gwas/ldsc-master/ldsc_data/eur_w_ld_chr/w_hm3.snplist.sorted | LC_ALL=C sort -k1b,1 ) <(zcat sumstats/TotalPCL_MVP_eur.gz  | awk '{if (NR == 1) $1="SNP"; print}'   | LC_ALL=C sort -k1b,1 ) | gzip > results_filtered/TotalPCL_MVP_eur.gz.premunge.gz 
 python2 ~/trauma_gwas/ldsc-master/munge_sumstats.py --sumstats results_filtered/TotalPCL_MVP_eur.gz.premunge.gz  --N 186689  --out results_filtered/TotalPCL_MVP_eur.gz.munge.gz #add --N-col OBS_CT for the sample size
 
 
  #MVP C/C
 LC_ALL=C join -1 1 -2 1  <(awk '{print $1}' ~/trauma_gwas/ldsc-master/ldsc_data/eur_w_ld_chr/w_hm3.snplist.sorted | LC_ALL=C sort -k1b,1 ) <(zcat sumstats/dbGAP_CaseControl_eur.harmonized.gz  | awk '{if (NR == 1) {$1="SNP" ; $4="A1DL"; $5="A2DL";}; print}'   | LC_ALL=C sort -k1b,1 ) | gzip > results_filtered/dbGAP_CaseControl_eur.harmonized.gz.premunge.gz 
 python2 ~/trauma_gwas/ldsc-master/munge_sumstats.py --sumstats results_filtered/dbGAP_CaseControl_eur.harmonized.gz.premunge.gz --a1 EA --a2 NEA --signed-sumstats logOR,0  --out results_filtered/dbGAP_CaseControl_eur.harmonized.gz.munge.gz #add --N-col OBS_CT for the sample size
 
 TotalPCL_MVP_eur.gz
 
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

#PTSD C/C vs MDD C/C
    python2 /home/maihofer/trauma_gwas/ldsc-master/ldsc.py \
 --rg results_filtered/eur_ptsdcasecontrol_pcs_v5_jan4_2022.fuma.gz.munge.gz.sumstats.gz,results_filtered/pgc-mdd2022-no23andMe-eur-v3.49.24.09.pgc.gz2.munge.gz.sumstats.gz \
 --ref-ld-chr /home/maihofer/trauma_gwas/ldsc-master/ldsc_data/eur_w_ld_chr/ \
 --w-ld-chr  /home/maihofer/trauma_gwas/ldsc-master/ldsc_data/eur_w_ld_chr/ \
 --out results_filtered/rg_eur_ptsd_pcs_v4_aug3_2021.gz1.tbl.munge.gz.sumstats.gz_pgc-mdd2022-no23andMe-eur-v3.49.24.09.pgc.gz2.munge.gz

#PTSD QT vs MDD C/C
    python2 /home/maihofer/trauma_gwas/ldsc-master/ldsc.py \
 --rg results_filtered/eur_ptsd_pcs_v4_aug3_2021.gz1.tbl.munge.gz.sumstats.gz,results_filtered/pgc-mdd2022-no23andMe-eur-v3.49.24.09.pgc.gz2.munge.gz.sumstats.gz\
 --ref-ld-chr /home/maihofer/trauma_gwas/ldsc-master/ldsc_data/eur_w_ld_chr/ \
 --w-ld-chr  /home/maihofer/trauma_gwas/ldsc-master/ldsc_data/eur_w_ld_chr/ \
 --out results_filtered/rg_eur_ptsd_pcs_v4_aug3_2021.gz1.tbl.munge.gz.sumstats.gz_eur_ptsdcasecontrol_pcs_v5_jan4_2022.fuma.gz.tbl.ldsc2

#PTSD QT 2.5 vs MDD C/C
    python2 /home/maihofer/trauma_gwas/ldsc-master/ldsc.py \
 --rg results_filtered/eur_ptsd_pgconly_pcs_v4_aug3_2021.fuma.gz.munge.gz.sumstats.gz,results_filtered/pgc-mdd2022-no23andMe-eur-v3.49.24.09.pgc.gz2.munge.gz.sumstats.gz\
 --ref-ld-chr /home/maihofer/trauma_gwas/ldsc-master/ldsc_data/eur_w_ld_chr/ \
 --w-ld-chr  /home/maihofer/trauma_gwas/ldsc-master/ldsc_data/eur_w_ld_chr/ \
 --out results_filtered/rg_eur_ptsd_pgconly_pcs_v4_aug3_2021.fuma.gz.munge.gz_eur_ptsdcasecontrol_pcs_v5_jan4_2022.fuma.gz.tbl.ldsc2


#PTSD EHR vs MDD C/C
    python2 /home/maihofer/trauma_gwas/ldsc-master/ldsc.py \
 --rg results_filtered/eur_ptsd_ehronly_pcs_v4_aug3_2021.fuma.gz.munge.gz.sumstats.gz,results_filtered/pgc-mdd2022-no23andMe-eur-v3.49.24.09.pgc.gz2.munge.gz.sumstats.gz\
 --ref-ld-chr /home/maihofer/trauma_gwas/ldsc-master/ldsc_data/eur_w_ld_chr/ \
 --w-ld-chr  /home/maihofer/trauma_gwas/ldsc-master/ldsc_data/eur_w_ld_chr/ \
 --out results_filtered/rg_eur_ptsd_ehronly_pcs_v4_aug3_2021.fuma.gz.munge.gz_eur_ptsdcasecontrol_pcs_v5_jan4_2022.fuma.gz.tbl.ldsc2

#PTSD MVP PCL vs MDD C/C
    python2 /home/maihofer/trauma_gwas/ldsc-master/ldsc.py \
 --rg results_filtered/TotalPCL_MVP_eur.gz.munge.gz.sumstats.gz,results_filtered/pgc-mdd2022-no23andMe-eur-v3.49.24.09.pgc.gz2.munge.gz.sumstats.gz\
 --ref-ld-chr /home/maihofer/trauma_gwas/ldsc-master/ldsc_data/eur_w_ld_chr/ \
 --w-ld-chr  /home/maihofer/trauma_gwas/ldsc-master/ldsc_data/eur_w_ld_chr/ \
 --out results_filtered/rg_TotalPCL_MVP_eur.gz_eur_ptsdcasecontrol_pcs_v5_jan4_2022.fuma.gz.tbl.ldsc2

#PTSD MVP C/C vs MDD C/C
    python2 /home/maihofer/trauma_gwas/ldsc-master/ldsc.py \
 --rg results_filtered/dbGAP_CaseControl_eur.harmonized.gz.munge.gz.sumstats.gz,results_filtered/pgc-mdd2022-no23andMe-eur-v3.49.24.09.pgc.gz2.munge.gz.sumstats.gz\
 --ref-ld-chr /home/maihofer/trauma_gwas/ldsc-master/ldsc_data/eur_w_ld_chr/ \
 --w-ld-chr  /home/maihofer/trauma_gwas/ldsc-master/ldsc_data/eur_w_ld_chr/ \
 --out results_filtered/rg_dbGAP_CaseControl_eur.harmonized.gz.munge.gz_eur_ptsdcasecontrol_pcs_v5_jan4_2022.fuma.gz.tbl.ldsc2


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


