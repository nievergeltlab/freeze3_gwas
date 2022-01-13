#Compare all decently sized studies to MVP

#MVP 
LC_ALL=C join -1 1 -2 1  <(awk '{print $1}' w_hm3.snplist.sorted | LC_ALL=C sort -k1b,1 ) <(zcat sumstats/TotalPCL_MVP_eur.gz  | awk '{if (NR==1){$1 ="SNP"; $11="zdir"}{print}}' |  LC_ALL=C sort -u -k1b,1 ) > rgtemp/TotalPCL_MVP_eur.gz.premunge
python2 /home/maihofer/trauma_gwas/ldsc-master/munge_sumstats.py --chunksize 500000  --sumstats rgtemp/TotalPCL_MVP_eur.gz.premunge --N 186689 --out rgtemp/TotalPCL_MVP_eur.gz.munge.gz
zcat sumstats/TotalPCL_MVP_eur.gz | awk '{if(NR==1) print "SNP","A1","A2"; else print $1,$4,$5}' > rgtemp/TotalPCL_MVP_eur.gz.alleles

#74 DAI2 
LC_ALL=C join -1 1 -2 2  <(awk '{print $1}' w_hm3.snplist.sorted | LC_ALL=C sort -k1b,1 ) <(zcat sumstats/daner_iPSYCH2015_PTSDbroad_chr1to22_HRC_MAF01.gz  | awk '{if (NR==1){$7="FRQ"}{print}}' |  LC_ALL=C sort -u -k2b,2 ) > rgtemp/daner_iPSYCH2015_PTSDbroad_chr1to22_HRC_MAF01.gz.premunge
python2 /home/maihofer/trauma_gwas/ldsc-master/munge_sumstats.py --chunksize 500000  --sumstats rgtemp/daner_iPSYCH2015_PTSDbroad_chr1to22_HRC_MAF01.gz.premunge --merge-alleles rgtemp/TotalPCL_MVP_eur.gz.alleles --N 55804 --out rgtemp/daner_iPSYCH2015_PTSDbroad_chr1to22_HRC_MAF01.gz.munge.gz

#79 BIOV 
LC_ALL=C join -1 1 -2 1  <(awk '{print $1}' w_hm3.snplist.sorted | LC_ALL=C sort -k1b,1 ) <(zcat sumstats/ptsd_gwas_dosage_chr.all.PHENO2.glm.firth.additive.with.header.and.freqs.good.gz  | awk '{if (NR==1){$1="SNP"}{print}}' |  LC_ALL=C sort -u -k1b,1 ) > rgtemp/ptsd_gwas_dosage_chr.all.PHENO2.glm.firth.additive.with.header.and.freqs.good.gz.premunge
python2 /home/maihofer/trauma_gwas/ldsc-master/munge_sumstats.py --chunksize 500000  --sumstats rgtemp/ptsd_gwas_dosage_chr.all.PHENO2.glm.firth.additive.with.header.and.freqs.good.gz.premunge --merge-alleles rgtemp/TotalPCL_MVP_eur.gz.alleles --N 24266 --out rgtemp/ptsd_gwas_dosage_chr.all.PHENO2.glm.firth.additive.with.header.and.freqs.good.gz.munge.gz

#81 MGBB  
LC_ALL=C join -1 1 -2 1  <(awk '{print $1}' w_hm3.snplist.sorted | LC_ALL=C sort -k1b,1 ) <(zcat sumstats/pbk_eur_ptsd_gwas_broad_share.txt.gz  | awk '{if (NR==1){$1="SNP"}{print}}' |  LC_ALL=C sort -u -k1b,1 ) > rgtemp/pbk_eur_ptsd_gwas_broad_share.txt.gz.premunge
python2 /home/maihofer/trauma_gwas/ldsc-master/munge_sumstats.py --chunksize 500000  --sumstats rgtemp/pbk_eur_ptsd_gwas_broad_share.txt.gz.premunge --merge-alleles rgtemp/TotalPCL_MVP_eur.gz.alleles  --N 16113 --out rgtemp/pbk_eur_ptsd_gwas_broad_share.txt.gz.munge.gz

#85 HUNT - Reverse coded
LC_ALL=C join -1 1 -2 1  <(awk '{print $1}' w_hm3.snplist.sorted | LC_ALL=C sort -k1b,1 ) <(zcat sumstats/HUNT_ptsd_pgc_allchr_filter_info_results.txt.gz.fuma.gz  | awk '{if (NR==1){$1="SNP"}{print}}' |  LC_ALL=C sort -u -k1b,1 ) > rgtemp/HUNT_ptsd_pgc_allchr_filter_info_results.txt.gz.fuma.gz.premunge
python2 /home/maihofer/trauma_gwas/ldsc-master/munge_sumstats.py --chunksize 500000  --sumstats rgtemp/HUNT_ptsd_pgc_allchr_filter_info_results.txt.gz.fuma.gz.premunge --merge-alleles rgtemp/TotalPCL_MVP_eur.gz.alleles --a1 Allele2 --a2 Allele1 --N 11938 --out rgtemp/HUNT_ptsd_pgc_allchr_filter_info_results.txt.gz.fuma.gz.munge.gz

#86 AGDS - Reverse coded
LC_ALL=C join -1 1 -2 1  <(awk '{print $1}' w_hm3.snplist.sorted | LC_ALL=C sort -k1b,1 ) <(zcat sumstats/PTSDsum_AGDS_full_19052021.QIMRB.txt.gz  | awk '{SNP=$16;if(NR==1){SNP="SNP";$3="NS"};print SNP,$0}' |  LC_ALL=C sort -u -k1b,1 ) > rgtemp/PTSDsum_AGDS_full_19052021.QIMRB.txt.gz.premunge
python2   /home/maihofer/trauma_gwas/ldsc-master/munge_sumstats.py  --chunksize 500000  --sumstats rgtemp/PTSDsum_AGDS_full_19052021.QIMRB.txt.gz.premunge --merge-alleles rgtemp/TotalPCL_MVP_eur.gz.alleles --a1 Allele2 --a2 Allele1  --N 10687 --a1 Allele2 --a2 Allele1 --out rgtemp/PTSDsum_AGDS_full_19052021.QIMRB.txt.gz.munge.gz
 
#90 FING  - Reverse coded
LC_ALL=C join -1 1 -2 1  <(awk '{print $1}' w_hm3.snplist.sorted | LC_ALL=C sort -k1b,1 ) <(zcat sumstats/file_download_22122020_chia_yen_PTSD_broad_bothsex_hg19.gz.fuma.gz  | awk '{if (NR==1){$1="SNP"}{print}}' |  LC_ALL=C sort -u -k1b,1 ) > rgtemp/file_download_22122020_chia_yen_PTSD_broad_bothsex_hg19.gz.fuma.gz.premunge
python2 /home/maihofer/trauma_gwas/ldsc-master/munge_sumstats.py --chunksize 500000  --sumstats rgtemp/file_download_22122020_chia_yen_PTSD_broad_bothsex_hg19.gz.fuma.gz.premunge --merge-alleles rgtemp/TotalPCL_MVP_eur.gz.alleles --a1 Allele2 --a2 Allele1  --N 37725 --out rgtemp/file_download_22122020_chia_yen_PTSD_broad_bothsex_hg19.gz.fuma.gz.munge.gz

#91 UKB2 
LC_ALL=C join -1 1 -2 3  <(awk '{print $1}' w_hm3.snplist.sorted | LC_ALL=C sort -k1b,1 ) <(zcat sumstats/PTSD_EHR_Broad_2020_12_10_ReGENIE_Step2_WG.txt.regenie.gz  | LC_ALL=C sort -u -k3b,3 ) > rgtemp/PTSD_EHR_Broad_2020_12_10_ReGENIE_Step2_WG.txt.regenie.gz.premunge
python2 /home/maihofer/trauma_gwas/ldsc-master/munge_sumstats.py --chunksize 500000  --sumstats rgtemp/PTSD_EHR_Broad_2020_12_10_ReGENIE_Step2_WG.txt.regenie.gz.premunge --merge-alleles rgtemp/TotalPCL_MVP_eur.gz.alleles --N 36541 --out rgtemp/PTSD_EHR_Broad_2020_12_10_ReGENIE_Step2_WG.txt.regenie.gz.munge.gz

#98 ESBB  - Reverse coded
LC_ALL=C join -1 1 -2 3  <(awk '{print $1}' w_hm3.snplist.sorted | LC_ALL=C sort -k1b,1 ) <(zcat sumstats/PTSD_broad_EstBB_GWAS_results.txt.noinfo.gz  | awk '{if (NR==1){$3="SNP"}{print}}' |  LC_ALL=C sort -u -k3b,3 ) > rgtemp/PTSD_broad_EstBB_GWAS_results.txt.noinfo.gz.premunge
python2 /home/maihofer/trauma_gwas/ldsc-master/munge_sumstats.py --chunksize 500000  --sumstats rgtemp/PTSD_broad_EstBB_GWAS_results.txt.noinfo.gz.premunge --merge-alleles rgtemp/TotalPCL_MVP_eur.gz.alleles --a1 Allele2 --a2 Allele1 --N 71022 --out rgtemp/PTSD_broad_EstBB_GWAS_results.txt.noinfo.gz.munge.gz


for analysis in daner_iPSYCH2015_PTSDbroad_chr1to22_HRC_MAF01.gz.munge.gz ptsd_gwas_dosage_chr.all.PHENO2.glm.firth.additive.with.header.and.freqs.good.gz.munge.gz pbk_eur_ptsd_gwas_broad_share.txt.gz.munge.gz HUNT_ptsd_pgc_allchr_filter_info_results.txt.gz.fuma.gz.munge.gz PTSDsum_AGDS_full_19052021.QIMRB.txt.gz.munge.gz file_download_22122020_chia_yen_PTSD_broad_bothsex_hg19.gz.fuma.gz.munge.gz PTSD_EHR_Broad_2020_12_10_ReGENIE_Step2_WG.txt.regenie.gz.munge.gz PTSD_broad_EstBB_GWAS_results.txt.noinfo.gz.munge.gz
do
 outname=$(echo $analysis)
 echo $outname
  python2 /home/maihofer/trauma_gwas/ldsc-master/ldsc.py \
 --rg  rgtemp/TotalPCL_MVP_eur.gz.munge.gz.sumstats.gz,rgtemp/"$analysis".sumstats.gz \
 --ref-ld-chr /home/maihofer/trauma_gwas/ldsc-master/ldsc_data/eur_w_ld_chr/ \
 --w-ld-chr  /home/maihofer/trauma_gwas/ldsc-master/ldsc_data/eur_w_ld_chr/ \
 --out rgs/mvppcl_$outname
 done
  