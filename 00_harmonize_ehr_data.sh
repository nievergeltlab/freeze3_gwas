#52 RCOG: Add rsids

    #the columns are: CHROM POS REF ALT ID
    cat vcf_header1b.txt vcf_header2.txt <(zcat PTSDdx.trauma_PTSDdx_oct2017EA.fuma.tbl.gz | awk '{gsub(/\./,":",$2); print}' | awk '{OFS="\t"}{if (NR>1) print $1,$3, ".", $12,$13 , "100", "PASS", "MVP="$2}' |  sort    -n -k1r,1 -k2,2  ) > PTSDdx.trauma_PTSDdx_oct2017EA.fuma.tbl.gz.vcf

    bcftools view  PTSDdx.trauma_PTSDdx_oct2017EA.fuma.tbl.gz.vcf -Oz -o PTSDdx.trauma_PTSDdx_oct2017EA.fuma.tbl.gz.vcf.gz

    tabix -f PTSDdx.trauma_PTSDdx_oct2017EA.fuma.tbl.gz.vcf.gz

    #Some will not be annotated, need to identify these for my merge steps..

    bcftools annotate  -a /mnt/ukbb/adam/tinnitus_gwas/All_20180423_hg19.vcf.gz -c INFO PTSDdx.trauma_PTSDdx_oct2017EA.fuma.tbl.gz.vcf.gz -o PTSDdx.trauma_PTSDdx_oct2017EA.fuma.tbl.gz.vcf.gz.annotated

    echo "ID SNP"  | awk 'BEGIN {OFS="\t"} {print $1,$2}' > snpheader.txt


    tail -n+82 PTSDdx.trauma_PTSDdx_oct2017EA.fuma.tbl.gz.vcf.gz.annotated  | grep RS | awk '{print $8}' | awk 'BEGIN {FS=";"}{OFS="\t"} {print $1,$2}' | sed 's/RS=/rs/g' | sed 's/MVP=//g'  | cat snpheader.txt -  | LC_ALL=C sort -k1b,1 > PTSDdx.trauma_PTSDdx_oct2017EA.fuma.tbl.gz.vcf.gz.annotated.success
    tail -n+82 PTSDdx.trauma_PTSDdx_oct2017EA.fuma.tbl.gz.vcf.gz.annotated | grep -v RS | awk '{print $8}' | awk 'BEGIN {FS=";"}{OFS="\t"} {print $1,$1}'                  | sed 's/MVP=//g'  | cat snpheader.txt - | LC_ALL=C sort -k1b,1   > PTSDdx.trauma_PTSDdx_oct2017EA.fuma.tbl.gz.vcf.gz.annotated.failed

    wc -l PTSDdx.trauma_PTSDdx_oct2017EA.fuma.tbl.gz
    wc -l PTSDdx.trauma_PTSDdx_oct2017EA.fuma.tbl.gz.vcf.gz.annotated
    wc -l PTSDdx.trauma_PTSDdx_oct2017EA.fuma.tbl.gz.vcf.gz.annotated.success
    wc -l PTSDdx.trauma_PTSDdx_oct2017EA.fuma.tbl.gz.vcf.gz.annotated.failed

    LC_ALL=C join -1 1 -2 2  PTSDdx.trauma_PTSDdx_oct2017EA.fuma.tbl.gz.vcf.gz.annotated.success  <(zcat PTSDdx.trauma_PTSDdx_oct2017EA.fuma.tbl.gz   | awk '{gsub(/\./,":",$2); print}' | LC_ALL=C sort -k2b,2 ) > PTSDdx.trauma_PTSDdx_oct2017EA.fuma.tbl.gz.success
    LC_ALL=C join -1 1 -2 2  PTSDdx.trauma_PTSDdx_oct2017EA.fuma.tbl.gz.vcf.gz.annotated.failed  <(zcat PTSDdx.trauma_PTSDdx_oct2017EA.fuma.tbl.gz  | awk '{gsub(/\./,":",$2); print}'  | LC_ALL=C sort -k2b,2) >  PTSDdx.trauma_PTSDdx_oct2017EA.fuma.tbl.gz.failed

    wc -l PTSDdx.trauma_PTSDdx_oct2017EA.fuma.tbl.gz.success
    wc -l  PTSDdx.trauma_PTSDdx_oct2017EA.fuma.tbl.gz.failed

    cat PTSDdx.trauma_PTSDdx_oct2017EA.fuma.tbl.gz.success PTSDdx.trauma_PTSDdx_oct2017EA.fuma.tbl.gz.failed | sort -g -k 10 | cut -d " " -f 2- | gzip > PTSDdx.trauma_PTSDdx_oct2017EA.fuma.tbl.gz.fuma.gz

  #the columns are: CHROM POS REF ALT ID
    cat /mnt/ukbb/adam/tinnitus_gwas/vcf_header1b.txt /mnt/ukbb/adam/tinnitus_gwas/vcf_header2.txt <(zcat PTSDdx.trauma_PTSDdx_oct2017AA.fuma.tbl.gz | awk '{gsub(/\./,":",$2); print}' | awk '{OFS="\t"}{if (NR>1) print $1,$3, ".", $12,$13 , "100", "PASS", "MVP="$2}' |  sort    -n -k1r,1 -k2,2  ) > PTSDdx.trauma_PTSDdx_oct2017AA.fuma.tbl.gz.vcf

    bcftools view  PTSDdx.trauma_PTSDdx_oct2017AA.fuma.tbl.gz.vcf -Oz -o PTSDdx.trauma_PTSDdx_oct2017AA.fuma.tbl.gz.vcf.gz

    tabix -f PTSDdx.trauma_PTSDdx_oct2017AA.fuma.tbl.gz.vcf.gz

    #Some will not be annotated, need to identify these for my merge steps..

    bcftools annotate  -a /mnt/ukbb/adam/tinnitus_gwas/All_20180423_hg19.vcf.gz -c INFO PTSDdx.trauma_PTSDdx_oct2017AA.fuma.tbl.gz.vcf.gz -o PTSDdx.trauma_PTSDdx_oct2017AA.fuma.tbl.gz.vcf.gz.annotated

    echo "ID SNP"  | awk 'BEGIN {OFS="\t"} {print $1,$2}' > snpheader.txt


    tail -n+82 PTSDdx.trauma_PTSDdx_oct2017AA.fuma.tbl.gz.vcf.gz.annotated  | grep RS | awk '{print $8}' | awk 'BEGIN {FS=";"}{OFS="\t"} {print $1,$2}' | sed 's/RS=/rs/g' | sed 's/MVP=//g'  | cat snpheader.txt -  | LC_ALL=C sort -k1b,1 > PTSDdx.trauma_PTSDdx_oct2017AA.fuma.tbl.gz.vcf.gz.annotated.success
    tail -n+82 PTSDdx.trauma_PTSDdx_oct2017AA.fuma.tbl.gz.vcf.gz.annotated | grep -v RS | awk '{print $8}' | awk 'BEGIN {FS=";"}{OFS="\t"} {print $1,$1}'                  | sed 's/MVP=//g'  | cat snpheader.txt - | LC_ALL=C sort -k1b,1   > PTSDdx.trauma_PTSDdx_oct2017AA.fuma.tbl.gz.vcf.gz.annotated.failed

    wc -l PTSDdx.trauma_PTSDdx_oct2017AA.fuma.tbl.gz
    wc -l PTSDdx.trauma_PTSDdx_oct2017AA.fuma.tbl.gz.vcf.gz.annotated
    wc -l PTSDdx.trauma_PTSDdx_oct2017AA.fuma.tbl.gz.vcf.gz.annotated.success
    wc -l PTSDdx.trauma_PTSDdx_oct2017AA.fuma.tbl.gz.vcf.gz.annotated.failed

    LC_ALL=C join -1 1 -2 2  PTSDdx.trauma_PTSDdx_oct2017AA.fuma.tbl.gz.vcf.gz.annotated.success  <(zcat PTSDdx.trauma_PTSDdx_oct2017AA.fuma.tbl.gz   | awk '{gsub(/\./,":",$2); print}' | sed 's/SNP/ID/g' | LC_ALL=C sort -k2b,2 ) > PTSDdx.trauma_PTSDdx_oct2017AA.fuma.tbl.gz.success
    LC_ALL=C join -1 1 -2 2  PTSDdx.trauma_PTSDdx_oct2017AA.fuma.tbl.gz.vcf.gz.annotated.failed  <(zcat PTSDdx.trauma_PTSDdx_oct2017AA.fuma.tbl.gz  | awk '{gsub(/\./,":",$2); print}'  | sed 's/SNP/ID/g'  | LC_ALL=C sort -k2b,2) >  PTSDdx.trauma_PTSDdx_oct2017AA.fuma.tbl.gz.failed

    wc -l PTSDdx.trauma_PTSDdx_oct2017AA.fuma.tbl.gz.success
    wc -l  PTSDdx.trauma_PTSDdx_oct2017AA.fuma.tbl.gz.failed

    cat PTSDdx.trauma_PTSDdx_oct2017AA.fuma.tbl.gz.success  | sort -g -k 10 | cut -d " " -f 2- | gzip > PTSDdx.trauma_PTSDdx_oct2017AA.fuma.tbl.gz.fuma.gz



#75 WTCS : Annotate RSids
    LC_ALL=C join -1 2 -2 3 test_wtc_dosage_HDS_SamClean_MAF.afreq test_wtc_dosage_HDS_SamClean_MAF.PHENO1.glm.linear | sed 's/#//g' | awk '{if(NR>1){ FRQ=$5; A1=$11;} if(A1==$10) { A2=$9 }  else if(A1 != $10){ A2=$10; FRQ=1-$5}; if(NR==1) {A2="A2"; FRQ="FRQ";} ; print $0,A2,FRQ}'   > test_wtc_dosage_HDS_SamClean_MAF.PHENO1.glm.linear.frq

    #the columns are: CHROM POS REF ALT ID
    #important to use ref and alt1 for this!
    cat vcf_header1b.txt vcf_header2.txt <(awk '{OFS="\t"}{if (NR>1) print $7,$8, ".", $3,$4 , "100", "PASS", "MVP="$1}' test_wtc_dosage_HDS_SamClean_MAF.PHENO1.glm.linear.frq |  sort    -n -k1r,1 -k2,2  ) > test_wtc_dosage_HDS_SamClean_MAF.PHENO1.glm.linear.frq.vcf

    bcftools view  test_wtc_dosage_HDS_SamClean_MAF.PHENO1.glm.linear.frq.vcf -Oz -o test_wtc_dosage_HDS_SamClean_MAF.PHENO1.glm.linear.frq.vcf.gz

    tabix -f test_wtc_dosage_HDS_SamClean_MAF.PHENO1.glm.linear.frq.vcf.gz

    #Some will not be annotated, need to identify these for my merge steps..

    bcftools annotate  -a /mnt/ukbb/adam/tinnitus_gwas/All_20180423_hg19.vcf.gz -c INFO test_wtc_dosage_HDS_SamClean_MAF.PHENO1.glm.linear.frq.vcf.gz -o test_wtc_dosage_HDS_SamClean_MAF.PHENO1.glm.linear.frq.vcf.gz.annotated

    echo "ID SNP"  | awk 'BEGIN {OFS="\t"} {print $1,$2}' > snpheader.txt


    tail -n+82 test_wtc_dosage_HDS_SamClean_MAF.PHENO1.glm.linear.frq.vcf.gz.annotated  | grep RS | awk '{print $8}' | awk 'BEGIN {FS=";"}{OFS="\t"} {print $1,$2}' | sed 's/RS=/rs/g' | sed 's/MVP=//g'  | cat snpheader.txt -  | LC_ALL=C sort -k1b,1 > test_wtc_dosage_HDS_SamClean_MAF.PHENO1.glm.linear.frq.vcf.gz.annotated.success
    tail -n+82 test_wtc_dosage_HDS_SamClean_MAF.PHENO1.glm.linear.frq.vcf.gz.annotated | grep -v RS | awk '{print $8}' | awk 'BEGIN {FS=";"}{OFS="\t"} {print $1,$1}'                  | sed 's/MVP=//g'  | cat snpheader.txt - | LC_ALL=C sort -k1b,1   > test_wtc_dosage_HDS_SamClean_MAF.PHENO1.glm.linear.frq.vcf.gz.annotated.failed

    wc -l test_wtc_dosage_HDS_SamClean_MAF.PHENO1.glm.linear.frq
    wc -l test_wtc_dosage_HDS_SamClean_MAF.PHENO1.glm.linear.frq.vcf.gz.annotated
    wc -l test_wtc_dosage_HDS_SamClean_MAF.PHENO1.glm.linear.frq.vcf.gz.annotated.success
    wc -l test_wtc_dosage_HDS_SamClean_MAF.PHENO1.glm.linear.frq.vcf.gz.annotated.failed

    LC_ALL=C join -1 1 -2 1  test_wtc_dosage_HDS_SamClean_MAF.PHENO1.glm.linear.frq.vcf.gz.annotated.success  <(LC_ALL=C sort -k1b,1 test_wtc_dosage_HDS_SamClean_MAF.PHENO1.glm.linear.frq) > test_wtc_dosage_HDS_SamClean_MAF.PHENO1.glm.linear.frq.success
    LC_ALL=C join -1 1 -2 1  test_wtc_dosage_HDS_SamClean_MAF.PHENO1.glm.linear.frq.vcf.gz.annotated.failed  <(LC_ALL=C sort -k1b,1 test_wtc_dosage_HDS_SamClean_MAF.PHENO1.glm.linear.frq) >  test_wtc_dosage_HDS_SamClean_MAF.PHENO1.glm.linear.frq.failed

    wc -l test_wtc_dosage_HDS_SamClean_MAF.PHENO1.glm.linear.frq.success
    wc -l  test_wtc_dosage_HDS_SamClean_MAF.PHENO1.glm.linear.frq.failed

    cat test_wtc_dosage_HDS_SamClean_MAF.PHENO1.glm.linear.frq.success test_wtc_dosage_HDS_SamClean_MAF.PHENO1.glm.linear.frq.failed | sort -g -k 18  | awk '{if(NR == 2  || ($20 >= 0.01 && $20 <= 0.99)) print $2,$1,$8,$9,$13,$12,$19,$20,$14,$15,$16,$17,$18}' | grep -v NA | gzip > test_wtc_dosage_HDS_SamClean_MAF.PHENO1.glm.linear.frq.fuma.gz
    echo "SNP ID CHROM POS TEST A1 A2 FRQ OBS_CT BETA SE T_STAT P" > header_x.txt

    zcat test_wtc_dosage_HDS_SamClean_MAF.PHENO1.glm.linear.frq.fuma.gz | cat header_x.txt - | gzip > test_wtc_dosage_HDS_SamClean_MAF.PHENO1.glm.linear.frq.fuma2.gz





#75A WTCS , as case control: Annotate RSids
    LC_ALL=C join -1 2 -2 3 test_wtc_dosage_HDS_SamClean_MAF.afreq test_wtc_dosage_HDS_SamClean_MAF.PHENO2.glm.logistic | sed 's/#//g' | awk '{if(NR>1){ FRQ=$5; A1=$11;} if(A1==$10) { A2=$9 }  else if(A1 != $10){ A2=$10; FRQ=1-$5}; if(NR==1) {A2="A2"; FRQ="FRQ";} ; print $0,A2,FRQ}'   > test_wtc_dosage_HDS_SamClean_MAF.PHENO2.glm.logistic.frq

    #the columns are: CHROM POS REF ALT ID
    #important to use ref and alt1 for this!
    cat vcf_header1b.txt vcf_header2.txt <(awk '{OFS="\t"}{if (NR>1) print $7,$8, ".", $3,$4 , "100", "PASS", "MVP="$1}' test_wtc_dosage_HDS_SamClean_MAF.PHENO2.glm.logistic.frq |  sort    -n -k1r,1 -k2,2  ) > test_wtc_dosage_HDS_SamClean_MAF.PHENO2.glm.logistic.frq.vcf

    bcftools view  test_wtc_dosage_HDS_SamClean_MAF.PHENO2.glm.logistic.frq.vcf -Oz -o test_wtc_dosage_HDS_SamClean_MAF.PHENO2.glm.logistic.frq.vcf.gz

    tabix -f test_wtc_dosage_HDS_SamClean_MAF.PHENO2.glm.logistic.frq.vcf.gz

    #Some will not be annotated, need to identify these for my merge steps..

    bcftools annotate  -a /mnt/ukbb/adam/tinnitus_gwas/All_20180423_hg19.vcf.gz -c INFO test_wtc_dosage_HDS_SamClean_MAF.PHENO2.glm.logistic.frq.vcf.gz -o test_wtc_dosage_HDS_SamClean_MAF.PHENO2.glm.logistic.frq.vcf.gz.annotated

    echo "ID SNP"  | awk 'BEGIN {OFS="\t"} {print $1,$2}' > snpheader.txt


    tail -n+82 test_wtc_dosage_HDS_SamClean_MAF.PHENO2.glm.logistic.frq.vcf.gz.annotated  | grep RS | awk '{print $8}' | awk 'BEGIN {FS=";"}{OFS="\t"} {print $1,$2}' | sed 's/RS=/rs/g' | sed 's/MVP=//g'  | cat snpheader.txt -  | LC_ALL=C sort -k1b,1 > test_wtc_dosage_HDS_SamClean_MAF.PHENO2.glm.logistic.frq.vcf.gz.annotated.success
    tail -n+82 test_wtc_dosage_HDS_SamClean_MAF.PHENO2.glm.logistic.frq.vcf.gz.annotated | grep -v RS | awk '{print $8}' | awk 'BEGIN {FS=";"}{OFS="\t"} {print $1,$1}'                  | sed 's/MVP=//g'  | cat snpheader.txt - | LC_ALL=C sort -k1b,1   > test_wtc_dosage_HDS_SamClean_MAF.PHENO2.glm.logistic.frq.vcf.gz.annotated.failed

    wc -l test_wtc_dosage_HDS_SamClean_MAF.PHENO2.glm.logistic.frq
    wc -l test_wtc_dosage_HDS_SamClean_MAF.PHENO2.glm.logistic.frq.vcf.gz.annotated
    wc -l test_wtc_dosage_HDS_SamClean_MAF.PHENO2.glm.logistic.frq.vcf.gz.annotated.success
    wc -l test_wtc_dosage_HDS_SamClean_MAF.PHENO2.glm.logistic.frq.vcf.gz.annotated.failed

    LC_ALL=C join -1 1 -2 1  test_wtc_dosage_HDS_SamClean_MAF.PHENO2.glm.logistic.frq.vcf.gz.annotated.success  <(LC_ALL=C sort -k1b,1 test_wtc_dosage_HDS_SamClean_MAF.PHENO2.glm.logistic.frq) > test_wtc_dosage_HDS_SamClean_MAF.PHENO2.glm.logistic.frq.success
    LC_ALL=C join -1 1 -2 1  test_wtc_dosage_HDS_SamClean_MAF.PHENO2.glm.logistic.frq.vcf.gz.annotated.failed  <(LC_ALL=C sort -k1b,1 test_wtc_dosage_HDS_SamClean_MAF.PHENO2.glm.logistic.frq) >  test_wtc_dosage_HDS_SamClean_MAF.PHENO2.glm.logistic.frq.failed

    wc -l test_wtc_dosage_HDS_SamClean_MAF.PHENO2.glm.logistic.frq.success
    wc -l  test_wtc_dosage_HDS_SamClean_MAF.PHENO2.glm.logistic.frq.failed

    cat test_wtc_dosage_HDS_SamClean_MAF.PHENO2.glm.logistic.frq.success test_wtc_dosage_HDS_SamClean_MAF.PHENO2.glm.logistic.frq.failed | sort -g -k 18  | awk '{if(NR == 2  || ($20 >= 0.01 && $20 <= 0.99)) print $2,$1,$8,$9,$13,$12,$19,$20,$14,$15,$16,$17,$18}' | grep -v NA | gzip > test_wtc_dosage_HDS_SamClean_MAF.PHENO2.glm.logistic.frq.fuma.gz
    echo "SNP ID CHROM POS TEST A1 A2 FRQ OBS_CT OR LOG(OR)_SE Z_STAT P" > header_x2.txt

    zcat test_wtc_dosage_HDS_SamClean_MAF.PHENO2.glm.logistic.frq.fuma.gz | cat header_x2.txt - | gzip > test_wtc_dosage_HDS_SamClean_MAF.PHENO2.glm.logistic.frq.fuma2.gz





# 79 BIOV: PLINK coding adjustment: If A1 is equal to ALT, use REF as A2. IF A1 is equal to REF, use ALT as A2. MUST also change the AF
    awk '{if (NR == 1) print $0, "A2", "A1_Freq"; else if ($6==$5) print $0,$4,$13; else if ($6==$4) print $0,$5,1-$13}' 79_biov/ptsd_gwas_dosage_chr.all.PHENO2.glm.firth.additive.with.header.and.freqs > 79_biov/ptsd_gwas_dosage_chr.all.PHENO2.glm.firth.additive.with.header.and.freqs.good
  
    awk '{if (NR == 1) print $0, "A2", "A1_Freq"; else if ($6==$5) print $0,$4,$13; else if ($6==$4) print $0,$5,1-$13}' 79_biov/ptsd_gwas_dosage_chr.all.PHENO2.glm.firth.additive.with.header.and.freqs > 79_biov/ptsd_gwas_dosage_chr.all.PHENO2.glm.firth.additive.with.header.and.freqs.good
    
  #Chr X: 
  
    #Need to annotate the rs-ids.
  
    #the columns are: CHROM POS REF ALT ID
    #important to use ref and alt1 for this!
    cat /mnt/ukbb/adam/tinnitus_gwas/vcf_header1b.txt /mnt/ukbb/adam/tinnitus_gwas/vcf_header2.txt <(awk '{OFS="\t"}{if (NR>1) print $1,$2, ".", $4,$5 , "100", "PASS", "MVP="$3}' ptsd_chrX_all_broad.PHENO1.glm.logistic |  sort    -n -k1r,1 -k2,2  ) > ptsd_chrX_all_broad.PHENO1.glm.logistic.vcf

    bcftools view  ptsd_chrX_all_broad.PHENO1.glm.logistic.vcf -Oz -o ptsd_chrX_all_broad.PHENO1.glm.logistic.vcf.gz

    tabix -f ptsd_chrX_all_broad.PHENO1.glm.logistic.vcf.gz

    #Some will not be annotated, need to identify these for my merge steps..

    bcftools annotate  -a /mnt/ukbb/adam/tinnitus_gwas/All_20180423_hg19.vcf.gz -c INFO ptsd_chrX_all_broad.PHENO1.glm.logistic.vcf.gz -o ptsd_chrX_all_broad.PHENO1.glm.logistic.vcf.gz.annotated

    echo "ID SNP"  | awk 'BEGIN {OFS="\t"} {print $1,$2}' > snpheader.txt


    tail -n+82 ptsd_chrX_all_broad.PHENO1.glm.logistic.vcf.gz.annotated  | grep RS | awk '{print $8}' | awk 'BEGIN {FS=";"}{OFS="\t"} {print $1,$2}' | sed 's/RS=/rs/g' | sed 's/MVP=//g'  | cat snpheader.txt -  | LC_ALL=C sort -k1b,1 > ptsd_chrX_all_broad.PHENO1.glm.logistic.vcf.gz.annotated.success
    tail -n+82 ptsd_chrX_all_broad.PHENO1.glm.logistic.vcf.gz.annotated | grep -v RS | awk '{print $8}' | awk 'BEGIN {FS=";"}{OFS="\t"} {print $1,$1}'                  | sed 's/MVP=//g'  | cat snpheader.txt - | LC_ALL=C sort -k1b,1   > ptsd_chrX_all_broad.PHENO1.glm.logistic.vcf.gz.annotated.failed

    wc -l ptsd_chrX_all_broad.PHENO1.glm.logistic
    wc -l ptsd_chrX_all_broad.PHENO1.glm.logistic.vcf.gz.annotated
    wc -l ptsd_chrX_all_broad.PHENO1.glm.logistic.vcf.gz.annotated.success
    wc -l ptsd_chrX_all_broad.PHENO1.glm.logistic.vcf.gz.annotated.failed

    LC_ALL=C join -1 1 -2 3  ptsd_chrX_all_broad.PHENO1.glm.logistic.vcf.gz.annotated.success  <(LC_ALL=C sort -k3b,3 ptsd_chrX_all_broad.PHENO1.glm.logistic) > ptsd_chrX_all_broad.PHENO1.glm.logistic.success
  
  
   #Get a1 and a2 
   awk '{if (NR == 1) print $0, "A2"; else if ($7==$6) print $0,$5; else if ($7==$5) print $0,$6}' ptsd_chrX_all_broad.PHENO1.glm.logistic.success  >  ptsd_chrX_all_broad.PHENO1.glm.logistic.success.a2
   
   #get frequencies from a reference file.
   zcat daner_iPSYCH2015_PTSDbroad_chrX_HRC_MAF01.gz  |  awk '{print $2,$7,$4,$5}'    >  temp/daip_AF.txt
   
   #Join frequencies.. do flip if 
   LC_ALL=C join -1 2 -2 1  <(LC_ALL=C sort -k2b,2 ptsd_chrX_all_broad.PHENO1.glm.logistic.success.a2) <(LC_ALL=C sort -k1b,1 temp/daip_AF.txt) \
   | awk '{FRQNEW=$15; if(NR>1){if($16!=$7){FRQNEW=(1-FRQNEW)}} print $0, FRQNEW}' \
   | cut -d " " -f 1-15,19  | sed 's/#//g'   > ptsd_chrX_all_broad.PHENO1.glm.logistic.success.a2.afs1
  
  #Determine A2
  
  
   awk '{if (NR == 1) print $0, "A2", "A1_Freq"; else if ($6==$5) print $0,$4,$13; else if ($6==$4) print $0,$5,1-$13}' 
   #VCF tools annotate the the rsids...then merge in a frequency...
   /home/pgcdac/DWFV2CJb8Piv_0116_pgc_data/pts/wave3/summary_stats/79_biov_v2//ptsd_chrX_all_broad.PHENO1.glm.logistic
   


 #For this, need rs-id annotated, along with frequency
  ptsd_chrX_all_broad.PHENO1.glm.logistic
#85 hunt: SNP annotation?? filtering? Find what I did

#90 FING : Annotate RSids and port to hg19. 

     #the columns are: CHROM POS REF ALT ID
     #important to use ref and alt1 for this!

 zcat file_download_22122020_chia_yen_PTSD_broad_bothsex.gz | awk '{SNPID=$3; SNPID2=$3; gsub ("_","\ ",SNPID) ;gsub ("_",":",SNPID2) ;  if(NR==1) print $0,"SNPID", "CHR","POS","A1","A2"; else print $0, SNPID2,SNPID}' > file_download_22122020_chia_yen_PTSD_broad_bothsex.gz2
     
 cat /mnt/ukbb/adam/tinnitus_gwas/vcf_header1b.txt /mnt/ukbb/adam/tinnitus_gwas/vcf_header2.txt <(awk '{OFS="\t"}{if (NR>1) print $26,$27, ".", $28,$29 , "100", "PASS", "MVP="$25}' file_download_22122020_chia_yen_PTSD_broad_bothsex.gz2 | sed 's/chr//g'  |  sort    -n -k1r,1 -k2,2  ) > file_download_22122020_chia_yen_PTSD_broad_bothsex.gz.vcf

 bcftools view  file_download_22122020_chia_yen_PTSD_broad_bothsex.gz.vcf -Oz -o file_download_22122020_chia_yen_PTSD_broad_bothsex.gz.vcf.gz

 tabix -f file_download_22122020_chia_yen_PTSD_broad_bothsex.gz.vcf.gz


 #Some will not be annotated, need to identify these for my merge steps..

 bcftools annotate  -a /mnt/ukbb/adam/tinnitus_gwas/All_20180418_hg38.vcf.gz -c INFO file_download_22122020_chia_yen_PTSD_broad_bothsex.gz.vcf.gz -o file_download_22122020_chia_yen_PTSD_broad_bothsex.gz.vcf.gz.annotated

 echo "SNPID SNP"  | awk 'BEGIN {OFS="\t"} {print $1,$2}' > snpheader.txt


 tail -n+82 file_download_22122020_chia_yen_PTSD_broad_bothsex.gz.vcf.gz.annotated  | grep RS | awk '{print $8}' | awk 'BEGIN {FS=";"}{OFS="\t"} {print $1,$2}' | sed 's/RS=/rs/g' | sed 's/MVP=//g'  | cat snpheader.txt -  | LC_ALL=C sort -k1b,1 > file_download_22122020_chia_yen_PTSD_broad_bothsex.gz.vcf.gz.annotated.success 
 
 tail -n+82 file_download_22122020_chia_yen_PTSD_broad_bothsex.gz.vcf.gz.annotated | grep -v RS | awk '{print $8}' | awk 'BEGIN {FS=";"}{OFS="\t"} {print $1,$1}'                  | sed 's/MVP=//g'  | cat snpheader.txt - | LC_ALL=C sort -k1b,1   > file_download_22122020_chia_yen_PTSD_broad_bothsex.gz.vcf.gz.annotated.failed

 wc -l file_download_22122020_chia_yen_PTSD_broad_bothsex.gz2
 wc -l file_download_22122020_chia_yen_PTSD_broad_bothsex.gz.vcf.gz.annotated
 wc -l file_download_22122020_chia_yen_PTSD_broad_bothsex.gz.vcf.gz.annotated.success
 wc -l file_download_22122020_chia_yen_PTSD_broad_bothsex.gz.vcf.gz.annotated.failed 


 LC_ALL=C join -1 1 -2 25  <(awk '{print "chr"$1,$2}' file_download_22122020_chia_yen_PTSD_broad_bothsex.gz.vcf.gz.annotated.success | cat snpheader.txt -)  <(LC_ALL=C sort -k25b,25 file_download_22122020_chia_yen_PTSD_broad_bothsex.gz2) > file_download_22122020_chia_yen_PTSD_broad_bothsex.gz.success
 #LC_ALL=C join -1 1 -2 25  file_download_22122020_chia_yen_PTSD_broad_bothsex.gz.vcf.gz.annotated.failed  <(LC_ALL=C sort -k25b,25 file_download_22122020_chia_yen_PTSD_broad_bothsex.gz2) >  file_download_22122020_chia_yen_PTSD_broad_bothsex.gz.failed


 #We are at the point where all markers are annotated, time to lift over. 
 #Convert to UCSC bed format
 #I am only interseted in the rs-id markers
 #Note that I do -1 for the starting position of the marker, this is what was done in  # BEGIN{OFS="\t"} {if(NR==1) print "chrom","chromStart","chromEnd";  #chrID is some junk that i accidentally merged
 tail -n+83 file_download_22122020_chia_yen_PTSD_broad_bothsex.gz.vcf.gz.annotated.success  | sed 's/:/ /g' | awk '{print "chr"$1,$2-1,$2,$5}' | grep -v chrID > file_download_22122020_chia_yen_PTSD_broad_bothsex.gz.vcf.gz.annotated.bed

 #Liftover positions
  ./liftOver file_download_22122020_chia_yen_PTSD_broad_bothsex.gz.vcf.gz.annotated.bed  hg38ToHg19.over.chain.gz file_download_22122020_chia_yen_PTSD_broad_bothsex.gz.vcf.gz.annotated.bed.lift file_download_22122020_chia_yen_PTSD_broad_bothsex.gz.vcf.gz.annotated.bed.unmapped
 #liftover rsid - not going to lift back - if an rs-id was deleted, its not real - rsids seem to be updated within builds also , so no telling which is used
 
 #merge in the new positions. Just 
  echo "CHRNEW POSNEWMIN1 BPNEW SNP" > liftheader.txt 
  LC_ALL=C join -1 4 -2 2 <(cat liftheader.txt file_download_22122020_chia_yen_PTSD_broad_bothsex.gz.vcf.gz.annotated.bed.lift | LC_ALL=C sort -k 4b,4) <(LC_ALL=C sort -k2b,2 file_download_22122020_chia_yen_PTSD_broad_bothsex.gz.success ) > file_download_22122020_chia_yen_PTSD_broad_bothsex.gz.success_lifted
 wc -l file_download_22122020_chia_yen_PTSD_broad_bothsex.gz.success_lifted

 #Get rid of 'chr' prefix..
 cat file_download_22122020_chia_yen_PTSD_broad_bothsex.gz.success_lifted |  sort -g -k 23  | awk '{if(NR == 1 || ($13 >= 0.01 && $13 <= 0.99)) print $1,$2,$4,$10,$11,$13,$18,$19,$20,$21,$22,$23}' | sed 's/chr//g' | gzip > file_download_22122020_chia_yen_PTSD_broad_bothsex_july122021.gz.fuma.gz

#Now can be lifted over using tool..




#Give all markers rs-ids


Liftover the coordinates first.

All_20180418_hg38.vcf.gz


#95. Biom: Get missing A2. allele freqs

95_biom/PTSD_broad_EA_clean.head.assoc.logistic
#lookup the A1/A2 for each Rs then apply an algorith to find the missing allele
#Don't use the reference genome for this, just use another dataset imputed using the same source strand, otherwise the strand info will be fucked.
#It may be a big assumption that they're on the same strand though
zcat /mnt/ukbb/adam/ptsd/ehr_reformat/74_dai2/daner_iPSYCH2015_PTSDbroad_chr1to22_HRC_MAF01.gz | awk '{print $2,$4,$5,$7}' | LC_ALL=C sort -k1b,1 > /mnt/ukbb/adam/ptsd/ehr_reformat/74_dai2/daner_iPSYCH2015_PTSDbroad_chr1to22_HRC_MAF01.afs
LC_ALL=C join -1 1 -2 2 <(cat /mnt/ukbb/adam/ptsd/ehr_reformat/74_dai2/daner_iPSYCH2015_PTSDbroad_chr1to22_HRC_MAF01.afs) <(LC_ALL=C sort -k2b,2 95_biom/PTSD_broad_EA_clean.head.assoc.logistic) > 95_biom/PTSD_broad_EA_clean.head.assoc.logistic.merged
awk '{FRQ=$4; A1=$7; if(A1==$2) { A2=$3}  else if(A1 != $2){ A2=$2; FRQ=1-$4}; print $0,A2,FRQ}'   95_biom/PTSD_broad_EA_clean.head.assoc.logistic.merged |  grep -v -E "A T$|T A$|C G$|G C$" > 95_biom/PTSD_broad_EA_clean.head.assoc.logistic.merged.gwas1
awk '{print $1,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15}' 95_biom/PTSD_broad_EA_clean.head.assoc.logistic.merged.gwas1 > 95_biom/PTSD_broad_EA_clean.head.assoc.logistic.merged.gwas2

zcat /mnt/ukbb/adam/ptsd/ehr_reformat/74_dai2/daner_iPSYCH2015_PTSDbroad_chr1to22_HRC_MAF01.gz | awk '{print $2,$4,$5,$7}' | LC_ALL=C sort -k1b,1 > /mnt/ukbb/adam/ptsd/ehr_reformat/74_dai2/daner_iPSYCH2015_PTSDbroad_chr1to22_HRC_MAF01.afs
LC_ALL=C join -1 1 -2 2 <(cat /mnt/ukbb/adam/ptsd/ehr_reformat/74_dai2/daner_iPSYCH2015_PTSDbroad_chr1to22_HRC_MAF01.afs) <(LC_ALL=C sort -k2b,2 95_biom/PTSD_broad_AA_clean.head.assoc.logistic) > 95_biom/PTSD_broad_AA_clean.head.assoc.logistic.merged
awk '{FRQ=$4; A1=$7; if(A1==$2) { A2=$3}  else if(A1 != $2){ A2=$2; FRQ=1-$4}; print $0,A2,FRQ}'   95_biom/PTSD_broad_AA_clean.head.assoc.logistic.merged |  grep -v -E "A T$|T A$|C G$|G C$" > 95_biom/PTSD_broad_AA_clean.head.assoc.logistic.merged.gwas1
awk '{print $1,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15}' 95_biom/PTSD_broad_AA_clean.head.assoc.logistic.merged.gwas1 > 95_biom/PTSD_broad_AA_clean.head.assoc.logistic.merged.gwas2




#96 WTCM: Get missing A2, allele freqs. god damn it, the markers are only on by position... jesus christ...
zcat daner_iPSYCH2015_PTSDbroad_chr1to22_HRC_MAF01.gz | awk '{CHRBP=$1":"$3 ; if(NR==1){CHRBP="SNP"} ; print CHRBP,$2,$4,$5,$7}' | LC_ALL=C sort -k1b,1 > daner_iPSYCH2015_PTSDbroad_chr1to22_HRC_MAF01.afs
LC_ALL=C join -1 1 -2 2 <(cat daner_iPSYCH2015_PTSDbroad_chr1to22_HRC_MAF01.afs) <(zcat PGC_WTC.W.assoc.logistic.gz | awk '{if (NR ==1 || ($5 == "ADD" && $9!= "NA") ) print}' | LC_ALL=C sort -k2b,2 )  | sort -g -k 13 > PGC_WTC.W.assoc.logistic.merged
awk '{FRQ=$5; A1=$8; if(A1==$3) { A2=$4}  else if(A1 != $3){ A2=$3; FRQ=1-$5}; print $0,A2,FRQ}'   PGC_WTC.W.assoc.logistic.merged |  grep -v -E "A T$|T A$|C G$|G C$" > PGC_WTC.W.assoc.logistic.merged.gwas1
awk '{print $2,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15}' PGC_WTC.W.assoc.logistic.merged.gwas1 | gzip > PGC_WTC.W.assoc.logistic.merged.gwas2.gz

zcat /mnt/ukbb/adam/ptsd/ehr_reformat/74_dai2/daner_iPSYCH2015_PTSDbroad_chr1to22_HRC_MAF01.gz | awk '{CHRBP=$1":"$3 ; if(NR==1){CHRBP="SNP"} ; print CHRBP,$2,$4,$5,$7}' | LC_ALL=C sort -k1b,1 > daner_iPSYCH2015_PTSDbroad_chr1to22_HRC_MAF01.afs
LC_ALL=C join -1 1 -2 2 <(cat daner_iPSYCH2015_PTSDbroad_chr1to22_HRC_MAF01.afs) <(zcat PGC_WTC.AA.assoc.logistic.gz | awk '{if (NR ==1 || ($5 == "ADD" && $9!= "NA") ) print}' | LC_ALL=C sort -k2b,2 )  | sort -g -k 13 > PGC_WTC.AA.assoc.logistic.merged
awk '{FRQ=$5; A1=$8; if(A1==$3) { A2=$4}  else if(A1 != $3){ A2=$3; FRQ=1-$5}; print $0,A2,FRQ}'   PGC_WTC.AA.assoc.logistic.merged |  grep -v -E "A T$|T A$|C G$|G C$" > PGC_WTC.AA.assoc.logistic.merged.gwas1
awk '{print $2,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15}' PGC_WTC.AA.assoc.logistic.merged.gwas1 | gzip > PGC_WTC.AA.assoc.logistic.merged.gwas2.gz





#98 ESBB: Imputation info column is NA , need to remove it
zcat 98_esbb/PTSD_broad_EstBB_GWAS_results.txt.gz  | awk '{if (NR==1){$8="infNA"}; print}' | gzip > 98_esbb/PTSD_broad_EstBB_GWAS_results.txt.noinfo.gz
cat PTSD_broad_EstBB_GWAS_saige_chr23.txt | awk '{if (NR==1){$8="infNA"}; print}' | gzip > PTSD_broad_EstBB_GWAS_saige_chr23.txt.noinfo.gz


#99 MAYO: Adjusted Mayo file to remove duplicate column names
#This is the correct one, per email with brandon on 5/4/21
zcat  99_mayo/PLINK/ptsd_def2.ALL.pcAdj.MAF0.01.glm.logistic.gz | awk '{if(NR==1){A2="A2"}; if(NR > 1 ){ if ($6==$5){A2=$4}; if($6==$4) {A2=$5}}; print $2,$3,$26,$13,$6,A2,$16,$17,$18,$19,$20,$21}' > 99_mayo/PLINK/ptsd_def2.ALL.pcAdj.MAF0.01.glm.logistic4

awk '{if($1=="23") print}' 99_mayo/PLINK/ptsd_def2.ALL.pcAdj.MAF0.01.glm.logistic4 > 99_mayo/PLINK/ptsd_def2.ALL.pcAdj.MAF0.01.glm.logistic4.chr23
#need to annotate rs-ids to this

#need to add rsids to x chromsoome... fuck..

    #the columns are: CHROM POS REF ALT ID
    #important to use ref and alt1 for this!
    cat /mnt/ukbb/adam/tinnitus_gwas/vcf_header1b.txt /mnt/ukbb/adam/tinnitus_gwas/vcf_header2.txt <(zcat ptsd_def2.ALL.pcAdj.MAF0.01.glm.logistic4.gz_23.gz | awk '{OFS="\t"}{if (NR>1) print $1,$2, ".", $5,$6 , "100", "PASS", "MVP="$1":"$2}'  |  sort    -n -k1r,1 -k2,2  ) > ptsd_def2.ALL.pcAdj.MAF0.01.glm.logistic4.gz_23.vcf

    bcftools view  ptsd_def2.ALL.pcAdj.MAF0.01.glm.logistic4.gz_23.vcf -Oz -o ptsd_def2.ALL.pcAdj.MAF0.01.glm.logistic4.gz_23.vcf.gz

    tabix -f ptsd_def2.ALL.pcAdj.MAF0.01.glm.logistic4.gz_23.vcf.gz

    #Some will not be annotated, need to identify these for my merge steps..

    bcftools annotate  -a /mnt/ukbb/adam/tinnitus_gwas/All_20180423_hg19.vcf.gz -c INFO ptsd_def2.ALL.pcAdj.MAF0.01.glm.logistic4.gz_23.vcf.gz -o ptsd_def2.ALL.pcAdj.MAF0.01.glm.logistic4.gz_23.vcf.gz.annotated

    echo "ID SNP"  | awk 'BEGIN {OFS="\t"} {print $1,$2}' > snpheader.txt


    tail -n+82 ptsd_def2.ALL.pcAdj.MAF0.01.glm.logistic4.gz_23.vcf.gz.annotated  | grep RS | awk '{print $8}' | awk 'BEGIN {FS=";"}{OFS="\t"} {print $1,$2}' | sed 's/RS=/rs/g' | sed 's/MVP=//g'  | cat snpheader.txt -  | LC_ALL=C sort -k1b,1 > ptsd_def2.ALL.pcAdj.MAF0.01.glm.logistic4.gz_23.vcf.gz.annotated.success
    tail -n+82 ptsd_def2.ALL.pcAdj.MAF0.01.glm.logistic4.gz_23.vcf.gz.annotated | grep -v RS | awk '{print $8}' | awk 'BEGIN {FS=";"}{OFS="\t"} {print $1,$1}'                  | sed 's/MVP=//g'  | cat snpheader.txt - | LC_ALL=C sort -k1b,1   > ptsd_def2.ALL.pcAdj.MAF0.01.glm.logistic4.gz_23.vcf.gz.annotated.failed

    wc -l ptsd_def2.ALL.pcAdj.MAF0.01.glm.logistic4.gz_23
    wc -l ptsd_def2.ALL.pcAdj.MAF0.01.glm.logistic4.gz_23.vcf.gz.annotated
    wc -l ptsd_def2.ALL.pcAdj.MAF0.01.glm.logistic4.gz_23.vcf.gz.annotated.success
    wc -l ptsd_def2.ALL.pcAdj.MAF0.01.glm.logistic4.gz_23.vcf.gz.annotated.failed

    LC_ALL=C join -1 1 -2 3  ptsd_def2.ALL.pcAdj.MAF0.01.glm.logistic4.gz_23.vcf.gz.annotated.success  <(LC_ALL=C sort -k3b,3 ptsd_def2.ALL.pcAdj.MAF0.01.glm.logistic4.gz_23) > ptsd_def2.ALL.pcAdj.MAF0.01.glm.logistic4.gz_23.success
  

