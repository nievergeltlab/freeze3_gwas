
#List all A/T G/C alleles
#Filter to Freq < 40% or freq >60% only
#Compare A1Freq CHIA to other studies. 
#Find all markers with difference > 10%. If 1 - A1 CHIA matches better (within 10%), flip the allele and genotype coding.
#estonia might be a good bet - same strand as other data, most markers included, genetically similar


library(data.table)

finngen <- fread('90_fing/file_download_22122020_chia_yen_PTSD_broad_bothsex_hg19.gz.fuma.gz',data.table=F)
all_nofinngen <- fread('98_esbb/PTSD_broad_EstBB_GWAS_results.txt.noinfo.gz',data.table=F)

finngen$ambiguous <- 0
finngen[which(finngen$Allele2 == "A" & finngen$Allele2 == "T"),]$ambiguous <- 1
finngen[which(finngen$Allele2 == "T" & finngen$Allele2 == "A"),]$ambiguous <- 1
finngen[which(finngen$Allele2 == "C" & finngen$Allele2 == "G"),]$ambiguous <- 1
finngen[which(finngen$Allele2 == "G" & finngen$Allele2 == "C"),]$ambiguous <- 1

finn1 <- subset(finngen,ambiguous == 1)
esbb1 <- subset(esbb,SNPID %in% finn1$SNPID)

dm <- merge(finn1,esbb1,by="SNPID",suffixes=c("_finn","_esbb"))


  ALLELELABELS  Allele2 Allele1 
  PVALUELABEL   p.value
  EFFECT   BETA
  FREQLABEL AF_Allele2
#Load all genotype data

grep -P "A\tT" pgbd_unfiltered.bim > ambiguous_snps.txt
grep -P "T\tA" pgbd_unfiltered.bim >> ambiguous_snps.txt
grep -P "C\tG" pgbd_unfiltered.bim >> ambiguous_snps.txt
grep -P "G\tC" Hpgbd_unfiltered.bim >> ambiguous_snps.txt
