# freeze3_gwas
Code base for freeze 3 GWAS

00_harmonize_ehr_data.sh  
Data manipulation scripts for reformatting specific summary data sets, including variant annotation using BCFtools+UCSC SNP lists, lifting data to hg19, dealing with duplicates, etc. This does not encompass file type differences (e.g. A1 versus Allele1), which are dealt with in the meta-analysis scripts  

00_pariwise_rgs_ehr.sh 
Calculate pairwise genetic correlations between EHR datasets and the MVP
