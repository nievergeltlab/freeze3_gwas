# freeze3_gwas
Code base for freeze 3 GWAS

00_harmonize_ehr_data.sh  
Data manipulation scripts for reformatting specific summary data sets, including variant annotation using BCFtools+UCSC SNP lists, lifting data to hg19, dealing with duplicates, etc. This does not encompass file type differences (e.g. A1 versus Allele1), which are dealt with in the meta-analysis scripts  

00_pariwise_rgs_ehr.sh  
Calculate pairwise genetic correlations between EHR datasets and the MVP

##MiXeR pipeline
00_prepare_files.sh  
Reformat files for MiXeR analysis

01_mixer_univariate.sh  
Univariate MiXeR

02_mixer_univariate_plots.sh 
Summarize univariate MiXeR runs and produce output plots/tables

03_mixer_bivariate.sh  
Bivariate MiXeR analysis

04_mixer_bivariate_plots.sh  
Summarize and plot bivariate run results
