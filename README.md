# freeze3_gwas
Code base for freeze 3 GWAS

## Data harmonization steps
00_harmonize_ehr_data.sh  
Data manipulation scripts for reformatting specific summary data sets, including variant annotation using BCFtools+UCSC SNP lists, lifting data to hg19, dealing with duplicates, etc. This does not encompass file type differences (e.g. A1 versus Allele1), which are dealt with in the meta-analysis scripts  

00_pariwise_rgs_ehr.sh  
Calculate pairwise genetic correlations between EHR datasets and the MVP

## Accounting and book keeping
00_sample_size_get.txt  
Using the .fam files and .pheno/.cov files, get exact counts of analyzed subjects - per the analyzed phenotype and per the case/control phenotype

05_count_sample_sizes.txt  
Crummier script that just takes the N listed in the output file. Check against results of previous, total N must match.


## MiXeR pipeline
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


## Meta-analysis
reformat.slurm  
Reformat summary stats to be split by chromosome, only take markers wiht >1% MAF

00_vetsa_gwas.sh  
VETSA twins cohort GWAS in Bolt LMM

01_f3_gwas_v1.sh  
Run individual GWAS (where I have data), meta-analyze results, preliminary LDSC analysis of results

  run_trauma_gwas_v2_freeze3.slurm,run_trauma_gwas_v2_freeze3_case_control.slurm  
  Job scripts for running GWAS
  
  run_meta_v2_loo_v2.slurm  
  Job script for meta-analysis
  
00_translate_meta_pluses_to_samplesize.sh  
Translate the +/- outputs from metal into exact sample sizes, to tabulate actual N/N cases/N controls for every marker.

02_examine_tophits.sh  
Grep out a list of top hits from the summar results 


## Forest plots
0_forest_Reformater.r  
Reformat study level summary statistics (for a select marker) to be used in forest plots

0_forest_plot.r  
Generate the forest plot for a select marker

06_forest_plot.sh  
User script to reformat study level sumstats and generate forest plots

## PRS 
07_prscs.sh  
Calculate PRS in individual studies, combine results, test association with PTSD

## Experimental
01_rescale_qt_for_ivw.sh  
For continuous phenotypes, rescale outputs to 0-1 range based on, for use in IVW meta-analysis
