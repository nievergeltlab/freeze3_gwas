#!/bin/bash

cd /home/maihofer/freeze3_gwas


 #conda create --name prscs  python=3.7.3
 
 #conda activate prscs
 #pip install scipy --user 
 #pip install  h5py --user 
 
#Process the training dataset into PRS-cs format 
 zcat  /home/maihofer/freeze3_gwas/results_filtered/eur_ptsd_pcs_v4_aug3_2021_no25.fuma.gz  | awk '{if (NR==1) {$3="SNP";$4="A1";$5="A2";$8="BETA";$9="P"}; print $3,toupper($4),toupper($5),$8,$9}' > eur_ptsd_pcs_v4_aug3_2021_no25.fuma.gz.prscs


#For each $study, run PRS-cs using weights from $infile
 infile=/home/maihofer/freeze3_gwas/prscs/eur_ptsd_pcs_v4_aug3_2021_no25.fuma.gz.prscs 
 n_gwas=409843 #number of participants
#Initial run with PRScs
for study in vets  auro ongb psy5 # $(cat prs_studylist.txt | grep -v betr)
do
 sbatch --array=1-22  --time=16:05:00 --error errandout/prscs_"$study"_%a.e --output errandout/prscs_"$study"_%a.o --export=ALL,study=$study,infile=$infile,n_gwas=$n_gwas prscs_command.slurm -D /home/maihofer/freeze3_gwas/prscs
done

#Combine all files per study
for study in vets  auro ongb psy5  #  $(cat prs_studylist.txt | grep auro)
do

 cat "$study"/"$study"_prscs_pst_eff_a1_b0.5_phi1e-02_chr*.txt | awk '{print $2,$4,$6}' > "$study"/"$study"_prscs_pst_eff_a1_b0.5_phi1e-02_allchr.txt
 awk '{print $1}' "$study"/"$study"_prscs_pst_eff_a1_b0.5_phi1e-02_allchr.txt  >  "$study"/"$study"_prscs_pst_eff_a1_b0.5_phi1e-02_allchr.txt.snplist

 for chr in {1..22}
 do
 cat "$study"/"$study"_prscs_pst_eff_a1_b0.5_phi1e-02_chr"$chr".txt | awk '{print $2,$4,$6}' > "$study"/"$study"_prscs_pst_eff_a1_b0.5_phi1e-02_chr"$chr".txt.plink
 done


done

#Calculate PRS
for study in  vets  auro ongb psy5  #  $(cat prs_studylist.txt  )
do
 sbatch  --time=16:05:00 --error errandout/prscal_"$study"_%a.e --output errandout/prscal_"$study"_%a.o --export=ALL,study=$study prscal_command.slurm -D /home/maihofer/freeze3_gwas/prscs
done

#Run association analysis
IFS=$'\n'

for study in  $(cat prs_studylist_substudy.txt  ) #grep ftca
do
 study1=$(echo $study | awk '{print $1}')
 study2=$(echo $study | awk '{print $2}')
 echo $study
 Rscript prstest_command.r $study1 $study2
done

cat */*.results > results_prs.txt

cat */*.results.beta > results_beta.txt
cat */*.results.se > results_se.txt

#add uKBB and MVP results to this file

R
library(metafor)
library(data.table)
d1 <- fread('results_prs.txt',data.table=F)

meta<- rma(yi=d1[,2],sei=d1[,3],slab=d1[,1],method="FE")

pdf('metaanalysis.pdf',7,7)
plot(meta)
dev.off()
pdf('forestplot.pdf',7,7)
forest(meta)
dev.off()


meta_noukb <- rma(yi=d1[-nrow(d1),2],sei=d1[-nrow(d1),3],slab=d1[-nrow(d1),1],method="FE")


#SE weighted average nagelkerke R2... 
h2=sum(d1$V6* 1/d1$V3^2 )/ sum(1/d1$V3^2)

h2=0.055
#To user: Put in these values:
K=0.10 #Population prevalence
P=0.12 #Sample prevalence
#h2=0.05 #h2 on the observed scale
#seh2=0.01 #standard error of h2 on the observed scale

#Do calculations
zv <- dnorm(qnorm(K))
h2_liab <- h2 * K^2 * ( 1 - K)^2 / P / (1-P) / zv^2
var_h2_liab <- ( seh2 * K^2 * ( 1 - K)^2 / P / (1-P) / zv^2) ^2

#Report liability scale h2snp and se 
h2_liab #0.0797
sqrt(var_h2_liab) #0.0159


#Quintile meta analysis

R
library(metafor)
library(data.table)
betas <- fread('results_beta.txt',data.table=F)
ses <- fread('results_se.txt',data.table=F)

meta2 <- rma(yi=betas$Q2,sei=ses$Q2,slab=betas$study,method="FE")
meta3 <- rma(yi=betas$Q3,sei=ses$Q3,slab=betas$study,method="FE")
meta4 <- rma(yi=betas$Q4,sei=ses$Q4,slab=betas$study,method="FE")
meta5 <- rma(yi=betas$Q5,sei=ses$Q5,slab=betas$study,method="FE")

ms <- as.data.frame(matrix(ncol=2,nrow=5))
ms[1,] <- c(0,0.00001)

 ms[2,1] <- meta2$beta
 ms[2,2] <- meta2$se

 ms[3,1] <- meta3$beta
 ms[3,2] <- meta3$se

 ms[4,1] <- meta4$beta
 ms[4,2] <- meta4$se

 ms[5,1] <- meta5$beta
 ms[5,2] <- meta5$se

names(ms) <- c("beta","se")


ms$OR <- exp(ms$beta)
ms$LCI = exp(ms$beta -1.96*ms$se)
ms$UCI = exp(ms$beta +1.96*ms$se)
ms$Decile <-c(1:5)

petrolblue <- rgb(62,81,192,maxColorValue=255 )

ms$color <- petrolblue

res2 <- ms

library(Hmisc)
library(fmsb)
library(data.table)
library(plotrix)
library(lmtest)

res2$color="blue"
res2 <- subset(res2,select=c(Decile,OR,LCI,UCI,color))


pdf('prs_decile_f3_meta.pdf',7,7)

plotCI(x=res2$Decile,y=res2$OR,li=res2$LCI,ui=res2$UCI,lwd=2,ylim=c(1,2.25),pch=19,cex.axis=1.25,xlab="PRS Decile",ylab="Quintile Odds Ratio (95% CI)",main="",cex.lab=1.45,col=res2$color,scol=alpha(res2$color,.4),sfrac=0,xaxt='n',yaxt='n',cex=1.8)
legend("topleft",col=c("white","blue","red"),legend=c("EHR + MVP -> PGC2.5"),bty="n",pch=19,cex=1.5)
axis(1,at=c(1:10),cex.axis=1.25)
axis(2,at=c(1,1.25,1.5,1.75,2,2.25,2.5,2.75), labels=c("1","1.25","1.5","1.75","2","2.25","2.5","2.75"),cex.axis=1.25)

dev.off()


#Meta analysis of all 

#Unit increases may not be on the same scale?

#Do association analyses

# R Script on phenotypes?

#PRS-csx not practical because most AAM training data are in 2.5?


#Make a file with the .bim name next to the study abbr
 
# python PRScsx-master/PRScsx.py --ref_dir=/home/maihofer/freeze3_gwas/ --bim_prefix=/home/pgcdac/DWFV2CJb8Piv_0116_pgc_data/pts/wave3/v1/gtpc/cobg_dir_genome_wide/pts_gtpc_mix_am-qc.hg19.ch.fl.bgs \
# --sst_file=eur_ptsd_pcs_v4_aug3_2021.fuma.gz.prscx,aam_ptsd_pcs_v5_jan4_2022_nogtp.fuma.gz.prscx --n_gwas=641533,37201 --pop=EUR,AFR \
# --out_dir=/home/maihofer/freeze3_gwas/prscspred/   --chrom=${SLURM_ARRAY_TASK_ID} --out_name=gtppred_${SLURM_ARRAY_TASK_ID} --phi=1e-2  --meta=True


# sbatch --array=1-21   --time=16:05:00 --error errandout/prscs_gtp_%a.e --output errandout/prscs_gtp_%a.o --export=ALL,study1=$study1,study2=$study2 07_prscs.sh 


#1-22