#Most files have different names for the header column. We want to make sure we get the header column.
#The tophits file should have all header names included!
#From FUMA, get list of leading SNPs. 

#Grep list of top snps from the meta-results to have on hand.
zgrep -w -f f3_tophits.txt results_filtered/eur_ptsd_pcs_v3_jul6_2021.fuma.gz > toploci.results

#1) Pull out top hits for all studies (make a file where each row is an input filename)
# f3_inputfiles.txt is a list of all contributing study result filenames
#Don't let there be a blank row in the top hits file or it will copy everything!!

for files in $(cat f3_inputfiles.txt )
do
  if [ -e sumstats/$files ]
  then
    echo trying $files
    zgrep -w -f  f3_tophits.txt sumstats/$files > tophits_meta/$files.tophits
  fi
  if [ -e results_cat/$files ]
  then
   echo trying $files
    zgrep -w -f  f3_tophits.txt results_cat/$files > tophits_meta/$files.tophits
  fi
done

#2) Do detailed meta-analysis to check for heterogeneity, perform meta-regression
 #Take a list of top hits to make plots
 #eur_descriptors_july2_2021.csv gives study information, including description of the header for each file, so that they can be harmonized
for fp in  rs13161130 #rs13161130 rs295017 #rs762897 # $(cat f3_tophits.txt)
do

 echo "Extracting data for $fp"

Rscript AA_forest_plot_v2.r $fp eur_descriptors_july2_2021.csv forest_plots/"$fp"

done


#3) For each study, perform a sign test to see if effect directions for the top hits deviate from the average direction
  zgrep -w -f f3_tophits.txt results_filtered/eur_ptsd_pcs_v2_jun30_2021.fuma.gz > tophits_meta/meta_tophits.txt
  awk '{print $10}' tophits_meta/meta_tophits.txt |  awk '{gsub("+","1 "); print}' | awk '{gsub("-","-1 "); print}' | awk '{gsub("?","NA "); print}' > tophits_meta/meta_tophits.signs.txt

  module load R
  R
   library(data.table)
   library(EnvStats)

  #Read in the effect directions for each marker and study (row = marker, column = study)
   d1 <- fread('tophits_meta/meta_tophits.signs.txt',data.table=F)

  #Number of studies is the number of columns
   nstudy <- dim(d1)[2]

  #For each marker, calculate the typical sign (i.e. the mode)
    getmode <- function(v) {
       uniqv <- unique(v)
       uniqv[which.max(tabulate(match(v, uniqv)))]
    }
    avgsign <- apply(d1,1,getmode)

  #For each study, see if they deviate from the typical sign at each marker. Do a sign test on these deviations
  deviations <- matrix(ncol=2,nrow=nstudy)
   for (i in 1:dim(d1)[2])
   {
    #'sign' outputs values 1 and 0. Convert 0 into -1 for the sign test to work
    testvals <- sign(d1[,i] == avgsign)
    testvals[which(testvals==0)] <- -1
    #test in the direction that the sign DOES NOT match (less than expected matching)
    #Otherwise studies that often match the typical sign (i.e. studies with high weights) will be significant!
    res <- signTest(testvals,alternative="less")
    deviations[i,] <- c(res$statistic,res$p.value)
    rm(res)
   }
  #Output the test results into a file
   write.table(deviations,file="tophits_meta/signtest_tophits.txt",quote=F,row.names=F)

# 4)Take the genetic correlation between each dataset and the meta-analysis? Possibly leave one out.
