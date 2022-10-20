#!/bin/bash
#SBATCH -n 1

#Each file gets encoded with a linking column

 # #Make preliminary linking data. THe header columns are known, but the order is not. This is to obtain the order.
 IFS=$'\n'
for study in $(cat f3_studies_forest.txt | tail -n+2 | grep biov)
do
  echo $study
  filename=$(echo $study | awk '{print $1}' | sed s'/XXX/1/g' ) #Just look at chr 1 to obtain the header order
  outfilename=$(echo $study | awk '{print $2}')
  zcat "$filename" | head -n1  | awk '{OFS="\t"}{$1=$1; print}'> results_filtered/forest_plots/"$outfilename".header  #Gives the sample header (order of columns)
  head -n1 f3_studies_forest.txt | cat - <(echo $study) > results_filtered/forest_plots/"$outfilename".coding #gives the coding meanings
done


#Loop over SNPs to reformat all
IFS=$'\n'
for fp in $(cat f3_hits_eur.txt | tail -n+2  )
do

  forestsnplocus=$(echo $fp | awk ' {print $1}') 
  forestsnp=$(echo $fp | awk ' {print $2}') 
  forestsnpchr=$(echo $fp | awk '{print $3}')
  echo $forestsnp
  #Loop over studies
  for study in $(cat f3_studies_forest.txt | tail -n+2 | grep biov )
  do

  
   filename=$(echo $study | awk '{print $1}' | sed s/XXX/$forestsnpchr/g ) #replace the XXX with the chromosome
   outfilename=$(echo $study | awk '{print $2}')
    echo $outfilename 
    cat results_filtered/forest_plots/"$outfilename".header <(zgrep -w -m1 $forestsnp $filename | awk '{OFS="\t"}{$1=$1; print}' ) >  results_filtered/forest_plots/"$outfilename"."$forestsnp" 
   
   #reformat data
   Rscript 0_forest_Reformater.r results_filtered/forest_plots/"$outfilename".coding results_filtered/forest_plots/"$outfilename"."$forestsnp" results_filtered/forest_plots/"$outfilename"."$forestsnp".formatted $forestsnp

 done
 
 done
 

 #WorkingDir=$(pwd)
 
 # sbatch  --time=2:05:00 --error errandout/snpformatter.e --output errandout/snpformatter.o --export=ALL 06_forest_plot.sh -d WorkingDir
  
  
 # # #Now do the plotting step
 for fp in $(cat f3_hits_eur.txt | tail -n+2  )
 do
   forestsnp=$(echo $fp | awk ' {print $2}') 
  cat  results_filtered/forest_plots/*."$forestsnp".formatted | awk '{if (NR==1 || ( $3 != "SNP")) print}' > results_filtered/forest_plots/"$forestsnp".allstudies
  Rscript 0_forest_plot.r results_filtered/forest_plots/"$forestsnp".allstudies f3_studies_forest_descriptors.txt results_filtered/forest_plots/"$forestsnp" eur
done
