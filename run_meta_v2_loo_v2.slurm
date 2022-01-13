#!/bin/bash
nodeuse=15

cd /home/maihofer/freeze3_gwas

counter=1
for files in $(cat $metafile)
do
/home/maihofer/trauma_gwas/METAL-master/build/bin/metal $files > errandout/"$files"_errorlogs 2>&1 &
counter=$((counter+1))
if [ $counter  -eq 8 ]
then
 wait
 counter=1
fi

done

wait

 