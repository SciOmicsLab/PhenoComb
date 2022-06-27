#!/bin/bash

# export DISPLAY=:0

for i in {6..10}
do

echo "Job (${i} markers) started at `date +"%T %a %d %b %Y"`" > ./outputs/flowType_bm_$i.log
# /Library/Frameworks/R.framework/Versions/4.0/Resources/bin/
CMD="time Rscript run_flowType_benchmark2.R $i >> ./outputs/flowType_bm_$i.log 2>&1"
echo $CMD
eval "${CMD}"
#Rscript run_flowType_benchmark.R $"{i}" > ./outputs/flowType_bm_${i}.log 
echo "Job (${i} markers) finished at `date +"%T %a %d %b %Y"`" >> ./outputs/flowType_bm_$i.log

done
