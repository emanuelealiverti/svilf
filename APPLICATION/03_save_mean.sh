#!/bin/bash
data=$(ls DATA)
methods=(svilf svilf_ada)
rs=/home/emanuele/bin/R-3.6.2/bin/Rscript
#++++++
for ww in ${data[@]}; do
	for m in ${methods[@]}; do
		$rs save_postMeans.R $ww $m
		echo "done $m on $ww"
	done
done
