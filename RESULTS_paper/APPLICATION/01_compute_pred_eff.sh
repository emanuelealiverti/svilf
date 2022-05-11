#!/bin/bash
data=$(ls DATA)
methods=(svilf svilf_ada)
rs=$(which Rscript)
#++++++
for ww in ${data[@]}; do
	for m in ${methods[@]}; do
		$rs run_postMean.R $ww $m
		echo "done $m on $ww"
	done
done
