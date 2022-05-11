#!/bin/bash
data=$(dir DATA/)
methods=(svilf svilf_ada)
rs=$(which Rscript)

for ww in ${data[@]}; do
	for m in ${methods[@]}; do
		#$rs spline_roc.R $ww $m
		$rs run_eval_ROC.R $ww $m
		$rs run_eval_predictions.R $ww $m
		echo "done $m on $ww"
	done
done

