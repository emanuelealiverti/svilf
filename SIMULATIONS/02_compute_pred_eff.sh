#!/bin/bash
scen=(1 2 3)
nodes=(3000 5000)
rs=$(which Rscript)
methods=(vbg)
#++++++
# SVILF
#++++++
for f in ${nodes[@]}; do
	for n in ${scen[@]}; do
		for m in ${methods[@]}; do
			$rs run_predictions_eff.R $f $n $m
			echo "done $m, scenario $n with $f nodes"
		done
	done
done
