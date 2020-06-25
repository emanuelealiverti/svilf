#!/bin/bash
scen=(1 2 3)
nodes=(100 200 300 500 1000 2000)
rs=$(which Rscript)
methods=(vbg)
#methods=(svilf svilf_ada vbg)
#++++++
# SVILF
#++++++
for f in ${nodes[@]}; do
	for n in ${scen[@]}; do
		for m in ${methods[@]}; do
			$rs run_predictions.R $f $n $m
			echo "done $m, scenario $n with $f nodes"
		done
	done
done
