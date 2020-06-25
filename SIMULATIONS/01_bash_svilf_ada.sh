#!/bin/bash
scen=(1 2 3)
nodes=(100 200 300 500 1000 2000 3000 5000)

rs=$(which Rscript)
#++++++
# SVILF
#++++++
for f in ${nodes[@]}; do
	for n in ${scen[@]}; do
		touch V$f/RESULTS/svilf$n.txt
		psrecord "$rs run_svilf_ada.R $f $n" --log V$f/RESULTS/svilf_ada$n.txt --interval 0.5
		echo "done scenario $n with $f nodes"
	done
done
