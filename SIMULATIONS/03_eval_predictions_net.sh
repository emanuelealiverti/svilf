#!/bin/bash
scen=(1 2 3)
nodes=(100 200 300 500 1000 2000 3000 5000)
rs=$(which Rscript)

#++++++
# SVILF
#++++++
for f in ${nodes[@]}; do
	for n in ${scen[@]}; do
		$rs run_eval_graph.R $f $n
		echo "done scenario $n with $f nodes"
	done
done
