#!/bin/bash
rs=$(which Rscript)
dd=$(dir DATA)
methods=(svilf svilf_ada)

#++++++
# SVILF
#++++++
for f in ${dd[@]}; do
	for n in ${methods[@]}; do 
		$rs run_svilf.R $f $n
	done
done
