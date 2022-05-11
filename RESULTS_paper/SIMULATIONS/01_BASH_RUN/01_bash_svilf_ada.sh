#!/bin/bash
Nscen=100
seeds=( $(seq 2 $Nscen) )
scen=(1 2 3)
nodes=(100 200 300 500 1000 2000 3000)
#nodes=(100 200 300 500 1000 2000 3000 5000)
rs=$(which Rscript)

#++++++
# SVILF
#++++++
for ss in ${seeds[@]}; do
	for f in ${nodes[@]}; do
		for n in ${scen[@]}; do
			#touch V$f/seed$ss/RESULTS/svilf_ada$n.txt
			#psrecord "$rs 00_MODELS_RUN/00_run_svilf_ada.R -seed $ss -nodes $f -scen $n" --log V$f/seed$ss/RESULTS/svilf_ada$n.txt --interval 0.5
			$rs 00_MODELS_RUN/00_run_svilf_ada.R -seed $ss -nodes $f -scen $n
			$rs 02_COMPUTE_PREDICTIONS/00_compute_predictions.r -seed $ss -method 'svilf_ada' -nodes $f -scen $n
			echo "done scenario $n with $f nodes"
		done
	done
done

