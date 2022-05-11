#!/bin/bash
Nscen=$(ls -l V100/ | wc -l)
seeds=( $(seq 1 $Nscen) )
scen=(1 2 3)
nodes=(100 200 300 500 1000 2000 3000 5000)
rs=$(which Rscript)

for f in ${nodes[@]}; do
	for ss in ${seeds[@]}; do
		for n in ${scen[@]}; do
			#touch V$f/seed$ss/RESULTS/lvm$n.txt
			#psrecord "$rs 00_MODELS_RUN/00_run_lvm4net.r -seed $ss -nodes $f -scen $n" --log V$f/seed$ss/RESULTS/lvm$n.txt --interval 0.5
			$rs 00_MODELS_RUN/00_run_lvm4net.r -seed $ss -nodes $f -scen $n
			echo "done scenario $n with $f nodes"
		done
	done
done

