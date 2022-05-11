#!/bin/bash

# 1] prerequirements
./00_COMPILATION/00_compile_and_check.sh

# 3] Run across all scenarios, methods and seeds [most consuming, probably you might want to run them at least in parallel across methods, for example using tmux or slurm].
#./00_MODELS_RUN/01_bash_svilf.sh
#./00_MODELS_RUN/01_bash_svilf_ada.sh
#./00_MODELS_RUN/01_bash_vblpcm.sh
#./00_MODELS_RUN/01_bash_vgb_init.sh
#./00_MODELS_RUN/01_bash_lvm4net.sh

# e.g.
tmux new-session -d -s "RunSvilfSIMS" 
for f in $(ls ./01_BASH_RUN/01*); do
	echo $f
	tmux new-window -n $f
	tmux send-keys $f 'C-m'
done

