# Application
This folder focuses on applying `svilf` to three different high-dimensional networks, available at [http://snap.stanford.edu/data/index.html](http://snap.stanford.edu/data/index.html). See Section 5 of the paper for further details.
Dataset are placed inside the `DATA` folder, each in its own separate folder, stored with an edge-list format.

# Run posterior inference 
Posterior inference can be performed launching the script `00_bash_svilf.sh` or, alternatively, the `slurm` job `slurm_svilf.job`. Both methods call the `R` script `run_svilf.R` multiple times, either with the uniform or adaptive method. Analysis with the probit link can be conduced passing `link = "probit"` to the function `svilf`.

# Evaluate predictions
Predictions are evaluated in terms of quality in recovering the network connectivity pattern, evaluated in terms of AUC and ROC curve; these quantities can be computed calling `01_compute_pred_eff.sh` and `02_eval_predictions.sh`.
