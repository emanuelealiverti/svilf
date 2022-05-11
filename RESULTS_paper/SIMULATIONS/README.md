# Setup
Simulations require an installed version of `svilf`, provided in the root folder of this repository. Competitors require a working installation of the library `rstan`.

- `cd COMPILE`, run `./compile.R` to save a dynamic shared stan model (avoid recompiling every time)
- Install `svilf` package. We recommend installing the package locally, using the default location which is inside the `lib` folder at the root of this repostiory. Then, the library can be called as `library(svilf, lib.loc = "lib/")` or its relative path.

All simulations and post-processing operations are evaluation are implemented in `R`. We provide bash scripts calling `R` scripts multiple times. This format allows us to measure the computational time and memory requirements of the different calls separately.
All bash scripts should be executable (e.g., `chmod u+x *.sh`).
If you're using a high-performance cluster, the different scripts can be called within a process manager such as `slurm`; see `slurm_svilf.job` for an example.


# Generate artificial data
- Run `00_generate_data.R n1 n2 n3 ...` to generate syntethic data according to the models describes in Section 4 of the paper, where `n1 n2...` denotes the number of nodes. 
Each call generates a folder `Vnx`, a sub-folder `Vnx/RESULTS`,  and $3$ artificial networks, respectively generated from a latent factor model (s1 in the paper), a latent distance model (s2) and a stochastic block model (s3). Each network is saved twice, using a dense matrix and sparse edge-list representation, leveraging the `igraph` package.
The analysis of the paper focus on `n=100,200,500,1000,2000,5000`. 


# Run simulations
- run `01_bash_**.sh`. Each script is associated with a different model, called through the associated `run_**.R` file. 
Specifically, `01_bash_svilf.sh` (and `run_svilf.R`) refer to `svilf` with uniform weights, `01_bash_svilf_ada.sh` to `svilf` with adaptive weights and `01_bash_vgb_init.sh` to `advi` variational approach, implemented in the `stan` library.
Each bash script define a set of local variables, defined in the first lines of each file. `scen` and `nodes` specify the simulation scenarios and number of nodes over which simulations are conducted, while `rs` should point to your local `Rscript` installation (`rs=$(which Rscript)` should work oob, but in some settings it might be preferable to call a different `R` version.). 
The focus of the bash script is to evaluate RAM usage and CPU time for each run (results for the approximate posterior are saved at the end of the `R` script, and used later). This measurements are taken with the python utility `psrecord`, which should be available on your `$PATH` (`pip install psrecord` or  `python3.5 -m pip install psrecord --user` should do the job).
Results with the probit link function (Appendix A.2 of the paper) can be obtained passing `link = "probit"` to the function `svilf`.

# Compute predictions
- Posterior predictive probabilities are computed via Monte Carlo integration using the samples stored in the previous phase.
This step is likely to be time and memory consuming in large $n$ cases. Therefore, we compute it *after* vb routines, in order to make time and memory comparison consistent.
For $n<=2000$, we recommend using `02_compute_pred.sh`, which calls  `run_predictions.R`. For larger $n$ `02_compute_pred_eff.sh` is likely to perform better in terms of memory, even though its execution time might be longer.

# Evaluate predictions
The last step involves the evaluation of the quality of the computed predictions, calling the bash script `03_eval_predictions.sh`. See `run_eval_predictions.R` for details on the metrics considered (auc, precision and recall).
