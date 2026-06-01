# CoV-immune

NOTE 6/1/26: this repo is in progress. All code and data to generate results are present, but container / workflow testing is ongoing. Workflow has been tested through "QTL mapping."

Environment prep
-----------------------
This project uses a Docker container to produce an environment similar to that used in the original analysis (i.e. R v4.4.1 and R package versions available in March, 2026). In order to run this container you will need [Docker](https://docs.docker.com/get-docker/) installed. 

Build the docker container:

```
docker build . -t imm 
```

Run the docker container, opening a terminal session within the container:

```
docker run -e PASSWORD=pw123 --rm -v $(pwd):/home/rstudio/work -p 8787:8787 -it imm /bin/bash
```

Navigate to the working directory: 

```
cd home/rstudio/work 
```

Data prep 
-----------------------

Run the lymphoid integration protocol: 

```
Rscript lymphoid_integration.R
```

This script integrates the three datasets in `source_data/lymphoid-combined.xlsx` into one dataset, stored in `derived_data/lymphoid-harmonization-3.csv`. At the moment, this data is integrated with the other immune and phenotype data manually in Excel, resulting in `source_data/F2_phenos_lymph3.xlsx`.

Prepare the Rqtl file: 

```
Rscript rqtl_file_prep.R
```

This integrates genotype data (`source_data/Cr_RB05_miniMUGA-013024_paddedIDs.csv`) and phenotype data (`source_data/F2_phenos_lymph3.xlsx`) into a .csv file formatted for analysis with `R/qtl` ([Broman et al., 2006](https://doi.org/10.1093/bioinformatics/btg112)), stored in `derived_data/RqtlCC006xCC044_ctrlAndSARS.csv`. This file and the `R/qtl` data structures are often used even when the `R/qtl` package is not being used for mapping. 

Prepare dataset with transformed and centered/scaled phenotypes (`derived_data/cross_data.csv`), as well as an `Rqtl` cross object with imputed genotypes (`derived_data/cross_imputed.csv`):

```
Rscript data_processing.R
```

Produce a version of the dataset with missing immune trait data imputed (`derived_data/cross_data_flow_imp.csv`):  

```
Rscript imputation.R
```

<!--
Set `--compareCV=TRUE` to run cross-validation (CV) on the following imputation models:
1. Bayesian imputation based on immune trait variance-covariance matrix $\Sigma$ only 
2. Bayesian imputation based on $\Sigma$ and fixed effect covariates (sex, infection)
3. Bayesian imputation based on $\Sigma$, fixed effect covariates and random polygenic effect 
4. Random Forest (RF) imputation with genotypes 
5. RF imputation without genotypes 

Set `--runCV=TRUE` to run cross-validation on only the chosen imputation method, method (2).
-->

# Data analysis 

Perform variance component analysis on a high-performance computing (HPC) cluster using SLURM (requests enough CPUs to run in parallel by trait):

```
sbatch var_comp.sh
```

or locally (runs in parallel using available CPUs):

```
Rscript var_comp.R
```

Results are stored in `results/var_comp_res.csv`, and visualized (in the following step) in Supp. Fig. 1E. 

*Note that variance component analysis is computationally intensive, with each trait potentially taking 1-3 hours, depending on if bootstrapping is necessary to define confidence intervals. It will be much faster to perform on an HPC cluster than locally so that analysis for all traits can be run completely in parallel. UNC's HPC cluster does not allow Docker containers - if this is the case for your cluster, you can either 1) FTP the repo to the cluster, manually ensure the necessary packages are installed, run the script, and FTP back the results file; or 2) create an Apptainer (fka Singularity) container from the Docker container (we do not provide any guidance on performing this analysis with Apptainer).*

Perform miscellaneous analyses/visualizations on trait data:

```
Rscript trait_eda.R
```

This script performs:
* Levene's test for variance heterogeneity (results printed to `results/eda.txt`)
* Classification of infection group by immune trait data with `glmnet` and random forest (results printed to `results/eda.txt`) 
* Frequentist variable importance analysis with `glmnet` (saved to `results/var_selection/frequentist/glmnet_var_importance.csv`)
* Create Figure 1 and Supp. Figure 1 (saved to `figures/Figure1.png` and `figures/supplemental/SuppFig1.png`)

## Variable selection 

Perform variable selection:

```
Rscript var_selection.R 
```

This is set up to run 4 chains of a Gibbs sampler with 5000 each, followed by 5 repeats of 10-fold cross-validation (CV). The Gibbs sampler is run in parallel by chain (4 chains), and CV is run in parallel by repeat (5 repeats). This script creates the following:
1. Log file in `results/var_selection/bayesian_log.txt`
2. R object with MCMC samples in `results/var_selection/res_dnm.RDS`. If this file is created, subsequent runs will skip the MCMC sampler and load the results from file. To re-run the sampler, you must delete this file. 
3. R object with CV results in `results/var_selection/cv_res_list.rds`. If this file is created, subsequent runs will skip the CV code and load the results from file. To re-run CV, you must delete this file. 
4. List of important phenotypes with stability information in `results/var_selection/important_phenos.csv`
5. Figure 2 in `figures/Figure2.png`
6. Supp. Figure 2 in `figures/supplemental/heatmap_labeled.pdf` and Supp. Figure 3 in `figures/supplemental/gibbs_performance.png`
7. Supp. Table 1 in `results/var_selection/STable1.xlsx`, which will be updated in later analyses 

## QTL mapping  
Perform QTL mapping on a high-performance computing cluster using SLURM (requests enough CPUs to run in parallel by trait): 

```
sbatch map_qtl.sh
```

*You must edit map_qtl.sh and set `hpc_dir` to the appropriate directory on your cluster.*

or locally (runs in parallel using available CPUs):

```
Rscript map_qtl.R
```

This code will generate: 
1. Log file in `results/qtl_mapping/mapping_log.txt`
2. QTL scan R objects in `results/qtl_mapping/modRDS` 
3. Permutation R objects in `results/qtl_mapping/permRDS`
4. QTL scan plots in `figures/QTL_scans/`.

*Note that due to permutation testing in combined QTL mapping, this is computationally intensive. With 1000 permutations, the set of analyses (stratified in control, SARS-CoV, and SARS-CoV-2, and combined) can take ~15 hours per trait. It will be much faster to run on an HPC cluster.*

Summarize QTL mapping results:

```
Rscript summarizeQTL.R
```

This code will generate:
1. QTL mapping results in `results/qtl_mapping/condensed_qtl.xlsx` (equivalent to Supp. Table 2, with formatting changes)
2. Figure 3 in `figures/Figure3.png`
3. Supp. Table 1 in `results/var_selection/STable1.xlsx` (this is primarily variable selection results, but is finalized with this script because it has a column `order` indicating the order of phenotypes on the y-axis in Fig. 3A)

## Mediation analysis  
Run mediation analysis and create Figure 4:

```
Rscript mediation.R
```

