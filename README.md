# CoV-immune

NOTE 4/1/26: this repo is in progress. All code and data to generate results are present, but container / workflow testing is ongoing. 

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

Prepare Rqtl file: 

```
Rscript rqtl_file_prep.R
```

This integrates genotype data (`source_data/Cr_RB05_miniMUGA-013024_paddedIDs.csv`) and phenotype data (`source_data/F2_phenos_lymph3.xlsx`) into a .csv file formatted for analysis with `R/qtl` ([Broman et al., 2006](https://doi.org/10.1093/bioinformatics/btg112)), stored in `derived_data/RqtlCC006xCC044_ctrlAndSARS.csv`. This file and the `R/qtl` data structures are often used even when the `R/qtl` package is not being used for mapping. 

Prepare dataset with transformed and centered/scaled phenotypes (`derived_data/cross_data.csv`), as well as an `Rqtl` cross object with imputed genotypes:

```
Rscript data_processing.R
```

Produce a version of the dataset with missing immune trait data imputed:  

<!--
Set `--compareCV=TRUE` to run cross-validation (CV) on the following imputation models:
1. Bayesian imputation based on immune trait variance-covariance matrix $\Sigma$ only 
2. Bayesian imputation based on $\Sigma$ and fixed effect covariates (sex, infection)
3. Bayesian imputation based on $\Sigma$, fixed effect covariates and random polygenic effect 
4. Random Forest (RF) imputation with genotypes 
5. RF imputation without genotypes 

Set `--runCV=TRUE` to run cross-validation on only the chosen imputation method, method (2).
-->

```
Rscript imputation.R
```

This imputed dataset will be stored in `derived_data/cross_data_flow_imp.csv`.

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

*Note that variance component analysis is computationally intensive, with each trait potentially taking 1-3 hours, depending on if bootstrapping is necessary. It will be much faster to perform on an HPC cluster than locally so that analysis for all traits can be run completely in parallel. UNC's HPC cluster does not allow Docker containers - if this is the case for your cluster, you can either 1) FTP the repo to the cluster, manually ensure the necessary packages are installed, run the script, and FTP back the results file; or 2) create an Apptainer (fka Singularity) container from the Docker container (we do not provide any guidance on performing this analysis with Apptainer).*

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
Perform variable selection and generate Figure 2. This is set up to run 5000 iterations of a Gibbs sampler, in addition to 5 repeats of 10-fold cross-validation, so is computationally intensive. 

```
Rscript var_selection.R 
```

## QTL mapping  
Perform QTL mapping on a high-performance computing cluster using SLURM (requests enough CPUs to run in parallel by trait): 

```
sbatch map_qtl.sh
```

or locally (runs in parallel using available CPUs):

```
Rscript map_qtl.R
```

Summarize QTL mapping results and create Figure 3:

```
Rscript summarizeQTL.R
```

## Mediation analysis  
Run mediation analysis and create Figure 4:

```
Rscript mediation.R
```

