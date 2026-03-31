# CoV-immune

NOTE 3/20/26: this repo is in progress

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

Run the lymphoid integration protocol, which integrates the three datasets in `source_data/lymphoid-combined.xlsx` into one dataset, stored in `derived_data/lymphoid-harmonization-3.csv`.

```
Rscript lymphoid_integration.R
```

At the moment, this data is integrated with the other immune and phenotype data manually in Excel, resulting in `source_data/F2_phenos_lymph3.xlsx`.

Prepare Rqtl file from genotype data (`source_data/Cr_RB05_miniMUGA-013024_paddedIDs.csv`) and phenotype data (`source_data/F2_phenos_lymph3.xlsx`). This creates the main file for QTL analyses, `derived_data/RqtlCC006xCC044_ctrlAndSARS.csv`. This file and the `Rqtl` data structures are often used even when the `Rqtl` package is not being used for mapping. 

```
Rscript rqtl_file_prep.R
```

Prepare dataset with transformed and centered/scaled phenotypes (`derived_data/cross_data.csv`), as well as an `Rqtl` cross object with imputed genotypes. 

```
Rscript data_processing.R
```

Produce a version of the dataset with missing immune trait data imputed, which will be stored in `derived_data/cross_data_flow_imp.csv`.  

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

# Data analysis 

Perform variance component analysis on a high-performance computing cluster using SLURM (requests 106 CPUs to run in parallel by trait):

```
sbatch vc.sh
```

or locally (runs in parallel using available CPUs):

```
Rscript var_comp.R
```

Perform miscellaneous analyses/visualizations on trait data. 
* Levene's test for variance heterogeneity 
* Classification of infection group by immune trait data with `glmnet` and random forest 
* Create Figure 1 and Supp. Figure 1 

```
Rscript trait_eda.R
```

## Variable selection 
Perform variable selection and generate Figure 2. This is set up to run 5000 iterations of a Gibbs sampler, in addition to 5 repeats of 10-fold cross-validation, so is computationally intensive. 

```
Rscript var_selection.R 
```

## QTL mapping  


## Mediation analysis  


