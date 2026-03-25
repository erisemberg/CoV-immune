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

Prepare dataset for non-Rqtl analyses (with data transformed, centered/scaled, etc.)

```
Rscript data_processing.R
```
