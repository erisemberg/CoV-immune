# CoV-immune

Environment prep
-----------------------

# update R version
This project uses a Docker container to produce an environment similar to that used in the original analysis (e.g. R v4.2.1 and R package versions available on August 1, 2023). In order to run this container you will need [Docker](https://docs.docker.com/get-docker/) installed. 

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

At the moment, this data is integrated with the other immune and phenotype data in Excel, resulting in `source_data/F2_phenos_lymph3.xlsx`.

Prepare Rqtl file:

```
Rscript rqtl_file_prep.R
```

Prepare dataset for non-Rqtl analyses (with data transformed, centered/scaled, etc.)

```
Rscript data_processing.R
```
