FROM rocker/verse:4.4.1

# Install R packages
RUN R -e "install.packages("qtl")"
RUN R -e "install.packages("car")"
RUN R -e "install.packages("lme4")"
RUN R -e "install.packages("tidyverse")"
RUN R -e "install.packages("stringi")"

# packages for lymphoid-integration.R, rqtl_file_prep.R, and data_processing.R included
# ADD qtl_functions.R PACKAGES	
RUN R -e "install.packages("qtl")"
RUN R -e "install.packages("qtl")"
RUN R -e "install.packages("qtl")"
RUN R -e "install.packages("qtl")"


# RUN R -e "devtools::install_github('variani/lme4qtl')"
# RUN R -e "BiocManager::install('rtracklayer')"