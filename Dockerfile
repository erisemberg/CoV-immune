FROM rocker/verse:4.4.1

# Install R packages
RUN R -e "install.packages('tidyverse')"
RUN R -e "install.packages('plyr')"
RUN R -e "install.packages('lme4')"
RUN R -e "install.packages('MESS')"
RUN R -e "install.packages('shades')"
RUN R -e "install.packages('ranger')"
RUN R -e "install.packages('scales')"
RUN R -e "install.packages('extRemes')"
RUN R -e "install.packages('RColorBrewer')"
RUN R -e "install.packages('parallel')"
RUN R -e "install.packages('snow')"
RUN R -e "install.packages('car')"
RUN R -e "install.packages('pzfx')"
RUN R -e "install.packages('ggh4x')"
RUN R -e "install.packages('qtl')"
RUN R -e "install.packages('AGHmatrix')"
RUN R -e "devtools::install_github('variani/lme4qtl')"
RUN R -e "install.packages('doParallel')"
RUN R -e "install.packages('mvtnorm')"

# packages for lymphoid-integration.R, rqtl_file_prep.R, and data_processing.R included
# ADD qtl_functions.R PACKAGES	

# RUN R -e "BiocManager::install('rtracklayer')"