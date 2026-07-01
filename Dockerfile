FROM rocker/verse:4.4.0

# Install R packages
RUN R -e "install.packages('tidyverse')"
RUN R -e "install.packages('plyr')"
RUN R -e "install.packages('lme4')"
RUN R -e "install.packages('MESS')"
RUN R -e "install.packages('ranger')"
RUN R -e "install.packages('extRemes')"
RUN R -e "install.packages('car')"
RUN R -e "install.packages('pzfx')"
RUN R -e "install.packages('qtl')"
RUN R -e "install.packages('AGHmatrix')"
RUN R -e "devtools::install_github('variani/lme4qtl')"
RUN R -e "install.packages('mvtnorm')"
RUN R -e "install.packages('r2glmm')"
RUN R -e "install.packages('ordinal')"
RUN R -e "install.packages('caret')"
RUN R -e "install.packages('glmnet')"
RUN R -e "install.packages('randomForest')"
RUN R -e "install.packages('MLmetrics')"

# parallel computing 
RUN R -e "install.packages('snow')"
RUN R -e "install.packages('parallel')"
RUN R -e "install.packages('doParallel')"

# plotting 
RUN R -e "install.packages('shades')"
#RUN R -e "install.packages('scales')" # should be automatically included in tidyverse
RUN R -e "install.packages('ggh4x')"
RUN R -e "install.packages('factoextra')"
RUN R -e "install.packages('cowplot')"
RUN R -e "install.packages('RColorBrewer')"
RUN R -e "install.packages('heatmaply')"
#RUN R -e "install.packages('BiocManager')"
RUN R -e "BiocManager::install('ComplexHeatmap')"
RUN R -e "install.packages('ComplexHeatmap')"
RUN R -e "install.packages('patchwork')"
RUN R -e "install.packages('ggridges')"
RUN R -e "install.packages('gridExtra')"
RUN R -e "install.packages('ggplotify')"

# Bayesian sampling 
RUN R -e "install.packages('matrixStats')"
RUN R -e "install.packages('coda')"
RUN R -e "install.packages('posterior')"

RUN R -e "install.packages('readxl')"
RUN R -e "install.packages('writexl')"
RUN R -e "install.packages('UpSetR')"
# RUN R -e "install.packages('ggplotify')"


# RUN R -e "BiocManager::install('rtracklayer')"