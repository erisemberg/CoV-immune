# This script processes genotype and phenotype data for R/qtl analysis. Performs
# marker filtering and recoding (ATCG -> AA/AB/BB) and integrates phenotype data 
# into a csv file in the appropriate format for analysis with R/qtl. 
#
# Input:
#   Genotype file - first three columns must be marker name, chr, pos; followed 
#       any number of marker metadata columns; followed by parent/F1 genotypes; 
#       followed by F2 genotypes 
#       ### F2 genotypes must be the last columns 
#       ### Change Position column name to "Position_b38"
#       ### column names (mouse IDs) MUST match Geno_ID in phenotype file - verify
#           that UNC IDs are padded in mouse IDs from miniMUGA. The ONLY checking
#           this script does is to ensure all IDs are uppercase. 
#       ### If there are duplicates, read.csv will rename the second instance 
#           so the column names are unique, and this code will use the first instance
#           only. 
#   Phenotype file - each row is an F2 mouse, each column is a phenotype. Must 
#       have a "Geno_ID" column in the format "Cr_RB05_<sex>_<padded UNC ID>"
#
# Output: csv file in R/qtl format
#
# Note: this script does not yet include any filtering/handling of Y/MT/PAR markers

library(tidyverse)
library(readxl) 
source("code-dependencies/cov_qtl_functions.R")

ensure_directory("derived_data")
ensure_directory("derived_data/data_processing")
log <- make_logger("derived_data/data_processing/file_processing_notes.md")

#---------------------------------Parameters-----------------------------------#
geno_file <- "source_data/Cr_RB05_miniMUGA-013024_paddedIDs.csv"
pheno_file <- "source_data/F2_phenos_lymph3.xlsx"
ensure_directory("derived_data")
out_file <- "derived_data/RqtlCC006xCC044_ctrlAndSARS.csv"
Y_MT_out_file <- "derived_data/data_processing/CC006xCC044_YandMTandPARgenos.csv"

num_F2s <- 1203 # number of genotyped F2 mice
crosstype <- "f2" # cross type (f2 or bc)
# Phenotypes to include - add "G_Eosinophils_pctCD69", "G_Eosinophils_pctMHCII" when phenotyped 
pheno.names <- read_delim(file = "source_data/phenotype_names.txt", delim = "\n", col_names = FALSE)[[1]]

# infection type: SARSCoV, SARS2CoV, PBS, HKU3CoV, or NULL
# Can also be a vector of >1 infection types
infection = c("PBS", "SARSCoV", "SARS2CoV")  

A.ref = "Cr_RB05_CC006.UncB6_F_7082" # parent mouse representing A genotype 
B.ref = "Cr_RB05_CC044.UncB6._F_2397" # parent mouse representing B genotype 

#------------------------------------------------------------------------------#
#-------------------------------Genotype data----------------------------------#
#------------------------------------------------------------------------------#
geno <- read.csv(geno_file)
log(paste("Genotype data loaded from ", geno_file, ".", sep = ""))

F2_start = ncol(geno)-num_F2s+1 # column index of first F2 mouse column
F2_end = ncol(geno) # column index of last F2 mouse column 
log(paste("Columns ", F2_start, " to ", F2_end, " are F2 mice."))

log(paste("Starting with ", nrow(geno), " markers.", sep=""))
geno <- geno[geno$Chromosome != 0, ] # Remove chr0 rows
log(paste("After removing chr0 markers: ", nrow(geno), sep=""))

geno$Position_b38 <- geno$Position_b38/1e6 # convert bp position to cM

#----------------------------Metric calculation--------------------------------#
# Get ref, alt, het and N calls in the F2 mice at each marker
geno$num_ref <- rowSums(geno[,F2_start:F2_end] == geno$reference)
geno$num_alt <- rowSums(geno[,F2_start:F2_end] == geno$alternate)
geno$num_N <- rowSums(geno[,F2_start:F2_end] == "N")
geno$pct_N <- geno$num_N/num_F2s
geno$num_H <- rowSums(geno[,F2_start:F2_end] == "H")

geno$ref_alt <- geno$num_ref/(geno$num_ref + geno$num_alt) # Calculate ref as proportion of ref+alt
geno$het_all <- geno$num_H/(geno$num_ref + geno$num_alt + geno$num_H) # Calculate het relative to good (non-N) calls

#--------------------------------Filtering-------------------------------------#
### Separate out MT and Y markers // what to do with PAR? 
geno.Y.MT <- geno[which((geno$Chromosome == "Y") | geno$Chromosome == "MT"), ]
### Filter Y chromosome markers - all females should have N calls for each Y marker
#geno.Y.MT <- geno[which((geno$Chromosome == "Y") & (geno$pct_N == 0))]
geno <- geno[which((geno$Chromosome != "Y") & (geno$Chromosome != "MT")),]
log(paste("After removing Y/MT markers: ", nrow(geno), sep=""))

# Filter out markers with failure rate > 5% N  
geno <- geno[which(geno$pct_N <= 0.05),]
log(paste("After removing markers with > 5 %% rate of N calls: ", nrow(geno), sep=""))

# Filter out bad and uninformative markers
# Anything with almost all ref or alt is uninformative - remove markers where ref/alt is > 95% 
geno <- geno[which((geno$num_alt/num_F2s) <= 0.95),] 
log(paste("After removing markers with > 95 %% alt calls: ", nrow(geno), sep=""))
geno <- geno[which((geno$num_ref/num_F2s) <= 0.95),] 
log(paste("After removing markers with > 95 %% ref calls: ", nrow(geno), sep=""))
# Remove markers with almost all het calls (or all H/N calls)
geno <- geno[which(geno$het_all != 1),] 
log(paste("After removing markers with all het calls: ", nrow(geno), sep=""))

# Plot x = ref/(ref+alt) and y = het/(ref+alt+het)
png("derived_data/data_processing/marker_qc_plot.png")
plot(x=geno$het_all, y=geno$ref_alt, main="Autosome/X markers", 
     ylab="Ref/(Ref + Alt)", xlab="Het/(Het + Ref + Alt)")
dev.off()

if (crosstype == "f2"){
  # For all chromosomes (autosomes and X): looking for ref:alt of ~1:1
  geno <- geno[which((geno$ref_alt >= 0.4) & (geno$ref_alt <= 0.6)),] # Looking for ref:alt of ~1:1
  # For autosomes, looking for het_all to be ~0.5. For X-chr, het_all should be ~0.25
  geno <- geno[which(((geno$Chromosome != "X") & (geno$het_all >= 0.4) & (geno$het_all <= 0.6)) | 
                     ((geno$Chromosome == "X") & (geno$het_all >= 0.15) & (geno$het_all <= 0.35))),]
} else if (crosstype == "bc") {
  geno <- geno[which((geno$het_all >= 0.4) & (geno$het_all <= 0.6)),] # Looking for het:hom ratio of ~1:1
  geno <- geno[which((geno$ref_alt == 0) | (geno$ref_alt == 1)),] # Need hard cutoff to satisfy R/qtl 
}

log(paste("After filtering for informative ", ifelse(crosstype=="bc", "backcross", "F2"), " markers: ", nrow(geno), sep=""))

log("Marker filtering complete.")

#---------------------------------Recoding-------------------------------------#
recoded.genos <- matrix(NA, nrow=dim(geno[,F2_start:F2_end])[1], 
                        ncol=dim(geno[,F2_start:F2_end])[2])
for (i in 1:ncol(recoded.genos)){ # For each F2 column 
  recoded.genos[(geno[,i+F2_start-1] == geno[[A.ref]]), i] <- "AA"
  recoded.genos[(geno[,i+F2_start-1] == geno[[B.ref]]), i] <- "BB"
  recoded.genos[(geno[,i+F2_start-1] == "H"),i] <- "AB"
}
geno[,F2_start:F2_end] <- recoded.genos 

geno[geno=="N"] <- "-" # Replace all "N" with "-"

log("Genotype data recoded from A/T/C/G/H to AA/AB/BB.")

#------------------------------Prep for R/qtl----------------------------------#
# Remove marker stats
rqtl <- geno[,-c((F2_end+1):ncol(geno))]  
# Remove marker metadata (other than name, chr, and pos) and marker data for parent and F1 mice
rqtl <- rqtl[,-c(4:(F2_start-1))]  
rqtl <- t(rqtl) # transpose

# Replace numbered colnames with marker names 
colnames(rqtl) <- rqtl[1,]
rqtl <- rqtl[-1,]

rqtl <- as.data.frame(rqtl)
rqtl <- rownames_to_column(rqtl, var="mouse_ID") %>% as_tibble

#------------------------------------------------------------------------------#
#-------------------------------Phenotype data---------------------------------#
#------------------------------------------------------------------------------#
pheno <- read_xlsx(pheno_file)
log(paste("Phenotype data loaded from ", pheno_file, ".", sep=""))

# Convert pheno$Geno_ID and geno$mouse_ID to uppercase for matching
pheno$Geno_ID <- toupper(pheno$Geno_ID) 
rqtl$mouse_ID <- c(rqtl$mouse_ID[1:2], toupper(rqtl$mouse_ID[3:nrow(rqtl)]))

# If specified, select mice with noted infection type 
if (!is.null(infection)){
  pheno <- pheno[which(pheno$infection %in% infection),]
}

# Before adding phenotype data to genotype data, filter down to mice that we have 
# phenotype data for (keep first two rows, chromosome and position for marker cols)
rqtl <- rqtl %>% filter(mouse_ID == "Chromosome" | 
                        mouse_ID == "Position_b38" | 
                        mouse_ID %in% pheno$Geno_ID)

# Empty phenotype column (with first two rows empty, as required for R/qtl)
new.col = c(rep("",2), rep(NA, nrow(rqtl)-2))

# Add phenotype data to genotype data
col.pos = 2
for (k in 1:length(pheno.names)){
  pheno.name = pheno.names[k]
  rqtl <- add_pheno(pheno.name, new.col, col.pos)
  col.pos = col.pos+1
}

log("Genotype and phenotype data integrated.")
#------------------------------------------------------------------------------#
#--------------------------------Output data-----------------------------------#
#------------------------------------------------------------------------------#
rqtl[1:2,"mouse_ID"] <- "" # Remove chr/pos row labels

write.csv(rqtl, out_file, row.names=F)
write.csv(geno.Y.MT, Y_MT_out_file)

log(paste("Output to csv file: ", out_file, ".", sep=""))
