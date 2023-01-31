library(gdata)
library(affy)
library(Biobase)
library(xml2)
library(abind)
library("rat2302rnensgcdf")

args <- commandArgs(trailingOnly = TRUE)
download_dir <- args[1]
process_dir <- args[2]
# download_dir <- "/Users/minoru/Code/bhklab/DataProcessing/ToxicoSet/getDrugMatrix/download"
# process_dir <- "/Users/minoru/Code/bhklab/DataProcessing/ToxicoSet/getDrugMatrix/processed"

cdf <- "rat2302rnensgcdf"

celfn <- list.files(download_dir, full.names = T, "\\.CEL$")

eset <- just.rma(filenames = celfn, verbose = TRUE, cdfname = cdf)

saveRDS(eset, file.path(process_dir, "eset_DM.rds"))
