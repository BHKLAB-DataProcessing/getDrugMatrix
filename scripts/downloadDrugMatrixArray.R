require(downloader)
library(curl)

options(stringsAsFactors = FALSE)

args <- commandArgs(trailingOnly = TRUE)
download_dir <- args[1]
# download_dir <- "/Users/minoru/Code/bhklab/DataProcessing/ToxicoSet/getDrugMatrix/download"

dir1 <- "ftp://ftp.ebi.ac.uk/pub/databases/microarray/data/dixa/DrugMatrix/archive/hepatocyte/"

tt <- read.csv("https://orcestradata.blob.core.windows.net/toxico/DrugMatrix_array_samples.txt") # CEL files to download (only require 939 for our TSet)

samples <- tt$x

tt <- split(samples, ceiling(seq_along(samples) / 100)) # split into chunks to avoid time-out

for (i in 1:length(tt)) {
  print(i)
  samples <- tt[[i]]

  lapply(samples, function(filename) {
    curl_download(paste(dir1, filename, sep = ""), destfile = file.path(download_dir, filename))
  })
}
