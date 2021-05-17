library(ToxicoGx)
library(Biobase)
library(affy)
library(affyio)
library(BiocManager)
library(dplyr)
library(biomaRt)
library(ggplot2)
library(data.table)
library(SummarizedExperiment)
library(biocompute)

#reading the metadata downloaded from diXa
s_Hepatocyte <- read.delim("/pfs/getDrugMatrix/s_Hepatocyte.csv", stringsAsFactors = F, sep = "\t")
s_Hepatocyte$Factor.Value.Compound. <- gsub("17\\?-ethynylestradiol", "17a-ethynylestradiol", s_Hepatocyte$Factor.Value.Compound.)

s_Hepatocyte_fil <- s_Hepatocyte[, colSums(is.na(s_Hepatocyte)) != nrow(s_Hepatocyte)]
colnames(s_Hepatocyte_fil) <- c("source_name", "samplename", "STRAIN_TYPE", "species", "org_id_abbr", "Term.Accession.Number", "sex_type", "cellid",
                                "Term.Source.REF.4", "Term.Accession.Number.4", "test_type", "Biological.Replicate", "Technical.Replicate", "dataset_drugid", 
                                "Term.Source.REF.5", "CHEBI.ID", "inchikey", "Comment.chEMBL.ID", "Control", "Sample.Match", "concentration", "Dose_Unit", 
                                "Term.Source.REF.6", "Term.Accession.Number.6", "Dose.Duration.inDays", "Dose.Duration.Unit", "Term.Source.REF.7", "Accession.Number.7", 
                                "Vehicle", "Treatment.Group.", "Protocol.REF", "Sample.Name")

s_Hepatocyte_fil$species <- "R.norvegicus"
s_Hepatocyte_fil$individual_id <- s_Hepatocyte_fil$Biological.Replicate
s_Hepatocyte_fil$batchid <- as.integer(0)
s_Hepatocyte_fil$chiptype <- "Rat_Genome_230_2.0"
s_Hepatocyte_fil$organ_id <- "Liver"
s_Hepatocyte_fil$Dose_Unit <- "μM"
s_Hepatocyte_fil$celfilename <- paste(s_Hepatocyte_fil$samplename, ".CEL", sep = "")
s_Hepatocyte_fil$duration <- floor(s_Hepatocyte_fil$Dose.Duration.inDays * 24)
s_Hepatocyte_fil$dose_level <- NA


s_Hepatocyte_fil$xptype <- NA
#add perturbation labels
for(xp in 1:nrow(s_Hepatocyte_fil)){
  
  if(s_Hepatocyte_fil$Treatment.Group.[xp] == "treated"){
    s_Hepatocyte_fil$xptype[xp] <- "perturbation" 
  }
  else if (s_Hepatocyte_fil$Treatment.Group.[xp] == "untreated"){
    s_Hepatocyte_fil$xptype[xp] <- "control" 
  }
}
#correcting wrong labeling of loratadine as desloratadine in the raw file
for(lr in 1:nrow(s_Hepatocyte_fil)){
  if(s_Hepatocyte_fil$Comment.chEMBL.ID[lr] == "CHEMBL998"){
    s_Hepatocyte_fil$dataset_drugid[lr] <- "loratadine"
    s_Hepatocyte_fil$CHEBI.ID[lr] <- as.integer(6538)
  }
}
#add lab curated drug names
curationDrug <- readRDS("/pfs/getDrugMatrix/curationDrug.rds")
s_Hepatocyte_fil$drugid <- NA

for(dt in 1:nrow(s_Hepatocyte_fil)){
  for(cd in 1:nrow(curationDrug)){
    if(s_Hepatocyte_fil$dataset_drugid[dt] == curationDrug$dataset_drugid[cd]){
      s_Hepatocyte_fil$drugid[dt] <- curationDrug$unique.drugid[cd]
    }
  }
}

#label control rows as DMSO
for(xp in 1:nrow(s_Hepatocyte_fil)){
    if(s_Hepatocyte_fil$xptype[xp] == "control"){
      s_Hepatocyte_fil$drugid[xp] <- "DMSO"
      s_Hepatocyte_fil$dataset_drugid[xp] <- "DMSO"
    
  }
}

#dose unit of TGF beta is in ng/ml
for(du in 1:nrow(s_Hepatocyte_fil)){
  if(s_Hepatocyte_fil$dataset_drugid[du] == "tgf beta-1, human recombinant"){
    s_Hepatocyte_fil$concentration[du] <- gsub(" ng/ml","",s_Hepatocyte_fil$concentration[du])
    s_Hepatocyte_fil$Dose_Unit[du] <- "ng/ml"
  }
}

#manually mapped dose levels
doselevel_mapping <- readRDS("/pfs/getDrugMatrix/doselevel_mapping.rds")
s_Hepatocyte_fil$concentration <- as.numeric(s_Hepatocyte_fil$concentration)

for(ds in 1:nrow(s_Hepatocyte_fil)){
  for(dl in 1:nrow(doselevel_mapping)){
    if((s_Hepatocyte_fil$drugid[ds] == doselevel_mapping$drugid[dl]) && (s_Hepatocyte_fil$concentration[ds] == doselevel_mapping$concentration[dl])){
      s_Hepatocyte_fil$dose_level[ds] <- doselevel_mapping$dose_level[dl]
    }
  }
}

rownames(s_Hepatocyte_fil) <- s_Hepatocyte_fil$samplename

dropcols <- c("Term.Accession.Number","Term.Source.REF.4","Term.Accession.Number.4","Term.Source.REF.5","Term.Source.REF.6","Term.Accession.Number.6",
              "Term.Source.REF.7","Protocol.REF","Sample.Name" )

phenodata_DM <- s_Hepatocyte_fil[,!colnames(s_Hepatocyte_fil) %in% dropcols]

#reorder columns
phenodata_DM <- phenodata_DM[,c("source_name","samplename","STRAIN_TYPE","species","org_id_abbr","sex_type","cellid","test_type",          
                  "Biological.Replicate","individual_id", "Technical.Replicate","drugid","dataset_drugid","CHEBI.ID","inchikey","Comment.chEMBL.ID","Control","Sample.Match",        
                  "concentration","dose_level", "Dose_Unit","Dose.Duration.inDays","Dose.Duration.Unit","duration","Accession.Number.7","Vehicle","Treatment.Group.",       
                  "batchid","chiptype","organ_id","celfilename","xptype")]

create_exprsdata_DM <- function(species=c("Rat"), phenodata_DM, verbose = TRUE) {
  if (verbose) {message("Creating eset object")}
  if(species == "Rat"){
  }
  
  eset <- readRDS("/pfs/processDrugMatrixArray/eset_DM.rds")
  
  storageMode(eset)<-"environment"
  eset <-subset(eset, substr(rownames(eset@assayData$exprs), 0, 4) != "AFFX")
  storageMode(eset)<-"lockedEnvironment"
  annotation(eset)<-"rna"
  if (verbose) {message("eset object created!")}
  return(eset)
}
eset <- create_exprsdata_DM("Rat", phenodata_DM)
phenodata_DM <- phenodata_DM[match(rownames(eset@phenoData@data), phenodata_DM$celfilename),] #added since eset had different ordering of CEL names as phenoData, which raises an error
colnames(eset) <- gsub(".CEL", "", colnames(eset))

create_featuredata_DM <- function(species=c("Rat"), eset, verbose = TRUE){
  if (verbose) {message("Creating featuredata_DM object...")}
  if (species == "Rat"){
    ensembl <- useMart("ensembl")
    datasets <- listDatasets(ensembl)
    ensembl = useMart("ensembl", dataset = "rnorvegicus_gene_ensembl", host="uswest.ensembl.org")
    storageMode(eset) <- "environment"
    affxrows <- rownames(eset@assayData$exprs)
    rownames(eset@assayData$exprs) <- substr(rownames(eset@assayData$exprs), 1, nchar(affxrows)-3)
    CELgenes <- affxrows
    CELgenes1 <- gsub(".at", " ", CELgenes)
    results <-getBM(attributes=c("external_gene_name","ensembl_gene_id","gene_biotype","entrezgene_id"
                                 ,"external_transcript_name","ensembl_transcript_id"), filters = "ensembl_gene_id"
                    ,values=CELgenes1, mart=ensembl,checkFilters = TRUE)
    uniqueB <- results[!duplicated(results$ensembl_gene_id),]
    CELnotB <- unique(CELgenes1) [!unique(CELgenes1) %in% uniqueB$ensembl_gene_id]
    names(uniqueB) <- c("gene_name", "gene_id", "gene_biotype", "EntrezGene.ID", "transcript_name", "transcript_id")
    finalFeature <- uniqueB
    
    finalFeature$BEST <- NA
    names(finalFeature) <- c("Symbol", "gene_id", "gene_biotype", "EntrezGene.ID", "transcript_name", "transcript_id", "BEST")
    rownames(finalFeature) <- finalFeature$gene_id
    finalFeature$gene_id
    geneid1 <- finalFeature$gene_id
    
    for (i in 1:length(geneid1)) {
      geneid1[i] = paste(geneid1[i], "at", sep="_")
    }
    geneid1
    finalFeature$gene_id <- geneid1
    finalFeature$gene_id
    finalFeature[,1]
    rownames(finalFeature) = finalFeature$gene_id
    
    if(verbose) {message("featuredata_DM object created!")}
    return(finalFeature)
    
  }
}
featuredata_DM <- create_featuredata_DM("Rat", eset)
create_Expressionset <- function(species=c("Rat"), eset, verbose = TRUE){
  if (verbose) {message("Creating expressionset...")}
  if (species == "Rat"){
    stopifnot(all(rownames(pData(eset)) == rownames(phenodata_DM)))
   
    pData(eset) <- phenodata_DM
    fData(eset) <- featuredata_DM
    #sorting rownames to maintain feature data mapping that m=is otherwise shuffled after converting to SE
    fData(eset) <- fData(eset)[sort(rownames(fData(eset))),]
    stopifnot(all(rownames(fData(eset)) == rownames(exprs(eset))))
    stopifnot(all(rownames(pData(eset)) == colnames(exprs(eset))))
    
    storageMode(eset) <- "lockedEnvironment"
    return(eset)
  }
}

ExpressionSet <- create_Expressionset("Rat", eset)

#the conversion function might be incorporated to the package later. This step needs to be updated then
print("Creating summarized experiment object...")
#new_SE_DM <-  as(ExpressionSet, value="SummarizedExperiment")
new_SE_DM <- SummarizedExperiment::makeSummarizedExperimentFromExpressionSet(ExpressionSet)
stopifnot(all(rownames(colData(new_SE_DM)) == rownames(pData(ExpressionSet))))
stopifnot(all(rownames(rowData(new_SE_DM)) == rownames(fData(ExpressionSet))))

print("Done!")
######################toxicoset constructor function ######################

create_curationCell <- function(phenodata_DM, verbose = TRUE){
  curationCell <- unique(subset(phenodata_DM, select=c(cellid)))
  curationCell$dataset_cellid <- curationCell$cellid
  names(curationCell) <- c("unique.cellid", "dataset_cellid")
  rownames(curationCell) <- curationCell$unique.cellid
  
  return(curationCell)
}
curationCell <- create_curationCell(phenodata_DM)


create_curationTissue <- function(phenodata_DM, verbose = TRUE){
  curationTissue <- unique(subset(phenodata_DM, select=c(organ_id)))
  curationTissue$dataset_tissueid <- "Liver"
  names(curationTissue)[1] <- "unique.tissueid"
  rownames(curationTissue) <- "Hepatocyte"
  
  return(curationTissue)
}
curationTissue <- create_curationTissue(phenodata_DM)


#reading the metadata downloaded from diXa

create_drug <- function(phenodata_DM){
  
  sub_hepa <- subset(phenodata_DM, subset = !duplicated(phenodata_DM$dataset_drugid),select = c("dataset_drugid","CHEBI.ID","inchikey","Comment.chEMBL.ID", "drugid" ), drop = F)
  colnames(sub_hepa) <- c("dataset_drugid", "CHEBI.ID", "inchikey", "CHEMBL.ID","drugid")
  rownames(sub_hepa) <- sub_hepa$drugid
  #reorder based on rows of curationDrug 
  sub_hepa <- sub_hepa[rownames(curationDrug),]
  stopifnot(all(rownames(sub_hepa) == rownames(curationDrug)))
  return(sub_hepa)
}

drug <- create_drug(phenodata_DM)

create_cell <- function(phenodata_DM, verbose = TRUE){
  cell <- unique(subset(phenodata_DM,select=c(cellid, organ_id, species, test_type)))
  names(cell)<-c("cellid","tissueid", "species","testType")
  cell$tissueid<-"Liver"
  rownames(cell) <- cell$cellid
  
  return(cell)
}
cell <- create_cell(phenodata_DM)

drugMatrix <- ToxicoSet("drugMatrix_rat",
                        molecularProfiles=list("rna"= new_SE_DM),
                        cell=cell,
                        drug=drug,
                        sensitivityInfo=NULL,
                        sensitivityRaw=NULL,
                        sensitivityProfiles=NULL,
                        curationDrug=curationDrug,
                        curationCell=curationCell,
                        curationTissue=curationTissue,
                        datasetType = c("perturbation"),
                        verify = TRUE)

saveRDS(drugMatrix, "/pfs/out/drugMatrix.rds")

####CREATE BIOCOMPUTE OBJECT####

library(biocompute)
###########################
#####Provenance Domain#####
###########################

#Created and modified dates
#Sys.setenv(TZ = "EST")
created <- as.POSIXct(Sys.time(), format = "%Y-%m-%dT%H:%M:%S", tz = "EST")
modified <- as.POSIXct(Sys.time(), format = "%Y-%m-%dT%H:%M:%S", tz = "EST")

#Contributions
contributors <- data.frame(
  "name" = c("Anthony Mammoliti", "Petr Smirnov", "Benjamin Haibe-Kains"),
  "affiliation" = c(rep("University Health Network", 3)),
  "email" = c("anthony.mammoliti@uhnresearch.ca", "petr.smirnov@utoronto.ca", "Benjamin.Haibe-Kains@uhnresearch.ca"),
  "contribution" = c("createdBy","createdBy","authoredBy"),
  "orcid" = c(NA,NA,"https://orcid.org/0000-0002-7684-0079"),
  stringsAsFactors = FALSE
)

#License
license <- "https://opensource.org/licenses/Apache-2.0"

#Name of biocompute object
name <- "DrugMatrix Rat"

#Version of biocompute object
bio_version <- "1.0.0"

#Embargo (none)
embargo <- c()

#Derived from and obsolete after (none)
derived_from <- c()
obsolete_after <- c()

#reviewers (none)
review <- c()

#compile domain
provenance <- compose_provenance_v1.3.0(
  name, bio_version, review, derived_from, obsolete_after,
  embargo, created, modified, contributors, license
)
provenance %>% convert_json()


############################
#####Description Domain#####
############################
times_rnaseq <- as.POSIXct("2020-01-22T2:40:38", format = "%Y-%m-%dT%H:%M:%S", tz = "EST")
#Keywords and platform info
keywords <- c("Biomedical", "Toxicogenomics")
platform <- c("Pachyderm", "ORCESTRA (orcestra.ca)", "Linux/Ubuntu")

#Metadata for each pipeline step
pipeline_meta <- data.frame(
  "step_number" = c("1","2","3"),
  "name" = c("Microarray processing",
             "Curated Sample and treatment identifier compilation",
             "Build data object"),
  "description" = c("Normalization of microarray data", 
                    "Download of appropriate sample and treatment identifiers from GitHub (curations performed by BHK Lab - http://bhklab.ca)",
                    "Building of ORCESTRA data object"),
  "version" = c(1.0,1.0,1.0),
  stringsAsFactors = FALSE
)

#Inputs for each pipeline step
pipeline_input <- data.frame(
  "step_number" = c("1","2","2","3"),
  "filename" = c("Raw microarray data",
                 "Sample annotation data",
                 "Treatment annotations",
                 "Script for data object generation"),
  "uri" = c(
    "ftp://ftp.ebi.ac.uk/pub/databases/microarray/data/dixa/DrugMatrix/archive/hepatocyte/",
    "/pfs/getDrugMatrix/s_Hepatocyte.csv",
    "/pfs/getDrugMatrix/curationDrug.rds",
    "https://github.com/BHKLAB-Pachyderm/getDrugMatrix/blob/main/getDM.R"
  ),
  "access_time" = c(times_rnaseq,created,created,created),
  stringsAsFactors = FALSE
)


#Outputs for each pipeline step
pipeline_output <- data.frame(
  "step_number" = c("1","2","2","3"),
  "filename" = c("Processed microarray data", 
                 "Downloaded sample annotations",
                 "Downloaded treatment annotations",
                 "Data object"),
  "uri" = c(
    "/pfs/processDrugMatrixArray/eset_DM.rds",
    "/pfs/getDrugMatrix/s_Hepatocyte.csv",
    "/pfs/getDrugMatrix/curationDrug.rds",
    "/pfs/out/drugMatrix.rds"
  ),
  "access_time" = c(times_rnaseq,created,created,created),
  stringsAsFactors = FALSE
)

#xref (none)
xref <- c()

#pipeline prereq (none)
pipeline_prerequisite <- c()

#compile domain
description <- compose_description_v1.3.0(
  keywords, xref, platform,
  pipeline_meta, pipeline_prerequisite, pipeline_input, pipeline_output
)
description %>% convert_json()


############################
######Execution Domain######
############################

script <- c()
script_driver <- c()

#software/tools and its versions used for data object creation
software_prerequisites <- data.frame(
  "name" = c("Pachyderm", "Docker Image"),
  "version" = c("1.9.3", "v3"),
  "uri" = c(
    "https://www.pachyderm.com", "https://hub.docker.com/r/bhklab/toxicogx"
  ),
  stringsAsFactors = FALSE
)

software_prerequisites[,"access_time"] <- rep(NA, length(software_prerequisites$name))
software_prerequisites[,"sha1_chksum"] <- rep(NA, length(software_prerequisites$name))

external_data_endpoints <- c()
environment_variables <- c()

execution <- compose_execution_v1.3.0(
  script, script_driver, software_prerequisites, external_data_endpoints, environment_variables
)
execution %>% convert_json()


############################
######Extension Domain######
############################

#repo of scripts/data used
scm_repository <- data.frame("extension_schema"= c("https://github.com/BHKLAB-Pachyderm"))
scm_type <- "git"
scm_commit <- c()
scm_path <- c()
scm_preview <- c()

scm <- compose_scm(scm_repository, scm_type, scm_commit, scm_path, scm_preview)
scm %>% convert_json()

extension <- compose_extension_v1.3.0(scm)
extension %>% convert_json()

############################
######Parametric Domain#####
############################

df_parametric <- data.frame(
  "param" = c(),
  "value" = c(),
  "step" = c(),
  stringsAsFactors = FALSE
)

parametric <- compose_parametric_v1.3.0(df_parametric)
parametric %>% convert_json()



############################
######Usability Domain######
############################

#usability of our data objects
text <- c(
  "Pipeline for creating DrugMatrix Rat data object through ORCESTRA (orcestra.ca), a platform for the reproducible and transparent processing, sharing, and analysis of biomedical data."
)

usability <- compose_usability_v1.3.0(text)
usability %>% convert_json()


######################
######I/O Domain######
######################

input_subdomain <- data.frame(
  "filename" = c("Raw microarray data",
                 "Sample annotation data",
                 "Treatment annotations",
                 "Script for data object generation"),
  "uri" = c(
    "ftp://ftp.ebi.ac.uk/pub/databases/microarray/data/dixa/DrugMatrix/archive/hepatocyte/",
    "/pfs/getDrugMatrix/s_Hepatocyte.csv",
    "/pfs/getDrugMatrix/curationDrug.rds",
    "https://github.com/BHKLAB-Pachyderm/getDrugMatrix/blob/main/getDM.R"
  ),
  "access_time" = c(times_rnaseq,created,created,created),
  stringsAsFactors = FALSE
)

output_subdomain <- data.frame(
  "mediatype" = c("RDS", "csv", "RDS", "RDS"),
  "uri" = c(
    "/pfs/processDrugMatrixArray/eset_DM.rds",
    "/pfs/getDrugMatrix/s_Hepatocyte.csv",
    "/pfs/getDrugMatrix/curationDrug.rds",
    "/pfs/out/drugMatrix.rds"
  ),
  "access_time" = c(created,created,created,created),
  stringsAsFactors = FALSE
)

io <- compose_io_v1.3.0(input_subdomain, output_subdomain)
io %>% convert_json()


########################
######Error Domain######
########################

empirical <- c()
algorithmic <- c()

error <- compose_error(empirical, algorithmic)
error %>% convert_json()


####Retrieve Top Level Fields####
tlf <- compose_tlf_v1.3.0(
  provenance, usability, extension, description,
  execution, parametric, io, error
)
tlf %>% convert_json()


####Complete BCO####

bco <- biocompute::compose_v1.3.0(
  tlf, provenance, usability, extension, description,
  execution, parametric, io, error
)
bco %>% convert_json() %>% export_json("/pfs/out/DrugMatrix_Rat_BCO.json") %>% validate_checksum()

