####################### Powrie MuSiC Tool #######################

# Metadata
## Author: Amit Halkhoree & Nicholas Ilott- Powrie Lab, Kenendy Institute of Rheumatology, University of Oxford
## Correspondance Email: amit.halkhoree@kennedy.ox.ac.uk
## Date: 04-Jan-2023
## Licence: MIT
## R Version: 4.2.1
## Version 0.06

# Description: 
## This tool is designed to provide a command line interface for cell type proportion estimations using the MuSiC algorithm 
## implemented by Xuran Wang https://doi.org/10.1038/s41467-018-08023-x.
## The tool takes bulk RNA Seq inputs from Kalisto outputs and references them versus a single cell reference set, to generate cell proportion estimations.
## The single cell reference set is known as the Powrie Album, and can be found on the shared drive of the BMRC
## For brevity, all functionality will be contained in this script.

# Usage Instructions:
## To run the script: Rscript powrieMuSiC.r -s SINGLE_CELL_EXPERIMENT_OBJECT -b BULK_DATA_FOLDER -t SPECIES_TYPE


################################################################

# Loading and installing the required packages


## Checking the existance of packages, and installing if not found
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.16")

if (!require("SingleCellExperiment", quietly = TRUE)) {
  biocmanager::install("SingleCellExperiment")
}

if (!require("tximport", quietly = TRUE)) {
  biocmanager::install("tximport")
}

if (!require("biomaRt", quietly = TRUE)) {
  biocmanager::install("biomaRt")
}

if (!require("Biobase", quietly = TRUE)) {
  biocmanager::install("Biobase")
}
if (!require("devtools", quietly = TRUE)) {
  install.packages("devtools")
}
devtools::install_github("xuranw/MuSiC")

if (!require("argparse", quietly = TRUE)) {
  install.packages("argparse")
}

if (!require("here", quietly = TRUE)) {
  install.packages("here")
}

## Loading in the required packages
library(SingleCellExperiment)
library(tximport)
library(biomaRt)
library(Biobase)
library(MuSiC)
library(argparse)
library(here)

####################### Utility Functions #######################

## Species handling

#' Handles species-specific data for mouse or human.
#'
#' @param species A character string indicating the species, either "mouse" or "human".
#' @return A list with two elements:
#'   ensembl_dataset: A character string with the Ensembl dataset name for the species.
#'   t2g_path: A character string with the file path to the transcript-to-gene mapping file for the species.
#' @examples
#' species_handler("mouse")
#' species_handler("human")
species_handler <- function(species){
  if(species == 'mouse'){
    ensembl_dataset <- 'mmusculus_gene_ensembl'
    t2g_path <- t2g_path <- file.path(here::here(), "data/t2g", "mouse_t2g.tsv")
  } else if(species == 'human'){
    ensembl_dataset <- 'hsapiens_gene_ensembl'
    t2g_path <- t2g_path <- file.path(here::here(), "data", "human_t2g.tsv")
  } else {
    stop("Error: Invalid species. Must be mouse or human. Please open a Github issue if you require a different species!")
  }
  return(list(ensembl_dataset, t2g))
}

## Bulk RNA Seq

#' Convert TPMs to counts using Tximport from the output of Kalisto
#'
#' This function reads in a transcript mapping file and a set of input files,
#' converts the TPM values in the input files to counts using Tximport, and
#' writes the resulting counts table to a file.
#'
#' @param transcript_mapping_file Path to the transcript mapping file.
#' @param input_directory Path to the input directory containing the input files.
#' @return A counts table
convert_to_counts <- function(transcript_mapping_file, input_directory) {
  tx2gene <- read.csv(transcript_mapping_file, header=TRUE, stringsAsFactors=FALSE, sep="\t")
  currwd <- getwd()
  setwd(input_directory)
  files <- list.files()[grep("*abundance.tsv", list.files())]
  names(files) <- gsub("_abundance.tsv", "", files)
  txi <- tximport(files, type = "kallisto", tx2gene = tx2gene, countsFromAbundance="lengthScaledTPM")
  counts <- txi$counts
  setwd(currwd)
  row.names(counts) <- gsub("\\..*", "", row.names(counts))
  write.table(counts, "counts.tsv", sep = '\t')
  return(counts)
}

#' Aggregate counts by external gene name
#'
#' This function aggregates a counts table by external gene name, using
#' Ensembl gene ID to map the counts to external gene names.
#'
#' @param counts A counts table.
#' @return The aggregated counts table.
aggregate_counts_by_gene_name <- function(counts, ensembl_dataset) {
  ensembl <- useMart('ensembl', dataset = ensembl_dataset)
  attributes = c("ensembl_gene_id", "external_gene_name")
  gene_info = getBM(attributes = attributes,
                    filters = 'ensembl_gene_id', 
                    values = row.names(counts),
                    mart = ensembl)
  row.names(gene_info) <- gene_info$ensembl_gene_id
  counts <- counts[gene_info$ensembl_gene_id,]
  counts <- aggregate(counts, by=list(gene_info$external_gene_name), FUN='sum')
  row.names(counts) <- counts$Group.1
  counts <- counts[,2:ncol(counts)]
  counts <- as.matrix(counts)
  counts <- counts[-1,]
  return(counts)
}



#' Convert counts table to ExpressionSet object
#'
#' This function takes a counts table and returns an ExpressionSet object.
#'
#' @param counts A counts table.
#' @return An ExpressionSet object.
bulkdata <- function(counts){
  bulkexp_set <- ExpressionSet(assayData = counts)
  bulk.mtx = exprs(bulkexp_set)
  return(bulk.mtx)
}


## Single Cell RNA Seq

#' Read SingleCellExperiment object from file
#'
#' This function reads a SingleCellExperiment (SCE) object from a file and
#' returns the object.
#'
#' @param file Path to the file containing the SCE object.
#' @return The SCE object.
readSCE <- function(file){
  sce <- readRDS(file)
  return(sce)
}

## MuSiC

#' Estimate cell type proportions using MuSiC
#'
#' This function estimates cell type proportions using the MuSiC package. It takes
#' an expression matrix and a SingleCellExperiment object and returns the
#' MuSiC estimations object.
#'
#' @param bulk.mtx An expression matrix.
#' @param sc.sce A SingleCellExperiment object.
#' @return The MuSiC estimations object.
estimate_cell_type_proportions <- function(bulk.mtx, sc.sce) {
  estimations <- music_prop(bulk.mtx = bulk.mtx, 
                            sc.sce = sc.sce, 
                            clusters = 'celltypes',
                            samples = 'sampleIDs', 
                            select.ct = NULL , 
                            verbose = T)
  est_cell_types <- t(estimations$Est.prop.weighted)
  return(est_cell_types)
}

## Input Validation

#rds_type <- function(string) {
#  if (!grepl("\\.RDS$", string)) {
#    stop(paste0(string, " is not an RDS file"), call.=FALSE)
#  }
#  return(string)
#}
#
#dir_type <- function(string) {
#  if (!dir.exists(string)) {
#    stop(paste0(string, " is not a directory"), call.=FALSE)
#  }
#  return(string)
#}

####################### Pipeline #######################

# Parse command line arguments
arg_parser <- ArgumentParser()
arg_parser$add_argument("-s", "--sce", dest="sce_file",
                        required=TRUE,
                        help="path to the Single Cell Expression Object within the Powrie Album")

arg_parser$add_argument("-b", "--bulk", dest="bulk_dir",
                        required=TRUE,
                        help="path to the directory containing the bulk Kalisto outputs")

arg_parser$add_argument("-t", "--type", dest="species",
                        required=TRUE,
                        help="specify if the data is mouse or human derived")

arguments <- arg_parser$parse_args()

# Detect and load the appropriate species

ensembl_dataset <- list(species_handler(tolower(arguments$species)))[[1]]
t2g <- list(species_handler(tolower(arguments$species)))[[2]]

# Load the single cell experiment object
sce <- readSCE(arguments$sce_file)

# Convert the kalisto output to a bulk mtx object
counts <- convert_to_counts(t2g, arguments$bulk_dir)

counts <- aggregate_counts_by_gene_name(counts, ensembl_dataset)

bulk.mtx <- bulkdata(counts = counts)

# Run MuSiC and isolate the MuSiC estimations
est_cell_types <- estimate_cell_type_proportions(bulk.mtx, sce)
write.table(est_cell_types, "est_cell_types.tsv", sep="\t")

print('Done!')





