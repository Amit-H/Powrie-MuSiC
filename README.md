# Powrie-MuSiC


## Description: 
This tool is designed to provide a command line interface for cell type proportion estimations using the MuSiC algorithm implemented by Xuran Wang https://doi.org/10.1038/s41467-018-08023-x. The tool takes bulk RNA Seq inputs from Kalisto outputs and references them versus a single cell reference set to generate cell proportion estimations. The single cell reference set is known as the Powrie Album, and can be found on the shared drive of the BMRC. For brevity, all functionality will be contained in this script.

## Usage Instructions:
To run the script: `Rscript powrieMuSiC.r -s SINGLE_CELL_EXPERIMENT_OBJECT -b BULK_DATA_FOLDER -t SPECIES_TYPE`

- -s The single cell experiment object, which you can find in data/PowrieAlbum/*SPECIES*/
- -b The folder where your bulk RNA seq data is located. This currently only works with Kallisto outputs.
- -t Mouse or Human. Other species are not yet supported.

Mouse single cell reference data is generated using data from C57BL/6 mice. You can find the sources of the reference data in the markdown file PowrieAlbum.md 

## To do list:
1. Add input validation

## Branch structure

For development purposes:
1. main - latest working version, works fine for deployment
2. dev - developement, for internal use, may not work as expected or at all

