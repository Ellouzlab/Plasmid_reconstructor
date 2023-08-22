# Plasmid_reconstructor
A program to reconstruct plasmids using homology to plasmid clusters.

## Overview

**Plasmid_recon.py**  
Main reconstruction program.

**Build_plas_db.py**  
Builds the plasmid database for Plasmid_recon.

**initialize_db.py**  
Downloads ORI and plasmid databases.

**Performance_test.py**  
Tests performance using complete genomes. (Currently under work) 

**genome_frag.py**  
Randomly fragments genomes for testing.

**utils.py**  
Functions to help abstract Plasmid_recon.py, Build_plas_db.py, Performance_test.py and download_db.py

## Installation
1. Clone this repository and enter directory:    
`cd Plasmid_reconstructor`
2. Create environment and activate it:    
`conda env create -f environment.yaml -n plas`    
`conda activate plas`    

## Usage (_Pantoea agglomerans_ only)    
1. Begin by downloading the databases using the initialize_db.py script.    
`python initialize_db.py -t <your_num_threads>`    
2. Use Plasmid_recon.py on your fasta file.  
`python Plasmid_recon.py -i INPUT_FASTA   -o OUTPUT_DIR -t NUM_OF_THREADS (optional) `

Note: if you are not searching for Pantoea agglomerans plasmids, you need to provide a chromosome from your species of interest, and run  
`makeblastdb -dbtype nucl -in your_chromosome_fasta -out your_location`  

Use the `--chromdb` flags for usage of Plasmid_reconstructor.py

## Testing Performance
1. Create artificially fragmented genomes (all fragments from chromosome will be named chromosome, and plasmid with a number depending on the plasmid they are from).  
`python genome_frag.py -i Test_dataset -n NUM_FRAGS -o OUTDIR`

2. Run Plasmid_recon.py individually on fragmented genomes and see performance.

Note: An automated method to assess performance is being developed and will be available soon. 

## Options
To see option available for any program, please run the following:  
`python Program_name.py -h` 
