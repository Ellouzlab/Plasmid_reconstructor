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
Tests performance using complete genomes

**utils.py**  
Functions to help abstract Plasmid_recon.py, Build_plas_db.py, Performance_test.py and download_db.py

## Installation

1. Clone this repository and enter directory:    
`git clone https://github.com/Ellouzlab/Plasmid_reconstructor`    
`cd Plasmid_reconstructor`
2. Create environment and activate it:    
`conda env create -f environment.yml -n plas`    
`conda activate plas`    

## Usage    
1. Begin by downloading the databases using the initialize_db.py script    
`python initialize_db.py`
Note: if you are not searching for Pantoea agglomerans plasmids, you need to provide a chromosome from your species of interest.


