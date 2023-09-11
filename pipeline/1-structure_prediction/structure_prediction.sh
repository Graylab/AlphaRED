#!/bin/bash
module unload python
module load cuda/11.1.0
module load anaconda

#Activate the conda environment
conda activate colabfold-conda

# Output dir is the directory where the structure will be stored. 
# fasta_dir is the directory that has the fasta files. You can also provide a single fasta file, i.e. the input from the user must be written to a .fasta file and should be provided to this command below. 
OUTPUT_DIR=results/
FASTA_DIR=fastas/
mkdir -p $OUTPUT_DIR

# We run colabfold for generating 5 predictions using the latest version of AlphaFold_v3
colabfold_batch $FASTA_DIR $OUTPUT_DIR --num-models 1 --model-type alphafold2_multimer_v3