#!/bin/bash

echo "Running on $(hostname)"
echo "Input file: $1"

##activating conda environment
source ~/miniconda3/etc/profile.d/conda.sh
conda activate /cms/kaur/coffea_latest

# Run your processor
python processor.py "$1" 
#python processor.py "$1" --outdir /cms/kaur/output/

# deactivate conda on the node
conda deactivate

