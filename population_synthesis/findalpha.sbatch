#!/bin/bash
#SBATCH --time=04:00:00
#SBATCH --partition=short.q
#SBATCH --mem=0
#SBATCH --cpus-per-task=1

# Get input name
name="$1"

if [ -z "$name" ]; then
    echo "Error: No name provided."
    echo "Usage: sbatch findalpha.sbatch NAME"
    exit 1
fi

# Load modules
module load gcc/10
module load anaconda/2/2019.03

# Run Python script
python2 find_alpha.py "$name"

