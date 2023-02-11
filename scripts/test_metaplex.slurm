#!/bin/bash

#SBATCH --partition=batch
#SBATCH --job-name=metaplex-benchmark
#SBATCH --output=res.txt

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=64G
#SBATCH --time=00:30:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=cedric.simal@unamur.be

ml load Julia/1.8.2-linux-x86_64 

#julia install.jl

julia metaplex_test.jl

