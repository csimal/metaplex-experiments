#!/bin/bash

#SBATCH --partition=batch
#SBATCH --job-name=metaplex-benchmark
#SBATCH --output=res.txt

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4G
#SBATCH --time=00:10:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=cedric.simal@unamur.be

julia metaplex_test.jl

