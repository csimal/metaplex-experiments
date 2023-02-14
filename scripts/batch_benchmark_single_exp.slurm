#!/bin/bash

#SBATCH --partition=batch
#SBATCH --job-name=single_exp-benchmark
#SBATCH --output=output/res_%x_%j.txt
#SBATCH --output=output/error_%x_%j.txt

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=32G
#SBATCH --time=00:45:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=cedric.simal@unamur.be

ml load Julia/1.8.2-linux-x86_64 


julia benchmark_single_experiment.jl

