#!/bin/bash

#SBATCH -t 40:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=32
#SBATCH --account=evolgen
#SBATCH -D ./

module load gcc bcftools

bcftools mpileup --skip-indels -f /project/evolgen/vshastry/wolbachia/genome/norm.SuperWolbachia_sequence.fasta -b /project/evolgen/vshastry/wolbachia/matches/bamfiles/Super/sortedbamfiles/listoffiles.txt -o wolbachia_super_16Dec20.bcf -O b --threads 32
