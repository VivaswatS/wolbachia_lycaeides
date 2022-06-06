#!/bin/bash

#SBATCH -t 48:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=32
#SBATCH --account=evolgen

module load gcc bcftools

bcftools call -m --ploidy 1 --variants-only --format-fields GQ,GP --skip-variants indels wolbachia_super_16Dec20.bcf | bcftools filter --set-GTs . | bcftools view -m2 -M2 -v snps --apply-filter "PASS" --output-type v --output-file rawvariants_super_16Dec20.vcf &> log_bcftools_super_16Dec20.txt

