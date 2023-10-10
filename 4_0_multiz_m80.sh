#!/bin/bash -e
#SBATCH --account=uoo02820 
#SBATCH --job-name=multiz_m80
#SBATCH --time=00:15:00
#SBATCH --mem=6GB
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mail-user=anna.clark@postgrad.otago.ac.nz
#SBATCH --mail-type=ALL
#SBATCH --qos=debug

export PATH=/nesi/nobackup/uoo02820/gene_conservation_multiz/multiz:${PATH}

#female candidates 
cat /nesi/nobackup/uoo02820/gene_conservation_multiz/assembly_databases/genomes/candidate_genes.txt | while read -r gene; do
	cd /nesi/nobackup/uoo02820/gene_conservation_multiz/assembly_databases/genomes/final_fasta_files/$gene/
	mkdir -p results_5_20_75
	multiz M=80 $gene\_lastz_rodent_5_20_75.maf $gene\_lastz_nonrodent_5_20_75.maf 1 >results_5_20_75/"$gene"_multiz_align_5_20_75_m80.maf;
done
