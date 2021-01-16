#!/bin/bash -e
#SBATCH --account=uoo02820
#SBATCH --job-name=extract_CDS
#SBATCH --time=00:15:00
#SBATCH --mem=6GB
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mail-user=anna.clark@postgrad.otago.ac.nz
#SBATCH --mail-type=ALL
#SBATCH --qos=debug
#SBATCH --chdir=/nesi/nobackup/uoo02820/gene_conservation_multiz/assembly_databases/genomes

#load modules
module load BEDTools/2.29.2-GCC-9.2.0
module load BEDOPS/2.4.30-gimkl-2017a

#get regions from GFF files
gffs=(/nesi/nobackup/uoo02820/gene_conservation_multiz/assembly_databases/GFF_files/*.gff)
for gff in "${gffs[@]}"; do
	gff_name="$(basename -- ${gff%.*})"
	cat /nesi/nobackup/uoo02820/gene_conservation_multiz/assembly_databases/genomes/candidate_genes.txt | while read -r gene; do
		gene_name=$(echo "$gene" | sed -e 's/gene=//g'| sed -e 's/;//g')
		mkdir -p /nesi/nobackup/uoo02820/gene_conservation_multiz/assembly_databases/genomes/final_fasta_files/$gene_name
		grep -i $gene $gff | \
		awk '{if ($3=="CDS") print $0}' | \
		convert2bed --input=gff | sed -e 's/ /_/g' | bedtools merge -i - -s -c 6 -o distinct | \
		#select for large regions (>80nt = 4x20nt gRNA)
		awk 'BEGIN{OFS="\t";} {if (($3 - $2)>80) print $0}' \
		>/nesi/nobackup/uoo02820/gene_conservation_multiz/assembly_databases/genomes/final_fasta_files/$gene_name/$gff_name\_CDS.bed;
	done;
done

#get seqs
genomes=(/nesi/nobackup/uoo02820/gene_conservation_multiz/assembly_databases/genomes/*.fna)
for genome in "${genomes[@]}"; do
        species="$(basename -- ${genome%.*})"
	cat /nesi/nobackup/uoo02820/gene_conservation_multiz/assembly_databases/genomes/candidate_genes.txt | while read -r gene; do
		gene_name=$(echo "$gene" | sed -e 's/gene=//g'| sed -e 's/;//g')
		bedtools getfasta -s -fi $genome -bed /nesi/nobackup/uoo02820/gene_conservation_multiz/assembly_databases/genomes/final_fasta_files/$gene_name/$species\_CDS.bed \
		-fo /nesi/nobackup/uoo02820/gene_conservation_multiz/assembly_databases/genomes/final_fasta_files/$gene_name/$species\_CDS.fasta;
	done;

done
