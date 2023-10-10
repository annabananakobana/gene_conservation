#!/bin/bash -e

cat /nesi/nobackup/uoo02820/gene_conservation_multiz/assembly_databases/genomes/candidate_genes.txt | while read -r gene; do
	cat /nesi/nobackup/uoo02820/gene_conservation_multiz/assembly_databases/genomes/fasta_file_list.txt | while read -r file; do
	gene_name=$(echo "$gene" | sed -e 's/gene=//g'| sed -e 's/;//g')
	spp=$(echo "$file" | sed -e 's/_genomic.*$//g'| sed -e 's/^GCF_[0-9]*\.//g' | sed -e 's/[0-9]*_//g')
	assembly=$(basename -- "${file%.*}")
	echo "$spp"
	cat /nesi/nobackup/uoo02820/gene_conservation_multiz/assembly_databases/genomes/final_fasta_files/"$gene_name"/"$file" | sed -e "s/>/>$spp\_/g" >/nesi/nobackup/uoo02820/gene_conservation_multiz/assembly_databases/genomes/final_fasta_files/"$gene_name"/"$assembly"\_labelled.fasta;
	done;
done
