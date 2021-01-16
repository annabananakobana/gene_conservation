#!/bin/bash -e
#script by Alana Alexander, adapted by Anna C Clark

module load LASTZ/1.04.03-GCC-9.2.0

cat /nesi/nobackup/uoo02820/gene_conservation_multiz/assembly_databases/genomes/candidate_genes_no_Il11ra.txt | while read -r gene; do
                cd /nesi/nobackup/uoo02820/gene_conservation_multiz/assembly_databases/genomes/final_fasta_files/$gene/
        # Counting the number of files for when we step through things
        numfiles=`wc -l /nesi/nobackup/uoo02820/gene_conservation_multiz/assembly_databases/genomes/rodent_labelled_files.txt | awk '{print $1}'`;
        # Stepping through the 2nd through last filenames to make the query file
        for j in `seq 2 $numfiles`;
        # Getting a tempfilename from filename list
                do tempfilename=`head -n $j /nesi/nobackup/uoo02820/gene_conservation_multiz/assembly_databases/genomes/rodent_labelled_files.txt | tail -n 1`
                # Catting all these into one file
                cat $tempfilename >> temp_query;
        done;
        # Getting the target sequence
        tempfilename=`head -n 1 /nesi/nobackup/uoo02820/gene_conservation_multiz/assembly_databases/genomes/rodent_labelled_files.txt | tail -n 1`;
        cat $tempfilename >> temp_target;
        # Now we've got all the bits to run lastz
        lastz temp_target[multiple] temp_query --rdotplot=$gene\_rdotplot_rodent_5_20_75.txt --census=$gene\_census_stats_rodent_5_20_75.txt --mismatch=5,20 --identity=75 --output=$gene\_lastz_rodent_5_20_75.maf --format=maf;
        # Finally need to remove temp_target and temp_query for next time round
        rm temp_target;
        rm temp_query;
done
