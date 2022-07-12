#!/usr/bin/env bash

#Note: This script was used to search files for Escherichia genomes downloaded from NCBI for virulence and resistance genes using ABRicate and  AMRFinder.
#Input directory should be named after the accession or strain of the genome of interest as the output will be named after the input directory.
#Script should be called using "bash database_amr.sh /file/path"

    #Path to annotation
    current_dir=$1
    echo "Input directory is $current_dir"

    #Output path
    out_dir="/home/amy/Documents/Genomics_Project/genomes/"
    echo "Output directory is located at $out_dir"

    #Path to image files
    image_path="/home/amy/Documents/Genomics_Project/scripts"
    echo "Images are found are $(ls $image_path/*.img)"

    #Finding the name of the strain
    current_dir=$1
    strain="${current_dir%"${current_dir##*[!/]}"}"
    strain="${strain##*/}"
    fasta=$1/$strain*.fna
    echo "Strain is $strain"
    echo "File to use: $fasta"

    #Running abricate against CARD, E.coli VF and plasmidfinder databases using the contig nucleotide file
    singularity run $image_path/abricate.img \
    	abricate --db card $fasta > $out_dir/card/"$strain"_card.tab
    singularity run $image_path/abricate.img \
    	abricate --db ecoli_vf $fasta > $out_dir/ecolivf/"$strain"_ecoli_vf.tab
    singularity run $image_path/abricate.img \
       abricate --db vfdb $fasta > $out_dir/vfdb/"$strain"_vfdb.tab
    singularity run $image_path/abricate.img \
    	abricate --db plasmidfinder $fasta > $out_dir/plasmidfinder/"$strain"_plasmidfinder.tab

    #Path to amr_finder directory
    amr_dir="/home/amy/Documents/Genomics_Project/genomes/amrfinder"
    echo "Output directory is $amr_dir"
    prokka_dir="$1/prokka_output"

    #Path to amr_finder directory
    amr_dir="/home/amy/Documents/Genomics_Project/genomes/amrfinder"
    echo "Output directory is $amr_dir"

    #Running AMR-finder for Escherichia
    singularity run $image_path/amrfinder.img \
	amrfinder --plus \
 	--protein $1/*.faa \
 	--nucleotide $1/"$strain"*.fna \
        --gff $1/*.gff \
        --organism "Escherichia" \
        --name $strain \
        --output $amr_dir/"$strain"_amrfinder.txt \
        --mutation_all $amr_dir/"$strain"_muations.txt \
        --protein_output $amr_dir/"$strain"_amr_protein.fasta \
        --nucleotide_output $amr_dir/"$strain"_amr_nucleotide.fasta \
        --log $amr_dir/"$strain"_amrfinder_log.txt
