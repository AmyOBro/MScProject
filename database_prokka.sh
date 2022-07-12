#!/usr/bin/env bash

#Note: This script requires comma separated metadata file where the 5th column corresponds to the strain name, on the same line as the name of the input folder. The strain name will be used to name outputs.
#This script was used to reannotate genomes downloaded from NCBI using the NCBI datasets tool to obtain .gff files of the format accepted by prokka and the .fna files as input for parsnp.
#Script should be called using "bash database_prokka.sh /file/path"

#Setting inputs

    #Set current directory
    current_dir=$1

    #Path to image files
    image_path="/home/amy/Documents/Genomics_Project/scripts"
    echo "Images are found are $(ls $image_path/prokka.img)"

    #Finding the name of the sample, assuming the set folder has the sample name
    samp="${current_dir%"${current_dir##*[!/]}"}"
    samp="${samp##*/}"
    strain="$( grep "$samp" /home/amy/Documents/Genomics_Project/genomes/E.mar_dataset/genomes_info.csv | awk -F "," '{ print $5 }' )"
    if [[ $strain == "NA" ]]; then strain=$samp
    else strain=$strain
    fi
    echo "Sample is $strain"

    #Setting the genus and species names
    genus="Escherichia"
    species="marmotae"
    echo "Species is $genus $species"

#Running Prokka

    #Finding the assembly file
    assembly=$1/$samp*.fna
    mkdir $1/prokka
    cp $assembly $1/prokka/
    mv $1/prokka/* $1/prokka/$samp.fa
    assembly=$1/prokka/$samp.fa
    echo "assembly path is $assembly"

    #Running the prokka pipeline
    singularity run $image_path/prokka.img \
        prokka --outdir $1/prokka_output \
        --cpus 4 \
        --genus "$genus" \
        --species "$species" \
        --prefix "$strain"_annotated \
	--locustag "$samp" \
	--centre X \
	--compliant \
        $assembly
