#!/usr/bin/env bash

#Note: This script requires a main directory containing fastq read files and a file named "species_file.txt" with the genus and species on separate lines.
#It also assumes that the main directory is named after the strain, and uses that name as a label for outputs.
#Software should be available in the folder at the image_path in the form of singularity images. See supplementary information in the manuscript for details of where images were sourced.
#Reads should be paired end, and contained within separate gzipped fastq files.
#Script should be called using "bash assembly_pipeline.sh /folder/path"

#Setting inputs

    #Moving fastq files to trimmomatic folder, if they weren't there already
    current_dir=$1
    mkdir $current_dir/trimmomatic
    mv $current_dir/*.fastq.gz $current_dir/trimmomatic

    #Setting path to fastq files
    fastq_dir=$current_dir/trimmomatic
    pair_1_raw=$(find $fastq_dir -name "*R1*.fastq.gz")
    echo "Read 1 is $pair_1_raw"
    pair_2_raw=$(find $fastq_dir -name "*R2*.fastq.gz")
    echo "Read 2 is $pair_2_raw"
    echo "fastq directory is $fastq_dir"

    #Path to image files
    image_path="/home/amy/Documents/Genomics_Project/scripts"
    echo "Images are found are $(ls $image_path/*.img)"

    #Finding the name of the sample, assuming the set folder has the sample name
    strain="${current_dir%"${current_dir##*[!/]}"}"
    strain="${strain##*/}"
    echo "Sample is $strain"

    #Setting the genus and species names
    mapfile -t lines <$1/species_file.txt
    echo "Species is ${lines[0]} ${lines[1]}"

    #Shovill and kraken2 images
    shovill_img=$image_path/shovill.img
    kraken_img=$image_path/kraken2.img

#Running first fastQC on raw reads
    mkdir $1/trimmomatic/fastqc_run1
    singularity run $image_path/fastqc.img fastqc -o $1/trimmomatic/fastqc_run1 \
        --extract --nogroup --format fastq --threads 6 --dir $1/trimmomatic \
        $pair_1_raw $pair_2_raw

#Making Kraken2 directory and running initiall check on reads
    mkdir $1/kraken2
    mkdir $1/kraken2/run1
    singularity run $kraken_img kraken2 \
        --db "/home/amy/Documents/Genomics_Project/kraken2-db/standard_16GB" \
        --memory-mapping \
        --threads 6 \
        --report $1/kraken2/run1/"$strain"_kraken2.txt \
        --use-names \
        --output $1/kraken2/run1 \
        --paired $pair_1_raw $pair_2_raw

#Running trimmomatic on reads, using the trimmomatic executable from the shovill image
    singularity run $shovill_img trimmomatic PE -threads 6 -phred33 \
        $pair_1_raw $pair_2_raw $1/trimmomatic/R1.fq.gz /dev/null $1/trimmomatic/R2.fq.gz /dev/null CROP:300 \
        ILLUMINACLIP:/shovill/shovill-1.1.0/db/trimmomatic.fa:1:30:11 LEADING:3 TRAILING:3 SLIDINGWINDOW:10:20 MINLEN:30 TOPHRED33 2>&1 | sed 's/^/[trimmomatic] /' | tee -a shovill.log

#Setting reads names
    pair_1=$1/trimmomatic/R1.fq.gz
    echo "Read 1 is $pair_1"
    pair_2=$1/trimmomatic/R2.fq.gz
    echo "Read 2 is $pair_2"

#Running first fastQC on raw reads
    mkdir $1/trimmomatic/fastqc_run2
    singularity run $image_path/fastqc.img fastqc -o $1/trimmomatic/fastqc_run2 \
        --extract --nogroup --format fastq --threads 6 --dir $1/trimmomatic \
        $pair_1 $pair_2

#Running shovill
    singularity run $shovill_img shovill \
	--cpus 6 \
	--ram 15 \
	--outdir $1/shovill \
	--R1 $pair_1 --R2 $pair_2

#Running Kraken2 a second time on the assembly
    mkdir $1/kraken2/run2
    singularity run $kraken_img kraken2 \
	--db "/home/amy/Documents/Genomics_Project/kraken2-db/standard_16GB" \
	--memory-mapping \
	--threads 6 \
	--report $1/kraken2/run2/"$strain"_kraken2.txt \
	--use-names \
	--output $1/kraken2/run2 \
	$1/shovill/contigs.fa

#Running Prokka

    #Copying assembly to prokka folder
    assembly=$1/shovill/contigs.fa
    echo "assembly path is $assembly"

    #Running the prokka pipeline
    singularity run $image_path/prokka.img \
        prokka --outdir $1/prokka \
        --cpus 4 \
        --genus "${lines[0]}" \
        --species "${lines[1]}" \
        --prefix "$strain"_annotated \
	--centre X \
	--compliant \
        $assembly

#Running Quast

    #Finding the prokka genome feature file output and each paired end read
    genome_feat=$(find $current_dir/prokka -name "*annotated.gff")
    echo "feature file is $genome_feat"

    #Running quast
    singularity run $image_path/quast.img python /usr/local/bin/quast.py $assembly \
	-g $genome_feat \
	--pe1 $pair_1 \
	--pe2 $pair_2 \
        -o $1/quast \
	-t 4 \
	-m 0

#Running abricate

    #Path to annotation
    prokka_dir=$1/prokka
    echo "Input directory is $prokka_dir"
    echo "File to use: $prokka_dir/"$strain"_annotated.fna"

    #Output path
    mkdir $1/abricate
    abr_dir=$1/abricate
    echo "Output will be located at $abr_dir"

    #Running the abricate against ncbi, CARD, resfinder and plasmidfinder databases using the contig nucleotide file
    singularity run $image_path/abricate.img \
        abricate --db card $prokka_dir/"$strain"_annotated.fna > $abr_dir/"$strain"_card.tab
    singularity run $image_path/abricate.img \
        abricate --db plasmidfinder $prokka_dir/"$strain"_annotated.fna > $abr_dir/"$strain"_plasmidfinder.tab
    singularity run $image_path/abricate.img \
        abricate --db vfdb $prokka_dir/"$strain"_annotated.fna > $abr_dir/"$strain"_vfdb.tab
    singularity run $image_path/abricate.img \
        abricate --db ecoli_vf $prokka_dir/"$strain"_annotated.fna > $abr_dir/"$strain"_ecoli_vf.tab

    #Creating a summary file
    singularity run $image_path/abricate.img \
        abricate --summary $abr_dir/"$strain"_card.tab $abr_dir/"$strain"_plasmidfinder.tab $abr_dir/"$strain"_vfdb.tab $abr_dir/"$strain"_ecoli_vf.tab > $abr_dir/"$strain"_summary.tab

#Running amrfinder

    #Path to amr_finder directory
    mkdir $1/amrfinder
    amr_dir=$1/amrfinder
    echo "Output directory is $amr_dir"

    #Setting genus and species variables
    genus=${lines[0]}
    species=${lines[1]}

    #Reformatting the gff file to make it compatible with amrfinder
    perl -pe '/^##FASTA/ && exit; s/(\W)Name=/$1OldName=/i; s/ID=([^;]+)/ID=$1;Name=$1/' $prokka_dir/"$strain"_annotated.gff  > $prokka_dir/"$strain"_amrfinder.gff
    echo "Created reformatted gff file $prokka_dir/"$strain"_amrfinder.gff"

    #Species list and genus list
    valid_species="Acinetobacter_baumannii Clostridioides_difficile Enterococcus_faecalis Enterococcus_faecium Pseudomonas_aeruginosa Staphylococcus_aureus Staphylococcus_pseudintermedius Streptococcus_agalactiae Streptococcus_pneumoniae Streptococcus_pyogenes Vibrio_cholerae"
    valid_genus="Escherichia Campylobacter Klebsiella Neisseria Salmonella"

    #Changing the command based on whether a valid species is present
    if [[ $valid_species =~ (^|[[:space:]])""$genus"_"$species""($|[[:space:]]) ]]; then
	      echo "Valid species "$genus"_"$species" found"
	      singularity run $image_path/amrfinder.img \
		 amrfinder --plus \
         	 --protein $prokka_dir/"$strain"_annotated.faa \
         	 --nucleotide $prokka_dir/"$strain"_annotated.fna \
         	 --gff $prokka_dir/"$strain"_amrfinder.gff \
         	 --organism "$genus"_"$species" \
         	 --name $strain \
         	 --output $amr_dir/"$strain"_amrfinder.txt \
         	 --mutation_all $amr_dir/"$strain"_muations.txt \
         	 --protein_output $amr_dir/"$strain"_amr_protein.fasta \
        	 --nucleotide_output $amr_dir/"$strain"_amr_nucleotide.fasta \
	         --log $amr_dir/"$strain"_amrfinder_log.txt

    #Changing command based on if a valid genus is present
    elif [[ $valid_genus =~ (^|[[:space:]])"$genus"($|[[:space:]]) ]]; then
              echo "Valid genus for "$genus"_"$species" found"
              singularity run $image_path/amrfinder.img \
              	 amrfinder --plus \
                 --protein $prokka_dir/"$strain"_annotated.faa \
                 --nucleotide $prokka_dir/"$strain"_annotated.fna \
                 --gff $prokka_dir/"$strain"_amrfinder.gff \
                 --organism "$genus" \
                 --name $strain \
                 --output $amr_dir/"$strain"_amrfinder.txt \
                 --mutation_all $amr_dir/"$strain"_muations.txt \
                 --protein_output $amr_dir/"$strain"_amr_protein.fasta \
                 --nucleotide_output $amr_dir/"$strain"_amr_nucleotide.fasta \
                 --log $amr_dir/"$strain"_amrfinder_log.txt

    #Command if no valid genus or species
    else
         echo "No valid species or genus for "$genus"_"$species" found"
         singularity run $image_path/amrfinder.img \
   	 	  amrfinder --plus \
        	  --protein $prokka_dir/"$strain"_annotated.faa \
	          --nucleotide $prokka_dir/"$strain"_annotated.fna \
	          --gff $prokka_dir/"$strain"_amrfinder.gff \
        	  --name $strain \
	          --output $amr_dir/"$strain"_amrfinder.txt \
	          --protein_output $amr_dir/"$strain"_amr_protein.fasta \
	          --nucleotide_output $amr_dir/"$strain"_amr_nucleotide.fasta \
	          --log $amr_dir/"$strain"_amrfinder_log.txt
    fi
