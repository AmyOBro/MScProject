# MScProject
## "An examination of genomic resistance and virulence of _Escherichia marmotae_ identified in deer faeces Portumna, Galway."

### Scripts

This repository holds scripts used for my masters project which involved the assembly of 11 strains of _Escherichia coli_ and _Escherichia marmotae_, and the investigation of their virulence and resistance gene content.

The scripts in this repository were used to do the following:

1. `assembly_pipeline.sh` - This script was used to trim raw paired-end reads in .fastq.gz files and assemble trimmed reads using the shovill pipeline. It also runs prokka for annotation of the assembly, quast for assemby quality statistics and ABRicate and AMRFinder for the identification of resistance, virulence and plasmid genes from assemblies.

2. `Project_Heatmaps.Rmd` This is an R markdown document containing the Copy_number_matrix() function used in R to format ABRicate and AMRFinder outputs into copy number matrices. It also contains the commands where these matrices were used as input for the pheatmap() function to create the heatmaps available as supplementary figures in the manuscript. It requires `names_total.csv` and `heatmaps_metadata.csv` to correctly name strains in the matrices and provide metadata for heatmaps.

3. `Plasmid_contig_genes.Rmd` This is an R markdown document containing the Plasmid_contents() function used in R to find virulence or resistance genes on the same sequence as a plasmid gene using the ABRicate PlasmidFinder data and another ABRicate or AMRFinder dataset. It also contains the commands used to get the output that was used for supplememtary table 2 of the manuscript. It requires `names_total.csv` to correctly name strains.

4. `database_amr.sh` was used to run abricate and AMRFinder for .gff, .fna and .faa files associated with genomes downloaded on NCBI, and `database_prokka.sh` was used to re-annotate the _E. marmotae_ genomes downloaded from NCBI to get the input .gff files used for Roary.

### Commands

Outside of these scripts, the following commands were used to run other software:

1. Running Roary: `singularity run ../scripts/roary.img roary -p 8 roary/*gff`

2. Running Parsnp: `singularity run ../scripts/parsnp.img parsnp -r ../phylogeny/assembled/HT073016-2_annotated.fna -d ../phylogeny/assembled/`
