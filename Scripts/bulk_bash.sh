#!/bin/bash

# Check if the correct number of arguments are provided
if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <input_directory> <output_directory>"
    exit 1
fi

# Assign arguments to variables
input_dir=$1
output_dir=$2

# Check if input directory exists
if [ ! -d "$input_dir" ]; then
    echo "Input directory not found."
    exit 1
fi

# Check if output directory exists, if not, create it
if [ ! -d "$output_dir" ]; then
    mkdir -p "$output_dir"
fi

# Run blastn on each FASTA file in the input directory
for file in "$input_dir"/*.fa; do
    filename=$(basename "$file")
    filename_noext="${filename%.*}"
    output_file="$output_dir/$filename_noext.blastn_out"

    # Run blastn command
    blastn -query "$file" -db /home/apaulson/working_dir/Local_blast/genomes/human -out "$output_file" -num_threads 4 -task blastn-short -evalue 0.8 -word_size 7 -gapopen 5 -gapextend 2 -perc_identity 100 -reward 1 -penalty -2 -outfmt '10 qseqid sseqid length qlen slen qstart qend sstart send pident mismatch gapopen evalue bitscore '
done

echo "Blastn run completed."
