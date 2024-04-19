#!/bin/bash

input_dir="/home/alumno10/ADO/ado_2/4"

output_dir="/home/alumno10/ADO/ado_2/analisis"

for file in "$input_dir"/*.fastq; do
    
    base=$(basename "$file" .fastq.gz)
    
   
    ./fastqc "$file" --outdir="$output_dir" --extract &
done


wait

echo "¡Análisis de calidad de los archivos FASTQ completado!"

