#!/bin/bash

# Directorio donde se encuentran los archivos FASTQ
input_dir="/home/alumno10/ADO/ado_2/4"

# Directorio donde se guardarán los resultados de FastQC
output_dir="/home/alumno10/ADO/ado_2/analisis"

# Bucle para analizar cada archivo .fastq.gz en el directorio de entrada
for file in "$input_dir"/*.fastq; do
    # Nombre base del archivo sin la extensión .fastq.gz
    base=$(basename "$file" .fastq.gz)
    
    # Ejecutar FastQC en modo no interactivo para el archivo actual
    ./fastqc "$file" --outdir="$output_dir" --extract &
done

# Esperar a que todos los procesos FastQC finalicen antes de continuar
wait

echo "¡Análisis de calidad de los archivos FASTQ completado!"

