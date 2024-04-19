#!/bin/bash

archivos_bam=(
    "alineamiento_78.bam"
    "alineamiento79.bam"
    "alineamiento80.bam"
    "alineamiento81.bam"
    "alineamiento82.bam"
    "alineamiento83.bam"
    "alineamiento84.bam"
    "alineamiento89.bam"
    "alineamiento89.bam"
    "alineamiento88.bam"
    "alineamiento87.bam"

)


for bam in "${archivos_bam[@]}"; do
    nombre_archivo="${bam%.bam}"
    samtools flagstat "$bam" > "flagstat_result/${nombre_archivo}.txt"
done

