#!/bin/bash

echo "Isolate,#bases,Suspected_species_#1,Suspected_species_#1_pct,Suspected_species_#2,Suspected_species_#2_pct,Suspected_species_#3,Suspected_species_#3_pct,MLST,Coverage,Contigs,N50,Assembly_length,ABRicate_NCBI"
for file in $@
do
	file=${file##*/}
	NUMBER_BASES=$(grep -A3 summary fastp_out/"$file"_fastp.json | awk -F ':' '/total_bases/ {print $2}' | tr -d ',')
	KRAKENSPECIES_1=$(awk '$4 == "S" {print $6" "$7}' kraken_out/"$file"_kraken2_report.txt | awk 'NR == 1 {print $0}')
	KRAKENSPECIES_1_PCT=$(awk '$4 == "S" {print $1}' kraken_out/"$file"_kraken2_report.txt | awk 'NR == 1 {print $0}')
	KRAKENSPECIES_2=$(awk '$4 == "S" {print $6" "$7}' kraken_out/"$file"_kraken2_report.txt | awk 'NR == 2 {print $0}')
	KRAKENSPECIES_2_PCT=$(awk '$4 == "S" {print $1}' kraken_out/"$file"_kraken2_report.txt | awk 'NR == 2 {print $0}')
	KRAKENSPECIES_3=$(awk '$4 == "S" {print $6" "$7}' kraken_out/"$file"_kraken2_report.txt | awk 'NR == 3 {print $0}')
	KRAKENSPECIES_3_PCT=$(awk '$4 == "S" {print $1}' kraken_out/"$file"_kraken2_report.txt | awk 'NR == 3 {print $0}')
	MLST=$(awk '{print $3}' mlst/"$file".tsv)
	COVERAGE=$(awk '{print $1}' coverage_out/"$file".txt)
	NUMBER_CONTIGS=$(grep '# contigs (>= 0 bp' quast_out/"$file"/report.tsv|awk '{print $6}')
	N50=$(grep '^N50' quast_out/"$file"/report.tsv|awk '{print $2}')
	ASSEMBLY_LENGTH=$(grep '^Total length (>= 0 bp' quast_out/"$file"/report.tsv|awk '{print $6}')
	ABRICATE_NCBI=$(awk 'NR>1 {print $5}' abricate_out/"$file".tsv | tr "\n" "|" | sed 's/.$//')
	echo "${file},${NUMBER_BASES},${KRAKENSPECIES_1},${KRAKENSPECIES_1_PCT},${KRAKENSPECIES_2},${KRAKENSPECIES_2_PCT},${KRAKENSPECIES_3},${KRAKENSPECIES_3_PCT},${MLST},${COVERAGE},${NUMBER_CONTIGS},${N50},${ASSEMBLY_LENGTH},${ABRICATE_NCBI}"
done
