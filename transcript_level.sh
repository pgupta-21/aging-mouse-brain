#!/bin/bash
#SBATCH -p fat
#SBATCH -n 24
#SBATCH -N 1
#SBATCH -C scratch
#SBATCH -t 1-10:00:00
#SBATCH —mail-type=ALL

module load salmon

# Downloading the reference transcriptome and genome
wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M23/gencode.vM31.transcripts.fa.gz
wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M23/GRCm39.primary_assembly.genome.fa.gz

# Creating metadata and decoy for indexing
grep "^>" <(gunzip -c GRCm39.primary_assembly.genome.fa.gz) | cut -d " " -f 1 > decoys.txt
sed -i.bak -e 's/>//g' decoys.txt
cat gencode.vM31.transcripts.fa.gz GRCm39.primary_assembly.genome.fa.gz > gentrome.fa.gz

# Building the index using salmon
salmon index -t ~/salmon_fa/gentrome.fa.gz -d decoys.txt -p 12 -i salmon_index --gencode

# Quantification using salmon
for file in ~/fastq_files/*.fastq
do
        salmon quant -i ~/salmon_fa/salmon_index --libType A -r ${file} -o ${file}_quant --writeMappings=${file}.sam
done
