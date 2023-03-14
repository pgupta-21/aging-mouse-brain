#!/bin/bash
#SBATCH -p fat
#SBATCH -n 24
#SBATCH -N 1
#SBATCH -C scratch
#SBATCH -t 1-10:00:00
#SBATCH —mail-type=ALL

module load fastqc
module load hisat2
module load samtools
module load subread

# Fasta file
wget http://ftp.ensembl.org/pub/release-103/fasta/mus_musculus/dna/Mus_musculus.GRCm39.dna.toplevel.fa.gz
gunzip Mus_musculus.GRCm39.dna.toplevel.fa.gz

# Gtf annotation file
wget http://ftp.ensembl.org/pub/release-103/gtf/mus_musculus/Mus_musculus.GRCm39.103.chr.gtf.gz
gunzip Mus_musculus.GRCm39.103.chr.gtf.gz

# Building index for HISAT2
hisat2-build /usr/users/pgupta2/Mus_musculus.GRCm39.dna.toplevel.fa /usr/users/pgupta2/mus_idx

# Alignment using HISAT2
for file in *.fastq; do
	hisat2 -p 8 --dta -x /usr/users/pgupta2/mus_index/mus_idx -U $file -S '${file}.bam'
done

# Sorting bam files
for file in *.bam; do
	samtools sort ~/align_files/bam_files/*.bam -o ~/align_files/bam_files/sorted_bam/${file}.sorted.bam
done

# Quantification using subread
featureCounts -T 1 -a ~/mus_anno/Mus_musculus.GRCm39.103.chr.gtf -o ~/align_files/fcount.txt ~/align_files/bam_files/sorted_bam/*.sorted.bam
