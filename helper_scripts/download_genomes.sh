#!/bin/bash

# BASH script to download genomes of interest into the current working
# directory from iGenomes:
# https://emea.support.illumina.com/sequencing/sequencing_software/igenome.html
# A website tool was used to assist with the downloading:
# https://ewels.github.io/AWS-iGenomes/
#
# This script requires AWS to be installed
#
# The script also needs to be run on a version of the Linux in which the command 
# "rename" takes regular expression.  Some versions of Linux "rename" do not do this.

# Make folder
mkdir Nextflow_iGenomes
cd Nextflow_iGenomes/

####################
#####################
# HUMAN
#####################
#####################

#####################
# GRCh38
#####################

# FASTA
################
aws s3 --no-sign-request --region eu-west-1 sync s3://ngi-igenomes/igenomes/Homo_sapiens/NCBI/GRCh38/Sequence/Chromosomes/ ./references/Homo_sapiens/NCBI/GRCh38/Sequence/Chromosomes/

#Move non-standard chromosomes
cd references/Homo_sapiens/NCBI/GRCh38/Sequence/Chromosomes
rename s/^chr/Homo_sapiens_GRCh38_chr/ chr*
mkdir extra_sequences
mv *.fa extra_sequences/
mv extra_sequences/*chr[a-z].* .
mv extra_sequences/Homo_sapiens_GRCh38_chr[0-9].* .
mv extra_sequences/Homo_sapiens_GRCh38_chr[0-9][0-9].* .
cd -

# Bowtie2
###############
# Needs to be made de novo, since iGenomes Bowtie2 index include non-regular
# sequences
mkdir references/Homo_sapiens/NCBI/GRCh38/Sequence/Chromosomes/Bowtie2Index

# GTF
###############
aws s3 --no-sign-request --region eu-west-1 sync s3://ngi-igenomes/igenomes/Homo_sapiens/NCBI/GRCh38/Annotation/Genes/ ./references/Homo_sapiens/NCBI/GRCh38/Annotation/Genes/ --exclude "*" --include "genes.gtf"
mv references/Homo_sapiens/NCBI/GRCh38/Annotation/Genes/genes.gtf references/Homo_sapiens/NCBI/GRCh38/Annotation/Genes/Homo_sapiens_GRCh38_genes.gtf


####################
#####################
# Mouse
#####################
#####################

#####################
# GRCm38
#####################

# FASTA
################
aws s3 --no-sign-request --region eu-west-1 sync s3://ngi-igenomes/igenomes/Mus_musculus/NCBI/GRCm38/Sequence/Chromosomes/ ./references/Mus_musculus/NCBI/GRCm38/Sequence/Chromosomes/
cd references/Mus_musculus/NCBI/GRCm38/Sequence/Chromosomes/
rename -e s/^/Mus_musculus_GRCm38_chr_/ *.fa
cd -

# Bowtie2
###############
aws s3 --no-sign-request --region eu-west-1 sync s3://ngi-igenomes/igenomes/Mus_musculus/NCBI/GRCm38/Sequence/Bowtie2Index/ ./references/Mus_musculus/NCBI/GRCm38/Sequence/Bowtie2Index/
cd references/Mus_musculus/NCBI/GRCm38/Sequence/Bowtie2Index/
rename -e s/genome/Mus_musculus_GRCm38/ *
cd -

# GTF
###############
aws s3 --no-sign-request --region eu-west-1 sync s3://ngi-igenomes/igenomes/Mus_musculus/NCBI/GRCm38/Annotation/Genes/ ./references/Mus_musculus/NCBI/GRCm38/Annotation/Genes/ --exclude "*" --include "genes.gtf"
mv references/Mus_musculus/NCBI/GRCm38/Annotation/Genes/genes.gtf references/Mus_musculus/NCBI/GRCm38/Annotation/Genes/Mus_musculus_GRCm38.gtf


####################
#####################
# Macaque
#####################
#####################

#####################
# GRCm38
#####################

# FASTA
################




echo "Done"
