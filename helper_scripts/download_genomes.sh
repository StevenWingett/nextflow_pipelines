#!/bin/bash

# BASH script to download genomes of interest into the current working
# directory from iGenomes:
#
# The script also needs to be run on a version of the Linux in which the command 
# "rename" takes regular expression.  Some versions of Linux "rename" do not do this.

# Make folder
mkdir Genome_References
cd Genome_References/

####################
####################
#####################
# Ensembl references
#####################
#####################
####################
mkdir Ensembl
cd Ensembl/


####################
#####################
# HUMAN
#####################
#####################
mkdir Homo_sapiens
cd Homo_sapiens/

#####################
# GRCh38
#####################
mkdir GRCh38
cd GRCh38

# FASTA
################
mkdir FASTA
cd FASTA/
mkdir primary_assembly
cd primary_assembly
lftp -e "mget Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz; bye" http://ftp.ensembl.org/pub/release-102/fasta/homo_sapiens/dna/
gunzip *.gz

# Bowtie2
################
genome=$(ls -m *.fa)
genome=$(echo $genome | sed 's/ //g')
bowtie2-build $genome Homo_sapiens.GRCh38.dna.primary_assembly > bowtie2-build.out
mkdir ../../Bowtie2
mv *.bt2  ../../Bowtie2/
mv bowtie2-build.out  ../../Bowtie2/

# HISAT2
################
hisat2-build $genome Homo_sapiens.GRCh38.dna.primary_assembly > hisat2-build.out
mkdir ../../HISAT2
mv *.ht2  ../../HISAT2/
mv hisat2-build.out  ../../HISAT2/
cd ../..

# GTF
#################
mkdir GTF
cd GTF/
wget -nv http://ftp.ensembl.org/pub/release-102/gtf/homo_sapiens/Homo_sapiens.GRCh38.102.gtf.gz
gunzip *.gz
gtf=$(ls *.gtf)
splice_sites=$(echo $gtf | sed 's/\.gtf$/.ss/g')
hisat2_extract_splice_sites.py $gtf > $splice_sites
exon=$(echo $gtf | sed 's/\.gtf$/.exon/g')
hisat2_extract_exons.py $gtf > $exon
cd ..

# Make NextFlow genome text
###########################
echo -e "name\tGRCh38" > GRCh38.genome
echo -e "species\tHomo sapiens" >> GRCh38.genome
echo -e "fasta\t$PWD/FASTA/primary_assembly/" >> GRCh38.genome
echo -e "bowtie2\t$PWD/Bowtie2/Homo_sapiens.GRCh38.dna.primary_assembly" >> GRCh38.genome
echo -e "gtf\t$PWD/GTF/Homo_sapiens.GRCh38.102.gtf" >> GRCh38.genome
echo -e "hisat2_splices\t$PWD/GTF/Homo_sapiens.GRCh38.102.ss" >> GRCh38.genome
echo -e "hisat2\t$PWD/HISAT2/Homo_sapiens.GRCh38.dna.primary_assembly" >> GRCh38.genome
echo -e "\n" >> GRCh38.genome
cd ../..


####################
#####################
# Mouse
#####################
#####################
mkdir Mus_musculus
cd Mus_musculus/

#####################
# GRCm38
#####################
mkdir GRCm38
cd GRCm38

# FASTA
################
mkdir FASTA
cd FASTA/
mkdir primary_assembly
cd primary_assembly

lftp -e "mget Mus_musculus.GRCm38.dna.primary_assembly.fa.gz; bye" http://ftp.ensembl.org/pub/release-100/fasta/mus_musculus/dna/
gunzip *.gz

# Bowtie2
################
genome=$(ls -m *.fa)
genome=$(echo $genome | sed 's/ //g')
bowtie2-build $genome Mus_musculus.GRCm38.dna.primary_assembly > bowtie2-build.out
mkdir ../../Bowtie2
mv *.bt2  ../../Bowtie2/
mv bowtie2-build.out  ../../Bowtie2/

# HISAT2
################
hisat2-build $genome Mus_musculus.GRCm38.dna.primary_assembly > hisat2-build.out
mkdir ../../HISAT2
mv *.ht2  ../../HISAT2/
mv hisat2-build.out  ../../HISAT2/
cd ../..

# GTF
#################
mkdir GTF
cd GTF/
wget -nv http://ftp.ensembl.org/pub/release-100/gtf/mus_musculus/Mus_musculus.GRCm38.100.gtf.gz
gunzip *.gz
gtf=$(ls *.gtf)
splice_sites=$(echo $gtf | sed 's/\.gtf$/.ss/g')
hisat2_extract_splice_sites.py $gtf > $splice_sites
exon=$(echo $gtf | sed 's/\.gtf$/.exon/g')
hisat2_extract_exons.py $gtf > $exon
cd ..

# Make NextFlow genome text
###########################
echo -e "name\tGRCm38" > GRCm38.genome
echo -e "species\tMus musculus" >> GRCm38.genome
echo -e "fasta\t$PWD/FASTA/primary_assembly/" >> GRCm38.genome
echo -e "bowtie2\t$PWD/Bowtie2/Mus_musculus.GRCm38.dna.primary_assembly" >> GRCm38.genome
echo -e "gtf\t$PWD/GTF/Mus_musculus.GRCm38.100.gtf" >> GRCm38.genome
echo -e "hisat2_splices\t$PWD/GTF/Mus_musculus.GRCm38.100.ss" >> GRCm38.genome
echo -e "hisat2\t$PWD/HISAT2/Mus_musculus.GRCm38.dna.primary_assembly" >> GRCm38.genome
echo -e "\n" >> GRCm38.genome
cd ..






echo "Done"
