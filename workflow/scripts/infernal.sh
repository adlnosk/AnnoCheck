#!/bin/bash

fasta=$1
gff_out=$2
wd=$3
n=$4
rfam_db=$5
cpus=$6


cd $wd

ln -s $fasta ./hap${n}.fasta


module load bioinfo/Infernal/1.1.4

# Download Rfam library if not already present
if [ ! -f $rfam_db/Rfam.cm ]; then
    mkdir -p $rfam_db; cd $rfam_db
    wget ftp://ftp.ebi.ac.uk/pub/databases/Rfam/CURRENT/Rfam.cm.gz
    gunzip Rfam.cm.gz
    wget ftp://ftp.ebi.ac.uk/pub/databases/Rfam/CURRENT/Rfam.clanin
fi


# Get total database size (multiply by 2 for both strands)
size=$(esl-seqstat ${fasta} | grep "Total # residues:" | awk '{print $4}' | awk '{print $1 * 2 / 1000000}')

# Run cmscan to annotate RNAs represented in Rfam
cmscan --cpu ${cpus} -Z $size --cut_ga --rfam --nohmmonly --tblout hap_${n}.tblout --fmt 2 --clanin $rfam_db/Rfam.clanin $rfam_db/Rfam.cm ${fasta} > hap_${n}.cmscan

# Remove overlapping hits with better scores
# !!!! if wish to remove tRNAs, then add in the pipe : grep -v RF00005 | grep -v RF01852 
cat hap_${n}.tblout | grep -v " = " > hap_${n}.deoverlapped.tblout

# Create GFF file from tblout
perl $wd/../../workflow/scripts/tblout2gff.pl --fmt2 --cmscan hap_${n}.deoverlapped.tblout > ${gff_out}

