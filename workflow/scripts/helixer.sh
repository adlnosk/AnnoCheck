#!/bin/bash

fasta=$1
hap=$2
wd=$PWD
echo $wd
lineage=$3

module load containers/singularity/3.9.9
module load devel/Miniconda/Miniconda3 bioinfo/Helixer/0.3.1
module load bioinfo/HelixerPost/4d4799b


if [ ! -f "$wd/../../../workflow/resources/helixer-docker_helixer_v0.3.3_cuda_11.8.0-cudnn8.sif" ]; then
    echo "Singularity image not found. Pulling the image..."
    singularity pull docker://gglyptodon/helixer-docker:helixer_v0.3.3_cuda_11.8.0-cudnn8
    mv helixer-docker_helixer_v0.3.3_cuda_11.8.0-cudnn8.sif $wd/../../workflow/resources/helixer-docker_helixer_v0.3.3_cuda_11.8.0-cudnn8.sif
else
    echo "Singularity image already exists."
fi


ln -s $wd/../../workflow/resources/helixer-docker_helixer_v0.3.3_cuda_11.8.0-cudnn8.sif helixer-docker_helixer_v0.3.3_cuda_11.8.0-cudnn8.sif
ln -s $fasta ./hap_$hap.fasta.gz

singularity run --nv helixer-docker_helixer_v0.3.3_cuda_11.8.0-cudnn8.sif Helixer.py --fasta-path hap_$hap.fasta.gz --lineage land_plant --gff-output-path hap_$hap.gff3

