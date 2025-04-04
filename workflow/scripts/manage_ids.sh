#!/bin/bash

# Arguments
inf2=$1
dict=$2
inf_out=$3
wd=$4
n=$5


inf3_gff="inf3.gff"
inf3g_gff="inf3g.gff"
inf4_gff3="inf4.gff3"

# Create header for GFF3 file
echo '##gff-version 3' > $inf3_gff

# Create feature-specific lists based on the dictionary file
cut -f 1 $dict | sort | uniq | while read feat; do
    grep "^${feat}" $dict | cut -f 2 > ${feat}rna_list.txt
done

# Process the features from the dictionary file and create the GFF entries
while read feat; do
    # Create the list of feature names for this feature type (e.g., mir-156, mir-160)
    feature_list="${feat}rna_list.txt"
    feature_type="${feat}"

    # Process ncRNA and gene features
    less $inf2 | grep -f $feature_list | awk -F "\t" \
    -v feat=$feature_type '{print $1,$2,"ncRNA",$4,$5,$6,$7,$8,"ID="feat"_"$1"_"$4"_"$5";Parent=gene_"$1"_"$4"_"$5";"$9";rfam="$3}' >> $inf3_gff

    less $inf2 | grep -f $feature_list | awk -F "\t" \
    -v feat=$feature_type '{print $1,$2,"gene",$4,$5,$6,$7,$8,"ID=gene_"$1"_"$4"_"$5";"$9";rfam="$3}' >> $inf3g_gff
done < <(cut -f 1 $dict | sort | uniq)

# Merge and sort the GFF entries
cat $inf3_gff $inf3g_gff | awk '$1 ~ /^#/ {print $0;next} {print $0 | "sort -k1,1 -k4,4n -k5,5n"}' | sed 's/ /\t/g' | awk '!seen[$9]++' > $inf4_gff3

# Check the validity of the GFF3 file
module load bioinfo/GenomeTools/1.6.5
gt gff3 -typecheck -sortlines yes -checkids yes -retainids yes -fixregionboundaries yes -tidy yes $inf4_gff3 > $inf_out

# Clean up temporary files
rm *rna_list*
rm $inf3_gff $inf3g_gff $inf4_gff3

