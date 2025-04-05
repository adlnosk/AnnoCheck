rule fasta_validation:
    input:
        fasta=FASTA_INPUT
    output:
        checked_fasta=results_dir + "/hap{n}.checked.fasta"
    threads: 5
    resources:
        time="00:30:00",
        mem_mb=10000
    params:
        resdir=results_dir
    shell:
        """
        cd {params.resdir}
      #fasta_validate - better to run manually, sometimes it cracs  
        #module load bioinfo/fasta_validator/a7cbc40;
        #fasta_validate {input.fasta}
        #EXIT_CODE=$?
        #if [ $EXIT_CODE -ne 0 ]; then
        #    echo "Error: FASTA validation failed with code $EXIT_CODE" >&2
        #    exit $EXIT_CODE
        #fi
        if [[ "{input.fasta}" == *.gz ]]; then
            module load bioinfo/bgzip/1.18
            bgzip -d -@ {threads} --stdout {input.fasta} | sed '/^>/ s/$/ {wildcards.species}/' > {output.checked_fasta}
        else
            cat {input.fasta} | sed '/^>/ s/$/ {wildcards.species}/' > {output.checked_fasta}
        fi
        """

rule helixer:
    input:
        fasta=rules.fasta_validation.output
    output:
        results_dir + "/helixer/hap_{n}.gff3"
    resources:
        time="96:00:00",
        mem_mb=50000
    params:
        resdir=results_dir,
        helixer_lineage=get_helixer_lineage
    shell:
        """
        cd {params.resdir}/helixer
        {PWD}/scripts/helixer.sh {input.fasta} {wildcards.n} {params.helixer_lineage}
        """

rule create_input_yaml:
    input:
        fasta=rules.fasta_validation.output
    output:
        results_dir + "/egapx/hap_{n}/input.yaml"
    params:
        taxid=get_egapx_taxid,
        reads=get_reads,
        resdir=results_dir
    shell:
        """
        cd {params.resdir}/egapx/hap_{wildcards.n}
        echo "genome: {input.fasta}" > {output}
        echo "taxid: {params.taxid}" >> {output}
        echo "reads:" >> {output}

        for read in {params.reads}; do
            echo "  - \\"$read\\"" >> {output}
        done
        """

rule egapx:
    input:
        fasta=rules.fasta_validation.output,
        yaml=rules.create_input_yaml.output
    output:
        results_dir + "/egapx/hap_{n}/output/complete.genomic.gtf",
        results_dir + "/egapx/hap_{n}/output/complete.proteins.faa",
        results_dir + "/egapx/hap_{n}/output/complete.transcripts.fna",
        results_dir + "/egapx/hap_{n}/output/complete.cds.fna"
    params:
        resdir=results_dir,
        config=PWD + "/resources/"
    shell:
        """
        cd {params.resdir}/egapx/hap_{wildcards.n};
        module load devel/python/Python-3.11.1 containers/singularity/3.9.9 bioinfo/Nextflow/24.10.0 devel/java/17.0.6 bioinfo/EGAPx/0.3.2-alpha
        egapx.py {input.yaml} -o output -e slurm --config-dir {params.config}
        """

rule repeat_modeler:
    input:
        fasta=results_dir + "/hap1.checked.fasta"
    output:
        families=results_dir + "/repeatmodeler/{species}-families.fa"
    threads: 20
    resources:
        time="96:00:00",
        mem_mb=50000
    params:
        resdir=results_dir
    shell:
        """
        cd {params.resdir}/repeatmodeler
        module load bioinfo/RepeatModeler/2.0.6
        BuildDatabase -name {wildcards.species} {input.fasta}
        RepeatModeler -database {wildcards.species} -threads {threads} -LTRStruct -quick
        """

rule repeat_masker:
    input:
        fasta=rules.fasta_validation.output,
        families=results_dir + "/repeatmodeler/{species}-families.fa"
    output:
        masked_fasta=results_dir + "/repeatmasker/hap{n}/hap{n}.checked.masked.fasta"
    threads: 16
    resources:
        time="24:00:00",
        mem_mb=50000
    params:
        resdir=results_dir
    shell:
        """
        mkdir -p {params.resdir}/repeatmasker/hap{wildcards.n}
        cd {params.resdir}/repeatmasker
        module load devel/python/Python-3.6.3 bioinfo/RepeatMasker/4.1.8
        RepeatMasker -xsmall -gff -pa 4 -q -e ncbi -dir hap{wildcards.n} -norna -lib {input.families} {input.fasta}
        """

rule infernal_rfam:
    input:
        fasta=rules.fasta_validation.output
    output:
        gff_out=results_dir + "/infernal_rfam/hap_{n}.deoverlapped.gff"
    threads: 16
    resources:
        time="96:00:00",
        mem_mb=50000
    params:
        rfam_db=PWD + "/resources/Rfam_lib/",
        resdir=results_dir
    shell:
        """
        cd {params.resdir}/infernal_rfam
        {PWD}/scripts/infernal.sh {input.fasta} {output.gff_out} {params.resdir}/infernal_rfam/ {wildcards.n} {params.rfam_db} {threads}
        """

rule manage_rfam_ids:
    input:
        inf = rules.infernal_rfam.output.gff_out,
        dict = PWD + "/resources/ncRNA_dictionary.txt"
    output:
        gff_out = results_dir + "/infernal_rfam/hap_{n}.deoverlapped.checked.gff"
    params:
        resdir=results_dir
    shell:
        """
        cd {params.resdir}
        {PWD}/scripts/manage_ids.sh {input.inf} {input.dict} {output.gff_out} {params.resdir}/infernal_rfam/ {wildcards.n}
        """

rule trnascan:
    input:
        fasta=rules.fasta_validation.output
    output:
        gff_out=results_dir + "/tRNAscan/hap_{n}.gff3"
    threads: 4
    params:
        prefix="hap{wildcards.n}",
        resdir=results_dir
    resources:
        time="24:00:00",
        mem_mb=10000
    shell:
        """
        cd {params.resdir}/tRNAscan/
        module load bioinfo/tRNAscan-SE/2.0.12
        tRNAscan-SE --gff {output.gff_out} --thread {threads} -o {params.prefix} -m {params.prefix}.stats -p {params.prefix} -d -Q -- {input.fasta}
        """

rule rnammer:
    input:
        fasta=rules.fasta_validation.output
    output:
        gff_out=results_dir + "/RNAmmer/hap_{n}.gff3"
    threads: 4
    resources:
        time="24:00:00",
        mem_mb=10000
    params:
        prefix="hap{wildcards.n}",
        resdir=results_dir
    shell:
        """
        cd {params.resdir}/RNAmmer/
        perl -MCPAN -e 'install XML::Simple'
        module load bioinfo/RNAmmer/1.2
        rnammer -gff {params.prefix}.gff -S euk -m lsu,ssu,tsu -f {params.prefix}.fa -h {params.prefix}.hmm -multi < {input.fasta}
        less {params.prefix}.gff | sed 's/[[:space:]]*$//' | sort -k1,1 -k4,4n -k5,5n > temp
        echo '##gff-version 3' > temp2
        less temp | awk -F "\t" '{{print $1,$2,"rRNA",$4,$5,$6,$7,$8,"ID=rRNA_"$1"_"$4"_"$5";Parent=gene_"$1"_"$4"_"$5";type="$9}}' OFS='\t' >> temp2
        less temp | awk -F "\t" '{{print $1,$2,"gene",$4,$5,$6,$7,$8,"ID=gene_"$1"_"$4"_"$5";type="$3}}' OFS='\t' >> temp2
        module load bioinfo/GenomeTools/1.6.5
        gt gff3 -typecheck -sortlines yes -checkids yes -retainids yes -fixregionboundaries yes -tidy yes temp  > {output.gff_out}
        rm temp*
        """
