# Snakefile to launch and combine quality checks

rule prepare_files:
    # Extract proteins from helixer for input to BUSCO and Egapx
    input:
        fasta=results_dir + "/hap{n}.checked.fasta",
        helixer=results_dir + "/helixer/hap_{n}.gff3"
    output:
        helixer_protein=results_dir + "/helixer/hap_{n}.helixer.protein.faa",
        helixer_cds=results_dir + "/helixer/hap_{n}.helixer.CDS.fna"
    params:
        resdir=results_dir
    shell:
        """
        cd {params.resdir}/helixer/
        module load bioinfo/AGAT/1.2.0
        agat_sp_extract_sequences.pl --gff {input.helixer} -f {input.fasta} -p -o {output.helixer_protein}
        module load bioinfo/gffread/0.12.6
        gffread -x {output.helixer_cds} -g {input.fasta} {input.helixer}
        """

rule compleasm_raw:
    input:
        protein_fasta=lambda wildcards: 
            results_dir + "/{dataset}/" + 
            ("hap_{n}/output/complete.proteins.faa" if wildcards.dataset == "egapx" else "hap_{n}.helixer.protein.faa")
    output:
        summary=results_dir + "/{dataset}/COMPL_hap{n}_summary.txt"
    params:
        lineage=get_busco_lineage,
        resdir=results_dir
    threads: 8
    resources:
        mem_mb=12000
    shell:
        """
        cd {params.resdir}/{wildcards.dataset}
        module load devel/python/Python-3.11.1 bioinfo/compleasm/0.2.5
        compleasm.py download {params.lineage}
        compleasm.py protein -a {input.protein_fasta} -o {output.summary} -l {params.lineage} -L /mb_downloads/{params.lineage}_odb10 --odb 10 -t {threads}
        """

rule psauron_raw:
    input:
        cds_fasta=lambda wildcards: 
            results_dir + "/{dataset}/" + 
            ("hap_{n}/output/complete.cds.fna" if wildcards.dataset == "egapx" else "hap_{n}.helixer.CDS.fna")
    output:
        summary=results_dir + "/{dataset}/PSAURON_hap{n}_summary.csv"
    resources:
        mem_mb=10000
    params:
        resdir=results_dir
    shell:
        """
        cd {params.resdir}/{wildcards.dataset}
        module load bioinfo/PSAURON/1.0.2
        psauron -i {input.cds_fasta} -c -o {output.summary}
        """    

rule combine_psauron_scores_raw:
    input:
        lambda wildcards: expand(
            results_dir + "/{dataset}/PSAURON_hap{n}_summary.csv",
            species=SPEC,
            dataset=["egapx", "helixer"],
            n=HAPS_DICT[wildcards.species]
        )
    output:
        results_dir + "/{dataset}/combined_psauron_scores.csv"
    run:
        import re
        scores = []
        for f in input:
            match = re.search(r'_hap(\d+)_', os.path.basename(f))
            hap = f"hap{match.group(1)}" if match else "unknown"
            with open(f, "r") as infile:
                for line in infile:
                    if line.startswith("psauron score:"):
                        score = line.strip().split(": ")[1]
                        scores.append(f"{hap},{score}")
                        break 
        with open(output[0], "w") as out:
            out.write("haplotype,score\n")
            out.write("\n".join(scores) + "\n")

rule omark_raw:
    input:
        protein_fasta=lambda wildcards: 
            results_dir + "/{dataset}/" + 
            ("hap_{n}/output/complete.proteins.faa" if wildcards.dataset == "egapx" else "hap_{n}.helixer.protein.faa")
    output:
        plot=results_dir + "/{dataset}/OMArk/hap{n}/hap{n}.png"
    params:
        db=results_dir + "/{dataset}/OMArk/hap{n}"
    resources:
        mem_mb=10000
    shell:
        """
        module load devel/Miniconda/Miniconda3 bioinfo/OMArk/0.3.0
        mkdir {params.db}
        omamer search --db /usr/local/bioinfo/src/OMArk/example_on_cluster/LUCA.h5 --query {input.protein_fasta} --out {params.db}
        omark -f {params.db}/search.omamer -d /usr/local/bioinfo/src/OMArk/example_on_cluster/LUCA.h -o {params.db}
        plot_all_results.py -i {params.db} -o {output.plot}
        """



