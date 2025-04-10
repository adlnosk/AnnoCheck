rule prepare_files:
    input:
        fasta = results_dir + "/hap{n}.checked.fasta",
        helixer = results_dir + "/helixer/hap_{n}.gff3"
    output:
        helixer_protein = results_dir + "/helixer/hap_{n}.helixer.protein.faa",
        helixer_cds = results_dir + "/helixer/hap_{n}.helixer.CDS.fna"
    resources:
        mem_mb=10000
    params:
        resdir = results_dir
    shell:
        """
        cd {params.resdir}/helixer/

        module load devel/Miniconda/Miniconda3 bioinfo/AGAT/1.2.0
        agat_sp_extract_sequences.pl --gff {input.helixer} -f {input.fasta} -p -o {output.helixer_protein}

        module load bioinfo/gffread/0.12.6
        gffread -x {output.helixer_cds} -g {input.fasta} {input.helixer}
        """


rule summarize_star_logs:
    input:
        trace = results_dir + "/egapx/hap_{n}/output/run.trace.txt"
    output:
        summary = results_dir + "/egapx/hap_{n}/STAR_logs/rna_summary.txt"
    params:
        resdir = results_dir
    threads: 1
    shell:
        """
        stardir="{params.resdir}/egapx/hap_{wildcards.n}/STAR_logs"
        workdir="{params.resdir}/egapx/hap_{wildcards.n}/work"
        mkdir -p "$stardir"
        grep "egapx:rnaseq_short_plane:star:run_star" {input.trace} | cut -f2 > "$stardir/wd.list"
        while read id; do
            grep -rl "Log.final.out" "$workdir/${{id}}"* | grep "Log.final.out" | xargs -I{{}} cp {{}} "$stardir"
        done < "$stardir/wd.list"
        declare -A sum
        fields=(
            "Number of input reads"
            "Uniquely mapped reads number"
            "Number of reads mapped to multiple loci"
            "Number of reads unmapped: too many mismatches"
            "Number of reads unmapped: other"
            "Number of chimeric reads"
        )
        for f in "$stardir"/*Log.final.out; do
            awk -F'|' '{{ gsub(/^ +| +$/, "", $1); gsub(/^ +| +$/, "", $2); print $1, $2 }}' "$f" |
            while read key val; do
                for fkey in "${{fields[@]}}"; do
                    [[ "$key" == "$fkey" ]] && ((sum["$key"] += val))
                done
            done
        done
        {{
            echo "RNA Summary (hap {wildcards.n}):"
            for k in "${{fields[@]}}"; do
                printf "%-50s %'15d\\n" "$k" "${{sum[$k]}}"
            done
        }} > {output.summary}
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
        compleasm.py protein -p {input.protein_fasta} -o {params.resdir}/{wildcards.dataset}/COMPL_hap{wildcards.n} -l {params.lineage} -L mb_downloads/ -t {threads}
        ln -s {params.resdir}/{wildcards.dataset}/COMPL_hap{wildcards.n}/summary.txt {output.summary} 
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
        mem_mb=50000
    shell:
        """
        module load devel/Miniconda/Miniconda3
        module load bioinfo/OMArk/0.3.0
        mkdir -p {params.db}
        cd {params.db}
        omamer search --db /usr/local/bioinfo/src/OMArk/example_on_cluster/LUCA.h5 --query {input.protein_fasta} --out {params.db}/search.omamer
        omark -f {params.db}/search.omamer -d /usr/local/bioinfo/src/OMArk/example_on_cluster/LUCA.h5 -o {params.db} -r family
        plot_all_results.py -i {params.db} -o {output.plot}
        """

rule multiqc_raw:
    input:
        star_logs = results_dir + "/egapx/hap_{n}/STAR_logs/",
        compl_summary = results_dir + "/{dataset}/COMPL_hap{n}_summary.txt",
        psa_summary = results_dir + "/{dataset}/PSAURON_hap{n}_summary.csv",
        omark_summary = results_dir + "/{dataset}/OMArk/hap{n}/
    output:
        multiqc = results_dir + "/MultiQC_raw/multiqc_report.html",
        star = results_dir + "/MultiQC_raw/STAR_logs/multiqc_report.html"
    params:
        resdir = results_dir,
        multiqc_inputs = directory(results_dir + "/MultiQC_raw/")
    threads: 1
    shell:
        """
        mkdir -p {params.multiqc_inputs}
        cd {params.multiqc_inputs}
        cp -r {input.star_logs} {output.multiqc_inputs}/STAR_logs
        module load bioinfo/MultiQC/1.27.1
        multiqc {output.multiqc_inputs}/STAR_logs
        cp {input.compl_summary} {output.multiqc_inputs}/
        cp {input.psa_summary} {output.multiqc_inputs}/
        cp {input.omark_summary} {output.multiqc_inputs}/
        multiqc .
        """


