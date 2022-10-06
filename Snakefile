"""
Snakefile for multiple fractions

Example command:
snakemake -c -s Snakefile -d out/snakemake_test --configfile config_ipsc_o18.yaml

Single fractions:
rule copy_files:
    input:
        mzml = lambda wildcards: config["data"][wildcards.timepoint]
    output:
        # dir="out/snakemake/{timepoint}",
        linked_mzml=temp("out/snakemake/{timepoint}/mzml/input.mzml.gz")
    shell:
        "cp {input.mzml} {output.linked_mzml}" # should change to symlink instead probably

rule comet:
    input: "out/snakemake/{timepoint}/mzml/input.mzml.gz" # lambda wildcards: config["data"][wildcards.timepoint]
    output:
        pin="out/snakemake/{timepoint}/comet.pin"
    log: "out/snakemake/{timepoint}/comet.log"
    threads: config["threads"]["comet"]
    benchmark: "out/snakemake/{timepoint}/comet.benchmark.txt"
    params:
        comet=config["paths"]["comet"],
        comet_params=config["paths"]["comet_params"],
        fasta=config["paths"]["fasta"]
    shell:
        "{params.comet} -P{params.comet_params} -D{params.fasta} {input} 1>> {log}; "
        "outname=$(basename {input} | cut -f 1 -d '.'); "
        "mv $(dirname {input})/${{outname}}.pin {output.pin}"
"""

rule all:
    input: "riana_fit_peptides.csv"

rule comet:
    input: mzml = lambda wildcards: config["data"][wildcards.timepoint]
    output:
        psms="{timepoint}/percolator/percolator.target.psms.txt"
    log: "{timepoint}/comet.log"
    threads: config["threads"]["comet"]
    benchmark: "{timepoint}/comet.benchmark.txt"
    params:
        comet=config["paths"]["comet"],
        percolator=config["paths"]["percolator"],
        comet_params=config["paths"]["comet_params"],
        fasta=config["paths"]["fasta"]
    shell: #
        """
        {params.comet} -P{params.comet_params} -D{params.fasta} {input}/*.mzML.gz 1>> {log}
        {params.percolator} --protein T {input}/*.pep.xml --decoy-prefix DECOY_ \\
         --overwrite T --output-dir {wildcards.timepoint}/percolator \\
         --spectral-counting-fdr 0.01 --maxiter 15 1>> {log} 2>>{log}
        """


rule riana_integrate:
    input:
        pin="{timepoint}/percolator/percolator.target.psms.txt",
        mzml= lambda wildcards: config["data"][wildcards.timepoint]
    output:
        riana="{timepoint}_riana.txt"
    params:
        iso=config["params"]["isotopomers"],
        mass_tol=config["params"]["mass_tol"]
    threads: config["threads"]["riana"]
    shell:
        "riana integrate {input.mzml} "
        "{input.pin} "
        "-i {params.iso} -q 0.01 -r 0.33 -m {params.mass_tol} -o {output} -s {wildcards.timepoint} "
        "-t {threads}"

rule riana_fit:
    input:
        integrated=expand("{timepoint}_riana.txt", timepoint=config["data"])
    output:
        riana="riana_fit_peptides.csv"
    params:
        ria=config["params"]["ria_max"],
        kp=config["params"]["kp"],
        kr=config["params"]["kr"],
        rp=config["params"]["rp"],
        depth=config["params"]["depth"],
        label_type=config["params"]["label_type"],
        model=config["params"]["model"]
    threads: config["threads"]["fitcurve"]
    shell:
        "riana fit {input.integrated} "
        "-q 0.01 -d {params.depth} -o . -m {params.model} --kp {params.kp} "
        "--kr {params.kr} --rp {params.rp} "
        "-t {threads} -r {params.ria} -l {params.label_type}"