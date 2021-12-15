"""
Snakefile for multiple fractions

Example command:
snakemake -c -s Snakefile -d out/snakemake_test --configfile config_multi.yaml


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
        pin="{timepoint}/comet.pin"
    log: "{timepoint}/comet.log"
    threads: config["threads"]["comet"]
    benchmark: "{timepoint}/comet.benchmark.txt"
    params:
        comet=config["paths"]["comet"],
        comet_params=config["paths"]["comet_params"],
        fasta=config["paths"]["fasta"]
    shell: #
        """
        {params.comet} -P{params.comet_params} -D{params.fasta} {input}/*.mzML.gz 1>> {log}
        tail -n +1 -q {input.mzml}/*.pin >> {output.pin}.temp
        head -1 {output.pin}.temp > {output.pin}.temp2
        tail -n +2 -q {input.mzml}/*.pin >> {output.pin}.temp2
        mv {output.pin}.temp2 {output.pin}
        rm {input.mzml}/*.pin 
        """
        # Comet does not allow output files to be specified and simply creates a pin or pepxml file of the same name
        # in the same directory. The standalone Percolator only takes in one pin file unlike the Crux distribution.
        # So the 4-line monstrosity above concatenates all the pin files while keeping only the first line
        # header in one of them, without using an open bracket which for some reason does not work with snakemake, e.g.,
        # head -1 {input.mzml}/*.pin([1]) >> {output.pin}.temp
        # tail -n +2 -q {input.mzml}/*.pin >> {output.pin}.temp
        # Note the {output.pin}.temp is done because I am not sure if snakemake will jump the gun and
        # start the percolator step as soon as the header is copied.

rule percolator:
    input:
        comet_pin= "{timepoint}/comet.pin"
    log: "{timepoint}/percolator.log"
    output:
        psms="{timepoint}/percolator/percolator.psms.txt"
    benchmark: "{timepoint}/percolator.benchmark.txt"
    params:
        percolator=config["paths"]["percolator"],
        fasta=config["paths"]["fasta"]
    shell:
        "percolator -Y -i 20 -P DECOY_ -f {params.fasta} {input.comet_pin} -m {output.psms} 2>> {log}"
    # 2021-11-11: might have to send stdout to log as well to stop it from spamming the shell

rule riana_integrate:
    input:
        pin="{timepoint}/percolator/percolator.psms.txt",
        mzml= lambda wildcards: config["data"][wildcards.timepoint]
    output:
        riana="{timepoint}_riana.txt"
    params:
        iso=config["params"]["isotopomers"]
    threads: config["threads"]["riana"]
    shell:
        "riana integrate {input.mzml} "
        "{input.pin} "
        "-i {params.iso} -q 0.01 -r 0.5 -m 25 -o {output} -s {wildcards.timepoint} "
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