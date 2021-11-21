configfile: "config_muscle.yaml"

rule all:
    input: "out/snakemake/riana_fit_peptides.csv"

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

rule percolator:
    input:
        comet_pin= "out/snakemake/{timepoint}/comet.pin"
    log: "out/snakemake/{timepoint}/percolator.log"
    output:
        psms="out/snakemake/{timepoint}/percolator/percolator.psms.txt"
    benchmark: "out/snakemake/{timepoint}/percolator.benchmark.txt"
    params:
        percolator=config["paths"]["percolator"],
        fasta=config["paths"]["fasta"]
    shell:
        "percolator -Y -i 20 -P DECOY_ -f {params.fasta} {input.comet_pin} -m {output.psms} 2>> {log}"
    # 2021-11-11: might have to send stdout to log as well to stop it from spamming the shell

rule riana_integrate:
    input:
        pin="out/snakemake/{timepoint}/percolator/percolator.psms.txt",
        mzml="out/snakemake/{timepoint}/mzml/input.mzml.gz"
    output:
        riana="out/snakemake/{timepoint}_riana.txt"
    threads: config["threads"]["riana"]
    shell:
        "python -m riana integrate out/snakemake/{wildcards.timepoint}/mzml "
        "out/snakemake/{wildcards.timepoint}/percolator/percolator.psms.txt"
        " -i 0,1,2,3,4,5 -q 0.01 -r 0.5 -m 25 -o {output} -s {wildcards.timepoint}"
        " -t {threads}"

rule riana_fit:
    input:
        integrated=expand("out/snakemake/{timepoint}_riana.txt", timepoint=config["data"])
    output:
        riana="out/snakemake/riana_fit_peptides.csv"
    threads: config["threads"]["fitcurve"]
    shell:
        "python -m riana fit {input.integrated} "
        "-q 0.01 -d 12 -o out/snakemake "
        "-t {threads} -r 0.045"