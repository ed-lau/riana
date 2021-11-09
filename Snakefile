configfile: "config.yaml"

rule all:
    input: expand("out/snakemake/{timepoint}_riana.txt", timepoint=config["data"])

rule copy_files:
    input:
        mzml = lambda wildcards: config["data"][wildcards.timepoint]
    output:
        # dir="out/snakemake/{timepoint}",
        linked_mzml="out/snakemake/{timepoint}/mzml/input.mzml.gz"
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

rule riana:
    input:
        pin="out/snakemake/{timepoint}/percolator/percolator.psms.txt",
        mzml="out/snakemake/{timepoint}/mzml/"
    output:
        riana="out/snakemake/{timepoint}_riana.txt"
    threads: config["threads"]["riana"]
    shell:
        "python -m riana integrate out/snakemake/{wildcards.timepoint}/mzml "
        "out/snakemake/{wildcards.timepoint}/percolator/percolator.psms.txt"
        " -i 0,1,2,3,4,5 -q 0.01 -r 0.5 -m 25 -o {output} -s {wildcards.timepoint}"
        " -t {threads}"

