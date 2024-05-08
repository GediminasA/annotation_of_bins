wdir = config["work_dir"]
bins = glob_wildcards(config["bin_dir"]+"/{bin}.fasta").bin

rule target:
    input:
        wdir +"/merged_proteins.fasta"
        # expand(
        #     wdir + "/prokka_annotations/{bin}",
        #     bin = bins
        # )

def get_protein_reference(wildcards):
    import os
    refd = config["references"]
    if wildcards.bin in refd.keys():
        out = os.getcwd()+"/"+refd[wildcards.bin]
        out = "--proteins " + out
    else:
        out = ""
    return(out)


rule prokka:
    input:
        fasta = config["bin_dir"]+"/{bin}.fasta"
    output:
        od = directory(wdir + "/prokka_annotations/{bin}"),
        prots = wdir + "/prokka_annotations/{bin}/{bin}.faa"
    conda:
        "prokka"
    threads: 96
    params:
        add = get_protein_reference
    shell:
        """
        prokka {params.add} --cpu  {threads} --prefix {wildcards.bin} --locustag {wildcards.bin} --addgenes  --force --outdir {output.od} {input}
        """
rule merge_fastas:
    input:
        expand(
            wdir + "/prokka_annotations/{bin}/{bin}.faa",
            bin = bins
        )
    output:
        wdir +"/merged_proteins.fasta"
    shell:
        """
        cat {input} > {output}
        """


        

    