wdir = config["work_dir"]
bins = glob_wildcards(config["bin_dir"]+"/{bin}.fasta").bin

#get rules
include: "transcription.smk"

### Preprocessing general

rule target:
    input:
        #wdir +"/merged_proteins.fasta"
        #wdir + "/proteins_splits"
        wdir + "/go_annots.tsv",
        wdir + "/de_annots.tsv",
        wdir + "/metagenome.fasta",

        # expand(
        #     wdir + "/prokka_annotations/{bin}/{bin}__sorted.gff",
        #     #wdir + "/prokka_annotations/{bin}",
        #     bin = bins
        # )

def get_protein_reference(wildcards):
    import os
    refd = config["references"]
    if wildcards.bin in refd.keys():
        out = os.getcwd()+"/"+refd[wildcards.bin]
        out = "--proteins " + out
        #TODO
    else:
        out = ""
    return(out)


rule prokka:
    input:
        fasta = config["bin_dir"]+"/{bin}.fasta"
    output:
        od = directory(wdir + "/prokka_annotations/{bin}"),
        prots = wdir + "/prokka_annotations/{bin}/{bin}.faa",
        gff = wdir + "/prokka_annotations/{bin}/{bin}.gff"
    conda:
        "prokka"
    threads: 96
    params:
        add = get_protein_reference
    shell:
        """
        prokka {params.add} --cpu  {threads} --prefix {wildcards.bin} --locustag {wildcards.bin} --addgenes  --force --outdir {output.od} {input}
        """


rule sortgff:
    input:
        "{stem}.gff"
    output:
        "{stem}__sorted.gff"
    conda:
        "genometools"
    shell:
        "gt gff3  -sort -retainids -tidy {input} > {output}"

rule get_gtf:
    input:
         wdir + "/metagenome.gff"
    output:
         wdir + "/metagenome.gtf"
    conda:
        "genometools"
    shell:
        """
        gt gff3_to_gtf {input} > {output}
        """

rule mergesortedgff:
    input:
        expand(
            wdir + "/prokka_annotations/{bin}/{bin}__sorted.gff",
            #wdir + "/prokka_annotations/{bin}",
            bin = bins
        )
    output:
         wdir + "/metagenome.gff"
    conda:
        "genometools"
    shell:
        """
        gt merge {input} > {output}
        """

rule get_contig_ids:
    input:
         wdir + "/metagenome.gff"
    output:
         wdir + "/metagenome.contigids"
    shell:
        """
        grep sequence-region {input}| cut -f4 -d " " > {output}
        """

rule get_contig_matchinggff:
    input:
        cids = wdir + "/metagenome.contigids",
        fastas = expand(
            config["bin_dir"]+"/{bin}.fasta",
            bin = bins
        )
    output:
         inter = temp(wdir + "/metagenome_pre.fasta"),
         final = wdir + "/metagenome.fasta"
    shell:
        """
        cat {input.fastas} | seqkit seq -w 0 -i | seqkit grep -w 0 -n -f {input.cids} > {output.inter}
        seqkit faidx  {output.inter} -l {input.cids} -w 0 -o {output.final}
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

checkpoint split:
    input:
        wdir +"/merged_proteins.fasta"
    output:
        directory(wdir + "/proteins_splits")
    shell:
        """
        seqkit split {input} -p 1000 --by-part-prefix split_  -O {output} 
        """
        
rule annotate:
    input:
        wdir + "/proteins_splits/split_{part}.fasta"
    output: 
        de = wdir + "/annots/split_{part}_DE.out",
        go = wdir + "/annots/split_{part}_GO.out",
        anno = wdir + "/annots/split_{part}_anno.out",
    log: 
        wdir + "/annots/split_{part}_anno.log",
    params:
        odir = wdir + "/annots",
        pannzer_bin = config["pannzer_bin"]
    threads: 8
    conda: "pannzer"
    shell:
        """
        python {params.pannzer_bin} -R -o ",{output.de},{output.go},{output.anno}" -i {input} &> {log} 
        """


def get_gos(wildcards):
    checkpoint_output = checkpoints.split.get(**wildcards).output[0]
    ids = glob_wildcards( checkpoint_output + "/split_{i}.fasta").i
    # print(checkpoint_output)
    # import sys
    # sys.exit()
    out = expand( wdir + "/annots/split_{part}_GO.out",part = ids)
    return(out)

def get_des(wildcards):
    checkpoint_output = checkpoints.split.get(**wildcards).output[0]
    ids = glob_wildcards( checkpoint_output + "/split_{i}.fasta").i
    # print(checkpoint_output)
    # import sys
    # sys.exit()
    out = expand( wdir + "/annots/split_{part}_DE.out",part = ids)
    return(out)


rule aggregate_gos:
    input:
        get_gos
    output:
        wdir + "/go_annots.tsv"
    shell:
        """
            echo "qpid    ontology    goid    desc    ARGOT_score ARGOT_PPV   ARGOT_rank  goclasscount" > {output}
            cat {input} | grep -v "ontology" >> {output}
        """        

rule get_prokka_defs:
    input:
        wdir +"/merged_proteins.fasta"
    output:
        wdir +"/merged_proteins.fasta_headers"
    shell:
        """
        grep ">" {input} > {output}
        """


rule aggregate_des:
    input:
        pro =  wdir +"/merged_proteins.fasta_headers",
        des = get_des
    output:
        wdir + "/de_annots.tsv"
    notebook:
        "notebooks/merge_des.r.ipynb"