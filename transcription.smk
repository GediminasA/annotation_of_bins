fids = glob_wildcards(config["sequencing_data"]+"/{fid}_1_sequence.fastq.gz").fid


rule target2:
    input:
        wdir + "/counts.tsv"
        # expand(
        #     #
        #     #wdir + "reads_trimmed/{fid}_R1.fast1.gz",
        #     #wdir + "/alignments/{fid}.bam",
        #     wdir + "/barcodes/{fid}.txt",
        #     fid = fids
        # )

rule trim:
    input:
        r1 = config["sequencing_data"]+"/{fid}_1_sequence.fastq.gz",
        r2 = config["sequencing_data"]+"/{fid}_2_sequence.fastq.gz"
    output:
        r1 = wdir + "/reads_trimmed/{fid}_R1.fast1.gz",
        r2 = wdir + "/reads_trimmed/{fid}_R2.fast1.gz",
    params:
        ref = config["bbduk_ref"]
    threads: 12
    shell:
        """
        bbduk.sh in={input.r1} in2={input.r2}  out={output.r1} out2={output.r2} qtrim=rl trimq=20 ktrim=r k=23 mink=11 hdist=1 ref={params.ref} tpe tbo threads={threads}

        """
rule get_ref:
    input:
        wdir + "/metagenome.fasta"
    output:
        wdir + "/metagenome.mmi"
    shell:
        """
        minimap2 -d {output} {input}
        """

rule extract_barcode:
    input:
        r1 = wdir + "/reads_trimmed/{fid}_R1.fast1.gz"
    output:
        b1 = wdir + "/barcodes/{fid}.txt",
    shell:
        """
        set +e
        seqkit seq -n {input} | head -n 10000 | cut -f 10 -d ":" | sort | uniq -c | sort -n -r | head -n 1 &> {output}
        exitcode=$?
        if [ $exitcode -eq 1 ]
        then
            exit 1
        else
            exit 0
        fi        
        """

rule mapping:
    input:
        r1 = wdir + "/reads_trimmed/{fid}_R1.fast1.gz",
        r2 = wdir + "/reads_trimmed/{fid}_R2.fast1.gz",
        ref = wdir + "/metagenome.mmi"
    output:
        wdir + "/alignments/{fid}.bam",
    log:
        wdir + "/alignments/{fid}.log"
    threads: 12
    shell:
        """
        minimap2 -t {threads} -ax sr {input.ref} {input.r1} {input.r2} | samtools sort -o {output} - &> {log}
        samtools index {output}
        """
rule get_only_csdsgff:
    input:
        gtf = wdir+"/metagenome.gff"
    output:
        gtf = wdir+"/metagenome_cds.gff"
    shell:
        """
        grep "CDS" {input} > {output}
        """
rule count_reads:
    input:
        bams = expand(wdir+"/alignments/{fid}.bam", fid=fids, wdir=wdir),
        gtf = wdir+"/metagenome_cds.gff"
    threads: 64
    output:
        counts = temp(wdir+"/counts_pre.tsv")
    shell:
        """
        featureCounts -p  -s 2 -F GFF -t CDS -g locus_tag -T {threads} -a {input.gtf} -o {output.counts} {input.bams} 

        """

rule clean_counts:
    input:
        annot = wdir + "/de_annots.tsv",
        cnts = wdir+"/counts_pre.tsv"
    output:
        wdir+"/counts.tsv",
        wdir+"/genes4de.tsv"
    notebook:
        "notrebook/cleancnt.r.ipynb"



