### Setup sample wildcard:
import os
from glob import glob

SAMPLE = [os.path.basename(fn).replace("_1.fastq.gz", "")
            for fn in glob(f"resources/reads/*_1.fastq.gz")]

print("Detected the following samples:")
print(SAMPLE)
################################################################################
### Setup the desired outputs
rule all:
    input:
        "results/coverm/gene_counts.tsv"
        
################################################################################
### Filter reads with fastp
rule fastp:
    input:
        r1i = "resources/reads/{sample}_1.fastq.gz",
        r2i = "resources/reads/{sample}_2.fastq.gz"
    output:
        r1o = temp("results/fastp/{sample}_1.fastq.gz"),
        r2o = temp("results/fastp/{sample}_2.fastq.gz"),
        fastp_html = "results/fastp/{sample}.html",
        fastp_json = "results/fastp/{sample}.json"
    conda:
        "Transcriptomics_conda.yaml"
    threads:
        10
    log:
        "logs/{sample}_fastp.log"
    message:
        "Using FASTP to trim adapters and low quality sequences for {wildcards.sample}"
    shell:
        """
        fastp \
            --in1 {input.r1i} --in2 {input.r2i} \
            --out1 {output.r1o} --out2 {output.r2o} \
            --trim_poly_g \
            --trim_poly_x \
            --n_base_limit 5 \
            --qualified_quality_phred 20 \
            --length_required 60 \
            --thread {threads} \
            --html {output.fastp_html} \
            --json {output.fastp_json} \
            --adapter_sequence CTGTCTCTTATACACATCT \
            --adapter_sequence_r2 CTGTCTCTTATACACATCT \
        &> {log}
        """

## Index host genomes:
# rule index_ref:
     input:
         "resources/reference/host"
     output:
         bt2_index = "resources/reference/host/CattedRefs.fna.gz.rev.2.bt2l",
         catted_ref = "resources/reference/host/CattedRefs.fna.gz"
     conda:
         "1_QC.yaml"
     threads:
         40
     log:
         "3_Outputs/0_Logs/host_genome_indexing.log"
     message:
         "Concatenating and indexing host genomes with Bowtie2"
     shell:
         """
         # Concatenate input reference genomes
         cat {input}/*.gz > {input}/CattedRefs.fna.gz

         # Index catted genomes
         bowtie2-build \
             --large-index \
             --threads {threads} \
             {output.catted_ref} {output.catted_ref} \
         &> {log}
         """

### Map to host reference genome using STAR
rule STAR_host_mapping:
    input:
        r1i = "results/fastp/{sample}_1.fastq.gz",
        r2i = "results/fastp/{sample}_2.fastq.gz"
    output:
        non_host_r1 = "results/star/{sample}_1.fastq.gz",
        non_host_r2 = "results/star/{sample}_2.fastq.gz",
        host_bam = "results/star/{sample}_host.bam"
    params:
        r1rn = "results/star/{sample}_1.fastq",
        r2rn = "results/star/{sample}_2.fastq",
        gene_counts = "results/star/{sample}_read_counts.tsv",
        sj = "results/star/{sample}_SJ.tsv",
        host_genome = "resources/reference/XXXXXXX",
    conda:
        "Transcriptomics_conda.yaml"
    threads:
        40
    log:
        "logs/{sample}_star.log"
    message:
        "Mapping {wildcards.sample} to host genome using STAR"
    shell:
        """
        # Map reads to host genome using STAR
        STAR \
            --runMode alignReads \
            --runThreadN {threads} \
            --genomeDir {params.host_genome} \
            --readFilesIn {input.r1i} {input.r2i} \
            --outFileNamePrefix {wildcards.sample} \
            --outSAMtype BAM Unsorted \
            --outReadsUnmapped Fastx \
            --readFilesCommand zcat \
            --quantMode GeneCounts \
        &> {log}

        # Rename files
        mv {wildcards.sample}Aligned.out.bam {output.host_bam}
        mv {wildcards.sample}ReadsPerGene.out.tab {params.gene_counts}
        mv {wildcards.sample}SJ.out.tab {params.sj}
        mv {wildcards.sample}Unmapped.out.mate1 {params.r1rn}
        mv {wildcards.sample}Unmapped.out.mate2 {params.r2rn}

        # Compress non-host reads
        pigz \
            -p {threads} \
            {params.r1rn}

        pigz \
            -p {threads} \
            {params.r2rn}

        # Clean up unwanted outputs
        rm -r {wildcards.sample}_STARtmp
        rm {wildcards.sample}*out*
        """

### Functionally annotate MAGs with DRAM
rule DRAM:
    input:
        "resources/reference/MAGs.fa"
    output:
        annotations = "results/dram/MAGs_annotations.tsv.gz",
        genes = "results/dram/MAGs_genes.fna.gz",
        genesfaa = "results/dram/MAGs_genes.faa.gz",
        genesgff = "results/dram/MAGs_genes.gff.gz",
        scaffolds = "results/dram/MAGs_scaffolds.fna.gz",
        gbk = "results/dram/MAGs.gbk.gz",
        distillate = directory("results/dram/MAGs_distillate")
    params:
        outdir = "results/dram/MAGs_annotate",
        mainout = "results/dramM",
        trnas = "results/dram/MAGs_trnas.tsv.gz",
        rrnas = "results/dram/MAGs_rrnas.tsv.gz",
    threads:
        2
    resources:
        mem_gb=24,
        time='04:00:00'
    log:
        "logs/dram.log"
    message:
        "Using DRAM to functionally annotate {wildcards.MAG}"
    shell:
        """
        source activate /projects/mjolnir1/people/ncl550/0_software/miniconda3/envs/DRAM_more_modules
        
        DRAM.py annotate \
            -i {input} \
            -o {params.outdir} \
            --threads {threads} \
            --min_contig_size 1500 

        if test ! -f {params.outdir}/trnas.tsv && test ! -f {params.outdir}/rrnas.tsv
        then
        echo "neither trnas nor rrnas found"
        DRAM.py distill \
            -i {params.outdir}/annotations.tsv \
            -o {output.distillate}

        elif test -f {params.outdir}/trnas.tsv && test -f {params.outdir}/rrnas.tsv
        then
        echo "both trnas and rrnas are available"
        DRAM.py distill \
            -i {params.outdir}/annotations.tsv \
            --rrna_path {params.outdir}/rrnas.tsv \
            --trna_path {params.outdir}/trnas.tsv \
            -o {output.distillate}
        mv {params.outdir}/rrnas.tsv {params.rrnas}
        mv {params.outdir}/trnas.tsv {params.trnas}

        elif test -f {params.outdir}/trnas.tsv && test ! -f {params.outdir}/rrnas.tsv
        then
        echo "only trnas found"
        DRAM.py distill \
            -i {params.outdir}/annotations.tsv \
            --trna_path {params.outdir}/trnas.tsv \
            -o {output.distillate}
        mv {params.outdir}/trnas.tsv {params.trnas}

        elif test ! -f {params.outdir}/trnas.tsv && test -f {params.outdir}/rrnas.tsv
        then
        echo "only rrnas found"
        DRAM.py distill \
            -i {params.outdir}/annotations.tsv \
            --rrna_path {params.outdir}/rrnas.tsv \
            -o {output.distillate}
        mv {params.outdir}/rrnas.tsv {params.rrnas}

        fi    

        pigz -p {threads} {params.outdir}/*.tsv
        pigz -p {threads} {params.outdir}/*.fna
        pigz -p {threads} {params.outdir}/*.faa
        pigz -p {threads} {params.outdir}/*.gff
        pigz -p {threads} {params.outdir}/genbank/*
        pigz -p {threads} {output.distillate}/*

        mv {params.outdir}/annotations.tsv.gz {output.annotations}
        mv {params.outdir}/scaffolds.fna.gz {output.scaffolds}
        mv {params.outdir}/genes.fna.gz {output.genes}
        mv {params.outdir}/*.faa.gz {output.genesfaa}
        mv {params.outdir}/*.gff.gz {output.genesgff}
        mv {params.outdir}/genbank/* {output.gbk}

        """

################################################################################
## Index MAGs:
rule index_MAGs:
    input:
        "resources/reference/MAGs_genes.fna.gz"
    output:
        bt2_index = "resources/reference/MAGs_genes.fna.gz.rev.2.bt2l",
        MAG_genes = "resources/reference/MAGs_genes.fna.gz"
    conda:
        "Transcriptomics_conda.yaml"
    threads:
        40
    log:
        "3_Outputs/0_Logs/MAG_genes_bowtie2_indexing.log"
    message:
        "Indexing MAG catalogue with Bowtie2"
    shell:
        """
        # Rename MAG gene catalogue
        mv {input} {output.MAG_genes}

        # Index MAG gene catalogue
        bowtie2-build \
            --large-index \
            --threads {threads} \
            {output.MAG_genes} {output.MAG_genes} \
        &> {log}
        """
        
################################################################################
### Map non-host reads to DRAM genes files using Bowtie2
rule bowtie2_MAG_mapping:
    input:
        non_host_r1 = "results/star/{sample}_1.fastq.gz",
        non_host_r2 = "results/star/{sample}_2.fastq.gz",
        bt2_index = "resources/reference/MAG_genes.fna.gz.rev.2.bt2l"
    output:
        bam = "results/bowtie/{sample}.bam"
    params:
        MAG_genes = "1_References/MAG_genes.fna.gz"
    conda:
        "Transcriptomics_conda.yaml"
    threads:
        20
    log:
        "logs/{sample}_bowtie.log"
    message:
        "Mapping {wildcards.sample} to MAG genes using Bowtie2"
    shell:
        """
        # Map reads to MAGs using Bowtie2
        bowtie2 \
            --time \
            --threads {threads} \
            -x {params.MAG_genes} \
            -1 {input.non_host_r1} \
            -2 {input.non_host_r2} \
            --seed 1337 \
        | samtools sort -@ {threads} -o {output.bam} \
        &> {log}
        """
################################################################################
### Calculate the number of reads that mapped to MAG catalogue genes with CoverM
rule coverM_MAG_genes:
    input:
        expand("results/bowtie/{sample}.bam", sample=SAMPLE),
    output:
        gene_counts = "results/coverm/gene_counts.tsv",
    params:
    conda:
        "Transcriptomics_conda.yaml"
    threads:
        40
    message:
        "Calculating MAG gene mapping rate using CoverM"
    shell:
        """
        coverm contig \
            -b {input} \
            -m count \
            -t {threads} \
            --proper-pairs-only \
            > {output.gene_counts}
        """
