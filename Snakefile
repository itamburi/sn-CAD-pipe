SAMPLES = [f"SRR{9130236 + i}" for i in range(1, 19)] # range doesnt include last number
SRA_DIR = "/share/crsp/lab/choljang/itamburi/sra_data/wirka_2019_GSE131778"
REF_FA = "/dfs6/pub/itamburi/refs/GRCh_ensembl/Homo_sapiens.GRCh38.dna.primary_assembly.fa"
REF_GTF = "/dfs6/pub/itamburi/refs/GRCh_ensembl/Homo_sapiens.GRCh38.114.gtf.gz"
WD = "/dfs6/pub/itamburi/CAD_snrnaseq"

rule all:
    input:
        #rule init_dirs outputs
        f"{WD}/logs/.init",
        # cellranger mkref outputs
        f"{WD}/cellranger_GRCh38/_mkref.done",
        #expand tells snakemake to expect all of these SRR fastq files/symlinks
        expand(f"{SRA_DIR}/{{sample}}_1.fastq.gz", sample=SAMPLES),
        expand(f"{SRA_DIR}/{{sample}}_2.fastq.gz", sample=SAMPLES),
        expand(f"{SRA_DIR}/{{sample}}_3.fastq.gz", sample=SAMPLES),
        expand(f"{WD}/fastq/{{sample}}_S1_L001_R1_001.fastq.gz", sample=SAMPLES),
        expand(f"{WD}/fastq/{{sample}}_S1_L001_R2_001.fastq.gz", sample=SAMPLES),
        expand(f"{WD}/fastq/{{sample}}_S1_L001_I1_001.fastq.gz", sample=SAMPLES)
        #expand(f"{WD}/counts/{{sample}}/outs/filtered_feature_bc_matrix/matrix.mtx.gz", sample=SAMPLES),
        #f"{WD}/GSE131780/outs/analysis/clustering/graphclust/clusters.csv"
        #expand(f"{WD}/counts/{{sample}}/outs/filtered_feature_bc_matrix/matrix.mtx.gz", sample=SAMPLES),
        #expand(f"{WD}/counts/{{sample}}/outs/molecule_info.h5", sample=SAMPLES)

# *** Initialize directories ***
# directory() tells snakemake to expect the dir as an output
# touch() checks that the rule is complete
rule init_dirs:
    output:
        touch(f"{WD}/logs/.init")
    shell:
        """
        mkdir -p {WD}/logs {WD}/fastq {WD}/counts
        touch {WD}/logs/.init
        """

# *** Build CellRanger reference ***
rule mkref:
    input:
        f"{WD}/logs/.init"
    output:
        f"{WD}/cellranger_GRCh38/_mkref.done" # make a "done" file after cellranger mkref is finished
    params:
        fasta=REF_FA,
        gtf=REF_GTF,
        builddir=WD,
        outdir="cellranger_GRCh38"
    log:
        f"{WD}/logs/mkref.log"
    threads: 8
    shell:
        """
        module load cellranger/8.0.1
        cd {params.builddir}
        cellranger mkref \
          --nthreads={threads} \
          --genome={params.outdir} \
          --fasta={params.fasta} \
          --genes={params.gtf} \
          &> {log}
          touch {output}
        """

#*** Download fastqs ***
rule download_fastq:
    input:
        f"{WD}/logs/.init",
        f"{SRA_DIR}"
    output:
        # SRA_DIR is a python variable, sample is a snakemake wildcard
        # f"{SRA_DIR}/{{sample}}_1.fastq.gz" uses double curly braces to escape {sample} inside f-string
        r1=f"{SRA_DIR}/{{sample}}_1.fastq.gz",
        r2=f"{SRA_DIR}/{{sample}}_2.fastq.gz",
        r3=f"{SRA_DIR}/{{sample}}_3.fastq.gz"
    params:
        # Use lambda to access {sample} in params; direct strings like "sra={sample}" won't expand
        sra = lambda wildcards: wildcards.sample,
        outdir=SRA_DIR
    log:
        f"{WD}/logs/download_fastq_{{sample}}.log"
    shell:
        """
        module load sra-tools/3.0.0
        fastq-dump --split-files --gzip {params.sra} --outdir {params.outdir}
        """

# *** Make symlinks to data in WD ***
# symlink names are formatted to match expected cellranger count input
rule make_symlink:
    input:
        r1=f"{SRA_DIR}/{{sample}}_1.fastq.gz",
        r2=f"{SRA_DIR}/{{sample}}_2.fastq.gz",
        r3=f"{SRA_DIR}/{{sample}}_3.fastq.gz"
    output:
        r1_out=f"{WD}/fastq/{{sample}}_S1_L001_R1_001.fastq.gz",
        r2_out=f"{WD}/fastq/{{sample}}_S1_L001_R2_001.fastq.gz",
        i1_out=f"{WD}/fastq/{{sample}}_S1_L001_I1_001.fastq.gz"
    shell:
        """
        ln -sf {input.r1} {output.r1_out}
        ln -sf {input.r2} {output.r2_out}
        ln -sf {input.r3} {output.i1_out}
        """


# CellRanger count
#rule cellranger_count:
#    input:
#        r1=f"{WD}/fastq/{{sample}}_S1_L001_R1_001.fastq.gz",
#        r2=f"{WD}/fastq/{{sample}}_S1_L001_R2_001.fastq.gz",
#        i1=f"{WD}/fastq/{{sample}}_S1_L001_I1_001.fastq.gz",
#        ref=f"{WD}/cellranger_GRCh38"
#    output:
#        f"{WD}/counts/{{sample}}/outs/filtered_feature_bc_matrix/matrix.mtx.gz",
#        f"{WD}/counts/{{sample}}/outs/molecule_info.h5"
#    params:
#        fastq=f"{WD}/fastq",
#        ref=f"{WD}/cellranger_GRCh38",
#        sample = lambda wildcards: wildcards.sample
#    log:
#        f"{WD}/logs/{{sample}}.cellranger_count.log"
#    threads: 8
#    shell:
#        """
#        set -e
#        module load cellranger/8.0.1
#        cd {WD}/counts
#        cellranger count \
#            --id={params.sample} \
#            --transcriptome={params.ref} \
#            --fastqs={params.fastq} \
#            --sample={params.sample} \
#            --localcores={threads} \
#            --localmem=24 \
#            --create-bam false
#        """

# Generate CSV for aggregation
#rule make_agg_csv:
#    input:
#        expand(f"{WD}/counts/{{sample}}/outs/molecule_info.h5", sample=SAMPLES)
#    output:
#        f"{WD}/GSE131780_libraries.csv"
#    run:
#        import os
#        with open(output[0], "w") as out:
#            out.write("sample_id,molecule_h5\n")
#            for sample in SAMPLES:
#                path = os.path.join(WD, "counts", sample, "outs", "molecule_info.h5")
#                out.write(f"{sample},{path}\n")




# Aggregate counts
# rule cellranger_aggr:
#    input:
#        csv=f"{WD}/GSE131780_libraries.csv",
#        f"{WD}/counts/{{sample}}/outs/molecule_info.h5"
#    output:
#        f"{WD}/GSE131780/outs/analysis/clustering/graphclust/clusters.csv"
#    params:
#        outdir=f"{WD}/GSE131780"
#    threads: 8
#    shell:
#        """
#        module load cellranger/8.0.1
#        cellranger aggr \
#            --id={params.outdir} \
#            --csv={input.csv} \
#            --normalize=none
#        """
