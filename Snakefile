SAMPLES = [f"SRR{9130236 + i}" for i in range(1, 19)]
SRA_DIR = "/share/crsp/lab/choljang/itamburi/sra_data/wirka_2019_GSE131778"
REF_FA = "/dfs6/pub/itamburi/refs/GRCh_ensembl/Homo_sapiens.GRCh38.dna.primary_assembly.fa"
REF_GTF = "/dfs6/pub/itamburi/refs/GRCh_ensembl/Homo_sapiens.GRCh38.114.gtf.gz"
WD = "/dfs6/pub/itamburi/pub/CAD_snrnaseq"

rule all:
    input:
        "logs/.init",
        # cellranger mkref outputs
        f"{WD}/cellranger_GRCh38/genes/genes.gtf*",
        f"{WD}/cellranger_GRCh38/fasta/genome.fa*",
        # expand tells snakemake to expect all of these SRR fastq files/symlinks
        expand(f"{WD}/fastq/{{sample}}_S1_L001_R1_001.fastq.gz", sample=SAMPLES),
        expand(f"{WD}/fastq/{{sample}}_S1_L001_R2_001.fastq.gz", sample=SAMPLES),
        expand(f"{WD}/fastq/{{sample}}_S1_L001_I1_001.fastq.gz", sample=SAMPLES),
        expand(f"{WD}/counts/{{sample}}/outs/filtered_feature_bc_matrix/matrix.mtx.gz", sample=SAMPLES),
        #f"{WD}/GSE131780/outs/analysis/clustering/graphclust/clusters.csv"
        f"{WD}/counts/{{sample}}/outs/filtered_feature_bc_matrix/matrix.mtx.gz",
        f"{WD}/counts/{{sample}}/outs/molecule_info.h5"

# Ensure all directories exist
# directory() tells snakemake to expect the dir as an output
# touch() checks that the rule is complete
rule init_dirs:
    output:
        touch(f"{WD}/logs/.init"),
        directory(f"{WD}/counts"),
        directory(f"{WD}/fastq"),
        directory(f"{SRA_DIR}")
    shell:
        """
        cd {WD}
        mkdir -p logs fastq counts
        touch {WD}/logs/.init
        """

# Build CellRanger reference
rule mkref:
    input:
        f"{WD}/logs/.init"
    output:
        f"{WD}/cellranger_GRCh38/genes/genes.gtf*",
        f"{WD}/cellranger_GRCh38/fasta/genome.fa*"
    params:
        fasta=REF_FA,
        gtf=REF_GTF,
        outdir=f"{WD}/cellranger_GRCh38"
    log:
        f"{WD}/logs/mkref.log"
    threads: 8
    shell:
        """
        module load cellranger/8.0.1
        cellranger mkref \
          --nthreads={threads} \
          --genome={params.outdir} \
          --fasta={params.fasta} \
          --genes={params.gtf}
        """

# Download fastqs
rule download_fastq:
    input:
        "logs/.init"
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
        f"{WD}/logs/{{sample}}/download_fastq.log"
    shell:
        """
        module load sra-tools/3.0.0
        fastq-dump --split-files --gzip {params.sra} --outdir {params.outdir}
        """

# Rename for CellRanger
rule rename_fastq:
    input:
        r1=f"{SRA_DIR}/{{sample}}_1.fastq.gz",
        r2=f"{SRA_DIR}/{{sample}}_2.fastq.gz",
        r3=f"{SRA_DIR}/{{sample}}_3.fastq.gz"
    output:
        r1_out=f"{SRA_DIR}/{{sample}}_S1_L001_R1_001.fastq.gz",
        r2_out=f"{SRA_DIR}/{{sample}}_S1_L001_R2_001.fastq.gz",
        i1_out=f"{SRA_DIR}/{{sample}}_S1_L001_I1_001.fastq.gz"
    shell:
        """
        mv {input.r1} {output.r1_out}
        mv {input.r2} {output.r2_out}
        mv {input.r3} {output.i1_out}
        """

# Make symbolic links in working directory
rule make_symlink:
    input:
        r1_out=f"{SRA_DIR}/{{sample}}_S1_L001_R1_001.fastq.gz",
        r2_out=f"{SRA_DIR}/{{sample}}_S1_L001_R2_001.fastq.gz",
        i1_out=f"{SRA_DIR}/{{sample}}_S1_L001_I1_001.fastq.gz"
    output:
        r1=f"{WD}/fastq/{{sample}}_S1_L001_R1_001.fastq.gz",
        r2=f"{WD}/fastq/{{sample}}_S1_L001_R2_001.fastq.gz",
        i1=f"{WD}/fastq/{{sample}}_S1_L001_I1_001.fastq.gz"
    shell:
        """
        ln -sf {input.r1_out} {output.r1},
        ln -sf {input.r2_out} {output.r2},
        ln -sf {input.i1_out} {output.i1}
        """



# CellRanger count
rule cellranger_count:
    input:
        r1=f"{WD}/fastq/{{sample}}_S1_L001_R1_001.fastq.gz",
        r2=f"{WD}/fastq/{{sample}}_S1_L001_R2_001.fastq.gz",
        i1=f"{WD}/fastq/{{sample}}_S1_L001_I1_001.fastq.gz",
        ref=f"{WD}/cellranger_GRCh38"
    output:
        f"{WD}/counts/{{sample}}/outs/filtered_feature_bc_matrix/matrix.mtx.gz",
        f"{WD}/counts/{{sample}}/outs/molecule_info.h5"
    params:
        fastq=f"{WD}/fastq",
        ref=f"{WD}/cellranger_GRCh38",
        sample = lambda wildcards: wildcards.sample
    log:
        f"{WD}/logs/{{sample}}.cellranger_count.log"
    threads: 8
    shell:
        """
        set -e
        module load cellranger/8.0.1
        cd {WD}/counts
        cellranger count \
            --id={params.sample} \
            --transcriptome={params.ref} \
            --fastqs={params.fastq} \
            --sample={params.sample} \
            --localcores={threads} \
            --localmem=24 \
            --create-bam false
        """

# Generate CSV for aggregation
rule make_agg_csv:
    input:
        expand(f"{WD}/counts/{{sample}}/outs/molecule_info.h5", sample=SAMPLES)
    output:
        f"{WD}/GSE131780_libraries.csv"
    run:
        import os
        with open(output[0], "w") as out:
            out.write("sample_id,molecule_h5\n")
            for sample in SAMPLES:
                path = os.path.join(WD, "counts", sample, "outs", "molecule_info.h5")
                out.write(f"{sample},{path}\n")




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
