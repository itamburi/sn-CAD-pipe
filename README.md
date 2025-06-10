# sn-CAD-pipe


This pipeline processes human single-nucleus RNA-seq (snRNA-seq) data using CellRanger and a Snakemake workflow. It is built for execution on the UCI RCIC cluster.

Current implementation processes SRA accessions from GSE131780, published by Wirka RC et al, Nat Med 2019 PMID: 31359001

---

## 🔬 Pipeline Overview

This pipeline:

- Initializes necessary folders
- Builds a CellRanger-compatible reference (from FASTA + GTF)
	- Human gtf and fasta files must be specified in snakefile header
- Downloads and renames SRA FASTQ files
	- SRA tools must be configured in advance on cluster
	- Specify download location in snakefile
	- Current setup to downloads to lab CRSP storage space
- Symbolically links renamed FASTQs
- Runs cellranger count per sample
- Optionally aggregates counts with cellranger aggr

---

## 📂 Input & Output Structure

**Input:**

- SRA data for 18 SRR samples: SRR9130237 to SRR9130254
- Reference genome FASTA and GTF

**Output:**

```
CAD_snrnaseq/
├── fastq/ # Linked FASTQs in CellRanger format
├── counts/ # Output folders per sample
│ ├── SRR9130237/outs/
│ │ ├── filtered_feature_bc_matrix/
│ │ │ └── matrix.mtx.gz
│ │ └── molecule_info.h5
├── cellranger_GRCh38/ # CellRanger reference folder
├── GSE131780_libraries.csv # Aggregation input CSV
├── GSE131780/outs/... # Aggregated output (optional)
├── logs/ # Log files per rule/sample
└── Snakefile # The Snakemake workflow
```

---

## ⚙️ Dependencies

Modules loaded in this pipeline:

- `sra-tools/3.0.0` (for fastq-dump)
- `cellranger/8.0.1`

Other requirements:

- Python 3
- Snakemake (v7+ recommended)

---

## 🧱 Rules Breakdown

### rule all

Top-level rule defining the expected outputs: initialized folders, reference files, FASTQs, counts, and (optionally) aggregated clustering results.

### init_dirs

Creates necessary directories (`logs/`, `counts/`, `fastq/`) and a `.init` file.

### mkref

Builds CellRanger-compatible reference from provided FASTA + GTF files.

Output: `cellranger_GRCh38/`

### fastq-dump (implied)

Downloads SRA reads and splits into 3 FASTQ files (R1, R2, I1).

*Note*: On RCIC HPC3, sra-tools cahce musy be preconfigured with vbd-config by: 
```
module load sra-tools/3.0.0
vdb-config -i
```
In the interactive windown navigate to CACHE window. Replace location of user repository (`/share/crsp/lab/choljang/itamburi/sra_data/cache`).
[See also](https://www.youtube.com/watch?v=ye4W6zTWtb4&ab_channel=Dr.Asif%E2%80%99sMol.Biology)

*Note*: vdb-config -i sets the VDB cache location (for .sra files), but fastq-dump writes FASTQ files to the current directory by default; use --outdir to specify a custom output location for FASTQ files.

*The actual FASTQs download location is specified in the snakefile header*

### rename_fastq

Renames SRA output files to CellRanger's naming convention:

```
SRRxxxxxx_1.fastq.gz -> SRRxxxxxx_S1_L001_R1_001.fastq.gz
```

*Explanation:*
`fastq-dump` generates FASTQ files in the format:
- `SRR#######_1.fastq.gz`
- `SRR#######_2.fastq.gz`
- `SRR#######_3.fastq.gz`

The following reference [here](https://davetang.org/muse/2018/06/06/10x-single-cell-bam-files/) gives us a breakdown of what each FASTQ contains
Using `zcat SRR##_1.fastq.gz | head`, we confirm:
- `_1`: 26 bp reads = cell barcode + UMI → `R1`
- `_2`: 98 bp reads = cDNA → `R2`
- `_3`: 8 bp reads = i7 sample index → `I1`

To make the files compatible with **Cell Ranger's naming convention**, we rename as follows:
- From: `SRR#######_1.fastq.gz`
- To: `<sample>_S1_L001_R1_001.fastq.gz` (and similarly for R2/I1)



### make_symlink
*Note* 

Creates symlinks to the renamed FASTQs in the working directory.

### cellranger_count

Runs CellRanger count per sample. Outputs `matrix.mtx.gz` and `molecule_info.h5` per sample.

### make_agg_csv

Generates a CSV listing all `molecule_info.h5` paths for use in CellRanger aggregation.

### cellranger_aggr (optional)

Aggregates counts from all samples into a unified dataset using the CSV.

---

## 🚀 How to Run

1. Load Conda & Activate Snakemake Environment (or module)

    ```bash
    module load anaconda
    conda activate snakemake-env  # or however you installed it
    ```

2. Submit with SLURM

    Use the provided `run_pipeline.sh`:

    ```bash
    sbatch run_pipeline.sh
    ```

3. Or Run Interactively (for testing)

    ```bash
    snakemake --cores 8 --use-envmodules --rerun-incomplete -p
    ```

---

## 📌 Notes for RCIC Cluster

- All rules use `module load` to avoid needing separate environments.
- Output logs are stored in `logs/`
- Rule logs can be viewed via:

    ```bash
    less logs/SRR9130237.cellranger_count.log
    ```

- For debugging, use flags:

    ```bash
    --printshellcmds --reason --keep-going
    ```

---

## 🧪 Testing

To test on a subset:

- Comment out some samples in `SAMPLES = [...]`
- Or pass a subset sample to Snakemake via `--config sample=SRR9130237`

---

## 👤 Author

Ian Tamburini – UCI RCIC cluster user

---

## 📎 References

- [10x Genomics Cell Ranger](https://support.10xgenomics.com/)
- [Snakemake Documentation](https://snakemake.readthedocs.io/en/stable/)
