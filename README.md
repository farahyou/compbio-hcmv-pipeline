HCMV Snakemake Pipeline Project

This project builds a reproducible Snakemake workflow to analyze RNA-seq data from cells infected with Human Cytomegalovirus (HCMV). The pipeline filters sequencing reads to retain only viral reads, assembles those reads into longer contigs, and compares the assembled sequences to known viral genomes to determine which strains they most closely match. The final output of the workflow is a file called PipelineReport.txt containing read statistics, assembly statistics, and BLAST results for each sample.


Step 1 :Manual Download of FASTQ Files


The following SRA runs were downloaded:

SRR5660030 (Donor 1, 2dpi)
SRR5660033 (Donor 1, 6dpi)
SRR5660044 (Donor 3, 2dpi)
SRR5660045 (Donor 3, 6dpi)

The FASTQ files were generated using fasterq-dump with paired-end splitting:

mkdir -p full_input_reads
cd full_input_reads

fasterq-dump SRR5660030 --split-files --threads 8
fasterq-dump SRR5660033 --split-files --threads 8
fasterq-dump SRR5660044 --split-files --threads 8
fasterq-dump SRR5660045 --split-files --threads 8

gzip *.fastq
cd ..

Small subsets of these FASTQ files (first 10,000 reads) are included in the folder:

sample_test_data/

These are used for testing the pipeline quickly.

Software that is required: 

The following software tools must be installed:

- Snakemake
- Bowtie2
- Samtools
- SPAdes
- BLAST+
- NCBI datasets / edirect
- Python 3

Running the Pipeline (Test Data): 

From the main project directory:

snakemake -j 8

This will run steps 2–5 using the files in sample_test_data/ and generate:
PipelineReport.txt

All intermediate output files are written to:

pipeline_outputs/


Running the Pipeline with Full Reads:

Place the full FASTQ files in:

full_input_reads/

Then update the READS_DIR variable in the Snakefile to:

READS_DIR = "full_input_reads"

Run:

snakemake -j 8
