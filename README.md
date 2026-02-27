# COMP 383 – HCMV Snakemake Pipeline by Farah Younis

This repository contains my Snakemake workflow for automating steps 2–5 of the HCMV transcriptome analysis project. The workflow does:

- Bowtie2 read filtering against the HCMV reference genome  
- SPAdes assembly of filtered reads (k=99)  
- Assembly statistics calculation (>1000 bp contigs and total bp)  
- Longest contig extraction  
- BLAST analysis restricted to the Betaherpesvirinae subfamily  

The pipeline is built so that the user can clone the repository and run the complete workflow on included sample test data with a single command.


## Repository Contents

- `Snakefile` – Complete Snakemake workflow 
- `sample_test_data/` – Small paired end FASTQ files intended for quick testing  
- `reference_genome/` – Reference genome and automatically generated BLAST database inputs  
- `Younis_PipelineReport.txt` – Final report generated from running the pipeline on the full dataset 

All intermediate outputs are written to `pipeline_outputs/` but are ignored via `.gitignore`.

---

## Software Requirements

The following programs have to be installed and available in your PATH:

- Snakemake  
- Bowtie2  
- SPAdes  
- BLAST+ (makeblastdb, blastn)  
- NCBI Datasets CLI (`datasets`)  
- unzip  
- Standard Unix utilities (awk, grep, find, xargs, gzip)



---

## Running the Pipeline on Sample Test Data

From the root of the repository:

```
snakemake --cores 4
```

This command will(without any manual editing):

1. Build a Bowtie2 index for the HCMV reference genome  
2. Filter reads that map to HCMV  
3. Assemble filtered reads using SPAdes (k-mer size 99)  
4. Compute assembly statistics (>1000 bp contigs and total bp)  
5. Download Betaherpesvirinae genomes from NCBI  
6. Build a local BLAST nucleotide database  
7. BLAST the longest contig from each assembly  
8. Generate the final report:

```
PipelineReport.txt
```



---

## Notes

- The `sample_test_data/` directory is included so the workflow can be tested quickly
- Betaherpesvirinae genomes are downloaded automatically using the NCBI `datasets` CLI
- Only the best HSP per subject is retained in the BLAST results
- Intermediate outputs are written to `pipeline_outputs/` and are not tracked in the repository