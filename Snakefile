import os
import glob

#Folder that contains the paired-end FASTQ files for testing purposes since it is sample_test_data 
reads_directory = "sample_test_data"

#Reference genome FASTA used for Bowtie2 mapping.
reference_fasta = (
    "reference_genome/ncbi_dataset/data/GCF_000845245.1/"
    "GCF_000845245.1_ViralProj14559_genomic.fna"
)

#Bowtie2 index prefix name (Bowtie2 creates multiple .bt2 files using this prefix).
reference_index_prefix = "reference_genome/hcmv_index"

#This FASTA will be created automatically from NCBI Datasets (Betaherpesvirinae genomes)
betaherpes_fasta = "reference_genome/betaherpesvirinae_genomes.fna"

#BLAST database prefix- makeblastdb creates multiple files using this prefix
betaherpes_blastdb_prefix = "reference_genome/betaherpes_blastdb"



# this is where the datasets download zip will be saved
BETAHERPES_ZIP = "reference_genome/betaherpes_datasets.zip"

#where the zip contents will be extracted.
BETAHERPES_DIR = "reference_genome/betaherpes_datasets"

#final report file
report_file = "PipelineReport.txt"

#SPAdes settings - k-mer size and number of threads to use 
spades_k = 99
spades_threads = 4

#Column order for BLAST output (sacc = subject accession, pident = percent identity, etc.
BLAST_COLS = "sacc pident length qstart qend sstart send bitscore evalue stitle"

#Output folders for each stage 
mapped_reads_folder = "pipeline_outputs/mapped_reads"
sam_folder = "pipeline_outputs/sam"
counts_folder = "pipeline_outputs/counts"
spades_folder = "pipeline_outputs/spades"
assembly_stats_folder = "pipeline_outputs/assembly_stats"
longest_contig_folder = "pipeline_outputs/longest_contigs"
blast_folder = "pipeline_outputs/blast"


def detect_samples(reads_dir):
    # This looks through the reads folder and finds sample IDs automatically.
    # It searches for files ending in _1.fastq or _1.fastq.gz and checks that the matching _2 file exists.
    patterns = [
        os.path.join(reads_dir, "*_1.fastq.gz"),
        os.path.join(reads_dir, "*_1.fastq"),
    ]

    #Collect every R1 file that matches either pattern.
    r1_files = []
    for p in patterns:
        r1_files.extend(glob.glob(p))

    samples = []
    for r1 in sorted(r1_files):
        base = os.path.basename(r1)

        #Strip the R1 suffix to get the sample name, then construct the expected R2 path.
        if base.endswith("_1.fastq.gz"):
            sample = base[:-len("_1.fastq.gz")]
            r2 = os.path.join(reads_dir, sample + "_2.fastq.gz")
        else:
            sample = base[:-len("_1.fastq")]
            r2 = os.path.join(reads_dir, sample + "_2.fastq")

        #Only accept the sample if R2 exists so we don’t accidentally run single-end data.
        if os.path.exists(r2):
            samples.append(sample)

    #If nothing was detected, stop early with an error message to not run the whole pipeline with no data
    if len(samples) == 0:
        raise ValueError(
            f"I couldn't find any paired FASTQ files in {reads_dir}. "
            "I expected <sample>_1.fastq(.gz) and <sample>_2.fastq(.gz)."
        )

    return samples


def r1_path(sample):
    # since some of the datasets are gzipped and some are not,  this checks both options.
    p1 = os.path.join(reads_directory, f"{sample}_1.fastq.gz")
    p2 = os.path.join(reads_directory, f"{sample}_1.fastq")
    return p1 if os.path.exists(p1) else p2


def r2_path(sample):
    p1 = os.path.join(reads_directory, f"{sample}_2.fastq.gz")
    p2 = os.path.join(reads_directory, f"{sample}_2.fastq")
    return p1 if os.path.exists(p1) else p2


#automatically detect samples from the reads folder.
sample_ids = detect_samples(reads_directory)


rule all:
    #This is the finish line for the workflow
    # When these files exist, all steps have run successfully.
    input:
        report_file,
        expand(spades_folder + "/{sample_id}/contigs.fasta", sample_id=sample_ids),
        expand(assembly_stats_folder + "/{sample_id}_assembly_stats.txt", sample_id=sample_ids),
        expand(longest_contig_folder + "/{sample_id}_longest_contig.fasta", sample_id=sample_ids),
        expand(blast_folder + "/{sample_id}_blast_top5.tsv", sample_id=sample_ids),
        betaherpes_blastdb_prefix + ".nin"


rule build_hcmv_index:
    #Builds Bowtie2 index files from the reference FASTA.
    #Bowtie2 produces several files like .1.bt2, .2.bt2
    input:
        reference_fasta
    output:
        expand(
            reference_index_prefix + ".{suffix}",
            suffix=["1.bt2", "2.bt2", "3.bt2", "4.bt2", "rev.1.bt2", "rev.2.bt2"]
        )
    shell:
        """
        # bowtie2-build reads the FASTA and creates the index files using the prefix.
        bowtie2-build {input} {reference_index_prefix}
        """


rule count_read_pairs_before:
    #Counts how many read pairs exist before filtering.
    #FASTQ format has 4 lines per read, so total reads = total lines / 4.
    #In paired-end data counting R1 reads gives the number of read pairs.
    input:
        lambda wc: r1_path(wc.sample_id)
    output:
        counts_folder + "/{sample_id}_before.txt"
    shell:
        """
        mkdir -p {counts_folder}

        # If the input is gzipped, use zcat; otherwise use cat.
        if echo {input} | grep -q ".gz$"; then
            zcat {input} | awk 'END{{print NR/4}}' > {output}
        else
            cat {input} | awk 'END{{print NR/4}}' > {output}
        fi
        """


rule bowtie2_map:
    #Maps paired reads to the reference genome.
    #The SAM output is kept as a record of alignments
    #The important outputs are the mapped paired reads written by --al-conc-gz, which are used for assembly
    input:
        #Snakemake waits for the index files before running Bowtie2.
        index_files = rules.build_hcmv_index.output,
        #These functions choose the right file paths for R1 and R2 based on the sample ID and whether they are gzipped or not.
        read1 = lambda wc: r1_path(wc.sample_id),
        read2 = lambda wc: r2_path(wc.sample_id)
    output:
        #SAM file with all alignments.
        sam = sam_folder + "/{sample_id}_map.sam",
        #These contain only the paired reads that aligned (filtered reads)
        mapped1 = mapped_reads_folder + "/{sample_id}_mapped_1.fq.gz",
        mapped2 = mapped_reads_folder + "/{sample_id}_mapped_2.fq.gz"
    shell:
        """
        mkdir -p {sam_folder} {mapped_reads_folder}

        # -x uses the index prefix (not the FASTA directly).
        # --al-conc-gz writes the aligned read pairs to two files:
        #   ..._mapped_1.fq.gz and ..._mapped_2.fq.gz
        bowtie2 --quiet \
          -x {reference_index_prefix} \
          -1 {input.read1} \
          -2 {input.read2} \
          -S {output.sam} \
          --al-conc-gz {mapped_reads_folder}/{wildcards.sample_id}_mapped_%.fq.gz
        """


rule count_read_pairs_after:
    #Counts how many read pairs remain after mapping/filtering.
    #The mapped output here is gzipped, so zcat is always correct.
    input:
        mapped1 = mapped_reads_folder + "/{sample_id}_mapped_1.fq.gz"
    output:
        counts_folder + "/{sample_id}_after.txt"
    shell:
        """
        mkdir -p {counts_folder}

        # Same FASTQ math: 4 lines per read.
        zcat {input.mapped1} | awk 'END{{print NR/4}}' > {output}
        """


rule spades_assembly:
    #Runs SPAdes on the filtered reads to assemble contigs
    #The output contigs.fasta file is what later steps use
    input:
        r1 = mapped_reads_folder + "/{sample_id}_mapped_1.fq.gz",
        r2 = mapped_reads_folder + "/{sample_id}_mapped_2.fq.gz"
    output:
        contigs = spades_folder + "/{sample_id}/contigs.fasta"
    threads: spades_threads
    shell:
        """
        mkdir -p {spades_folder}/{wildcards.sample_id}

        # -k sets the k-mer size.
        # -t uses the number of threads Snakemake allocates to this rule.
        # -o is the output directory for SPAdes.
        spades.py \
          -1 {input.r1} \
          -2 {input.r2} \
          -k {spades_k} \
          -t {threads} \
          -o {spades_folder}/{wildcards.sample_id}
        """


rule assembly_stats:
    # Calculates two things from the assembly:
    # 1) number of contigs longer than 1000 bp
    # 2) total base pairs across contigs longer than 1000 bp
    input:
        contigs = spades_folder + "/{sample_id}/contigs.fasta"
    output:
        stats = assembly_stats_folder + "/{sample_id}_assembly_stats.txt"
    shell:
        r"""
        mkdir -p {assembly_stats_folder}

        # This awk reads FASTA and tracks sequence length between headers.
        # When a new header shows up, it checks the previous contig length.
        gt1000_count=$(awk '/^>/ {{if (seqlen > 1000) c++; seqlen=0; next}} {{seqlen += length($0)}} END {{if (seqlen > 1000) c++; print c+0}}' {input.contigs})

        # Same logic as above, but sums lengths instead of counting contigs.
        total_bp=$(awk '/^>/ {{if (seqlen > 1000) t+=seqlen; seqlen=0; next}} {{seqlen += length($0)}} END {{if (seqlen > 1000) t+=seqlen; print t+0}}' {input.contigs})

        echo "In the assembly of sample {wildcards.sample_id}, there are $gt1000_count contigs > 1000 bp and $total_bp total bp." > {output.stats}
        """


rule longest_contig:
    #Extracts the longest contig sequence from contigs.fasta.
    #This walks through each contig keeps the longest one, and writes it as a new FASTA
    input:
        contigs = spades_folder + "/{sample_id}/contigs.fasta"
    output:
        longest = longest_contig_folder + "/{sample_id}_longest_contig.fasta"
    shell:
        r"""
        mkdir -p {longest_contig_folder}

        # awk logic:
        # - when a header starts (>), compare the previous contig to the current max
        # - store the sequence string and its length
        # - at the end, print the longest one
        awk '
        /^>/ {{
            if (seqlen > maxlen) {{maxlen=seqlen; maxseq=seq}}
            seq=""; seqlen=0; next
        }}
        {{
            seq=seq $0; seqlen+=length($0)
        }}
        END {{
            if (seqlen > maxlen) {{maxlen=seqlen; maxseq=seq}}
            print ">{wildcards.sample_id}_longest_contig"
            print maxseq
        }}' {input.contigs} > {output.longest}
        """


rule download_betaherpes_datasets:
    # Downloads Betaherpesvirinae virus genomes from NCBI using the datasets CLI.
    # This produces a single zip file that contains many genome FASTAs.
    output:
        BETAHERPES_ZIP
    shell:
        r"""
        mkdir -p $(dirname {output})

        datasets download virus genome taxon Betaherpesvirinae \
          --include genome \
          --filename {output}
        """


rule unzip_betaherpes_datasets:
    #Unzips the datasets download into a folder
    #This keeps the raw downloaded files separate from the final concatenated FASTA
    input:
        BETAHERPES_ZIP
    output:
        directory(BETAHERPES_DIR)
    shell:
        r"""
        rm -rf {output}
        mkdir -p {output}
        unzip -q {input} -d {output}
        """


rule build_betaherpes_fasta:
    #CThis is the file that makeblastdb uses.
    input:
        BETAHERPES_DIR
    output:
        betaherpes_fasta
    shell:
        r"""
        mkdir -p $(dirname {output})
        tmp="{output}.tmp"
        rm -f "$tmp"

        # find searches for FASTA-like files and concatenates them in sorted order.
        # sort makes the output deterministic across runs.
        find {input} -type f \( -name "*.fna" -o -name "*.fa" -o -name "*.fasta" \) -print0 \
          | sort -z \
          | xargs -0 cat > "$tmp"

        mv "$tmp" {output}
        """


rule make_betaherpes_blastdb:
    input:
        rules.build_betaherpes_fasta.output
    output:
        betaherpes_blastdb_prefix + ".nin"
    shell:
        r"""
        makeblastdb \
          -in {input} \
          -out {betaherpes_blastdb_prefix} \
          -title betaherpesvirinae \
          -dbtype nucl
        """


rule blast_longest_contig:
    #BLASTs the longest contig against the local Betaherpesvirinae database.
    #This saves only the top 5 subject hits and only the best HSP per subject.
    input:
        db = betaherpes_blastdb_prefix + ".nin",
        query = longest_contig_folder + "/{sample_id}_longest_contig.fasta"
    output:
        hits = blast_folder + "/{sample_id}_blast_top5.tsv"
    shell:
        r"""
        mkdir -p {blast_folder}

        blastn \
          -query {input.query} \
          -db {betaherpes_blastdb_prefix} \
          -out {output.hits} \
          -outfmt "6 {BLAST_COLS}" \
          -max_target_seqs 5 \
          -max_hsps 1
        """


rule make_report:
    # Writes one report file that includes:
    # - read counts before/after mapping
    # - assembly stats
    # - BLAST results for each sample
    input:
        before = expand(counts_folder + "/{sample_id}_before.txt", sample_id=sample_ids),
        after  = expand(counts_folder + "/{sample_id}_after.txt",  sample_id=sample_ids),
        stats  = expand(assembly_stats_folder + "/{sample_id}_assembly_stats.txt", sample_id=sample_ids),
        blast  = expand(blast_folder + "/{sample_id}_blast_top5.tsv", sample_id=sample_ids)
    output:
        report_file
    run:
        with open(output[0], "w") as out:
            for sample_id in sample_ids:
                with open(f"{counts_folder}/{sample_id}_before.txt") as f:
                    before_count = f.read().strip()
                with open(f"{counts_folder}/{sample_id}_after.txt") as f:
                    after_count = f.read().strip()

                out.write(
                    f"Sample {sample_id} had {before_count} read pairs before and {after_count} read pairs after Bowtie2 filtering.\n"
                )

            out.write("\n")

            for sample_id in sample_ids:
                with open(f"{assembly_stats_folder}/{sample_id}_assembly_stats.txt") as f:
                    out.write(f.read().strip() + "\n")

            out.write("\n")

            for sample_id in sample_ids:
                out.write(f"{sample_id}:\n")
                out.write(f"{BLAST_COLS}\n")

                blast_path = f"{blast_folder}/{sample_id}_blast_top5.tsv"

                if os.path.exists(blast_path) and os.path.getsize(blast_path) > 0:
                    with open(blast_path) as f:
                        out.write(f.read().rstrip() + "\n")
                else:
                    out.write("NO_HITS_FOUND\n")

                out.write("\n")


     