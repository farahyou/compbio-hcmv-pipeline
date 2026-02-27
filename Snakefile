import os

# step 2-5 pipeline settings
sample_ids = ["SRR5660030", "SRR5660033", "SRR5660044", "SRR5660045"]
reads_directory = "sample_test_data"

reference_fasta = "reference_genome/ncbi_dataset/data/GCF_000845245.1/GCF_000845245.1_ViralProj14559_genomic.fna"
reference_index_prefix = "reference_genome/hcmv_index"

mapped_reads_folder = "pipeline_outputs/mapped_reads"
sam_folder = "pipeline_outputs/sam"
counts_folder = "pipeline_outputs/counts"

spades_folder = "pipeline_outputs/spades"
assembly_stats_folder = "pipeline_outputs/assembly_stats"
longest_contig_folder = "pipeline_outputs/longest_contigs"
blast_folder = "pipeline_outputs/blast"

report_file = "Younis_PipelineReport.txt"

betaherpes_fasta = "reference_genome/betaherpesvirinae_genomes.fna"
betaherpes_blastdb_prefix = "reference_genome/betaherpes_blastdb"

BLAST_COLS = "sacc pident length qstart qend sstart send bitscore evalue stitle"


rule all:
    input:
        report_file,
        expand(spades_folder + "/{sample_id}/contigs.fasta", sample_id=sample_ids),
        expand(assembly_stats_folder + "/{sample_id}_assembly_stats.txt", sample_id=sample_ids),
        expand(longest_contig_folder + "/{sample_id}_longest_contig.fasta", sample_id=sample_ids),
        expand(blast_folder + "/{sample_id}_blast_top5.tsv", sample_id=sample_ids),
        betaherpes_blastdb_prefix + ".nin"


# step 2
rule build_hcmv_index:
    input:
        reference_fasta
    output:
        expand(reference_index_prefix + ".{suffix}",
               suffix=["1.bt2", "2.bt2", "3.bt2", "4.bt2", "rev.1.bt2", "rev.2.bt2"])
    shell:
        """
        bowtie2-build {input} {reference_index_prefix}
        """


# step 2
rule bowtie2_map:
    input:
        index_files = rules.build_hcmv_index.output,
        read1 = reads_directory + "/{sample_id}_1.fastq.gz",
        read2 = reads_directory + "/{sample_id}_2.fastq.gz"
    output:
        sam = sam_folder + "/{sample_id}_map.sam",
        mapped1 = mapped_reads_folder + "/{sample_id}_mapped_1.fq.gz",
        mapped2 = mapped_reads_folder + "/{sample_id}_mapped_2.fq.gz"
    shell:
        """
        mkdir -p {sam_folder} {mapped_reads_folder}

        bowtie2 --quiet \
          -x {reference_index_prefix} \
          -1 {input.read1} \
          -2 {input.read2} \
          -S {output.sam} \
          --al-conc-gz {mapped_reads_folder}/{wildcards.sample_id}_mapped_%.fq.gz
        """


# step 2
rule count_read_pairs_before:
    input:
        reads_directory + "/{sample_id}_1.fastq.gz"
    output:
        counts_folder + "/{sample_id}_before.txt"
    shell:
        """
        mkdir -p {counts_folder}
        zcat {input} | awk 'END{{print NR/4}}' > {output}
        """


# step 2
rule count_read_pairs_after:
    input:
        mapped1 = mapped_reads_folder + "/{sample_id}_mapped_1.fq.gz"
    output:
        counts_folder + "/{sample_id}_after.txt"
    shell:
        """
        mkdir -p {counts_folder}
        zcat {input.mapped1} | awk 'END{{print NR/4}}' > {output}
        """


# step 3
rule spades_assembly:
    input:
        r1 = mapped_reads_folder + "/{sample_id}_mapped_1.fq.gz",
        r2 = mapped_reads_folder + "/{sample_id}_mapped_2.fq.gz"
    output:
        contigs = spades_folder + "/{sample_id}/contigs.fasta"
    threads: 4
    shell:
        """
        mkdir -p {spades_folder}/{wildcards.sample_id}
        spades.py \
          -1 {input.r1} \
          -2 {input.r2} \
          -k 99 \
          -t {threads} \
          -o {spades_folder}/{wildcards.sample_id}
        """


# step 4
rule assembly_stats:
    input:
        contigs = spades_folder + "/{sample_id}/contigs.fasta"
    output:
        stats = assembly_stats_folder + "/{sample_id}_assembly_stats.txt"
    shell:
        r"""
        mkdir -p {assembly_stats_folder}

        gt1000_count=$(awk '/^>/ {{if (seqlen > 1000) c++; seqlen=0; next}} {{seqlen += length($0)}} END {{if (seqlen > 1000) c++; print c+0}}' {input.contigs})
        total_bp=$(awk '/^>/ {{if (seqlen > 1000) t+=seqlen; seqlen=0; next}} {{seqlen += length($0)}} END {{if (seqlen > 1000) t+=seqlen; print t+0}}' {input.contigs})

        echo "In the assembly of sample {wildcards.sample_id}, there are $gt1000_count contigs > 1000 bp and $total_bp total bp." > {output.stats}
        """


# step 5
rule longest_contig:
    input:
        contigs = spades_folder + "/{sample_id}/contigs.fasta"
    output:
        longest = longest_contig_folder + "/{sample_id}_longest_contig.fasta"
    shell:
        r"""
        mkdir -p {longest_contig_folder}

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


# step 5
rule make_betaherpes_blastdb:
    input:
        betaherpes_fasta
    output:
        betaherpes_blastdb_prefix + ".nin"
    shell:
        r"""
        makeblastdb -in {input} -out {betaherpes_blastdb_prefix} -title betaherpesvirinae -dbtype nucl
        """


# step 5
rule blast_longest_contig:
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


# steps 2, 4, 5 report
rule make_report:
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
                out.write(f"Sample {sample_id} had {before_count} read pairs before and {after_count} read pairs after Bowtie2 filtering.\n")

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