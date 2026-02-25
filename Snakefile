SAMPLES = ["SRR5660030", "SRR5660033", "SRR5660044", "SRR5660045"]

READS_DIR = "sample_test_data"

rule all:
    input:
        "PipelineReport.txt"

rule make_report:
    output:
        "PipelineReport.txt"
    shell:
        """
        echo "Pipeline is running successfully." > PipelineReport.txt
        """

