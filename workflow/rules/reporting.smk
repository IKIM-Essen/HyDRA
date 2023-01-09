rule multiqc_before_trim:
    input:
        nanostats = expand("results/reports/before_trim/{strain}/{strain}_before_trim_NanoStats.txt", strain=get_all_strain_ids()),
        fastqc_bt = expand("results/reports/before_trim/{strain}/{strain}_{lane}_{read}_fastqc.zip", strain=get_all_strain_ids(), lane=get_all_lanes(), read=get_all_read_ids())
    output:
        report("results/reports/multiqc/before_trim_multiqc.html", caption="../report/before_trim_multiqc.rst", category="before_trim"),
        ont = "results/reports/multiqc/before_trim_multiqc_data/multiqc_nanostat.txt",
        ill = "results/reports/multiqc/before_trim_multiqc_data/multiqc_fastqc.txt"
    params:
        extra=""  # Optional: extra parameters for multiqc.
    log:
        "logs/multiqc/before_trim.log"
    wrapper:
        "v1.21.0/bio/multiqc"

use rule multiqc_before_trim as multiqc_trimmed with:
    input:
        nanostats = expand("results/reports/trimmed/{strain}/{strain}_trimmed_NanoStats.txt", strain=get_all_strain_ids()),
        fastqc_t = expand("results/reports/trimmed/{strain}/{strain}_{read}_fastqc.zip", strain=get_all_strain_ids(), read=get_all_read_ids())
    output:
        "results/reports/multiqc/trimmed_multiqc.html",#, caption="../report/trimmed_multiqc.rst", category="trimmed")
        ont = "results/reports/multiqc/trimmed_multiqc_data/multiqc_nanostat.txt",
        ill = "results/reports/multiqc/trimmed_multiqc_data/multiqc_fastqc.txt"
    log:
        "logs/multiqc/trimmed.log"

rule coverages:
    input:
        ont = "results/reports/multiqc/{stage}_multiqc_data/multiqc_nanostat.txt",
        ill = "results/reports/multiqc/{stage}_multiqc_data/multiqc_fastqc.txt"
    output:
        "results/reports/{stage}_coverages.csv"
    params:
        genomesizes = get_genome_size_dict()
    log:
        "logs/coverage/{stage}_coverages.log"
    conda:
        "../envs/watstats.yaml"
    script:
        "../scripts/coverage.py"
    