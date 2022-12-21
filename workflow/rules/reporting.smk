rule multiqc_before_trim:
    input:
        nanostats = expand("results/reports/before_trim/{strain}_before_trim_NanoStats.txt", strain=get_all_strain_ids()),
        fastqc_bt = expand("results/reports/before_trim/{strain}_{lane}_{read}_fastqc.zip", strain=get_all_strain_ids(), lane=get_all_lanes(), read=get_all_read_ids())
    output:
        "results/reports/multiqc/before_trim_multiqc.html"
    params:
        extra=""  # Optional: extra parameters for multiqc.
    log:
        "logs/multiqc/before_trim.log"
    wrapper:
        "v1.21.0/bio/multiqc"

rule multiqc_trimmed:
    input:
        nanostats = expand("results/reports/trimmed/{strain}_trimmed_NanoStats.txt", strain=get_all_strain_ids()),
        fastqc_t = expand("results/reports/trimmed/{strain}_{read}_fastqc.zip", strain=get_all_strain_ids(), read=get_all_read_ids())
    output:
        "results/reports/multiqc/trimmed_multiqc.html"
    params:
        extra=""  # Optional: extra parameters for multiqc.
    log:
        "logs/multiqc/trimmed.log"
    wrapper:
        "v1.21.0/bio/multiqc"
