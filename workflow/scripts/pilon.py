#pilon ablauf automatiseren für mehrere Runden
#bwa index, mem
#samtools view, sort, index
import os

round = 2
pilon_folder = "results/pilon/"
pilon_1_assembly = "results/pilon/188/188.fasta"
ill_R1 = "results/cutadapt/188/188_R1_trimmed.fastq.gz"
ill_R2 = "results/cutadapt/188/188_R2_trimmed.fastq.gz"
strain = "188"
#os.system(f"bwa index {pilon_1_assembly}")
bwa_out = f"results/bwa/{strain}_p{round}/"
#os.system(f"mkdir {bwa_out}")
sam_out = f"{bwa_out}{strain}.sam"
#os.system(f"bwa mem {pilon_1_assembly} {ill_R1} {ill_R2} > {sam_out}")
bam_out = f"{bwa_out}{strain}.bam"
os.system(f"samtools view -Sb {sam_out} > {bam_out}")
bam_sorted_out = f"{bwa_out}{strain}.sorted.bam"
os.system(f"samtools sort {bam_out} -o {bam_sorted_out}")
os.system(f"samtools index {bam_sorted_out}")
os.system(f"pilon --genome {pilon_1_assembly} --frags {bam_sorted_out} --output {strain} --outdir {pilon_folder}{strain}_p{round}/")
