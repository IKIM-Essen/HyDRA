import pandas as pd
import json
import sys

sys.stderr = open(snakemake.log[0], "w")

prefilt_path = snakemake.input.prefilt #"results/report_prefilt/multiqc_data.json"
postfilt_path = snakemake.input.postfilt #"results/report/multiqc_data.json"

samples = snakemake.params.samples
samples.sort()

outfile_prefilt = snakemake.output.csv_prefilt #"results/report_prefilt/QC_all_prefilt.csv"
outfile_postfilt = snakemake.output.csv_postfilt #"results/report_postfilt/QC_all_postfilt.csv"

def save_QC_table(samples, json_data, outfile):
    all_metrics = {}
    for sample in samples:
        sample_data= json_data["report_saved_raw_data"]["multiqc_nanostat"][sample]
        all_metrics[sample]=sample_data

    df = pd.DataFrame.from_dict(all_metrics)
    df.drop("STDEV read length_fastq", inplace=True)
    df_T = df.T
    df_T.sort_index(inplace=True)
    df_T.index.name="sample"
    df_T.to_csv(outfile)

# output table pre filtering
f = open(prefilt_path)
json_data_prefilt = json.load(f)
save_QC_table(samples, json_data_prefilt, outfile_prefilt)

# output table after filtering
f = open(postfilt_path)
json_data_postfilt = json.load(f)
save_QC_table(samples, json_data_postfilt, outfile_postfilt)