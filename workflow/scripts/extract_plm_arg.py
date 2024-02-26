import pandas as pd
import sys

sys.stderr = open(snakemake.log[0], "w")

infile=snakemake.input.tsv

out_arg=snakemake.output.arg
out_non_arg=snakemake.output.non_arg


df=pd.read_table(infile)

seq_ids=df["seq_id"].tolist()
df["id"] = [elem.split(" ", maxsplit=1)[0] for elem in seq_ids]
func_list = [elem.split(" ", maxsplit=1)[1] for elem in seq_ids]
df.insert(0, "function", func_list)

df.drop(["seq_id"], axis=1 ,inplace=True)
df.set_index("id",inplace=True)

non_arg_df=df.loc[df['pred'].isna()]
non_arg_df= non_arg_df.drop(["pred"], axis=1)

arg_df=df.dropna()

non_arg_df.to_csv(out_non_arg)
arg_df.to_csv(out_arg)