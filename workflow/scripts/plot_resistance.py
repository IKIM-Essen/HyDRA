import pandas as pd
import altair as alt
import sys

sys.stderr = open(snakemake.log[0], "w")

infiles=snakemake.input.csv

out_html=snakemake.output.html


def get_agent_resistance(infile):
    sample=infile.split("/")[-2]

    arg_df=pd.read_csv(infile, index_col="id")
    #filtering out all ARG that are labeled as hypothetical protein
    arg_tab=arg_df[arg_df["function"] != "hypothetical protein"]

    # filter for the maximum probabilty in each agent column
    arg_max = pd.DataFrame(arg_tab.iloc[:,3:-1].max(), columns=["prob"])
    arg_max.astype({'prob': 'float'}).dtypes
    arg_max=arg_max.round(4)

    # add agent as column and replace it as index
    arg_max.index.name = "agent"
    arg_max.reset_index(inplace=True)

    arg_max["sample"] = sample

    return arg_max


def get_agent_resistance_all(sample_files):
    first=True

    for sample_file in sample_files:
        sample_max=get_agent_resistance(sample_file)
        if first:
            all_max = sample_max
            first=False
        else:
            all_max=pd.concat([all_max,sample_max], ignore_index=True)

    return all_max


def plot_resistance(sample_files, out_html):
    all_max=get_agent_resistance_all(sample_files)

    # number of samples -> determine height of plot
    no_samples=len(sample_files)

    chart=alt.Chart(all_max).mark_rect().encode(
        x=alt.X("agent:N", title="", axis=alt.Axis(labelAngle=-45)),
        y=alt.Y("sample:N", title="Sample"),
        color=alt.Color("prob", legend=None).scale(scheme="redyellowgreen",reverse=True).title(None),
        tooltip=[
        alt.Tooltip("sample", title="Sample:"),
        alt.Tooltip("agent", title="Antimicrobial agent:"),
        alt.Tooltip("prob", title="Probability:"),
    ],
    ).properties(
        width="container",
        height=(50 + 25*no_samples)
    )
    chart.save(out_html)

plot_resistance(infiles, out_html)

