import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from IPython.display import display
"""
abr_df = pd.read_table("ml_ABR/results/115_allcontigs.tab", header=0, index_col="ID")
#print(abr_df.head(5))
#abr_df.info()
#print(abr_df.describe().round(2))

#sns_plot = 
sns.histplot(abr_df["gc_content"])
#hist = abr_df.hist()
plt.show()
"""
penguins = sns.load_dataset("penguins")
sns.histplot(data=penguins, x="flipper_length_mm", hue="species", multiple="stack")