import os
import pandas as pd
import pickle
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

# specify the directory you want to use
directory = 'SIMULATIONS'

data = []

for filename in os.listdir(directory):
    if filename.endswith('.pkl'):
        try:
            split_filename = filename.split('_')
            lam_value = split_filename[1].replace("lam", "")  # extract lam value
            alpha_values = [x.replace("alpha", "") for x in split_filename[2:5]]  # extract alpha values

            # reading the pickle file
            with open(os.path.join(directory, filename), 'rb') as f:
                results_dict = pickle.load(f)

            # extract mean value
            mean_value = np.mean(results_dict["levenshtein_distances"])
            mean_value = results_dict["agg_df"].loc[2, "mode"]

            # add to our data
            data.append([lam_value, alpha_values, mean_value/4])
        except:
            continue

# create pandas DataFrame
df = pd.DataFrame(data, columns=['lam', 'alpha', 'mean'])
pd.set_option('display.max_rows', None)
print(df)
df['lam'] = pd.to_numeric(df['lam'], errors='coerce')
df = df.dropna(subset=['lam'])

# convert alpha values from lists to tuples so they are hashable
df['alpha'] = df['alpha'].apply(tuple)

# convert 'lam' column to numeric
df['lam'] = pd.to_numeric(df['lam'])

print(df)

# Set the grid layout: 5 columns and as many rows as needed
ncol = 5
nrow = int(np.ceil(len(df['alpha'].unique()) / ncol))

fig, axes = plt.subplots(nrow, ncol, figsize=(20, 4*nrow), sharey=True)
axes = axes.ravel()

alpha_values = df['alpha'].unique()

for i, alpha in enumerate(alpha_values):
    data = df[df['alpha'] == alpha]
    data.groupby('lam')['mean'].mean().plot(kind='bar', ax=axes[i], color='b', edgecolor='k')
    axes[i].set_title(f'Alpha: {alpha}')
    axes[i].set_xlabel('Lambda')
    axes[i].set_ylabel('Mean Value')

plt.tight_layout()
plt.savefig('pictures/barplot.png')

