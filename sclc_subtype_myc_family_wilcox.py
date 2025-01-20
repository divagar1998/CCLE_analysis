# python/3.12.5
import seaborn as sns
import pandas as pd 
import matplotlib.pyplot as plt 
from scipy.stats import wilcoxon
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('model', help='path for Model.csv')
parser.add_argument('tpm', help='path for CCLE_TPM.csv')
parser.add_argument('cn', help='path for CCLE_CN.csv')

args = parser.parse_args()

model_df = pd.read_csv(args.model)
tpm_df = pd.read_csv(args.tpm)
cn_df = pd.read_csv(args.cn)

#unique_values = model_df['OncotreeSubtype'].unique()
#print(unique_values)

sclc_model_df = model_df[model_df['OncotreeSubtype']=='Small Cell Lung Cancer']
sclc_model_df = sclc_model_df[['ModelID','StrippedCellLineName']]

sclc_genes_tpm_df = tpm_df[['Unnamed: 0','MYCL (4610)','MYC (4609)', 'MYCN (4613)', 'ASCL1 (429)', 'NEUROD1 (4760)', 'POU2F3 (25833)']]
sclc_genes_tpm_df.rename(columns=
    {'Unnamed: 0': 'ModelID',
    'MYCL (4610)':'MYCL TPM',
    'MYC (4609)': 'MYC TPM',
    'MYCN (4613)': 'MYCN TPM', 
    'ASCL1 (429)': 'ASCL1 TPM', 
    'NEUROD1 (4760)': 'NEUROD1 TPM', 
    'POU2F3 (25833)': 'POU2F3 TPM'
}, inplace=True)

sclc_genes_cn_df = cn_df[['Unnamed: 0','MYCL (4610)','MYC (4609)', 'MYCN (4613)', 'ASCL1 (429)', 'NEUROD1 (4760)', 'POU2F3 (25833)']]
sclc_genes_cn_df.rename(columns=
    {'Unnamed: 0': 'ModelID',
    'MYCL (4610)':'MYCL CN',
    'MYC (4609)': 'MYC CN',
    'MYCN (4613)': 'MYCN CN', 
    'ASCL1 (429)': 'ASCL1 CN', 
    'NEUROD1 (4760)': 'NEUROD1 CN', 
    'POU2F3 (25833)': 'POU2F3 CN'
}, inplace=True)

merge1_df = pd.merge(sclc_model_df,sclc_genes_tpm_df, on='ModelID')
sclc_merged_df = pd.merge(merge1_df, sclc_genes_cn_df, on='ModelID')

subtype_mapping = {
    'ASCL1 TPM': 'SCLC-A',
    'NEUROD1 TPM': 'SCLC-N',
    'POU2F3 TPM': 'SCLC-P'
}
# Step 1: Ensure the correct subtype assignment (based on max gene expression)
sclc_merged_df[['MYCL TPM', 'MYC TPM']] = (2**sclc_merged_df[['MYCL TPM', 'MYC TPM']])-1

sclc_merged_df['SCLC subtype'] = sclc_merged_df[['ASCL1 TPM', 'NEUROD1 TPM', 'POU2F3 TPM']].idxmax(axis=1).map(subtype_mapping)


sclc_a_df = sclc_merged_df[sclc_merged_df['SCLC subtype'] == 'SCLC-A']
stat, p_value = wilcoxon(sclc_a_df['MYCL TPM'], sclc_a_df['MYC TPM'])
print(f"Wilcoxon test statistic: {stat}")
print(f"P-value: {p_value}")

sclc_n_df = sclc_merged_df[sclc_merged_df['SCLC subtype'] == 'SCLC-N']
stat, p_value = wilcoxon(sclc_n_df['MYCL TPM'], sclc_n_df['MYC TPM'])
print(f"Wilcoxon test statistic: {stat}")
print(f"P-value: {p_value}")

sclc_p_df = sclc_merged_df[sclc_merged_df['SCLC subtype'] == 'SCLC-P']
stat, p_value = wilcoxon(sclc_p_df['MYCL TPM'], sclc_p_df['MYC TPM'])
print(f"Wilcoxon test statistic: {stat}")
print(f"P-value: {p_value}")

# Function to perform Wilcoxon signed-rank test and return p-value
def calculate_p_value(df, column1, column2):
    # Perform Wilcoxon signed-rank test
    stat, p_value = wilcoxon(df[column1], df[column2])
    return p_value

# Create a figure and two subplots (side by side)
fig, axes = plt.subplots(1, 2, figsize=(7, 4))  # Increased figsize

# SCLC-A Plot (Left)
df_long_a = sclc_a_df[['MYCL TPM', 'MYC TPM']].melt(var_name='Gene', value_name='TPM')
sns.boxplot(data=df_long_a, x='Gene', y='TPM', palette="Set2", width=0.3, showfliers=False, ax=axes[0])
sns.stripplot(data=df_long_a, x='Gene', y='TPM', color='black', alpha=0.7, jitter=True, dodge=False, ax=axes[0])
axes[0].set_title('SCLC-A', fontsize=12)
axes[0].set_ylabel('TPM', fontsize=10)

# Calculate p-value for SCLC-A
p_value_a = calculate_p_value(sclc_a_df, 'MYCL TPM', 'MYC TPM')
p_value_a_rounded = round(p_value_a, 2)

# Annotate p-value on the plot (SCLC-A)
axes[0].annotate(f'p-value = {p_value_a_rounded:.2f}', xy=(0.5, 0.9), ha='center', va='bottom', fontsize=8,
                 color='black', xycoords='axes fraction', fontweight='bold')

# SCLC-N Plot (Right)
df_long_n = sclc_n_df[['MYCL TPM', 'MYC TPM']].melt(var_name='Gene', value_name='TPM')
sns.boxplot(data=df_long_n, x='Gene', y='TPM', palette="Set2", width=0.3, showfliers=False, ax=axes[1])
sns.stripplot(data=df_long_n, x='Gene', y='TPM', color='black', alpha=0.7, jitter=True, dodge=False, ax=axes[1])
axes[1].set_title('SCLC-N', fontsize=12)
axes[1].set_ylabel('TPM', fontsize=10)

# Calculate p-value for SCLC-N
p_value_n = calculate_p_value(sclc_n_df, 'MYCL TPM', 'MYC TPM')
print(f"SCLC-N P-value: {p_value_n:.3e}")  # Debugging: print p-value

# Round p-value to 2 decimal places
p_value_n_rounded = round(p_value_n, 2)

# Annotate p-value on the plot (SCLC-N)
axes[1].annotate(f'p-value = {p_value_n_rounded:.2f}', xy=(0.5, 0.9), ha='center', va='bottom', fontsize=8,
                 color='black', xycoords='axes fraction', fontweight='bold')

# Save the figure
plt.savefig('sclc_subtype_myc_wilcox.png', dpi=300, bbox_inches='tight')
