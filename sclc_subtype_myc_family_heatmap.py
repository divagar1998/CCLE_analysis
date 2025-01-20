# python/3.12.5
import seaborn as sns
import pandas as pd 
import matplotlib.pyplot as plt 
import matplotlib.patches as mpatches
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
sclc_merged_df['SCLC subtype'] = sclc_merged_df[['ASCL1 TPM', 'NEUROD1 TPM', 'POU2F3 TPM']].idxmax(axis=1).map(subtype_mapping)

# Custom sort function to preserve your custom order based on gene expression
def custom_sort(group):
    if group['SCLC subtype'].iloc[0] == 'SCLC-A':
        return group.sort_values(by='ASCL1 TPM', ascending=False)
    elif group['SCLC subtype'].iloc[0] == 'SCLC-N':
        return group.sort_values(by='NEUROD1 TPM', ascending=False)
    elif group['SCLC subtype'].iloc[0] == 'SCLC-P':
        return group.sort_values(by='POU2F3 TPM', ascending=False)

# Apply the custom sorting
sclc_merged_df = sclc_merged_df.groupby('SCLC subtype', group_keys=False).apply(lambda group: custom_sort(group))

# Reset index after sorting
sclc_merged_df.reset_index(drop=True, inplace=True)

# Step 2: Prepare heatmap data (excluding the 'SCLC subtype' column)
plot_heatmap_tpm_df = sclc_merged_df[['StrippedCellLineName', 'MYCL TPM', 'MYC TPM', 'MYCN TPM', 'ASCL1 TPM', 'NEUROD1 TPM', 'POU2F3 TPM', 'SCLC subtype']]
plot_heatmap_tpm_df.set_index('StrippedCellLineName', inplace=True)

# Extract the expression data (excluding the 'SCLC subtype' column)
expression_df = plot_heatmap_tpm_df.drop(columns=['SCLC subtype'])

# Get the indices for the chunks
chunk1 = expression_df.loc[:'SBC5']
chunk2 = expression_df.loc['CORL279':'HCC33']  # From SBC5 to HCC33
chunk3 = expression_df.loc['NCIH526':]  # From HCC33 to the end

# Create subplots (3 subplots side by side)
fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(13, 8), gridspec_kw={'width_ratios': [39, 12, 10]})

# Plot the first chunk (up to SBC5) with y-ticks
sns.heatmap(chunk1.T, annot=False, cmap="viridis", cbar=False,
    linewidths=0.5, linecolor='black', xticklabels=True, yticklabels=True, ax=ax1)
ax1.set_title('SCLC-A', pad=10,fontweight='bold')
ax1.set_xlabel('')

# Plot the second chunk (from SBC5 to HCC33) without y-ticks and colorbar
sns.heatmap(chunk2.T, annot=False, cmap="viridis", cbar=False,
    linewidths=0.5, linecolor='black', xticklabels=True, yticklabels=False, ax=ax2)
ax2.set_title('SCLC-N', pad=10,fontweight='bold')
ax2.set_xlabel('')

# Plot the third chunk (from HCC33 to end) with colorbar and without y-ticks
sns.heatmap(chunk3.T, annot=False, cmap="viridis", cbar_kws={'label': 'Expression in log2(TPM+1)'},
    linewidths=0.5, linecolor='black', xticklabels=True, yticklabels=False, ax=ax3)
ax3.set_title('SCLC-P', pad=10,fontweight='bold')
ax3.set_xlabel('')

# Adjust layout to prevent overlap
plt.tight_layout()

# Save the figure
plt.savefig('sclc_subtype_myc_family_heatmap.png', dpi=300, bbox_inches='tight')