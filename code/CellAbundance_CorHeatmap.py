###############################################################################################
#Aim: Correlation between abundance of cells
#Description: To make heatmap showing correlation between abundance of cells (Figure S13F).
# Note that MCPcounter was used here due to its ability to estimate abundance of dendritic cells
#Run command: python CellAbundance_CorHeatmap.py
###############################################################################################


import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.stats import spearmanr
from statsmodels.sandbox.stats.multicomp import multipletests
import copy

fontSize  = 8
plt.rcParams['font.size'] = fontSize
plt.rcParams["font.family"] = "Arial"

print('Raw data read in ...')
data_survival_fn = '../02.Input/MSK/Morris_cell_abundance_MCPcounter.csv'
data_df = pd.read_csv(data_survival_fn, index_col=0)

all_features = ["B cells", "T cells", "NK cells", "Dendritic cells", "Macrophages/Monocytes", "Neutrophils", "Endothelial cells", "Fibroblasts"]
data_all = data_df[all_features]

corr_out = data_all.corr(method='spearman', min_periods=1)
# set column and row names
corr_out.columns = data_all.columns
corr_out.index = data_all.columns
# create a mask to only show the lower triangle
mask = np.zeros_like(corr_out)
mask[np.triu_indices_from(mask)] = True
# set heatmap color palette and breaks
palette_length = 400
my_color = sns.color_palette("RdBu_r", n_colors=palette_length)

corr_out.to_csv('../03.Results/Data/source_data_fig05g_scc.csv', index=False)

# plot correlation heatmap
#fig, ax = plt.subplots(figsize=(3.5*0.8, 3.2*0.8))
fig, ax = plt.subplots(figsize=(3.5, 3.2))
plt.subplots_adjust(left= 0.02, bottom=0.02, right=0.9, top=0.95)

heatmap = sns.heatmap(corr_out, mask=mask, cmap=my_color, center=0,
            vmin=-1, vmax=1, xticklabels=True, yticklabels=True,
            cbar=True, cbar_kws={"shrink": 0.5, "label": "SCC"}, cbar_ax=ax.inset_axes([0.9, 0.5, 0.04, 0.5]),
            linewidths=0.1, linecolor='white', square=True, ax=ax)

# Define significance levels
sig_levels = [(0.001, '****'), (0.01, '***'), (0.05, '**'), (0.1, '*')]
# calculate significance symbols
p_list = []
for i in range(corr_out.shape[1]):
    for j in range(corr_out.shape[1]-1,i,-1):
        if mask[i, j]:
            corr, pval = spearmanr(data_all.iloc[:,i], data_all.iloc[:,j], nan_policy='omit')
            p_list.append(pval)
# adjusting p-values with multiple tests
adjusted_p_values = multipletests(p_list, method='fdr_bh')[1]
# add significance symbols
adjusted_p_values_out = copy.deepcopy(corr_out)
adjusted_p_values_out.loc[:, :] = ""
count = 0
for i in range(corr_out.shape[1]):
    for j in range(corr_out.shape[1]-1,i,-1):
        if mask[i, j]:
            adjusted_p_values_out.iloc[i, j] = adjusted_p_values[count]
            anno_text = '%.2f' % corr_out.iloc[i,j]
            adj_pval = adjusted_p_values[count]
            count+=1
            for level in sig_levels:
                if adj_pval < level[0]:
                    anno_text = '%.2f\n%s' % (corr_out.iloc[i,j], level[1])
                    break
            ax.text(i + 0.5, j + 0.5, anno_text, ha='center', va='center', fontsize=7, color='k')
cbar = heatmap.collections[0].colorbar
cbar.ax.tick_params(labelsize=8)

adjusted_p_values_out.to_csv('../03.Results/Data/source_data_fig05g_adj_pval.csv', index=False)

#### display the column names at the diagonal
continuous_features_full = all_features
for i in range(len(corr_out.columns)):
    plt.text(i + 0.5, i + 0.5, continuous_features_full[i], ha='left', va='bottom', rotation=45, fontsize=8)
#### show the plot
plt.xticks([])
plt.yticks([])
output_fig = '../03.Results/Figures/corHeatmap_Morris.pdf'
plt.savefig(output_fig, transparent = True)
plt.close()
