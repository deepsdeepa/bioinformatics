import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import chi2_contingency

file_path = 'PupilBioTest_PMP_revA.csv'
df = pd.read_csv(file_path)

df.columns = df.columns.str.replace('`', '').str.strip()
print(df.columns)

#Coverage Analysis

#a) Calculate the median and coefficient of variation (CV) for single CpG coverage in each tissue
df['coverage'] = df[['000', '001', '010', '011', '100', '101', '110', '111']].sum(axis=1)
tissue_stats = df.groupby('Tissue')['coverage'].agg(['median', 'std']).reset_index()
tissue_stats['cv'] = tissue_stats['std'] / tissue_stats['median'] * 100 
print(tissue_stats) 

# Set up the plotting style
# sns.set(style="whitegrid")
# plt.figure(figsize=(18, 12))

#b) Generate plots summarizing the coverage statistics
plt.subplot(2, 3, 1)
sns.violinplot(x='Tissue', y='coverage', data=df, hue='Tissue', inner="quart", palette="muted")
plt.title('Violin Plot of Coverage by Tissue')

# Bar plot 
plt.subplot(2, 3, 2)
sns.barplot(x='Tissue', y='cv', data=tissue_stats, hue='Tissue', palette="coolwarm")
plt.title('Coefficient of Variation (CV) by Tissue')

# Pairplot 
plt.subplot(2, 3, 3)
sns.pairplot(df[['coverage', 'Tissue']], hue='Tissue', palette="Set1")
plt.title('Pairplot for Coverage by Tissue')

coverage_by_tissue = df.groupby(['Tissue', 'CpG_Coordinates'])[['000', '001', '010', '011', '100', '101', '110', '111']].mean().reset_index()
coverage_corr = coverage_by_tissue[['000', '001', '010', '011', '100', '101', '110', '111']].corr()

# Heatmap
# plt.figure(figsize=(10, 8))
sns.heatmap(coverage_corr, annot=True, cmap='coolwarm', fmt='.2f')
plt.title('Heatmap of Coverage Correlation between Methylation States')
plt.show()


#Biomarker Identification
#a) Identify PMPs with high specificity for tissue differentiation, minimizing false positives for Tissue #1 while allowing some false negatives. Use statistical or machine learning approaches to assign confidence (e.g., p-values) to each PMP
pmp_data = df[['Tissue', 'CpG_Coordinates', '000', '001', '010', '011', '100', '101', '110', '111']]
pmp_counts = pd.crosstab(pmp_data['Tissue'], pmp_data['CpG_Coordinates'])
chi2, p, dof, expected = chi2_contingency(pmp_counts)
print("Chi-Square p-value:", p)

# b) Calculate the mean variant read fraction (VRF) for each PMP in both tissues
pmp_data.loc[:, 'Variant_Read_Fraction'] = pmp_data[['000', '001', '010', '011', '100', '101', '110', '111']].mean(axis=1)

mean_vrf = pmp_data.groupby(['Tissue', 'CpG_Coordinates'])['Variant_Read_Fraction'].mean().reset_index()
print(f'VRF',mean_vrf)

def simulate_sequencing_depth(data, depth_factor):
    return data.sample(frac=depth_factor, random_state=42)

sampled_50 = simulate_sequencing_depth(df, 0.5)
sampled_20 = simulate_sequencing_depth(df, 0.2)

df['Cumulative_coverage'] = df.groupby('Tissue')['coverage'].cumsum()

threshold_tissue2 = df[df['Tissue'] == 'cfDNA'].iloc[0]['Cumulative_coverage']
print(f"Threshold for Tissue #2 (cfDNA) at 1 million reads: {threshold_tissue2}")

# Get the top 10 PMPs by mean VRF
top_10_pmp = mean_vrf.sort_values(by='Variant_Read_Fraction', ascending=False).head(10)

# Compare specificity using a chi-square test
individual_cpg_data = df.groupby(['Tissue', 'CpG_Coordinates']).size().reset_index(name='Coverage')
top_pmp_specificity = top_10_pmp.merge(individual_cpg_data, on=['Tissue', 'CpG_Coordinates'])

# Run chi-square test for specificity comparison
crosstab = pd.crosstab(top_pmp_specificity['Tissue'], top_pmp_specificity['Variant_Read_Fraction'])
chi2, p, dof, expected = chi2_contingency(crosstab)
print(f"Chi-square p-value for specificity comparison: {p}")

plt.tight_layout()
plt.show()
