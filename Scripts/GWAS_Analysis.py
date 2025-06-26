#!/usr/bin/env python
# coding: utf-8

#Importing Pandas library 
import pandas as pd


#Loading the GWAS catalog TSV file
df = pd.read_csv("gwas_catalog_v1.0-associations_e114_r2025-06-10.tsv", sep = "\t")
print("Total no. of rows and columns are: ", df.shape)
df.head()

#extracting first 25000 rows from the dataset
new_df = df.head(25000)
new_df.to_csv("new.tsv", sep="\t", index=False)


#Finding the unique genes in new.tsv
unique_gene = df['MAPPED_GENE'].dropna().unique()
unique_gene = pd.Series(unique_gene).str.split(',').explode()
unique_gene = unique_gene.str.strip()

print("No.of unique genes present : ", len(unique_gene))
print(unique_gene[:20])

genes_df=pd.DataFrame(unique_gene, columns=['Gene'])
genes_df.to_csv("unique_genes.csv", index=False)


#Finding unique traits/diseases
unique_traits=df['DISEASE/TRAIT'].dropna().unique()

print("No.of unique traits present : ", len(unique_traits))
print(unique_traits[:20])

traits_df=pd.DataFrame(unique_traits, columns=['Disease/Trait'])
traits_df.to_csv("unique_traits.csv", index=False)


#Top genes with highest frequency across traits
genes = df['MAPPED_GENE'].dropna()
genes_split = genes.str.split(',').explode()
genes_split = genes_split.str.strip()

gene_counts = genes_split.value_counts()
top_10_genes = gene_counts.head(10)

print("Top 10 genes with most associations (across traits):")
print(top_10_genes)

top_genes_df=top_10_genes.reset_index()
top_genes_df.columns=['Gene', 'Count']
top_genes_df.to_csv("Top_10_genes.csv", index=False)


#Finding Chr number and mutations for the top 10 genes
top_10_genes = gene_counts.head(10).index.tolist()

results=[]

for gene in top_10_genes:
    gene_rows=df[df['MAPPED_GENE'].str.contains(rf'\b{gene}\b', na=False)]

    chr_series = pd.Series(gene_rows['CHR_ID'].dropna().unique())
    chr_flat = chr_series.str.split(';').explode()
    chr_cleaned = pd.to_numeric(chr_flat.str.strip(), errors='coerce').dropna().astype(int).astype(str).unique()

    chr_id = chr_cleaned
    
    snp_allele=gene_rows['STRONGEST SNP-RISK ALLELE'].dropna().unique()
    snp_allele = pd.Series(snp_allele)
    snp_allele = snp_allele[~snp_allele.str.endswith('-?')].unique()

    results.append(
      {
        'Gene': gene,
        'Chromosome ID':', '.join(chr_id),
        'SNP-Risk Allele':', '.join(snp_allele)
      }
    )

result_df=pd.DataFrame(results)
print(result_df)

result_df.to_csv("ChrID_and_mutations.tsv", sep='\t', index=False)


#Population with most occurrences for top 10 genes
top_10_genes = gene_counts.head(10).index.tolist()
pop_keywords = ['European', 'African', 'Asian', 'East Asian', 'South Asian', 'Hispanic', 
                'Japanese', 'Chinese', 'British', 'Finnish']
pop_results=[]

for gene in top_10_genes:
    gene_rows=df[df['MAPPED_GENE'].str.contains(rf'\b{gene}\b', na=False)]

    pop=[]
    for entry in gene_rows['INITIAL SAMPLE SIZE'].dropna():
        for keyword in pop_keywords:
             if keyword.lower() in entry.lower():
                pop.append(keyword)
    if pop:
        most_common_pop = pd.Series(pop).value_counts().idxmax()
    else:
        most_common_pop = 'Unknown'

    pop_results.append({'Gene': gene, 'Most Common Population': most_common_pop})

pop_df = pd.DataFrame(pop_results)
print(pop_df)

pop_df.to_csv("populations_with_most_occurance_for_the_top_10_genes.csv", index=False)


#Genes having P-value > 0.5
df['P-VALUE'] = pd.to_numeric(df['P-VALUE'], errors='coerce')

high_pval_df = df[df['P-VALUE'] > 0.5]

high_pval_genes = high_pval_df['MAPPED_GENE'].dropna().unique()

print("Genes with P-value > 0.5:")
for gene in high_pval_genes:
    print(gene)
else:
    print("No genes found with P > 0.5")

#Confirming the above result
print("Min P-VALUE:", df['P-VALUE'].min())
print("Max P-VALUE:", df['P-VALUE'].max())


#Top 10 genes with highest Odds Ratio
import numpy as np

df['OR or BETA'] = pd.to_numeric(df['OR or BETA'], errors='coerce')

df_or = df.dropna(subset=['OR or BETA', 'MAPPED_GENE'])
df_or['MAPPED_GENE'] = df_or['MAPPED_GENE'].str.replace(' - ', ',')
df_or['MAPPED_GENE'] = df_or['MAPPED_GENE'].str.split(',')

df_genes = df_or.explode('MAPPED_GENE')
df_genes['MAPPED_GENE'] = df_genes['MAPPED_GENE'].str.strip()

df_clean = df_genes[np.isfinite(df_genes['OR or BETA'])]

df_top_or = df_clean.sort_values(by='OR or BETA', ascending=False).drop_duplicates('MAPPED_GENE').head(10)

top_10_or = df_top_or[['MAPPED_GENE', 'OR or BETA']].head(10)

top_10_or_df = pd.DataFrame(top_10_or)
print("Top 10 genes with highest valid Odds Ratio:")
print(top_10_or)

top_10_or_df.to_csv("top_10_OR.csv", index=False)





