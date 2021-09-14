#!/usr/bin/env python3
# -*- coding: utf-8 -*-


import os
import pandas as pd
from urllib.request import urlretrieve



# =============================================================================
# now = os.getcwd()
# 
# url = 'http://api.gdc.cancer.gov/data/9a4679c3-855d-4055-8be9-3577ce10f66e'
name = 'EBPlusPlusAdjustPANCAN_IlluminaHiSeq_RNASeqV2-v2.geneExp.tsv'
# exp_filepath = os.path.join(now, name)
# if not os.path.exists(now):
#     os.makedirs(now)
#     
# if not os.path.exists(exp_filepath):
#     urlretrieve(url, exp_filepath)
# else:
#     print('Downloaded data file already exists, skipping download')
# =============================================================================

# Only reading in the first 100 rows
# tcga_expr_df = pd.read_csv(name, index_col=0, sep='\t')

# Commit from https://github.com/cognoma/cancer-data/
sample_commit = 'da832c5edc1ca4d3f665b038d15b19fced724f4c'

# Reading in cancer types data
url = 'https://raw.githubusercontent.com/cognoma/cancer-data/{}/mapping/tcga_cancertype_codes.csv'.format(sample_commit)
cancer_types_df = pd.read_csv(url,
                              dtype='str',
                              keep_default_na=False)

cancertype_codes_dict = dict(zip(cancer_types_df['TSS Code'],
                                 cancer_types_df.acronym))

# Reading in sample types data
url_sample_types = 'https://raw.githubusercontent.com/cognoma/cancer-data/{}/mapping/tcga_sampletype_codes.csv'.format(sample_commit)
sample_types_df = pd.read_csv(url_sample_types, dtype='str')

sampletype_codes_dict = dict(zip(sample_types_df.Code,
                                 sample_types_df.Definition))


genes_commit = 'ad9631bb4e77e2cdc5413b0d77cb8f7e93fc5bee'
url_genes = 'https://raw.githubusercontent.com/cognoma/genes/{}/data/genes.tsv'.format(genes_commit)
gene_df = pd.read_csv(url_genes, sep='\t')

# Only consider protein-coding genes
gene_df = (
    gene_df.query("gene_type == 'protein-coding'")
)


# Load gene updater - define up to date Entrez gene identifiers where appropriate
url = 'https://raw.githubusercontent.com/cognoma/genes/{}/data/updater.tsv'.format(genes_commit)
# =============================================================================
# updater_df = pd.read_csv(url, sep='\t')
# 
# old_to_new_entrez = dict(zip(updater_df.old_entrez_gene_id,
#                              updater_df.new_entrez_gene_id))
# 
# tcga_expr_df.index = tcga_expr_df.index.map(lambda x: x.split('|')[1])
# 
# tcga_expr_df = (tcga_expr_df
#     .dropna(axis='rows')
#     .rename(index=old_to_new_entrez)
#     .groupby(level=0).mean()
#     .transpose()
#     .sort_index(axis='rows')
#     .sort_index(axis='columns')
# )
# 
# tcga_expr_df.index.rename('sample_id', inplace=True)
# =============================================================================


# Update sample IDs to remove multiple samples measured on the same tumor
# and to map with the clinical information
tcga_expr_df.index = tcga_expr_df.index.str.slice(start=0, stop=15)
tcga_expr_df = tcga_expr_df.loc[~tcga_expr_df.index.duplicated(), :]

# Filter for valid Entrez gene identifiers
tcga_expr_df = tcga_expr_df.loc[:, tcga_expr_df.columns.isin(gene_df.entrez_gene_id.astype(str))]

# Extract sample type in the order of the gene expression matrix
tcga_id = pd.DataFrame(tcga_expr_df.index)

# Extract the last two digits of the barcode and recode sample-type
tcga_id = tcga_id.assign(sample_type = tcga_id.sample_id.str[-2:])
tcga_id.sample_type = tcga_id.sample_type.replace(sampletype_codes_dict)

# Extract the first two ID numbers after `TCGA-` and recode cancer-type
tcga_id = tcga_id.assign(cancer_type = tcga_id.sample_id.str[5:7])
tcga_id.cancer_type = tcga_id.cancer_type.replace(cancertype_codes_dict)

# Append cancer-type with sample-type to generate stratification variable
tcga_id = tcga_id.assign(id_for_stratification = tcga_id.cancer_type.str.cat(tcga_id.sample_type))

# Get stratification counts - function cannot work with singleton strats
stratify_counts = tcga_id.id_for_stratification.value_counts().to_dict()

# Recode stratification variables if they are singletons
tcga_id = tcga_id.assign(stratify_samples_count = tcga_id.id_for_stratification)
tcga_id.stratify_samples_count = tcga_id.stratify_samples_count.replace(stratify_counts)
tcga_id.loc[tcga_id.stratify_samples_count == 1, "stratify_samples"] = "other"
tcga_id = tcga_id.set_index('sample_id')

# =============================================================================
# # Write out files for downstream use
# file = os.path.join(now, 'tcga_sample_identifiers.tsv')
# 
# (
#     tcga_id.drop(['stratify_samples', 'stratify_samples_count'], axis='columns')
#     .to_csv(file, sep='\t', index=False)
# )
# =============================================================================

#print(tcga_expr_df.head())

from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np

# Standardizing gene expression values for PCA
standardized = StandardScaler().fit_transform(tcga_expr_df.values)

# Putting standardized data into a dataframe.
standardized_df = pd.DataFrame(data = standardized,
                               index = tcga_expr_df.index,
                               columns = tcga_expr_df.columns)

# Fitting PCA object to data and determining how many components are necessary.
# =============================================================================
# pca = PCA()
# components = pca.fit_transform(standardized)
# 
# # Creating scree plot.
# comp_var = np.round(pca.explained_variance_ratio_.cumsum() * 100, decimals = 1)
# 
# plt.bar(x = range(1,len(comp_var) + 1), height = comp_var)
# plt.ylabel('Percentage of Explained Variance')
# plt.xlabel('Principal Component')
# plt.title('Scree Plot')
# plt.plot([0,65],[90,90], 'k-')
# 
# plt.show()
# # Scree plot shows we need about 50 principal components to cover 90% of the variance in the data
# 
# 
# 
# components_df = pd.DataFrame(data=components,
#                              index = tcga_expr_df.index,
#                              columns = ['PC' + str(x) for x in range(1, len(pca.components_) + 1)])
# 
# # Combining PCA components with cancer types and definitions.
# component_cancer_type_df = components_df.join(tcga_id[['cancer_type','id_for_stratification']])
# # =============================================================================
# # print(component_cancer_type['Principal Component 1'].max())
# # print(component_cancer_type['Principal Component 1'].min())
# # print(component_cancer_type['Principal Component 2'].max())
# # print(component_cancer_type['Principal Component 2'].min())
# # =============================================================================
# 
# # Creating a dataframe of only cancer types with 500 or more samples.
# reduced_df = component_cancer_type_df.loc[component_cancer_type_df
#                                           ['cancer_type'].isin(
#                                               ['BRCA','KIRC',
#                                                'LUAD','THCA',
#                                                'UCEC','HNSC',
#                                                'LUSC','PRAD','LGG'])]
#                                      
# 
# # Graphing PCA plot
# sns.scatterplot(data = reduced_df,
#                 x = 'PC1',
#                 y = 'PC2',
#                 hue = 'cancer_type',
#                 alpha = 0.5,
#                 )
# plt.title('PCA Plot')
# plt.show()
# =============================================================================




def tcga_pca_plotting(gene_data,cancer_types_iterable):
    """Takes data and returns scree plot of PCA components and PCA Plot"""
    standardized_data = StandardScaler().fit_transform(gene_data.values)
    
    pca_object = PCA()
    pca_components = pca_object.fit_transform(standardized_data)
    
    #Creating Scree Plot
    comp_var = np.round(pca_object.explained_variance_ratio_.cumsum() * 100, decimals = 1)

    plt.bar(x = range(1,len(comp_var) + 1), height = comp_var)
    plt.ylabel('Percentage of Explained Variance')
    plt.xlabel('Principal Components')
    plt.title('Scree Plot')
    plt.plot([0,65],[90,90], 'k-')
    
    plt.show()

    pca_components_df = pd.DataFrame(data=pca_components,
                                 index = gene_data.index,
                                 columns = ['PC' + str(x) for x in range(1, len(pca_object.components_) + 1)])

    combined_df = pca_components_df.join(tcga_id[['cancer_type','id_for_stratification']])

    final_df = combined_df.loc[combined_df['cancer_type'].isin(cancer_types_iterable)]
    
    # Graphing PCA Plot of selected cancer types
    sns.scatterplot(data = final_df,
                    x = 'PC1',
                    y = 'PC2',
                    hue = 'cancer_type',
                    alpha = 0.5,
                    )
    plt.title('PCA Plot')
    plt.show()    


import umap
import umap.plot

def tcga_umap(gene_data,cancer_types_iterable,n_neighbors = 15,min_dist = 0.1):
    standardized_data = StandardScaler().fit_transform(gene_data.values)

    reducer = umap.UMAP(n_neighbors=n_neighbors,min_dist=min_dist)

    umap_data = reducer.fit_transform(standardized_data)
    umap_df = pd.DataFrame(umap_data,
                       index = gene_data.index,
                       columns = ['UMAP1','UMAP2'])
    umap_cancer =  umap_df.join(tcga_id[['cancer_type','id_for_stratification']])

    final_df = umap_cancer.loc[umap_cancer['cancer_type'].isin(cancer_types_iterable)]
    
    sns.scatterplot(data = final_df,
                    x = 'UMAP1',
                    y = 'UMAP2',
                    hue = 'cancer_type',
                    alpha = 0.5,
                    )
    plt.title('UMAP Plot')




# Subsetting for cancers that only appear a certain a amount of times (500).
v = tcga_id['cancer_type'].value_counts()
cancer_list = tcga_id[tcga_id['cancer_type'].isin(v.index[v.gt(500)])]

# PCA Plot and Scree Plot
# tcga_pca_plotting(tcga_expr_df,cancer_list['cancer_type'].unique())

# UMAP Plot
# tcga_umap(tcga_expr_df,cancer_list['cancer_type'].unique(),n_neighbors = 15, min_dist=0.1)


gene_df_iso = gene_df[['entrez_gene_id','symbol']]

# Creating dictionary that will be used to map entrez ids to gene symols.
gene_dict = gene_df_iso.set_index('entrez_gene_id').to_dict()['symbol']
gene_dict = {str(k):str(v) for k,v in gene_dict.items()}

tcga_symbols = tcga_expr_df.rename(columns = gene_dict)



from tcga_functions import *

# Loading data set initially to find only most mutated genes.
# =============================================================================
# mutated_100 = find_most_mutated('pancan_mutation_freeze.tsv',100)    
# 
# 
# mutated_list = list(mutated_100.index)
# 
# mutated_list.append('SAMPLE_BARCODE')
# 
# 
# 
# # Filter for only valid genes and genes I have expression data for
# mutation_data = pd.read_csv('pancan_mutation_freeze.tsv',
#                             sep = '\t',
#                             usecols = mutated_list)
# 
# mutation_data = mutation_data.set_index('SAMPLE_BARCODE')
# 
# 
# mutation_data = mutation_data.loc[:, mutation_data.columns.isin(gene_df.symbol.astype(str))]
# =============================================================================


# mutation_data.columns = [str(x) + '_mutation' for x in mutation_data.columns]

# Finding only genes that are in both the mutation data and gene_symbol data.
genes_inters = set(mutation_data.columns).intersection(set(tcga_symbols.columns))
print('These are the genes you can choose from: {}'.format(genes_inters))

# Small function to that will create a scatterplot of gene expression data stratified by mutation status.
# plot_mutation(mutation_data,tcga_symbols,'DMD')

# Function that plots mutation data using UMAP
# UMAP_mutation_plot(tcga_symbols,mutation_data,'DMD')



# Function that filters for cancer type and gene mutation data, then creates a UMAP plot of the expression data.
def UMAP_cancer_plot_1(expression_data, cancer_data, cancer_type, gene_mutation_data, gene_name, n_neighbors = 15, min_dist = 0.1):
    """This function first filters for cancer type then trains UMAP only on that cancer type expression data."""
    single_cancer = cancer_data[cancer_data['cancer_type'] == cancer_type] # filter for the cancer type of interest
       
    expr_data_cancer = expression_data.merge(single_cancer['cancer_type'], left_index = True, right_index = True)
    expr_data_cancer = expr_data_cancer.drop('cancer_type', axis = 1)
    standardized_data = StandardScaler().fit_transform(expr_data_cancer.values) # standardizing gene expression data
    
    reducer = umap.UMAP(n_neighbors = n_neighbors, min_dist = min_dist) #initializing UMAP model
    
    umap_data = reducer.fit_transform(standardized_data) # Fitting UMAP model to expression data and putting back into dataframe.
    umap_df = pd.DataFrame(umap_data,
                           index = expr_data_cancer.index,
                           columns = ['UMAP1','UMAP2'])
    
    final_df = umap_df.merge(mutation_data[gene_name],left_index = True, right_index = True)
    
    sns.scatterplot(data = final_df,
                    x = 'UMAP1',
                    y = 'UMAP2',
                    hue = gene_name,
                    )
    plt.title(cancer_type + ' UMAP Plot')
    plt.show()


# Similar function that trains UMAP on ALL the expression data before filtering.
def UMAP_cancer_plot_2(expression_data, cancer_data, cancer_type, gene_mutation_data, gene_name, n_neighbors = 15, min_dist = 0.1):
    """This function first trains UMAP on ALL of the expression data, then does the filtering of cancer type."""
    single_cancer = cancer_data[cancer_data['cancer_type'] == cancer_type] # filter for the cancer type of interest
    
    standardized_data = StandardScaler().fit_transform(expression_data.values) # standardizing gene expression data
    reducer = umap.UMAP(n_neighbors = n_neighbors, min_dist = min_dist) #initializing UMAP model
    umap_data = reducer.fit_transform(standardized_data) # Fitting UMAP model to expression data and putting back into dataframe.

    umap_df = pd.DataFrame(umap_data,
                           index = expression_data.index,
                           columns = ['UMAP1','UMAP2'])
        
    expr_data_cancer = umap_df.merge(single_cancer['cancer_type'], left_index = True, right_index = True)
    expr_data_cancer = expr_data_cancer.drop('cancer_type', axis = 1)
    
    
    final_df = expr_data_cancer.merge(mutation_data[gene_name],left_index = True, right_index = True)
    
    sns.scatterplot(data = final_df,
                    x = 'UMAP1',
                    y = 'UMAP2',
                    hue = gene_name,
                    )
    plt.title(cancer_type + ' UMAP Plot')
    plt.show()                                              


# =============================================================================
# def UMAP_cancer_plot_3(expression_data, cancer_data, cancer_type, gene_mutation_data, gene_name, n_neighbors = 15, min_dist = 0.1):
#     """This function first filters for cancer type then trains UMAP only on that cancer type expression data."""
#     single_cancer = cancer_data[cancer_data['cancer_type'] == cancer_type] # filter for the cancer type of interest
#        
#     expr_data_cancer = expression_data.merge(single_cancer['cancer_type'], left_index = True, right_index = True)
#     expr_data_cancer = expr_data_cancer.drop('cancer_type', axis = 1)
#     standardized_data = StandardScaler().fit_transform(expr_data_cancer.values) # standardizing gene expression data
#     
#     reducer = umap.UMAP(n_neighbors = n_neighbors, min_dist = min_dist) #initializing UMAP model
#     
#     umap_data = reducer.fit(standardized_data) # Fitting UMAP model to expression data and putting back into dataframe.
#     umap_df = pd.DataFrame(umap_data,
#                            index = expr_data_cancer.index,
#                            columns = ['UMAP1','UMAP2'])
#     
#     final_df = umap_df.merge(mutation_data[gene_name],left_index = True, right_index = True)
#     
#     umap.plot.points(umap_data, labels = final_df[gene_name])
# =============================================================================
    

sing_can = cancer_list[cancer_list['cancer_type'] == 'LGG']

expr_data_c = tcga_symbols.merge(sing_can['cancer_type'], left_index = True, right_index = True) 
expr_data_c = expr_data_c.drop('cancer_type', axis = 1) 
standardized_d = StandardScaler().fit_transform(expr_data_c.values)

stand_df = pd.DataFrame(data = standardized_d,
                        index = expr_data_c.index,
                        columns = expr_data_c.columns)


mut_df = stand_df.merge(mutation_data['IDH1'],left_index = True, right_index = True)

reducer = umap.UMAP()

UMAP_data = reducer.fit(mut_df.loc[:, mut_df.columns != 'IDH1_y'])

umap.plot.points(UMAP_data, labels = mut_df['IDH1_y'], theme = 'fire')
plt.title('LGG UMAP Plot')



# UMAP_cancer_plot_1(tcga_symbols, cancer_list, 'LGG', mutation_data, 'IDH1')

# UMAP_cancer_plot_2(tcga_symbols, cancer_list, 'LGG', mutation_data, 'IDH1')

# UMAP_cancer_plot_3(tcga_symbols,cancer_list, 'LGG', mutation_data, 'IDH1')