#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 26 09:18:21 2021

@author: matthewayala
"""

import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
import pandas as pd
import umap


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


# Creates UMAP plot of 
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
    
    
# Loading data set initially to find only most mutated genes.
def find_most_mutated(data,amount_to_subset):
    """Loads full data set, and selects the most mutated genes and returns them as an index object"""
    mutation_data = pd.read_csv(data, sep = '\t')
   
    mutation_data = mutation_data.set_index('SAMPLE_BARCODE')
    
    mutations = mutation_data.sum()
    
    most_mutated = mutations.nlargest(amount_to_subset)
    
    
    return most_mutated    



# Small function to that will create a scatterplot of gene expression data stratified by mutation status.
def plot_mutation(mutation_df,expression_df,gene):
    mutation_expr_df = expression_df.join(mutation_df[gene],rsuffix = '_mutation')
    sns.scatterplot(data = mutation_expr_df,
                    x = range(1,len(mutation_expr_df.index) + 1),
                    y = gene,
                    hue = '{}_mutation'.format(gene))
    
    plt.show()


# Function that plot mutation data using UMAP
def UMAP_mutation_plot(expression_df,mutation_df,gene):
    
    stand_exp = StandardScaler().fit_transform(expression_df.values)
    
    reducer = umap.UMAP(n_neighbors=15,min_dist=0.1)
    
    umap_data = reducer.fit_transform(stand_exp)
    
    umap_df = pd.DataFrame(umap_data,
                           index = expression_df.index,
                           columns = ['UMAP1','UMAP2'])
        
    mut_expr = umap_df.join(mutation_df[gene])    
    sns.scatterplot(data = mut_expr,
                    x = 'UMAP1',
                    y = 'UMAP2',
                    hue = gene,
                    alpha = 0.5,)
    
    plt.show()

