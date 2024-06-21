# embeddingsAnalytics

This repository contains toolkits for analyzing embeddings data through dimension reduction techniques and multivariate analysis, facilitating insights into high-dimensional data structures.

## Introduction

This repository currently contains a script for analyzing embeddings and performing multivariate analysis. The script reads data, formats embeddings, and visualizes results using various functions from `scDiagnostics` packages.

## Required Packages

## Installation

To run the main R script, you will need some standard Bioconductor packages which you can install with the following command:

``` r
BiocManager::install(c("SingleCellExperiment", "scater"))
```

In addition, we will use the latest version of the `scDiagnostics` package which you can install with the following command:
``` r
BiocManager::install("ccb-hms/scDiagnostics")
```

## Data Preprocessing

### Read in Data

The current main script begins by reading in data from the file `bge-small-en-v1.5_embedding.csv` which contains the embeddings of corpus nodes from several scientific articles in the `references` folder.

### Format Embeddings

Next, it formats embeddings for both corpus and question data, and creates `SingleCellExperiment` objects using formatted embeddings for corpus and question data.

## Visualizations

### MDS Scatter Plot

Generates a multidimensional scaling (MDS) plot with file type coloring.

### PCA Plot

Runs PCA on both `SingleCellExperiment` objects and visualizes principal components.

### Discriminant Space Plot

Calculates discriminant space and plots the results.
