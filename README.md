# Repository for Glycogenomics Analysis in AD

This repository contains the R scripts and data for performing comprehensive analysis on AD glycogenomics, including data import, initial analysis, meta-analysis, forest plots, network analyses, and pathway analyses using PCA.

## Requirements

This project is built using R and several specific packages. Ensure that you have the following packages installed:

```R
# Load required libraries
library(data.table)
library(readxl)
library(ggplot2)
library(ggpie)
library(biomaRt)
library(metaviz)
library(fgsea)
library(igraph)
library(ggnet)
library(network)
library(sna)
library(ggnetwork)
library(DT)
library(Seurat)
library(CellChat)
library(nichenetr)
library(scCustomize)
library(RandomWalkRestartMH)
library(dorothea)
library(OmnipathR)
library(fedup)
library(RColorBrewer)
library(tidyverse)
library(gtools)
library(pheatmap)
library(reshape2)
library(hdWGCNA)
library(pathwayPCA)
library(metap)
```

## Project Structure

- **Data Import and Initial Analysis**: Scripts for loading and initial processing of glycogene data.
- **Meta Analysis**: Scripts for conducting meta-analyses on glycogenes.
- **Forest Plots for Pathway Activity**: Visualization of results from meta-analyses using forest plots.
- **Network Analyses**: Identification of core genes from leading edge of GSEA and network degree analyses.
- **PathwayPCA Analyses**: PCA analyses on different pathways using data from various studies.
- **Single Cell Analysis**: Single cell analysis using data from the Human Protein Atlas and GSE163577 dataset.
- **Enrichment Analysis**: Running enrichment analysis using fgsea and other methods.


## Data Sources

Data for this project can be downloaded from the following URL as described in our preprint publication: [10.1101/2023.12.25.573290v2](https://www.biorxiv.org/content/10.1101/2023.12.25.573290v2.full).

## Paper Summary
Glycosylation is increasingly recognized as a potential therapeutic target in Alzheimer’s disease. In recent years, evidence of Alzheimer’s disease-specific glycoproteins has been established. However, the mechanisms underlying their dysregulation, including tissue- and cell-type specificity, are not fully understood. We aimed to explore the upstream regulators of aberrant glycosylation by integrating multiple data sources using a glycogenomics approach. We identified dysregulation of the glycosyltransferase PLOD3 in oligodendrocytes as an upstream regulator of cerebral vessels and found that it is involved in COL4A5 synthesis, which is strongly correlated with amyloid fiber formation. Furthermore, COL4A5 has been suggested to interact with astrocytes via extracellular matrix receptors as a ligand. This study suggests directions for new therapeutic strategies for Alzheimer’s disease targeting glycosyltransferases.

## Usage

To run the analyses, open the R scripts associated with each section and execute them in your R environment. Make sure to adjust the paths to your data files as necessary.

## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details.
```
