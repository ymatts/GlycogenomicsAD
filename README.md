# Glycogene Analysis Repository

This repository contains the R scripts and data for performing comprehensive analysis on glycogenes, including data import, initial analysis, meta-analysis, forest plots, network analyses, and pathway analyses using PCA.

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

## Usage

To run the analyses, open the R scripts associated with each section and execute them in your R environment. Make sure to adjust the paths to your data files as necessary.

## Contributing

Contributions to this project are welcome. Please fork the repository and submit pull requests with your proposed changes.

## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details.
```

This README structure will help users and potential contributors understand and use your project effectively. Adjust the paths and any specific details as necessary for your setup.
