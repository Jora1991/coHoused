# Microbiome & Metabolome Analysis Pipeline

## Overview

This repository contains R scripts and source data for the integrated analysis of 16S rRNA sequencing and metabolomics data, including alpha/beta diversity, differential abundance analysis, and metabolic profiling across experimental groups.

## Directory Structure

```
├── 00.rawdata/                    # Source data for all analyses
├── 01.diversity.R                  # Alpha and beta diversity analysis
├── 02.top10.R                     # Top 10 taxa distribution (genus & species levels)
├── 03.differentialAnalysis.R      # Differential analysis of significant species and genera
├── 04.metabolome.R                # Metabolomics data analysis
└── 05.AfterCohoused.R             # Post-cohousing analysis
```

## File Descriptions

### Scripts

| Script | Description |
|---|---|
| `01.diversity.R` | Performs alpha diversity (within-sample richness/evenness) and beta diversity (between-sample dissimilarity) analyses based on 16S rRNA sequencing data. |
| `02.top10.R` | Generates bar plots of the top 10 most abundant taxa at genus- and species-levels. |
| `03.differentialAnalysis.R` | Identifies and visualizes differentially abundant species and genera across experimental conditions. |
| `04.metabolome.R` | Performs statistical analysis on metabolomics data. |
| `05.AfterCohoused.R` | Contains analysis code for samples collected after cohousing experiments. |

### Data Files

#### 16S rRNA Sequencing

| File | Description |
|---|---|
| `group.txt` | Group assignment information for 16S rRNA sequencing samples. |
| `alpha_diversity_index.xls` | Pre-computed alpha diversity indices for each sample. |
| `[pgs]_tax_rel.sort.screen.txt` | Relative abundance profiles at phylum, genus, and species levels (tab-separated). |

#### Metabolomics

| File | Description |
|---|---|
| `CK_B.vs.CT_B` | Blood metabolites comparing WT vs. KO groups. |
| `CK_F.vs.CT_F` | Fecal metabolites comparing WT vs. KO groups. |
| `SK_B.vs.ST_B` | Blood metabolites comparing co-housed WT (coWT) vs. co-housed KO (coKO) groups. |
| `SK_F.vs.ST_F` | Fecal metabolites comparing co-housed WT vs. co-housed KO groups. |

## Requirements

- R (>= 4.0.0)
- Key R packages: `vegan`, `phyloseq`, `ggplot2`, `pheatmap`, `dplyr`, `tidyr`
