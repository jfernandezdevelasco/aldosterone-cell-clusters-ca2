# aldosterone-cell-clusters-ca2

This repository contains the data processing pipeline and code used to evaluate intracellular calcium dynamics in primary human adrenocortical cells. The research and methodology are based on a university master thesis and a journal poster [1],[2].

This repository supports the methodologies detailed in:

* **Poster:** [Synchronized Ca²⁺ Dynamics in Glomeruli-Like Human Aldosterone-Producing Cell Clusters In Vitro.](https://www.ovid.com/jnls/jhypertension/abstract/10.1097/01.hjh.0001197000.34298.a2~synchronized-ca2-dynamics-in-glomeruli-like-human?redirectionsource=fulltextview)
* **University of Padua Masters Thesis:** [Temporal Dynamics of Calcium Signalling: Clustering Adrenal
Cells Producing Aldosterone through Ca2+ FURA-2 with the use of functional statistical analysis.](https://thesis.unipd.it/handle/20.500.12608/64485)

## Key Features and Methodology

* **Data Preprocessing:** The pipeline performs Z-score normalization on raw data. It applies a Savitzky-Golay filter using a third-order polynomial over a moving window of thirteen data points to reduce noise and photobleaching artifacts while maintaining peak features. Data extraction, normalization, and filtering steps are implemented in Python.

* **Clustering Analysis:** A correlation matrix is turned into a dissimilarity matrix based on Pearson correlation coefficients. Hierarchical agglomerative clustering is executed using the SciPy library's linkage function with complete linkage to identify synchronized subclusters.

* **Feature Extraction:** The code extracts features at various dendrogram distances (0.05 to 0.25) to explore optimal cluster numbers, excluding single-cell clusters to avoid redundant information.

* **Statistical Modeling:** Linear mixed effects models are implemented using the lme4 package. These models analyze variables such as the number of predicted subclusters, sample origin (primary normal vs. aldosterone producing adenoma), and patient level random effects to accommodate inter subject variation.

* **Functional ANOVA:** A non parametric bootstrap F-type test is implemented in R to perform functional ANOVA, evaluating if variability between labeled clusters is greater than the variability within each cluster. The modified heteroscedastic version of the F-type bootstrap test is found in the following repository [Fb-Heteroscedastic-Test](https://github.com/jfernandezdevelasco/Fb-Heteroscedastic-Test)


## System Requirements

* **Languages:** Python (for preprocessing and clustering) and R (for statistical modeling and FANOVA).


* **Core Libraries:** SciPy and Pandas for clustering and data management; lme4 package in RStudio for mixed-effects models.


**Input Data Format**
The dataset processed in this study includes Fura-2 Ca2+ imaging data from nine patients, comprising 1138 primary normal tissue (SN) cells and 1190 aldosterone-producing adenoma (APA) cells, totaling 2328 cells.

## Installation and Setup

1. **Clone the repository:**
```bash
git clone https://github.com/jfernandezdevelasco/aldosterone-cell-clusters-ca2.git
cd calcium-dynamics-pipeline

```

## Input Data Format

The pipeline expects exported time series data from microscopy software.

* **Rows:** Time (seconds)
* **Columns:** 340/380 nm fluorescence ratios for individual Regions of Interest (ROIs) corresponding to single cells.

## Source

1. H. Ajjour, G. Pallafacchina, J. F. de Velasco Biasiolo, M. G. Pedersen, B. Caroccia, L. Lenzini, T. M. Seccia, G. Rossi. (2026). Synchronized Ca²⁺ Dynamics in Glomeruli-Like Human Aldosterone-Producing Cell Clusters In Vitro. *Journal of Hypertension*. DOI: [10.1097/01.hjh.0001197000.34298.a2](https://doi.org/10.1097/01.hjh.0001197000.34298.a2)
2. J. F. de Velasco Biasiolo. (2024). Temporal Dynamics of Calcium Signalling: Clustering Adrenal Cells Producing Aldosterone through Ca2+ FURA-2 with the use of functional statistical analysis. *University of Padua*. [https://thesis.unipd.it/handle/20.500.12608/64485](https://thesis.unipd.it/handle/20.500.12608/64485)
