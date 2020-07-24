# Xeva: XEnograft Visualization & Analysis

## Integration of molecular and pharmacological profiles of patient-derived xenograft models.

### Abtsract
One of the key challenges in cancer precision medicine is finding robust biomarkers of drug response. Patient-derived tumor xenografts (PDXs) have emerged as reliable preclinical models since they better recapitulate tumor response to chemo- and targeted therapies. However, the lack of standard tools poses a challenge in the analysis of PDXs with molecular and pharmacological profiles. Efficient storage, access and analysis is key to the realization of the full potential of PDX pharmacogenomic data. We have developed Xeva (XEnograft Visualization & Analysis), an open-source software package for processing, visualization and integrative analysis of a compendium of in vivo pharmacogenomic datasets. The Xeva package follows the PDX minimum information (PDX-MI) standards and can handle both replicate-based and 1x1x1 experimental designs. We used Xeva to characterize the variability of gene expression and pathway activity across passages. We found that only a few genes and pathways have passage specific alterations (median intraclass correlation of 0.53 for genes and positive enrichment score for 92.5% pathways). For example, activity of the mRNA 3'-end processing and elongation arrest and recovery pathways were strongly affected by model passaging (gene set enrichment analysis false discovery rate [FDR] <5%). We then leveraged our platform to link the drug response and the pathways whose activity is consistent across passages by mining the Novartis PDX Encyclopedia (PDXE) data containing 1,075 PDXs spanning 5 tissue types and 62 anticancer drugs. We identified 87 pathways significantly associated with response to 51 drugs (FDR < 5%), including associations such as erlotinib response and signaling by EGFR in cancer pathways and MAP kinase activation in TLR cascade and binimetinib response. Among the significant pathway-drug associations, we found novel biomarkers based on gene expressions, Copy Number Aberrations (CNAs) and mutations predictive of drug response (concordance index > 0.60; FDR < 0.05). Xeva provides a flexible platform for integrative analysis of preclinical in vivo pharmacogenomics data to identify biomarkers predictive of drug response, a major step toward precision oncology.

## Citation

Mer, Arvind S., Wail Ba-Alawi, Petr Smirnov, Yi X. Wang, Ben Brew, Janosch Ortmann, Ming-Sound Tsao, David W. Cescon, Anna Goldenberg, and Benjamin Haibe-Kains. "Integrative pharmacogenomics analysis of patient-derived xenografts." Cancer research 79, no. 17 (2019): 4539-4550.
doi: http://dx.doi.org/10.1158/0008-5472.CAN-19-0349


## How to install

Note: This version of Xeva requires R>=3.6 . If you are looking for R 3.5 version of Xeva, 
download the source package from: https://figshare.com/s/bc4b892b38c04937b791

- Install latest version of Xeva directly from Github using `devtools`:
```
library(devtools)
devtools::install_github("bhklab/Xeva")
```
- To install from Bioconductor
```
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("Xeva", version = "3.9")
```
