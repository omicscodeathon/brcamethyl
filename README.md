<h1 align="center"> Comparative Transcriptomic Analysis of UHRF1 Knockout in Different Cancer Cell Lines Reveals Gene Expression and Pathway Dysregulation </h1>

## Table of Contents
- [Project Overview)(#Project Overview)
- [Objectives](#Objectives)
- [Methods](#Methods)
  - [Data Acquisition and Screening](#Data Acquisition and Screening)](#Data Acquisition and Screening)
  - Data Preprocessing
  - Differential Gene Expression Analysis Preprocessing
  - Correlation of  Gene expression changes across the different cancer cell lines
  - Identification of Protein-Protein Interactions among the Intersecting Genes
  - Selection of Hub Genes
- [Results](#Results)
  - Analysis of Differentially Expressed Genes in Four Cancer Types
  - Overlap analysis
  - Correlation of Gene Expression Profiles Between Cancer Cell Lines After UHRF1 Knockout
  - Gene Set Enrichment Analysis for Overlap DEGs
  - Protein-Protein Interaction Network
- [Acknowledgment](#Acknowledgment)
- [Team](#Team)
- [References](#References)
  
## Project Overview

The **UHRF1 (ubiquitin-like with PHD and Ring Finger domains 1)** gene is an important epigenetic regulator that plays a key role in modulating DNA methylation patterns and chromatin structure. 
The UHRF1 gene plays a critical role in cancer progression. Overexpression of UHRF1 has been studied to promote tumorigenesis through different mechanisms, while its knockdown is known to hinder cancer cell migration and induce apoptosis as described in the figure below. ​​The UHRF1 protein has gained considerable attention as a potential biomarker and key regulator in various cancers, though its precise role across different cancer types remains unclear. 
This study thus aims to determine whether UHRF1 downregulation can disrupt common tumor-promoting pathways and simultaneously activate pathways that suppress tumorigenesis in different cancer types. Identifying these mechanisms could reveal new biomarkers and therapeutic targets, expanding treatment options for cancer. 
To do this we analyzed publicly available datasets from Gene Expression Omnibus (GEO) containing RNA sequencing data on cancer cell lines with UHRF1 gene knockout (retinoblastoma (Y79), breast cancer (MCF-7), monocytic leukemia (Kasumi-1) and myeloid leukemia (THP). We performed differential gene expression and gene set enrichment analysis to assess differentially expressed genes (DEGs) and pathway dysregulation resulting from UHRF1 loss. Gene overlap analysis for the detection of shared DEGs across all cell lines was performed, and these shared DEGs were further investigated for functional enrichment and protein-protein interaction. 

## Objectives
### Main Objectives
To elucidate the impact of UHRF1 knockout on global gene expression patterns in leukemia, retinoblastoma, and breast cancer cell lines.

### Specific Objectives
- To perform differential expression analysis of UHRF-1 knockout data obtained from GEO database for leukemia, retinoblastoma, and breast cancer cell lines 
- To analyze and compare the differential gene expression profiles on the obtained cancer cell lines during UHR1 knockout
- To identify commonly expressed genes during UHRF-1 knockout in leukemia, retinoblastoma, and breast cancer cell lines
- To identify key biological processes, molecular pathways and protein-protein interactions(PPI) impacted by UHRF1 knockout

## Installation
1. Clone the repository:
   ```bash
   git clone https://github.com/omicscodeathon/brcamethyl/UHRF1_Knockout_Cancer_Analysis.git
2. 
  pip install -r requirements.txt

3. 


## Methods
**RNA-seq** data was used to perform differential gene expression analysis between UHRF1 knockout and control groups across leukemia, retinoblastoma, and breast cancer cell lines. **Quality control**, **genome alignment**, and **gene quantification** was done using the following tools; **FastQC**, **STAR**, and **featureCounts** respecitively. Differential expression will be analyzed with **DESeq2**, followed by **Gene Set Enrrichment Analysis** using Gene Ontology (GO) and KEGG databases. 

### Workflow

#### Gene Expression Pipeline
![image](https://github.com/user-attachments/assets/a3a0f5db-8c22-42b3-9d9c-fd0d3182c41c)

## Results
![pca-kas](https://github.com/user-attachments/assets/d23c9c86-8018-44f6-aab2-81d2aea18cbb)
#### Fig: PCA plot of top 500 most variable genes

![image](https://github.com/user-attachments/assets/83b9876d-a564-4584-96d1-9f85fe26f2dc)
#### Fig: MA plot

![image](https://github.com/user-attachments/assets/87af7723-d018-4c86-8d3d-bd7d0e5a1279)
#### Fig: Dispersion plot

![brca-heatmap](https://github.com/user-attachments/assets/78cdd1cf-b7ef-453a-b753-2c003b115061)
#### Fig: Heatmap of displaying 50 most significant DEGs of both control and treatment group

![volcano-plot-kasumi-1](https://github.com/user-attachments/assets/c6bcbdf0-d5c3-4390-88e6-89def56985cb)
#### Fig: Volcano plot of gene fold change and p-value

![y79-go](https://github.com/user-attachments/assets/5895fe1b-74da-41ee-8a4e-eb9c63bef64c)
#### Fig: Gene ontology of cellular components 

![y79-kegg](https://github.com/user-attachments/assets/92beff5a-03a8-4bce-aa03-cd6ccf6da906)
#### Fig: KEGG plot showing pathways that genes are enriched in

## Key findings
- Our analysis revealed distinct patterns of differentially expressed genes (DEGs) in each cancer type, with Y79, Kasumi-1, THP-1, and MCF-7 showing varying numbers of DEGs after UHRF1 knockout. 
- Gene overlap analysis identified 80 DEGs across all four cancer types and highlighted critical regulatory pathways 
- PPI network analysis with the 80 overlapping genes revealed 5 key functional genes GPI, SOD2, HSPD1, TXNRD1, and GLUL across the understudy cancer types during UHRF1 knockout.

## Conclusions 
The identification of commonly differentially expressed genes during UHRF1 downregulation is essential for understanding pathways affecting cancer cell survival in UHRF1-targeted therapies. Our study found that certain genes, including GPI, SOD2, HSPD1, TXNRD1, and GLUL, were activated across cancer types during UHRF1 knockout. These genes’ upregulation, linked to poor prognosis in prior studies, suggests that while UHRF1 knockdown may sometimes improve prognosis, this adaptive gene response might reduce its overall therapeutic benefit. These findings raise important questions about the complexity of UHRF1-targeted cancer treatments.

## Recommendations
**In vitro and in vivo studies** would be recommended to understand the molecular mechanisms by which UHRF1 influences these distinct processes across different cancer types. Specifically, it is important to look at how UHRF1-targeted therapies regulate these pathways 

## Usage

This documentation and tutorials outline how to use the pipeline to perform analysis in this study. 


## Acknowledgement
The authors thank the National Institutes of Health (NIH) Office of Data Science Strategy (ODSS) for their immense support before and during the October 2024 Omics codeathon organized in collaboration with the African Society for Bioinformatics and Computational Biology (ASBCB).

## Team
1. Jonathan Kalami
2. Miriam Eyram Lawson Gakpey
3. Sala Kotochi
4. Benthai Benjamin
5. Benson R. Kidenya
6. Olaitan I. Awe


## References
1. De Almeida, B. P., Apolónio, J. D., Binnie, A., & Castelo-Branco, P. (2019). Roadmap of DNA methylation in breast cancer identifies novel prognostic biomarkers. BMC Cancer, 19(1). https://doi.org/10.1186/S12885-019-5403-0
2. Geissler, F., Nesic, K., Kondrashova, O., Dobrovic, A., Swisher, E. M., Scott, C. L., & J. Wakefield, M. (2024). The role of aberrant DNA methylation in cancer initiation and clinical impacts. Therapeutic Advances in Medical Oncology, 16. https://doi.org/10.1177/17588359231220511
3. Kanwal, R., & Gupta, S. (2012). Epigenetic modifications in cancer. Clinical Genetics, 81(4), 303. https://doi.org/10.1111/J.1399-0004.2011.01809.X
4. Moore, L. D., Le, T., & Fan, G. (2012). DNA Methylation and Its Basic Function. Neuropsychopharmacology 2013 38:1, 38(1), 23–38. https://doi.org/10.1038/npp.2012.112



