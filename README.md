<h1 align="center"> Comparative Transcriptomic Analysis of UHRF1 Knockout in Different Cancer Cell Lines Reveals Gene Expression and Pathway Dysregulation </h1>

## Table of Contents
- [Background](#Background)
- [Objectives](#Objectives)
- [Methods](#Methods)
- [Results](#Results)
- [Acknowledgment](#Acknowledgment)
- [Team](#Team)
- [References](#References)

## Background
The **UHRF1** (ubiquitin-like with PHD and Ring Finger domains 1) protein is a critical epigenetic regulator involved in maintaining DNA methylation patterns and controlling chromatin dynamics. By interacting with both DNA methyltransferases and histone modifiers, UHRF1 plays a central role in the regulation of gene expression. Its dysregulation has been linked to cancer progression, contributing to uncontrolled cell proliferation, metastasis, and treatment resistance. However, the extent of UHRF1’s influence on global gene expression remains underexplored, particularly across different cancer types.

This study aims to bridge that gap by conducting a comparative transcriptomic analysis of UHRF1 knockout in four cancer cell lines: retinoblastoma (Y79), breast cancer (MCF-7), and leukemia (Kasumi-1 and THP-1). Through this multi-cancer approach, we seek to uncover the specific gene expression changes and pathway disruptions induced by UHRF1 loss, providing insights into its role as a master regulator of oncogenic processes across distinct cancer types. 

## Objectives
### Main Objectives
To elucidate the impact of UHRF1 knockout on global gene expression patterns in leukemia, retinoblastoma, and breast cancer cell lines.

### Specific Objectives
- To analyze and compare the differential gene expression profiles between UHRF1 knockout and control groups across leukemia, retinoblastoma, and breast cancer cell lines, identifying significantly upregulated and downregulated genes.
- To identify key biological processes and molecular pathways impacted by UHRF1 knockout, while identifying candidate genes or pathways that could serve as biomarkers for cancer prognosis or therapeutic intervention.


## Methods
**RNA-seq** data was used to perform differential gene expression analysis between UHRF1 knockout and control groups across leukemia, retinoblastoma, and breast cancer cell lines. **Quality control**, **genome alignment**, and **gene quantification** was done using the following tools; **FastQC**, **STAR**, and **featureCounts** respecitively. Differential expression will be analyzed with **DESeq2**, followed by **Gene Set Enrrichment Analysis** using Gene Ontology (GO) and KEGG databases. 

### Workflow

#### Gene Expression Pipeline
![image](https://github.com/user-attachments/assets/a3a0f5db-8c22-42b3-9d9c-fd0d3182c41c)

## Results

![image](https://github.com/user-attachments/assets/a5e386dc-1ce6-4d04-99f6-217ca3bb1b14)
#### Fig: PCA plot of top 500 most variable genes

![image](https://github.com/user-attachments/assets/83b9876d-a564-4584-96d1-9f85fe26f2dc)
#### Fig: MA plot

![image](https://github.com/user-attachments/assets/87af7723-d018-4c86-8d3d-bd7d0e5a1279)
#### Fig: Dispersion plot

![image](https://github.com/user-attachments/assets/2e12bd16-86ec-4807-bda3-37bfcb56cd1a)
#### Fig: Heatmap of displaying 50 most significant DEGs of both control and treatment group

![image](https://github.com/user-attachments/assets/55e400eb-c3ce-48ff-a767-74551e0f06fd)
#### Fig: Volcano plot of gene fold change and p-value

![GOCC_brca](https://github.com/user-attachments/assets/b64e9a6d-b7d7-48fe-8149-b4cb631a624e)
#### Fig: Gene ontology of cellular components 

![image](https://github.com/user-attachments/assets/f39ae78d-c747-4b8e-8097-6932272eeafd)
#### Fig: KEGG plot showing pathways that genes are enriched in

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



