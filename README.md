# Viewing_expressed_variants


## Contributors
-  Kevin Elaba, Ankita Murmu, Rajarshi Mondal, Anukrati Nigam, ChunHsuan LO (Jason) - **Sysadmin**
-  Ahmad Khleifat, Olaitan I. Awe, Varuna Chander - **Tech support**
-  Sara Carioscia, Yejie Yun - **Writer**
-  Kevin Elaba - **Slides construction**
-  Yejie Yun, Anukrati Nigam - **Results presenter & advertisements**
-  Yejie Yun, Rajarshi Mondal, Ankita Murmu, Olaitan I. Awe, ChunHsuan LO (Jason) - **Github maintenance**
-  ChunHsuan LO (Jason) - **Lead, Liaison**


## Goal
To visualize the expression profiles and pathways associated with pathogenic variants for suggesting clinical therapy target and drug usage for Colorectal Cancer. The tool is intended to generate visualizations with a clinical focus.


## Introduction
The development of next-generation sequencing revolutionized sequencing technology, allowing researchers to construct gene expression profiles and identify potential pathogenic variants. Although next-generation sequencing has led to a leap in sequencing technology, there remain significant challenges in translating the information to a clinical setting. Expression of variants has significant clinical relevance. Analysis of variant expression and its associated cellular pathways can be used to assess cancer risk, clinical outcomes, and possible treatment tagets. However, interpretation of variant expression in a clinical setting is an ongoing challenge. For those who are not genetics specialists, raw data regarding expression of genetic variants can be uninformative. Here we aim to translate variant expression data into clinically relevent infromation for health care professionals and patients identifying pathogenic variants, associated cellular pathways, and suggestions for clinical therapy targets.


## Idea Outlines
![](pictures/idea_outlines_00.png)


## Example Data
#### Option 1: Filtered .vcf files from github  <br/>
_(Provided by: https://github.com/collaborativebioinformatics/expression_and_SNPs_to_clinic)_  <br/>
-  testSample.cancer.tab <br/>
-  testSample.cancer.vcf <br/>
-  testv25.variants.HC_hard_cutoffs_applied.cancer.tab <br/>
-  testv25.variants.HC_hard_cutoff_applied.cancer.vcf <br/>
#### Option 2: Raw .vcf files & paired multi-omic data from TCGA  <br/>
_(Provided by: https://portal.gdc.cancer.gov/)_  <br/>
-  TCGA-44-6164 (case ID) <br/>
![](pictures/Example_data.png)
-- https://portal.gdc.cancer.gov/repository?facetTab=files&filters=%7B%22op%22%3A%22and%22%2C%22content%22%3A%5B%7B%22content%22%3A%7B%22field%22%3A%22cases.case_id%22%2C%22value%22%3A%5B%220c0b610e-fe4c-406d-a5ed-5cc3b11dabf5%22%5D%7D%2C%22op%22%3A%22in%22%7D%5D%7D&searchTableTab=files <br/>
-- https://portal.gdc.cancer.gov/repository?facetTab=files&filters=%7B%22op%22%3A%22and%22%2C%22content%22%3A%5B%7B%22content%22%3A%7B%22field%22%3A%22cases.case_id%22%2C%22value%22%3A%5B%220c0b610e-fe4c-406d-a5ed-5cc3b11dabf5%22%5D%7D%2C%22op%22%3A%22in%22%7D%5D%7D&searchTableTab=files


## Installation
**1.** Installing the Git Repository as a Package
```
devtools::install_github("collaborativebioinformatics/Viewing_expressed_variants")
```
**2.** Software Requirements

The following expression variants analysis tools have been installed in this container:

```
R version 4.0.4 (2021-02-15) -- "Lost Library Book"
Bioconductor: 3.12
ggplot2: 3.3.5
gridExtra: 2.3
dplyr: 1.0.7
tidyr: 1.1.4
magrittr: 2.0.2
rWikiPathways: 1.10.0
ggradar: 0.2
ggpubr: 0.4.0

```

**3.** Setting up the Environment

A singularity container was built to run the expression variants analysis and visualization pipeline. The recipe file (expressed_variants.def) will be available this Git repository.

To build the singularity container on your unix environment, do:
```
singularity build expressed_variants.sif expressed_variants.def
```

To run the container on your unix environment, do:
```
singularity run expressed_variants.sif
```

To run specific R packages by using the container, do:
```
singularity exec expressed_variants.sif R <path_to_Rscript>
```


## Methods

### Inputs:
_VCF file (sample-> online data base - 1000 genome, TCGA, or etc.) + RNAseq bam file_
1. Expressed variants (VCF files from RNA-seq data) -> NEED TO VISUALIZE THE COLUMNS IN THE FIELDS.

### Outputs:
1. Summary statistics of expressed variants and pathogenic variants.
2. Gene ontology and KEGG Pathway analysis for the expressed pathogenic variants.
3. Potential clinical targets.

### Detailed flow charts:
![](pictures/workflow_charts_00.png)


## Implementation (codes)

### (step 1) Preparing the sample files:<br/>
**1.**<br/>
```
(codes)
```

### (step 2) Data cleaning for Vcf files or tabulated files as input (Sara):<br/>
**1.**<br/>
```
(codes)
```

### (step 3) Focusing on pathogenic variants only (Sara & Varuna):<br/>
**1.**<br/>
```
(codes)
```

### (step 4) Gene ontology and pathway analysis for pathogenic variants & genes (Yejie, Varuna, Jason, Kevin, Olaitan & Sara):<br/>
**1.**<br/>
```
(codes)
```
**2. Molecular mechanisms:**<br/>
```
(codes)
```
**3. Identification of druggable target:**<br/>
```
(codes)
```

### (step 5) Results integration (Kevin):<br/>
_PS. Two results short list for clinician<br/>_

**1. Showing Top 5 important variants/ associated pathways for clinicians:**<br/>
```
(codes)
```
**2. Broader list for researchers:**<br/>
```
(codes)
```

### (step 6) Visualization (Ankita, Jason & Anukrati):<br/>
**1. Visualization of facts about each expressed variant: what genes/pathways it affects:**<br/>
```
(codes)
```
**2. Visualization of genome tracks where variants located:**<br/>
```
(codes)
```


## Example Data outputs & results
**1. Figures & Tables:**<br/>
<img width="323" alt="Screen Shot " src="https://XXX.png">
For use by https://github.com/collaborativebioinformatics/Viewing_expressed_variants


## References
- Data: https://XXX
- GATK Best Practices https://gatk.broadinstitute.org/hc/en-us/sections/360007226651-Best-Practices-Workflows
- DNANexus documentation https://documentation.dnanexus.com/developer/apps/execution-environment/connecting-to-jobs
- https://github.com/collaborativebioinformatics/omics_clinical_reporting
- https://github.com/collaborativebioinformatics/expression_and_SNPs_to_clinic
- https://github.com/collaborativebioinformatics/snpReportR
- https://github.com/collaborativebioinformatics/Differential_Expression_and_Variant_Association
- https://github.com/collaborativebioinformatics/DeepExpression
