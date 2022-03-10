# Viewing_expressed_variants

## Contributors 
-  Kevin Elaba, Ankita Murmu, Rajarshi Mondal, Anukrati Nigam, ChunHsuan LO (Jason) - **Sysadmin** 
-  Ahmad Khleifat, Olaitan Awe, Varuna Chander - **Tech support**
-  Sara Carioscia, Yejie Yun - **Writer**
-  Kevin Elaba - **Slides construction** 
-  Yejie Yun, Anukrati Nigam - **Results presenter & advertisements** 
-  Yejie Yun, Rajarshi Mondal, Ankita Murmu, Olaitan I. Awe, ChunHsuan LO (Jason) - **Github maintenance**
-  ChunHsuan LO (Jason) - **Lead, Liaison** 

## Goal 
To visualize the expression profiles and pathways associated with variants for suggesting clinical therapy target & drug usage.

## Introduction 
(Describe)

## Idea Outlines
![](pictures/idea_outlines.png)

## Example test Data 
(Describe)

## Installation & environments setup
**1.** Install the package
```
devtools::install_github("collaborativebioinformatics/Viewing_expressed_variants")
```
**2.** setting the environments
```
(codes)
```

## Methods

### Inputs:
**VCF file (sample-> online data base - 1000 genome or etc.) + RNAseq bam file**
1. Expressed variants (VCF files from RNA seq data) -> NEED TO VISUALIZE THE COLUMNS IN THE FIELDS.
2. Making technical framework for input of VCF and BAM files for both visualization of reports as well as looking at the input of the read coverage for isoforms and gene count distribution.
### Outputs:
1. figures. (Overlapping SNV site to Expression level - read coverage, distribution of pathogenic variants, which genes have the highest overlap for pathogenic variants, circular-omic plot of overlap of structural variants) 
2. Table for statistics. (Fot those annotation recorded by the VCF files.) 
3. Gene ontology analysis. (for the expressed variants.) & Pathway analysis & KEGG 
4. To suggest clinical therapy target & drug usage.
### Detailed flow charts:
![](pictures/workflow_charts.png)

## Implementation (codes)
I. Preparing the sample files:<br/>
**1.**<br/>
```
(codes)
```
----------//-------------//---------------------<br/>
II. Data cleaning for Vcf files or tabulated files as input (Sara):<br/>
**1.**<br/>
```
(codes)
```
----------//-------------//---------------------<br/>
III. Focusing on pathogenic variants only (Sara & Varuna):<br/>
**1.**<br/>
```
(codes)
```
----------//-------------//---------------------<br/>
IV. Pathway analysis for pathogenic variants & genes (Olaitan & Sara):<br/>
**1.**<br/>
```
(codes)
```
----------//-------------//---------------------<br/>
V. Gene ontology analysis - molecular mechanisms and identification of druggable target (Yejie, Varuna, Jason & Kevin):<br/>
**1.** Molecular mechanisms:<br/>
```
(codes)
```
**2.** Identification of druggable target:<br/>
```
(codes)
```
----------//-------------//---------------------<br/>
VI. Results integratiion (Kevin):<br/>
PS. Two results short list for clinician<br/>
**1.** Top 5 important variants/ associated pathways:<br/>
```
(codes)
```
**2.** Broader list for researchers:<br/>
```
(codes)
```
----------//-------------//---------------------
VII. Visualization (Ankita, Jason & Anukrati):<br/>
**1.** Visualization of facts about each expressed variant: what genes/pathways it affects:<br/>
```
(codes)
```
**2.** Visualization of genome tracks where variants located:<br/>
```
(codes)
```
----------//-------------//---------------------<br/>

##Example outputs & results
1. results pictures:
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
- https://github.com/collaborativebioinformatics/DeepExpression (might be useful for both projects)

