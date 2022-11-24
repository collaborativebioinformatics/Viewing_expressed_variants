Example (input files & output pictures) files can be downloaded from the following URLs: https://drive.google.com/drive/folders/1_tQ0vh9TqMWXvbBjVgIxEZ6KWOooqIOo?usp=sharing

(1) Folder structure:
The root directory is set as server="/Users/chunhsuanlojason/Desktop/CMU_Libraries_Hackathon" in Rscript, which should be changed accordingly when installing.  
The input files (Bigwigs, VCF, and methylation.txt) will be changed to be deposited in “~/data_raw/<TCGA_caseID.filetype>“.
The output picture will be saved in “~/pictures/<filename.png>“.
The output files will be changed to be saved in “~/data_output/<TCGA_caseID.filetype>“.

(2) Preparing:
1. Download HM450.hg38.manifest.tsv from "https://drive.google.com/file/d/1zfH4MqPsbC8vUP1jV0CQ6V4i9UXpJSo2/view?usp=share_link", and save it in “~/data_raw/" folder.
2. Download GRCh38_hg38_refFlat_annotation_primaryAssemblyOnly_NoXY.bed from "https://drive.google.com/file/d/1jaSkmuT22wFnLmxtixPfaMHcZQsue_45/view?usp=share_link", and save it in “~/data_raw/" folder

(3) Download the somatic mutations from TCGA dataportal "https://portal.gdc.cancer.gov/". Choose mutect2 vcf file for the targeted caseID.
Rename the downloaded vcf file as "<TCGA_case_ID>_WES_somaticvariants.vcf" and save it in “~/data_raw/" folder.
PS. Details can be referred to "https://drive.google.com/file/d/1xN3ka7ereNAT2JhpCdedRJhy3P6XC_BI/view?usp=share_link"

(4) Download epigenetic (ATAC-Seq) files from following TCGA database "https://gdc.cancer.gov/about-data/publications/ATACseq-AWG"
1. Download Methylation data from TCGA website for "tumor tissue", and save it in “~/data_raw/" folder with the name "<TCGA_case_ID>_Tumor_atacReads_raw.bw"
PS. Details can be referred to "https://drive.google.com/file/d/1wFKPDDpf7qwOQSEjGMWIrELYtdCUCI46/view?usp=share_link"

(5) Download epigenetic (Methylation) files from following TCGA database "https://portal.gdc.cancer.gov/"
1. Download Methylation data from TCGA website for "tumor tissue", and save it in “~/data_raw/" folder with the name "<TCGA_case_ID>_Tumor_HumanMethylation450array.txt"
2. Download Methylation data from TCGA website for "paired solid normal tissue", and save it in “~/data_raw/" folder with the name "<TCGA_case_ID>_Solidnormal_HumanMethylation450array.txt".
PS. pleace replace <TCGA_case_ID> according, and seperate by "_" rather than "-".
PS. Details can be referred to "https://drive.google.com/file/d/1DNfg6ysVh1XgmTpjLEs31xSeal3fPEIT/view?usp=share_link"

(6) Run ERVcaller to detect transposable element (TE) insertiioin sites, which belong to structural variiants (SVs). ERVcaller can be downloaded from "https://dllab.org/software/ERVcaller.html"
1. For "tumor tissue", ERVcaller will take WGSeq as input and save the output VCF files in “~/data_raw/" folder. Please rename the vcf file as "<TCGA_case_ID>_Tumor_WGS_TEinsertions.vcf"
1. For "blood tissue", ERVcaller will take WGSeq as input and save the output VCF files in “~/data_raw/" folder. Please rename the vcf file as "<TCGA_case_ID>_Blood_WGS_TEinsertions.vcf"
PS. pleace replace <TCGA_case_ID> according, and seperate by "_" rather than "-".

(7) Execute the R scripts "CMUHackathon_visualization_Genometrack_datapreparing.R" and then "CMUHackathon_visualization_Genometrack_plottingcore.R"


