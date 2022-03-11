#CMUHackathon_visualization_Genometrack_datapreparing.R
####＝＝＝＝＝＝＝＝＝＝Environments setting＝＝＝＝＝＝＝＝＝＝
target_case="TCGA_44_6146"

library("icesTAF")
server="/Users/chunhsuanlojason/Desktop/CMU_Libraries_Hackathon" 
dir_input_workingspace=paste(server, "workingspace", sep="/")
dir_input_epigenetics=paste(server, "TCGA_raw_data", sep="/")
dir_temp=paste(server, "temp_CMUHackathon_visualization_Genometrack", sep="/")
mkdir(dir_temp)
dir_output=paste(server, "output", sep="/")
mkdir(dir_output)


####＝＝＝＝＝＝＝＝＝＝Data importing & initializing: VariantsSites & ATAC-Seq & methylation_HM450 & RNAseq＝＝＝＝＝＝＝＝＝＝
##===data_import_VCF===
library("AnnotationHub")
library("VariantAnnotation")

input_targeted_vcf_file <- function(caseID){
  print(caseID)
  filelist_target <- as.data.frame(list.files(dir_input_workingspace, pattern=".vcf"))
  colnames(filelist_target) = "file"
  filename=paste(caseID, "_WES_somaticvariants",".vcf",sep="")
  
  targeted_vcf_file_path <- paste(dir_input_workingspace, filename, sep="/")
  targeted_vcf_file <- readVcf(targeted_vcf_file_path)
  
  return(targeted_vcf_file)
}

assign(paste(target_case, "_input_targeted_vcf_file", sep=""), input_targeted_vcf_file(paste(target_case, sep="")))
#View(`TCGA_44_6146_input_targeted_vcf_file`)
#View(get(paste(target_case, "_input_targeted_vcf_file", sep="")))
##=====================

##===data_import_CSV===
#(omit)
#=GenomicRanges_construction=
#GenomicRanges_construction <- function(caseID){
#  print(caseID)
#  if(!is.na(print(as.vector(as.character(get(paste(caseID, "targeted_genes_list", sep=""))[,4])))) && !identical(print(as.vector(as.character(get(paste(caseID, "targeted_genes_list", sep=""))[,4]))), character(0))){
#    BAEgenes_GenomicRanges <- GRanges(
#      #=====(necessary parameters)#
#      seqnames = as.vector(as.character(get(paste(caseID, "targeted_genes_list", sep=""))[,1])),
#      ranges = IRanges(start = as.vector(as.numeric(get(paste(caseID, "targeted_genes_list", sep=""))[,2])), end = as.vector(as.numeric(get(paste(caseID, "targeted_genes_list", sep=""))[,3]))),
#      strand = Rle("*",nrow(get(paste(caseID, "targeted_genes_list", sep="")))),
#      #=====(unnecessary parameters)#
#      symbol = as.vector(as.character(get(paste(caseID, "targeted_genes_list", sep=""))[,4])),
#      dominant = as.vector(as.character(get(paste(caseID, "targeted_genes_list", sep=""))[,5])),
#      snp_position = as.vector(as.character(get(paste(caseID, "targeted_genes_list", sep=""))[,10])),
#      allelic_fraction_dR_dA_rR_rA = as.vector(as.character(get(paste(caseID, "targeted_genes_list", sep=""))[,11])),
#      note = as.vector(as.character(get(paste(caseID, "targeted_genes_list", sep=""))[,14]))
#    )
#    return(BAEgenes_GenomicRanges)
#  }
#  return(0)
#}
#
#assign(paste(target_case, "_genes_GenomicRanges", sep=""), GenomicRanges_construction(target_case))
##View(`TCGA_44_6146_genes_GenomicRanges`)
##View(get(paste(target_case, "_genes_GenomicRanges", sep="")))
##=====================

##===building_genomicranges_target_variants===
#(omit)
#get(paste(target_case, "_input_targeted_vcf_file", sep=""))
#
#targeted_gene <- GRanges(
#  #=====(necessary parameters)#
#  seqnames = as.vector(seqnames(BAEgene[target_gene_rowcount, ])),
#  ranges = IRanges(start = start(ranges(BAEgene[target_gene_rowcount, ])), end = end(ranges(BAEgene[target_gene_rowcount, ]))),
#  #=====(unnecessary parameters)#
#  symbol = BAEgene[target_gene_rowcount, ]$symbol,
#  dominant = BAEgene[target_gene_rowcount, ]$dominant,
#  snp_position = BAEgene[target_gene_rowcount, ]$snp_position,
#  note = BAEgene[target_gene_rowcount, ]$note   
#)
##=====================

##===data_import_and_preprocess_(Srtructure_Variants_VCF)===
SV_vcf_files_input <- function(caseID){
  print(caseID)
  
  filelist_target <- as.data.frame(list.files(dir_input_workingspace, pattern=".vcf"))
  colnames(filelist_target) = "file"
  filename=paste(caseID, "_Tumor_WGS_TEinsertions",".vcf",sep="")
  targeted_vcf_file_path <- paste(dir_input_epigenetics, filename, sep="/")
  targeted_vcf_file_Tumor_SV <- readVcf(targeted_vcf_file_path)
  
  filelist_target <- as.data.frame(list.files(dir_input_workingspace, pattern=".vcf"))
  colnames(filelist_target) = "file"
  filename=paste(caseID, "_Blood_WGS_TEinsertions",".vcf",sep="")
  targeted_vcf_file_path <- paste(dir_input_epigenetics, filename, sep="/")
  targeted_vcf_file_Blood_SV <- readVcf(targeted_vcf_file_path)
  
  ERVcaller_WGS_vcf_files <- list("Tumor"=targeted_vcf_file_Tumor_SV, "Blood"=targeted_vcf_file_Blood_SV)
  return(ERVcaller_WGS_vcf_files)
}

print(target_case)
assign(paste(target_case, "SV_vcf_files", sep=""), SV_vcf_files_input(target_case))
#View(get(paste(target_case, "SV_vcf_files", sep="")))

##=====================

##===data_import_and_preprocess_(methylation_HM450)===
#$MethylationArray$(https://bioconductor.org/packages/release/workflows/vignettes/methylationArrayAnalysis/inst/doc/methylationArrayAnalysis.html)
#!PS: (SKIP:some files do not have DNAMethylation_HM450_solidnormal)
library("tidyr")
library("dplyr")

print("methylation_HM450")  
methylation_HM450_processing <- function(caseID){
  #===Data import===  
  print(caseID)
  filelist_target <- as.data.frame(list.files(dir_input_epigenetics, pattern="HumanMethylation450array.txt"))
  colnames(filelist_target) = "file"
  
  filename_Tumor=paste(caseID, "_Tumor_HumanMethylation450array",".txt",sep="")
  targeted_file_Tumor_path <- paste(dir_input_epigenetics, filename_Tumor, sep="/")
  targeted_case_HM450_tumor <- read.table(targeted_file_Tumor_path, sep = '\t', header= TRUE, na = "NA", stringsAsFactors = F)

  filename_Solidnormal=paste(caseID, "_Solidnormal_HumanMethylation450array",".txt",sep="")
  targeted_file_Solidnormal_path <- paste(dir_input_epigenetics, filename_Solidnormal, sep="/")
  targeted_case_HM450_solidnormal <- read.table(targeted_file_Solidnormal_path, sep = '\t', header= TRUE, na = "NA", stringsAsFactors = F)
  
  #==Density plot for Beta Value==
  targeted_case_HM450_tumor_d_Betavalue <- density(na.omit(targeted_case_HM450_tumor$Beta_value)) # returns the density data
  png(filename=paste(dir_temp,"/", caseID, "_HM450_tumor_Betavalue.png", sep=""))
  plot.new()
  plot(targeted_case_HM450_tumor_d_Betavalue)
  dev.off()
  
  targeted_case_HM450_solidnormal_d_Betavalue <- density(na.omit(targeted_case_HM450_solidnormal$Beta_value)) # returns the density data
  png(filename=paste(dir_temp,"/", caseID, "_HM450_solidnormal_Betavalue.png", sep=""))
  plot.new()
  plot(targeted_case_HM450_solidnormal_d_Betavalue)
  dev.off()
  
  
  targeted_case_HM450_tumor <- mutate(targeted_case_HM450_tumor, M_value=as.numeric(log2((Beta_value+0.0001)/((1-Beta_value)+0.0001))))
  targeted_case_HM450_tumor_d_Mvalue <- density(na.omit(targeted_case_HM450_tumor$M_value)) 
  png(filename=paste(dir_temp,"/", caseID, "_HM450_tumor_Mvalue.png", sep=""))
  plot.new()
  plot(targeted_case_HM450_tumor_d_Mvalue)
  dev.off()
  
  targeted_case_HM450_solidnormal <- mutate(targeted_case_HM450_solidnormal, M_value=as.numeric(log2((Beta_value+0.0001)/((1-Beta_value)+0.0001))))
  targeted_case_HM450_solidnormal_d_Mvalue <- density(na.omit(targeted_case_HM450_solidnormal$M_value)) 
  png(filename=paste(dir_temp,"/", caseID, "_HM450_solidnormal_Mvalue.png", sep=""))
  plot.new()
  plot(targeted_case_HM450_solidnormal_d_Mvalue)
  dev.off()
  
  
  targeted_case_HM450_data <- list("tumor_HM450_rawdata"=targeted_case_HM450_tumor, "tumor_d_Betavalue"=targeted_case_HM450_tumor_d_Betavalue, "tumor_d_Mvalue"=targeted_case_HM450_tumor_d_Mvalue, "solidnormal_HM450_rawdata"=targeted_case_HM450_solidnormal, "solidnormal_d_Betavalue"=targeted_case_HM450_solidnormal_d_Betavalue, "solidnormal_d_Mvalue"=targeted_case_HM450_solidnormal_d_Mvalue)
  return(targeted_case_HM450_data)
}

assign(paste(target_case, "_methylation_HM450_processing", sep=""), methylation_HM450_processing(target_case))
#View(`TCGA-44-6146_methylation_HM450_processing`)
#View(get(paste(target_case, "_methylation_HM450_processing", sep="")))



