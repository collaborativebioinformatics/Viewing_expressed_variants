#CMUHackathon_visualization_Genometrack_datapreparing.R

####＝＝＝＝＝＝＝＝＝＝Environments setting＝＝＝＝＝＝＝＝＝＝
target_case="TCGA_44_6146"
server="/Users/chunhsuanlojason/Desktop/CMU_Libraries_Hackathon" 

library("icesTAF")
library("AnnotationHub")
library("VariantAnnotation")
library("Rsamtools") 
library("GenomicAlignments") 
library("rtracklayer")
library("icesTAF")
library("magrittr")
library("Gviz")
library("GenomicRanges")
library("rtracklayer")
#library("liftOver")
library("tidyr")
library("dplyr")

dir_input=paste(server, "data_raw", sep="/")
mkdir(dir_input)
dir_output=paste(server, "data_output", sep="/")
mkdir(dir_output)
dir_pictures=paste(server, "pictures", sep="/")
mkdir(dir_pictures)

####＝＝＝＝＝＝＝＝＝＝Data importing & initializing: VariantsSites & ATAC-Seq & methylation_HM450 & RNAseq＝＝＝＝＝＝＝＝＝＝
##===data_import_VCF===
input_targeted_vcf_file <- function(caseID){
  print(caseID)
  filelist_target <- as.data.frame(list.files(dir_input, pattern=".vcf"))
  colnames(filelist_target) = "file"
  filename=paste(caseID, "_WES_somaticvariants",".vcf",sep="")
  
  targeted_vcf_file_path <- paste(dir_input, filename, sep="/")
  targeted_vcf_file <- readVcf(targeted_vcf_file_path)
  
  return(targeted_vcf_file)
}

assign(paste(target_case, "_input_targeted_vcf_file", sep=""), input_targeted_vcf_file(paste(target_case, sep="")))
#View(`TCGA_44_6146_input_targeted_vcf_file`)
#View(get(paste(target_case, "_input_targeted_vcf_file", sep="")))
##=====================


##===data_import_and_preprocess_(Srtructure_Variants_VCF)===
#TO detect TE insertiion sites by ERV_caller: https://github.com/xunchen85/ERVcaller

SV_vcf_files_input <- function(caseID){
  print(caseID)
  
  filelist_target <- as.data.frame(list.files(dir_input, pattern=".vcf"))
  colnames(filelist_target) = "file"
  filename=paste(caseID, "_Tumor_WGS_TEinsertions",".vcf",sep="")
  targeted_vcf_file_path <- paste(dir_input, filename, sep="/")
  targeted_vcf_file_Tumor_SV <- readVcf(targeted_vcf_file_path)
  
  filelist_target <- as.data.frame(list.files(dir_input, pattern=".vcf"))
  colnames(filelist_target) = "file"
  filename=paste(caseID, "_Blood_WGS_TEinsertions",".vcf",sep="")
  targeted_vcf_file_path <- paste(dir_input, filename, sep="/")
  targeted_vcf_file_Blood_SV <- readVcf(targeted_vcf_file_path)
  
  ERVcaller_WGS_vcf_files <- list("Tumor"=targeted_vcf_file_Tumor_SV, "Blood"=targeted_vcf_file_Blood_SV)
  return(ERVcaller_WGS_vcf_files)
}

print(target_case)
assign(paste(target_case, "SV_vcf_files", sep=""), SV_vcf_files_input(target_case))
#View(get(paste(target_case, "SV_vcf_files", sep="")))
##=====================


##===data_import_and_preprocess_(methylation_HM450)===
#TCGAportal for downloading methylatiion data (HM450_Beta_value): https://portal.gdc.cancer.gov/repository?facetTab=files&filters=%7B%22op%22%3A%22and%22%2C%22content%22%3A%5B%7B%22content%22%3A%7B%22field%22%3A%22cases.submitter_id%22%2C%22value%22%3A%5B%22TCGA-44-6146%22%5D%7D%2C%22op%22%3A%22in%22%7D%2C%7B%22op%22%3A%22in%22%2C%22content%22%3A%7B%22field%22%3A%22files.data_format%22%2C%22value%22%3A%5B%22txt%22%5D%7D%7D%2C%7B%22op%22%3A%22in%22%2C%22content%22%3A%7B%22field%22%3A%22files.experimental_strategy%22%2C%22value%22%3A%5B%22Methylation%20Array%22%5D%7D%7D%5D%7D
#$MethylationArray$(https://bioconductor.org/packages/release/workflows/vignettes/methylationArrayAnalysis/inst/doc/methylationArrayAnalysis.html)
#!PS: (SKIP:some files do not have DNAMethylation_HM450_solidnormal)

print("methylation_HM450")  
methylation_HM450_processing <- function(caseID){
  #===Data import===  
  print(caseID)
  filelist_target <- as.data.frame(list.files(dir_input, pattern="HumanMethylation450array.txt"))
  colnames(filelist_target) = "file"
  
  filename_Tumor=paste(caseID, "_Tumor_HumanMethylation450array",".txt",sep="")
  targeted_file_Tumor_path <- paste(dir_input, filename_Tumor, sep="/")
  targeted_case_HM450_tumor <- read.table(targeted_file_Tumor_path, sep = '\t', header= TRUE, na = "NA", stringsAsFactors = F)
  colnames(targeted_case_HM450_tumor) <- c("probeID","Beta_value")                                                                                   #nrow(targeted_case_HM450_tumor): 485576
  
  filename_Solidnormal=paste(caseID, "_Solidnormal_HumanMethylation450array",".txt",sep="")
  targeted_file_Solidnormal_path <- paste(dir_input, filename_Solidnormal, sep="/")
  targeted_case_HM450_solidnormal <- read.table(targeted_file_Solidnormal_path, sep = '\t', header= TRUE, na = "NA", stringsAsFactors = F)
  colnames(targeted_case_HM450_solidnormal) <- c("probeID","Beta_value")                                                                             #nrow(targeted_case_HM450_solidnormal): 485576
  
  #==Coverting probe ID into genome coordinates (hg38)==
  #To download "Basic manifest with mapping information - hg38" from :http://zwdzwd.github.io/InfiniumAnnotation
  HM450_hg38_manifest <- read.table(paste(dir_input, "/HM450.hg38.manifest.tsv", sep=""), sep = '\t', header= TRUE, na = "NA", stringsAsFactors = F) #nrow(HM450_hg38_manifest): 485577
  
  targeted_case_HM450_tumor <- left_join(targeted_case_HM450_tumor, HM450_hg38_manifest, by="probeID")                                               #nrow(temp): 485576
  targeted_case_HM450_solidnormal <- left_join(targeted_case_HM450_solidnormal, HM450_hg38_manifest, by="probeID")                                   #nrow(temp): 485576
  
  #==Density plot for Beta Value==
  targeted_case_HM450_tumor_d_Betavalue <- density(na.omit(targeted_case_HM450_tumor$Beta_value)) # returns the density data
  png(filename=paste(dir_pictures,"/", caseID, "_HM450_tumor_Betavalue.png", sep=""))
  plot.new()
  plot(targeted_case_HM450_tumor_d_Betavalue, main=paste(caseID, "_HM450_tumor_Betavalue", sep=""))
  dev.off()
  
  targeted_case_HM450_solidnormal_d_Betavalue <- density(na.omit(targeted_case_HM450_solidnormal$Beta_value)) # returns the density data
  png(filename=paste(dir_pictures,"/", caseID, "_HM450_solidnormal_Betavalue.png", sep=""))
  plot.new()
  plot(targeted_case_HM450_solidnormal_d_Betavalue, main=paste(caseID, "_HM450_solidnormal_Betavalue", sep=""))
  dev.off()
  
  
  targeted_case_HM450_tumor <- mutate(targeted_case_HM450_tumor, M_value=as.numeric(log2((Beta_value+0.0001)/((1-Beta_value)+0.0001))))
  targeted_case_HM450_tumor_d_Mvalue <- density(na.omit(targeted_case_HM450_tumor$M_value)) 
  png(filename=paste(dir_pictures,"/", caseID, "_HM450_tumor_Mvalue.png", sep=""))
  plot.new()
  plot(targeted_case_HM450_tumor_d_Mvalue, main=paste(caseID, "_HM450_tumor_Mvalue", sep=""))
  dev.off()
  
  targeted_case_HM450_solidnormal <- mutate(targeted_case_HM450_solidnormal, M_value=as.numeric(log2((Beta_value+0.0001)/((1-Beta_value)+0.0001))))
  targeted_case_HM450_solidnormal_d_Mvalue <- density(na.omit(targeted_case_HM450_solidnormal$M_value)) 
  png(filename=paste(dir_pictures,"/", caseID, "_HM450_solidnormal_Mvalue.png", sep=""))
  plot.new()
  plot(targeted_case_HM450_solidnormal_d_Mvalue, main=paste(caseID, "_HM450_solidnormal_Mvalue", sep=""))
  dev.off()
  
  
  targeted_case_HM450_data <- list("tumor_HM450_rawdata"=targeted_case_HM450_tumor, "tumor_d_Betavalue"=targeted_case_HM450_tumor_d_Betavalue, "tumor_d_Mvalue"=targeted_case_HM450_tumor_d_Mvalue, "solidnormal_HM450_rawdata"=targeted_case_HM450_solidnormal, "solidnormal_d_Betavalue"=targeted_case_HM450_solidnormal_d_Betavalue, "solidnormal_d_Mvalue"=targeted_case_HM450_solidnormal_d_Mvalue)
  return(targeted_case_HM450_data)
}

assign(paste(target_case, "_methylation_HM450_processing", sep=""), methylation_HM450_processing(target_case))
#View(`TCGA-44-6146_methylation_HM450_processing`)
#View(get(paste(target_case, "_methylation_HM450_processing", sep="")))
##=====================


