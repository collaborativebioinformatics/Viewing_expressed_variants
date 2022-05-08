#Utilities script
##Data input

## Hackathon
library(data.table) # for data.table functions
library(dplyr) # for pipe, filter, str_detect

process_vcf_file <- function(vcf_file_path){
  # read in full file, espite variable number of fields at first
  #all_content <- readLines("/Users/saracarioscia/Downloads/testSample.cancer.vcf")
  all_content <- readLines(vcf_file_path)
  # skip until we get to the fixed area header
  variant_rows <- all_content[-c(1:grep("#CHROM", all_content))]
  # read back in the fixed lines as a table
  variants <- read.table(textConnection(variant_rows))
  
  # get only pathogenic variants
  variants_dt <- variants %>% as.data.table()
  pathogenic_variants <- variants_dt %>%
    filter_all(any_vars(str_detect(., pattern = "PATHOGENIC")))
  # write it back as a csv
  write.csv(pathogenic_variants, "~/Downloads/pathogenic_variants.csv", row.names = FALSE, quote = FALSE)
  
  # Format is `|GENE=gene_name|`. Try to grab the content between GENE= and |
  for (i in 1:nrow(pathogenic_variants)) {
    # go through each row of the pathogenic variants
    row <- pathogenic_variants[i,]
    # grab the content after GENE=
    gene <- str_match(row, "GENE=\\s*(.*?)\\s*;")
    # other rows are grabbed; we want those that are not NA
    gene <- gene[,2][!is.na(gene[,2])]
    # assign the gene to the column in pathogenic_variants
    pathogenic_variants$genes[i] <- gene
  }
  
}
## Visualization of genes

library(gprofiler2)
library(ggplot2)
library(gridExtra)
library(data.table)
library(gridExtra)
# grabs from data output, check this path w group?
path_data <- read.csv("../data_output/pathogenic_variants_try.csv")


expressed_variants_viz <- function(path_data){
  #saves the table for top 5 variants
  fig_dir <- "../figures/"
  file_path1 <- paste0(fig_dir, "top_variants_table.png")
  png(file_path1, height = 50*nrow(head(path_data)), width = 200*ncol(path_data))
  grid.table(head(path_data))
  dev.off()
  
  gene_ids <- as.list(path_data[,5])
  
  # query for KEGG + REAC
  gostres <- gost(query = gene_ids, 
                  organism = "hsapiens", sources = c("KEGG","REAC"))
  
  ens_id <- data.frame(gostres$meta$genes_metadata$query)
  ens_id <- data.frame(t(ens_id))
  ens_id <- unique(ens_id)
  ens_id$geneInfo <- rownames(ens_id)
  ens_id_data <- gostres$result
  
  spl <- t(data.frame(strsplit(ens_id$geneInfo, split='.', fixed=TRUE)))
  colnames(spl) <- c("queryID","geneName")
  ens_id_df <- cbind(ens_id, spl)
  
  gene_pathway_association <- dplyr::left_join(ens_id_data, ens_id_df, by = c("query" = "queryID"))
  gene_pathway_association %>% ggplot2::ggplot(aes(x = term_name, y = geneName, fill =term_name)) +
    geom_tile(show.legend = F) +
    theme_bw() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.text.x = element_text(angle = 90))+
    labs(title = "Genes and Pathway association",x = "Pathways", y = "Genes")
  ggsave("../pictures/gene_pathway_association.png")
  
  
  
  pathdf <- gene_pathway_association %>% dplyr::group_by(term_name) %>% dplyr::tally()
  pathdf %>% ggplot2::ggplot(aes(x = term_name, y = n, fill = term_name, label = n)) +
    geom_bar(width = 1, stat = "identity", color = "white") +
    geom_text()+
    theme_bw()+
    labs(fill = "Pathways", title = "Common pathways")+
    theme(axis.title = element_blank(),
          axis.line = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks = element_blank(),
          axis.text.y = element_blank())+
    coord_polar()
  ggsave("../pictures/common_pathways.png")
  
  # query for GO
  gostres <- gost(query = gene_ids, 
                  organism = "hsapiens", sources = "GO")
  
  ens_id <- data.frame(gostres$meta$genes_metadata$query)
  ens_id <- data.frame(t(ens_id))
  ens_id <- unique(ens_id)
  ens_id$geneInfo <- rownames(ens_id)
  ens_id_data <- gostres$result
  
  spl <- t(data.frame(strsplit(ens_id$geneInfo, split='.', fixed=TRUE)))
  colnames(spl) <- c("queryID","geneName")
  ens_id_df <- cbind(ens_id, spl)
  
  GO_association <- dplyr::left_join(ens_id_data, ens_id_df, by = c("query" = "queryID"))
  GO_association %>% ggplot2::ggplot(aes(x = term_name, y = geneName, fill =term_name)) +
    geom_tile(show.legend = F) +
    theme_bw() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.text.x = element_text(angle = 90))+
    labs(title = "Genes and Gene Ontology association",x = "Gene Ontology", y = "Genes")
  ggsave("../pictures/gene_GO_association.png")
  
  
  ## Visualization of pathways
  
  # Read in df
  df <- read.csv("~/pathogenic_variants.csv", header = T)
  colnames(df) <- c("CHR","POS","REF","ALT","FILTER","VariantImpact","GENE")
  
  # Plot variant distribution figure
  
  df %>% ggplot2::ggplot(aes(x=VariantImpact, fill = VariantImpact))+
    geom_bar()+
    theme(axis.text.x = element_text(angle = 90)) +
    geom_text(aes(label = ..count..), stat="count",vjust = -0.5, colour = "black")+
    labs(title = "Pathogenic variants impact",
         x = "Variant Type",
         y = "Number of variants", 
         fill = "Variant Impact")
  ggsave("../pictures/pathogenic_variants_impact.png")
  
  df %>% ggplot2::ggplot(aes(x=VariantImpact, y = GENE, fill = VariantImpact, colour - "white"))+
    geom_tile(show.legend = F)+
    theme(axis.text.x = element_text(angle = 90)) +
    labs(title = "Genes with variants and their impact",
         x = "Variant Impact",
         y = "Genes") 
  ggsave("../pictures/gene-genetic_variant_association.png") 
  
}

###TCGA FILE #####
generate_sv_viz <- function(tcga_df){
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
  
  
  ### tcga genome track ###
  
  ####＝＝＝＝＝＝＝＝＝＝Data Visualization＝＝＝＝＝＝＝＝＝＝
  ##===(Plotting core)===
  library("Rsamtools")
  library("GenomicAlignments")
  library("rtracklayer")
  library("icesTAF")
  library("magrittr")
  library("Gviz")
  library("GenomicRanges")
  library("VariantAnnotation")
  
  #--Loading_additional_annotation_track_(additionaltrack)--#
  library("AnnotationHub")
  ah <- AnnotationHub()
  query(ah, c("Homo sapien", "CTCF", "hepG"))
  id <- names(query(ah, "wgEncodeUwTfbsHepg2CtcfStdPkRep2.narrowPeak.gz"))
  Hepg2Ctcf.gr <- ah[[tail(id, 1)]]
  
  library("rtracklayer")
  #library("liftOver")
  path = system.file(package="liftOver", "extdata", "hg38ToHg19.over.chain")
  ch = import.chain(path)
  seqlevelsStyle(Hepg2Ctcf.gr) = "UCSC"
  Hepg2Ctcf.gr_hg38 = liftOver(Hepg2Ctcf.gr, ch)
  Hepg2Ctcf.gr_hg38 <- unlist(Hepg2Ctcf.gr_hg38)
  #---------------------------------------------------------#
  
  gene_range <- read.table(paste(dir_input_epigenetics, "GRCh38_hg38_refFlat_annotation_primaryAssemblyOnly_NoXY.bed", sep="/"), header = FALSE)
  gene_range <- dplyr::rename(gene_range, chr = V1, start = V2, end = V3, ID = V4)
  gene_range_filter <- gene_range %>% group_by(ID) %>% filter(row_number() == 1) %>% ungroup()
  gene_range_filter <- as.data.frame(gene_range_filter)
  #---------------------------------------------------------#
  
  genomeTracksGrapgic_targetlocations <- function(target_case_in, temp_variants_filtered_in, count00_in, temp_targetgene_in){
    if(nrow(filter(gene_range_filter, ID==temp_targetgene_in))!=0){
      target_gene = filter(gene_range_filter, ID==temp_targetgene_in)[1,]
      target_gene_GRanges <- GRanges(
        #=====(necessary parameters)#
        seqnames = as.vector(target_gene$chr),
        ranges = IRanges(start = as.integer(target_gene$start), end = as.integer(target_gene$end)),
        #=====(unnecessary parameters)#
        symbol = as.character(target_gene$ID)
      )
      genome(target_gene_GRanges) = "hg38"
      
      TumorWGSTE <- get(paste(target_case, "SV_vcf_files", sep=""))$Tumor
      TumorWGSTE_rowRanges <- rowRanges(TumorWGSTE)
      genome(TumorWGSTE_rowRanges) = "hg38"
      
      BloodWGSTE <- get(paste(target_case, "SV_vcf_files", sep=""))$Blood
      BloodWGSTE_rowRanges <- rowRanges(BloodWGSTE)
      genome(BloodWGSTE_rowRanges) = "hg38"
      
      seqlevels(TumorWGSTE_rowRanges, pruning.mode="coarse") <- seqlevels(target_gene_GRanges)
      seqlevels(BloodWGSTE_rowRanges, pruning.mode="coarse") <- seqlevels(target_gene_GRanges)
      target_TumorWGSTE_rowRanges <- TumorWGSTE_rowRanges[(seqnames(TumorWGSTE_rowRanges) == as.character(seqnames(target_gene_GRanges))) & (start(TumorWGSTE_rowRanges) > (as.numeric(start(target_gene_GRanges))-1000000)) & (end(TumorWGSTE_rowRanges) < (as.numeric(end(target_gene_GRanges)) + 1000000))]
      target_BloodWGSTE_rowRanges <- BloodWGSTE_rowRanges[(seqnames(BloodWGSTE_rowRanges) == as.character(seqnames(target_gene_GRanges))) & (start(BloodWGSTE_rowRanges) > (as.numeric(start(target_gene_GRanges))-1000000)) & (end(BloodWGSTE_rowRanges) < (as.numeric(end(target_gene_GRanges)) + 1000000))]
      
      #=====visualization=====
      #target_CHR = as.character( seqnames(granges(temp_variants_filtered_in[count00_in])) )
      #target_START = as.numeric( start(granges(temp_variants_filtered_in[count00_in])) ) - 1000000
      #target_END = as.numeric( start(granges(temp_variants_filtered_in[count00_in])) ) + 1000000
      
      target_CHR = as.character( seqnames(target_gene_GRanges) )
      target_START = as.numeric( start(target_gene_GRanges) ) - 1000000
      target_END = as.numeric( end(target_gene_GRanges) ) + 1000000
      
      atrack <- AnnotationTrack(target_gene_GRanges, name=paste(target_gene_GRanges$symbol,"_somaticmutation",sep=""))
      gtrack <- GenomeAxisTrack()
      itrack <- IdeogramTrack(genome=genome(target_gene_GRanges)[1], chromosome=seqlevels(target_gene_GRanges)[1])
      grtrack_Tumor_TE <- GeneRegionTrack(target_TumorWGSTE_rowRanges, genome(target_gene_GRanges)[1], chromosome=seqlevels(target_gene_GRanges)[1], name="Tumor_TE")
      grtrack_Blood_TE <- GeneRegionTrack(target_BloodWGSTE_rowRanges, genome(target_gene_GRanges)[1], chromosome=seqlevels(target_gene_GRanges)[1], name="Blood_TE")
      additionaltrack <- AnnotationTrack(Hepg2Ctcf.gr_hg38[seqnames(Hepg2Ctcf.gr_hg38)==target_CHR], name="CTCF_sites")
      
      #ATAC_bigwig
      dtrack_Tumor_ATAC_open <- DataTrack(range = paste(dir_input_epigenetics, "/", target_case_in, "_Tumor_atacReads_Open",".bw", sep=""), name="Tumor_ATAC_open")
      
      #Methylation_tumor
      targeted_case_HM450_tumor <- get(paste(target_case, "_methylation_HM450_processing", sep=""))$tumor_HM450_rawdata
      targeted_case_HM450_tumor <- targeted_case_HM450_tumor[targeted_case_HM450_tumor$Chromosome==target_CHR, ]
      targeted_case_HM450_tumor <- targeted_case_HM450_tumor[targeted_case_HM450_tumor$Start >= target_START, ]
      targeted_case_HM450_tumor <- targeted_case_HM450_tumor[targeted_case_HM450_tumor$Start <= target_END, ]
      targeted_case_HM450_tumor_GRanges <- GRanges(
        #=====(necessary parameters)#
        seqnames = as.vector(targeted_case_HM450_tumor$Chromosome),
        ranges = IRanges(start = targeted_case_HM450_tumor$Start, end = targeted_case_HM450_tumor$End),
        #=====(unnecessary parameters)#
        Beta_value = targeted_case_HM450_tumor$Beta_value,  
        M_value = targeted_case_HM450_tumor$M_value
      )
      dtrack_tumor_methylation <- DataTrack(targeted_case_HM450_tumor_GRanges, type="histogram", name="tumor_met.")
      
      #Methylation_solidnormal
      targeted_case_HM450_solidnormal <- get(paste(target_case, "_methylation_HM450_processing", sep=""))$solidnormal_HM450_rawdata
      targeted_case_HM450_solidnormal <- targeted_case_HM450_solidnormal[targeted_case_HM450_solidnormal$Chromosome==target_CHR, ]
      targeted_case_HM450_solidnormal <- targeted_case_HM450_solidnormal[targeted_case_HM450_solidnormal$Start >= target_START, ]
      targeted_case_HM450_solidnormal <- targeted_case_HM450_solidnormal[targeted_case_HM450_solidnormal$Start <= target_END, ]
      targeted_case_HM450_solidnormal_GRanges <- GRanges(
        #=====(necessary parameters)#
        seqnames = as.vector(targeted_case_HM450_solidnormal$Chromosome),
        ranges = IRanges(start = targeted_case_HM450_solidnormal$Start, end = targeted_case_HM450_solidnormal$End),
        #=====(unnecessary parameters)#
        Beta_value = targeted_case_HM450_solidnormal$Beta_value,  
        M_value = targeted_case_HM450_solidnormal$M_value
      )
      dtrack_solidnormal_methylation <- DataTrack(targeted_case_HM450_solidnormal_GRanges, type="histogram", name="solidnormal_met.")
      
      #RNAseq_read_coverage
      #input_RNAseq_tumor_BAM_path = paste(dir_input_epigenetics, "/", target_case_in, "_Tumor_RNAseq_chrosome10&11&12",".bam", sep="")
      #alTrack_RNAseq_tumor <- AlignmentsTrack(input_RNAseq_tumor_BAM_path, genome=genome(target_gene_GRanges), chromosome=as.vector(target_CHR), start=target_START, end=target_END, name="RNAseq", isPaired=TRUE, mapq=20)
      
      #Plotting_core
      png(filename=paste(dir_output, "/", target_case_in, "_", target_gene_GRanges$symbol, "_Epigenetic_plotting.png", sep=""))
      plot.new()
      #plotTracks(list(itrack, gtrack, atrack, grtrack_Tumor_TE, grtrack_Blood_TE, alTrack_RNAseq_tumor, dtrack_Tumor_ATAC_open, dtrack_tumor_methylation, dtrack_solidnormal_methylation, additionaltrack), from=target_START, to=target_END) #(with_RNAseq_track)
      plotTracks(list(itrack, gtrack, atrack, grtrack_Tumor_TE, grtrack_Blood_TE, dtrack_Tumor_ATAC_open, dtrack_tumor_methylation, dtrack_solidnormal_methylation, additionaltrack), from=target_START, to=target_END)
      dev.off()
    }else{
      return(NULL)
    }
  }
  
  #==filtering target sites==
  print(target_case)
  temp_variants <- get(paste(target_case, "_input_targeted_vcf_file", sep=""))
  temp_variants_filtered <- temp_variants[(seqnames(temp_variants) == as.character("chr10"))]
  temp_variants_filtered <- temp_variants_filtered[(start(temp_variants_filtered) > 0) & (end(temp_variants_filtered) < 14000000)]
  #temp_targetgene <- strsplit(info(temp_variants_filtered)$CSQ[[1]][1], "\\|")[[1]][4]
  
  #==call plotting function==
  for(count00 in seq(1,nrow(temp_variants_filtered),1)){
    print("Data Visualization")  
    print(count00)
    temp_targetgene=strsplit(info(temp_variants_filtered)$CSQ[[count00]][1], "\\|")[[1]][4]
    print(temp_targetgene)
    if(temp_targetgene != ""){
      assign(paste(target_case, "_", temp_targetgene, "_genomeTracksGrapgic_targetlocations", sep=""), genomeTracksGrapgic_targetlocations(target_case, temp_variants_filtered, count00, temp_targetgene))
      print("plotting!!!")
    }
    #get(paste(target_case, "_", visualizing_gene, "_genomeTracksGrapgic_Epigenetic_plotting", sep=""))
  }
  
}


## Pathway graphs
library(pathfindR)
#data = read.csv("xyz.csv")
#head(data,15)
#tail(data,15)

pathway_viz <- function(pathway_data){
  output_df <- run_pathfindR(pathway_data)
  #knitr::kable(head(output_df, 2))
  fig_dir <- "../figures/"
  file_path1 <- paste0(fig_dir, "pathway_enrichment_chart.png")
  png(file_path1, width=6, height=5, units="in", res=1200)
  #Graphical summary of enrichment results for top 10 enriched terms
  enrichment_chart(result_df = output_df)
  dev.off()
  
  file_path2 <- paste0(fig_dir, "pathway_gene_graph.png")
  png(file_path2, width=6, height=5, units="in", res=1200)
  #Term-Gene Graph
  term_gene_graph(output_df, num_terms = 10, use_description = TRUE)
  dev.off()
}

## Create summary report

summary_report <- function(type){
  if(type == "TCGA"){
    ##Enter code for R markdown for the report for both epigenetics and expressed variants
    
    return(TCGA_report)
  } else if (type == "Expressed Variants only"){
    ##Enter code for R markdown for the report for expressed variants
    
    
    return(expressed_variants_report)
  }
}
