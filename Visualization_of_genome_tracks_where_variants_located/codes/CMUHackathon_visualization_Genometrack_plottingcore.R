#CMUHackathon_visualization_Genometrack_plottingcore.R
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