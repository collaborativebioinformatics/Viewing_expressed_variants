#CMUHackathon_visualization_Genometrack_plottingcore.R

####＝＝＝＝＝＝＝＝＝＝Data Visualization＝＝＝＝＝＝＝＝＝＝


##===Plotting_core===

#--Loading_additional_annotation_track_(additionaltrack)--#
ah <- AnnotationHub()
query(ah, c("Homo sapien", "CTCF", "hepG"))
id <- names(query(ah, "wgEncodeUwTfbsHepg2CtcfStdPkRep2.narrowPeak.gz")) 
Hepg2Ctcf.gr <- ah[[tail(id, 1)]]

path = system.file(package="liftOver", "extdata", "hg38ToHg19.over.chain")
ch = import.chain(path)
seqlevelsStyle(Hepg2Ctcf.gr) = "UCSC"
Hepg2Ctcf.gr_hg38 = liftOver(Hepg2Ctcf.gr, ch)
Hepg2Ctcf.gr_hg38 <- unlist(Hepg2Ctcf.gr_hg38)

gene_range <- read.table(paste(dir_input, "GRCh38_hg38_refFlat_annotation_primaryAssemblyOnly_NoXY.bed", sep="/"), header = FALSE)
gene_range <- dplyr::rename(gene_range, chr = V1, start = V2, end = V3, ID = V4)
gene_range_filter <- gene_range %>% group_by(ID) %>% filter(row_number() == 1) %>% ungroup()
gene_range_filter <- as.data.frame(gene_range_filter)
#---------------------------------------------------------#

genomeTracksGraphic_targetlocations <- function(target_case_in, temp_variants_in, count00_in, temp_targetgene_in){
  #target_case_in=target_case; temp_variants_in=temp_variants; count00_in=count00; temp_targetgene_in=temp_targetgene
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
    #target_CHR = as.character( seqnames(granges(temp_variants_in[count00_in])) )
    #target_START = as.numeric( start(granges(temp_variants_in[count00_in])) ) - 1000000
    #target_END = as.numeric( start(granges(temp_variants_in[count00_in])) ) + 1000000
    
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
    dtrack_Tumor_ATAC_raw <- DataTrack(range = paste(dir_input, "/", target_case_in, "_Tumor_atacReads_raw",".bw", sep=""), name="Tumor_ATAC_raw")
    
    #Methylation_tumor
    targeted_case_HM450_tumor <- get(paste(target_case, "_methylation_HM450_processing", sep=""))$tumor_HM450_rawdata
    targeted_case_HM450_tumor <- targeted_case_HM450_tumor[targeted_case_HM450_tumor$CpG_chrm==target_CHR, ]
    targeted_case_HM450_tumor <- targeted_case_HM450_tumor[targeted_case_HM450_tumor$CpG_beg >= target_START, ]
    targeted_case_HM450_tumor <- targeted_case_HM450_tumor[targeted_case_HM450_tumor$CpG_beg <= target_END, ]
    targeted_case_HM450_tumor = targeted_case_HM450_tumor %>% drop_na(Beta_value)
    targeted_case_HM450_tumor_GRanges <- GRanges(
      #=====(necessary parameters)#
      seqnames = as.vector(targeted_case_HM450_tumor$CpG_chrm),
      ranges = IRanges(start = targeted_case_HM450_tumor$CpG_beg, end = targeted_case_HM450_tumor$CpG_end),
      #=====(unnecessary parameters)#
      Beta_value = targeted_case_HM450_tumor$Beta_value,  
      M_value = targeted_case_HM450_tumor$M_value
    )
    dtrack_tumor_methylation <- DataTrack(targeted_case_HM450_tumor_GRanges, type="histogram", name="tumor_met.")
    
    #Methylation_solidnormal
    targeted_case_HM450_solidnormal <- get(paste(target_case, "_methylation_HM450_processing", sep=""))$solidnormal_HM450_rawdata
    targeted_case_HM450_solidnormal <- targeted_case_HM450_solidnormal[targeted_case_HM450_solidnormal$CpG_chrm==target_CHR, ]
    targeted_case_HM450_solidnormal <- targeted_case_HM450_solidnormal[targeted_case_HM450_solidnormal$CpG_beg >= target_START, ]
    targeted_case_HM450_solidnormal <- targeted_case_HM450_solidnormal[targeted_case_HM450_solidnormal$CpG_beg <= target_END, ]
    targeted_case_HM450_solidnormal = targeted_case_HM450_solidnormal %>% drop_na(Beta_value)
    targeted_case_HM450_solidnormal_GRanges <- GRanges(
      #=====(necessary parameters)#
      seqnames = as.vector(targeted_case_HM450_solidnormal$CpG_chrm),
      ranges = IRanges(start = targeted_case_HM450_solidnormal$CpG_beg, end = targeted_case_HM450_solidnormal$CpG_end),
      #=====(unnecessary parameters)#
      Beta_value = targeted_case_HM450_solidnormal$Beta_value,  
      M_value = targeted_case_HM450_solidnormal$M_value
    )
    dtrack_solidnormal_methylation <- DataTrack(targeted_case_HM450_solidnormal_GRanges, type="histogram", name="solidnormal_met.")
    
    #RNAseq_read_coverage
    #input_RNAseq_tumor_BAM_path = paste(dir_input, "/", target_case_in, "_Tumor_RNAseq",".bam", sep="")
    #alTrack_RNAseq_tumor <- AlignmentsTrack(input_RNAseq_tumor_BAM_path, genome=genome(target_gene_GRanges), chromosome=as.vector(target_CHR), start=target_START, end=target_END, name="RNAseq", isPaired=TRUE, mapq=20)
    
    #Plotting_core
    png(filename=paste(dir_pictures, "/", target_case_in, "_", target_gene_GRanges$symbol, "_Epigenetic_plotting.png", sep=""), width = 2000, height = 1200, res = 90)
    plot.new()
    #plotTracks(list(itrack, gtrack, atrack, grtrack_Tumor_TE, grtrack_Blood_TE, alTrack_RNAseq_tumor, dtrack_Tumor_ATAC_open, dtrack_tumor_methylation, dtrack_solidnormal_methylation, additionaltrack), from=target_START, to=target_END) #(with_RNAseq_track) 
    plotTracks(list(itrack, gtrack, atrack, grtrack_Tumor_TE, grtrack_Blood_TE, dtrack_Tumor_ATAC_raw, dtrack_tumor_methylation, dtrack_solidnormal_methylation, additionaltrack), from=target_START, to=target_END, background.title = "darkblue", title.width=NULL, main=paste(target_case_in, "_", temp_targetgene_in, sep=""))
    dev.off()
  }else{
    return(NULL)
  }
}
##=====================

##==call_plotting_function==
temp_variants <- get(paste(target_case, "_input_targeted_vcf_file", sep=""))

for(count00 in seq(1,nrow(temp_variants),1)){
  print("Data Visualization")  
  print(count00)
  temp_targetgene=strsplit(info(temp_variants)$CSQ[[count00]][1], "\\|")[[1]][4]
  print(temp_targetgene)
  if(temp_targetgene != ""){
    assign(paste(target_case, "_", temp_targetgene, "_genomeTracksGraphic_targetlocations", sep=""), genomeTracksGraphic_targetlocations(target_case, temp_variants, count00, temp_targetgene))
    print("plotting!!!")
  }
  #get(paste(target_case, "_", visualizing_gene, "_genomeTracksGraphic_Epigenetic_plotting", sep=""))
}
##=====================