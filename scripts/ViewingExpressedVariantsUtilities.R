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
  fig_dir <- "../pictures/"
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


generate_sv_viz <- function(tcga_vcf_df){
  #CMUHackathon_visualization_Genometrack_datapreparing.R
  source("../CMUHackathon_visualization_Genometrack_datapreparing.R")
  
  #CMUHackathon_visualization_Genometrack_plottingcore.R
  source("../CMUHackathon_visualization_Genometrack_plottingcore.R")
  
}




## Pathway graphs
library(pathfindR)
#data = read.csv("xyz.csv")
#head(data,15)
#tail(data,15)

pathway_viz <- function(pathway_data){
  output_df <- run_pathfindR(pathway_data)
  #knitr::kable(head(output_df, 2))
  fig_dir <- "../pictures/"
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

## Return html output
#Code obtained from stackoverflow
html_output <- function(Df,meta = NULL, cacheable = NA) {
  rmarkdown::render('./report.rmd',params=list(output_file = report.html))
} 

## Create summary report


summary_report <- function(type){
  if(type == "TCGA"){
    ##Enter code for R markdown for the report for both epigenetics and expressed variants
    #source("../markdown/report_SV.Rmd")
    TCGA_report <- rmarkdown::render('../markdown/report_SV.Rmd',params=list(output_file = report_SV.html))
    return(TCGA_report)
  } else if (type == "Expressed Variants only"){
    ##Enter code for R markdown for the report for expressed variants
    
    #source("../markdown/report.Rmd")
    expressed_variants_report <- rmarkdown::render('./report.Rmd',params=list(output_file = report.html))
    return(expressed_variants_report)
  }
}
