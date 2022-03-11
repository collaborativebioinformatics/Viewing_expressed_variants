library(gprofiler2)
# grabs from data output, check this path w group?
path_data <- read.csv("../data_output/pathogenic_variants_try.csv")
gene_ids <- as.list(path_data[,5])

gostres <- gost(query = gene_ids, 
                organism = "hsapiens", sources = c("KEGG","REAC"))

ens_id <- data.frame(gostres$meta$genes_metadata$query)
ens_id <- data.frame(t(ens_id))
ens_id <- unique(ens_id)
ens_id$geneInfo <- rownames(ens_id)
ens_id_data <- gostres$result

spl <- t(data.frame(strsplit(ens_id$geneInfo, split='.', fixed=TRUE)))
colnames(spl) <- c("queryID","geneName")
ens_id_df <- cbind(resTranspose, spl)

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



