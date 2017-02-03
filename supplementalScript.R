#probe statistics
library(tidyverse)
probe_list <- read_csv("PROBELIST.csv")
length(unique(probe_list %>% filter(group != 0) %>% .$TARGET_NAME))
length(unique(probe_list %>% filter(group == 0) %>% .$TARGET_NAME))
length(unique(probe_list %>% filter(group != 0) %>% .$UNIPROTKB))
length(unique(probe_list %>% filter(group != 0) %>% .$CHEMICAL_ID))
table(probe_list %>% filter(group != 0) %>% .$is_agonist)
table(probe_list %>% filter(group != 0) %>% .$orthogonal)
median(top_probes_onetarget$n_total)
median(top_probes_onetarget$n_select)

unique(intersect(probe_list %>% filter(group != 0) %>% .$TARGET_NAME, probe_list %>% filter(group == 0) %>%
                   .$TARGET_NAME))

# summary_table <- tibble(Source = c("CHEMBL", "chemicalprobes.org", "Total"), 
#                         Probes = c(717,110, 827), Targets = c(367, 117, 444))
pdf("summary_table.pdf", height=11, width=8.5)
grid.table(summary_table)
dev.off()

#-------------------------    GO ANALYSIS   --------------------------------------

#could find overlap between all gene go membership and probe gene go membership
#get all gene names----------
library(UniProt.ws)
library(clusterProfiler)
library(org.Hs.eg.db)
library(tibble)
library(dplyr)
up <- UniProt.ws(taxId=9606)
egs <- as.character(unique(keys(up, "ENTREZ_GENE")))
gene <- as.character(unique(probe_list$ENTREZ_GENE))

#enter either 3 or 4 for the level on the groupGO to see pathways
ggo_all <- groupGO(gene     = egs,
                   OrgDb    = "org.Hs.eg.db",
                   ont      = "MF",
                   level = 4,
                   readable = TRUE)

ggo_all_df <- ggo_all@result %>% dplyr::select(ID, Description, Count) %>% 
  dplyr::rename(count_all = Count)

ggo_probes <- groupGO(gene     = gene,
                      OrgDb    = "org.Hs.eg.db",
                      ont      = "MF",
                      level = 4,
                      readable = TRUE)

ggo_probes_df <- ggo_probes@result %>% dplyr::select(ID, Description, Count) %>% 
  dplyr::rename(count_probes = Count)

ggo_df <- left_join(ggo_all_df, ggo_probes_df)
ggo_df <- ggo_df %>% mutate(ratio = count_probes / count_all) %>% arrange(desc(ratio), desc(count_all))
ggo_df %>% filter(count_all > 50) %>% top_n(10)
head(ggo_df %>% filter(count_probes < 2) %>% arrange(desc(count_all)))
#### should make a table that shows the above two results for MF level 3, MF level 4, kegg
#### and shows the go description
analysis_table <- tibble(Source = as.character(c("GO MF Level 3", "GO MF Level 4", "KEGG")), 
                         `Most enriched terms` = as.character(c(
                                                        "hormone binding; 
                                                        amide binding;
                                                        receptor activity; 
                                                        channel regulator activity",
                                                       "peptide receptor activity; 
                                                        steroid hormone receptor activity;
                                                        receptor signaling protein serine/threonine kinase activity; 
                                                        p53 binding; 
                                                        peptide binding", 
                                                       "Neuroactive ligand receptor interaction; 
                                                        Calcium signaling pathway; 
                                                        cAMP signaling pathway; 
                                                        Pathways in cancer;
                                                        FoxO signaling pathway")), 
                         `Least enriched terms` = as.character(c(
                                                        "structural constituent of ribosome; 
                                                        lyase activity; 
                                                        structural constituent of cytoskeleton;
                                                        odorant binding;
                                                        carbohydrate binding", 
                                                        "hydrolase activity, acting on glycosyl bonds;
                                                        ubiquitinyl hydrolase activity;
                                                        phosphatase regulator activity;
                                                        Ras guanyl-nucleotide exchange factor activity", 
                                                        "Proteasome;
                                                        Insulin secretion;
                                                        RIG-I-like receptor signaling pathway;
                                                        Cytokine-cytokine receptor interaction")))
analysis_table$`Most enriched terms` <- stringr::str_wrap(analysis_table$`Most enriched terms`, width = 40)
analysis_table$`Least enriched terms` <- stringr::str_wrap(analysis_table$`Least enriched terms`, width = 40)


jpeg("analysis_table.jpeg", height = 400, width = 800, units = "px", quality = 100)
gridExtra::grid.arrange(gridExtra::tableGrob(analysis_table))
dev.off()

pdf("analysis_table.pdf")
gridExtra::grid.table(analysis_table)
dev.off()

#get ggo_df for above
x <- ggo_df %>% filter(count_all > 50) %>% slice(1:5)
analysis_table[2,2] <- str_c(as.character(x$Description), collapse = "; ")
x <- ggo_df %>% filter(count_all > 50 & count_probes == 0) %>% slice(1:5)
analysis_table[2,3] <- str_c(as.character(x$Description), collapse = "; ")


#### basic go enrichment plot
ggo_probes_df <- ggo_probes_df %>% dplyr::rename(Count = count_probes) %>% arrange(desc(Count))
ggo_probes_df$Description <- stringr::str_wrap(ggo_probes_df$Description, width = 20)
g <- ggplot(data = ggo_probes_df %>% slice(2:7), aes(x = reorder(Description, Count), y = Count, fill = Count), 
       colorRamp("blue")) + coord_flip() + xlab("") +
  geom_bar(stat="identity") + ggtitle("Top GO MF terms for targets") + 
  theme(legend.position="none", text = element_text(size=15), 
        title = element_text(hjust = .5))
ggsave("gomfbarallprobes.png", plot = g, dpi = 300, height = 4, width = 4)


tail(ggo_df %>% filter(count_all > 50))
#any correlation between number of probes and number of genes
cor(ggo_df$count_all, ggo_df$count_probes)



barplot(ggo_probes, drop=TRUE, showCategory=8)


kk <- enrichKEGG(gene         = gene,
                 organism     = 'hsa',
                 pvalueCutoff = 1)
barplot(kk, drop=TRUE, showCategory=15)
head(kk@result)
tail(kk@result)
analysis_table[3,2] <- str_c(as.character(kk@result$Description[1:5]), collapse = "; ")
analysis_table[3,3] <- str_c(as.character(kk@result$Description[(nrow(kk@result)-4):nrow(kk@result)]), collapse = "; ")

pdf("analysis_table.pdf", height=11, width=8.5)
grid.table(analysis_table)
dev.off()

####---------------- histograms of selectivity
hist_plots <- bind_rows(top_probes_onetarget %>% select(molregno, n_total, n_select), top_probes_twotarget %>% select(molregno, n_total, n_select))
hist_plots <- hist_plots %>% ungroup() %>% select(molregno, n_total, n_select) %>% 
  rename(`Selectivity measurements` = n_select, `Total activities` = n_total)
hist_plots <- hist_plots[!(duplicated(hist_plots)),]
hist_plots$`Total activities`[hist_plots$`Total activities`>100] <- 100
hist_plots$`Selectivity measurements`[hist_plots$`Selectivity measurements`>50] <- 50
ggplot(data = hist_plots, aes(`Total activities`), colors = "blue") + 
  #stat_count(color = "blue", fill = "light blue") +
  geom_histogram(breaks=c(seq(0, 100, by=5)), color = "blue", fill = "light blue") +
  scale_x_continuous(limits=c(0, 100), breaks=c(seq(0, 100, by=10)), labels=c(seq(0,90, by=10), "100+")) + 
  annotate("text", x = 65, y = 100, label = "Median = 11") + ggtitle("Total activities") + xlab("")
ggsave("totalactivitiesplot.png", dpi = 300)
ggplot(data = hist_plots, aes(`Selectivity measurements`), colors = "blue") + 
  #stat_count(color = "blue", fill = "light blue") +
  geom_histogram(breaks=c(seq(0, 50, by=1)), color = "blue", fill = "light blue") +
  scale_x_continuous(limits=c(0, 50), breaks=c(seq(0, 50, by=5)), labels=c(seq(0,45, by=5), "50+")) + 
  annotate("text", x = 35, y = 100, label = "Median = 4") + ggtitle("Selectivity measurements") + xlab("")
ggsave("selectivityplot.png", dpi = 300)
median(hist_plots$`Selectivity measurements`)
median(hist_plots$`Total activities`)

#figure showing pancreatic cancer kegg graph
#probes found
cat(unique(probe_list$CHEMICAL_ID[grepl("hsa04014", probe_list$KEGG)]), sep = "; ")
length(unique(probe_list$CHEMICAL_ID[grepl("hsa04014", probe_list$KEGG)]))

#targets found
unique(probe_list$TARGET_NAME[grepl("hsa04014", probe_list$KEGG)])
#36
unique(probe_list$TARGET_NAME[grepl("TP53 Activity", probe_list$REACTOME)])
#21















