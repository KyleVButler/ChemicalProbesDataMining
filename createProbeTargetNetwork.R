#This script makes the probe-target network dataset for use in data visualization functions

library(RefNet)
library(dplyr)
library(stringr)
library(tibble)
library(readr)
refnet <- RefNet()
probe_genes <- as.character(unique(probe_list$UNIPROTKB))

interaction_table <- tibble(A = as.character(), B = as.character(), aliasA = as.character(), aliasB = as.character(), 
                            publicationID = as.character(), type = as.character())

for(i in 1:length(probe_genes)) {
    network_table <- interactions(refnet, species="9606", id=probe_genes[i], provider= "IntAct")
    network_table <- dplyr::select(network_table, A, B, aliasA, aliasB, publicationID, type)
    interaction_table <- dplyr::bind_rows(interaction_table, network_table)
}

network_table <- interaction_table
network_table$A <- str_split(network_table$A, ":", simplify = TRUE)[, 2]
network_table$B <- str_split(network_table$B, ":", simplify = TRUE)[, 2]
probe_genes <- unique(c(network_table$A, network_table$B))

interaction_table <- tibble(A = as.character(), B = as.character(), aliasA = as.character(), aliasB = as.character(), 
                            publicationID = as.character(), type = as.character(), confidenceScore = as.character())


for(i in 1:length(probe_genes)) {
    network_table <- interactions(refnet, species="9606", id=probe_genes[i], provider= "IntAct")
    network_table <- dplyr::select(network_table, A, B, aliasA, aliasB, publicationID, type, confidenceScore)
    interaction_table <- dplyr::bind_rows(interaction_table, network_table)
    print(i)
}

network_table <- interaction_table
network_table$A <- str_split(network_table$A, ":", simplify = TRUE)[, 2]
network_table$B <- str_split(network_table$B, ":", simplify = TRUE)[, 2]
network_table$publicationID <- str_split(network_table$publicationID, "\\|", n = 2, simplify = TRUE)[, 1]
network_table$type <- str_split(network_table$type, "\\(", simplify = TRUE)[, 2]
network_table$type <- str_split(network_table$type, "\\)", simplify = TRUE)[, 1]
network_table$aliasA <- str_split(network_table$aliasA, "\\(gene", n = 2, simplify = TRUE)[, 1]
network_table$aliasB <- str_split(network_table$aliasB, "\\(gene", n = 2, simplify = TRUE)[, 1]
network_table$aliasA <- str_split(network_table$aliasA, "uniprotkb:", n = 2, simplify = TRUE)[, 2]
network_table$aliasB <- str_split(network_table$aliasB, "uniprotkb:", n = 2, simplify = TRUE)[, 2]
network_table$confidenceScore <- str_split(network_table$confidenceScore, "miscore\\:", simplify = TRUE)[, 2]
network_table <- network_table[!duplicated(network_table), ]
network_table <- network_table[network_table$A != network_table$B, ]
#network_table$A <- str_split(network_table$A, "-", n = 2, simplify = TRUE)[, 1]
#network_table$B <- str_split(network_table$B, "-", n = 2, simplify = TRUE)[, 1]
network_table <- filter(network_table, confidenceScore >= 0.1)
network_table <- network_table %>% group_by(A, B) %>% arrange(desc(confidenceScore)) %>% dplyr::slice(1)
network_table <- network_table %>% dplyr::select(A, B, aliasA, aliasB, confidenceScore, publicationID, type)
#change the aliases to the entry name
up <- UniProt.ws::UniProt.ws(taxId=9606)


res2 <- UniProt.ws::select(up,
                           keys = unique(network_table$A),
                           columns = c("ENTRY-NAME"),
                           keytype = "UNIPROTKB")
names(res2) <- c("A", "A_NAME")
network_table <- left_join(network_table, res2)
res2 <- UniProt.ws::select(up,
                           keys = unique(network_table$B),
                           columns = c("ENTRY-NAME"),
                           keytype = "UNIPROTKB")
names(res2) <- c("B", "B_NAME")
network_table <- left_join(network_table, res2)


network_table$A_NAME <- stringr::str_split(network_table$A_NAME, "_", simplify = TRUE)[,1]
network_table$B_NAME <- stringr::str_split(network_table$B_NAME, "_", simplify = TRUE)[,1]

network_table <- network_table[nchar(network_table$A) == 6, ]
network_table <- network_table[nchar(network_table$B) == 6, ]

network_table$aliasA <- network_table$A_NAME
network_table$aliasB <- network_table$B_NAME
network_table <- dplyr::select(network_table, -A_NAME, -B_NAME)
write_csv(network_table, "network_table.csv")


#below is old way to view networks
# library(igraph)
# 
# #network_table <- network_table_temp
# #network_table <- filter(network_table, confidenceScore >= 0.5)
# 
# 
# 
# 
# 
# #----------------------------------------------
# #get first and second degree connections to target, find all probes, and all nodes that connect target to probe
# network_table_temp <- read_csv("network_table.csv")
# network_table <- filter(network_table_temp, confidenceScore >= 0.25)
# target <- "P01116"
# probe_genes <- as.character(unique(probe_list$UNIPROTKB))
# #target_connections <- network_table[apply(network_table, 1, FUN = 
# #                                 function(x, y = target) {any(grepl(y, x, ignore.case = TRUE))}), ]
# target_connections <- network_table[network_table$A == target | network_table$B == target, ]
# target_primary <- unique(c(target_connections$A, target_connections$B))
# target_connections <- network_table[network_table$A %in% target_primary | network_table$B %in% target_primary, ]
# relevant_probes <- probe_genes[probe_genes %in% unique(c(target_connections$A, target_connections$B))]
# probe_primary <- network_table[network_table$A %in% relevant_probes | network_table$B %in% relevant_probes, ]
# probe_primary <- unique(c(probe_primary$A, probe_primary$B))
# middle_nodes <- probe_primary[probe_primary %in% target_primary]
# #build attribute list
# name_list <- tibble(node = c(target_connections$A, target_connections$B), alias = c(target_connections$A_NAME, target_connections$B_NAME))
# name_list <- name_list %>% group_by(node)  %>% slice(1)
# relevant_nodes <- tibble(node = unique(c(middle_nodes, target, relevant_probes)), node_type = "middle")
# relevant_nodes <- left_join(relevant_nodes, name_list)
# relevant_nodes[relevant_nodes$node == target, 2] <- "target"
# relevant_nodes[relevant_nodes$node %in% relevant_probes, 2] <- "probe"
# 
# 
# target_connections <- network_table[network_table$A %in% relevant_nodes$node & network_table$B %in% relevant_nodes$node, ]
# 
# g <- graph.data.frame(as.matrix(target_connections[, c(1, 2)]), vertices = relevant_nodes, directed = FALSE)
# 
# 
# V(g)[V(g)$node_type == "probe"]$color <- "red"
# V(g)[V(g)$node_type == "middle"]$color <- "yellow"
# V(g)[V(g)$node_type == "target"]$color <- "lightblue"
# V(g)$name <- V(g)$alias
# plot(g, layout = layout.fruchterman.reingold, vertex.size = 3, vertex.label.cex = 0.9, vertex.label.dist=0.3, edge.curved = FALSE)
# 
# 
