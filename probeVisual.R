
if(require("ReactomePA")){
} else {
  source("https://bioconductor.org/biocLite.R")
  biocLite("ReactomePA")
  require("ReactomePA")
}
if(require("org.Hs.eg.db")){
} else {
  source("https://bioconductor.org/biocLite.R")
  biocLite("org.Hs.eg.db")
  require("org.Hs.eg.db")
}
if(require("readr")){
} else {
  install.packages("readr")
  require("readr")
}



if(require("RCy3")){
} else {
  source("https://bioconductor.org/biocLite.R")
  biocLite("RCy3")
  require("RCy3")
}
if(require("plotly")){
} else {
  install.packages("plotly")
  require("plotly")
}
if(require("igraph")){
} else {
  install.packages("igraph")
  require("igraph")
}
if(require("dplyr")){
} else {
  install.packages("dplyr")
  require("dplyr")
}

if(require("pathview")){
} else {
  source("https://bioconductor.org/biocLite.R")
  biocLite("pathview")
  require("pathview")
}

probe_list <- read_csv("PROBELIST.csv")

reactome_probe_network <- function(pathway_name, data = probe_list){
  if(is.character(pathway_name) == FALSE) stop("Pathway name is not a character.")
  gene_names <- unique(data$ENTREZ_GENE)
  gene_names <- as.character(na.omit(gene_names))
  mygeneList <- rep(100, length(gene_names))
  names(mygeneList) <- gene_names
  tryCatch({
    print(viewPathway(pathway_name, organism = "human", fixed = TRUE, readable = TRUE, 
              foldChange=mygeneList, vertex.label.cex = 0.75))
  }, error=function(e){}, warning = function(w){})
}
#Regulation of TP53 Activity through Methylation
#png("tp53.png", width = 2000, height = 2000, res = 300)
#reactome_probe_network("Regulation of TP53 Activity through Methylation")
#dev.off()

kegg_probe_network <- function(pathway_name, data = probe_list){
  if(is.character(pathway_name) == FALSE) stop("Pathway name is not a character.")
  gene_names <- as.numeric(unique(data$ENTREZ_GENE))
  gene_names <- as.character(na.omit(gene_names))
  #keggview.native(gene.data = gene_names, res = 300, pathway.name = pathway_name, species = "hsa", col.key = FALSE)
  pathview(gene.data = gene_names, pathway.id = pathway_name, species = "hsa", col.key = FALSE, kegg.native = TRUE)
  #keggview.graph(gene.data = gene_names, res = 300, pathway.name = pathway_name, species = "hsa", col.key = FALSE)
  
  print("Pathway image saved in local directory")
}




highlight_nodes_cytoscape <- function(border_width = 6, data = probe_list){
  identifiers <- c(data$UNIPROTKB, data$ENTREZ_GENE)
  cyconn <- CytoscapeConnection()
  window_names <- getWindowList(cyconn)
  window_ids <- getWindowID(cyconn, window_names)
  cw <- existing.CytoscapeWindow(window_names, copy.graph.from.cytoscape.to.R = TRUE)
  nodes <- getAllNodes(cw)
  setNodeBorderColorDirect(cw, nodes[nodes %in% identifiers], "#ff000d")
  setNodeBorderOpacityDirect(cw, nodes[nodes %in% identifiers], 255)
  setNodeBorderWidthDirect(cw, nodes[nodes %in% identifiers], border_width)
}

network_table <- read_csv("network_table.csv")
#plot_probe_target_network("P45973", 0.5)
plot_probe_target_network <- function(target, cutoff = 0.4, data = probe_list, network_table_temp = network_table){
  ####
  ####add a stopifnot here for target and cutoff
  ####

  #----------------------------------------------
  #get first and second degree connections to target, find all probes, and all nodes that connect target to probe
  if(is.character(target) == FALSE) stop("Target name is not a character.")
  network_table_temp <- filter(network_table_temp, confidenceScore >= cutoff)
  probe_genes <- as.character(unique(data$UNIPROTKB))
  #target_connections <- network_table_temp[apply(network_table_temp, 1, FUN = 
  #                                 function(x, y = target) {any(grepl(y, x, ignore.case = TRUE))}), ]
  target_connections <- network_table_temp[network_table_temp$A == target | network_table_temp$B == target, ]
  if(nrow(target_connections) == 0) stop("Target not found.")
  target_primary <- unique(c(target_connections$A, target_connections$B))
  target_connections <- network_table_temp[network_table_temp$A %in% target_primary | network_table_temp$B %in% target_primary, ]
  relevant_probes <- probe_genes[probe_genes %in% unique(c(target_connections$A, target_connections$B))]
  probe_primary <- network_table_temp[network_table_temp$A %in% relevant_probes | network_table_temp$B %in% relevant_probes, ]
  probe_primary <- unique(c(probe_primary$A, probe_primary$B))
  middle_nodes <- probe_primary[probe_primary %in% target_primary]
  #build attribute list
  name_list <- tibble(node = c(target_connections$A, target_connections$B), alias = c(target_connections$aliasA, target_connections$aliasB))
  name_list <- name_list %>% group_by(node)  %>% dplyr::slice(1)
  relevant_nodes <- tibble(node = unique(c(middle_nodes, target, relevant_probes)), node_type = "middle")
  relevant_nodes <- left_join(relevant_nodes, name_list)
  relevant_nodes[relevant_nodes$node == target, 2] <- "target"
  relevant_nodes[relevant_nodes$node %in% relevant_probes, 2] <- "probe"
  # add the hoverinfo
  relevant_nodes$hoverinfo <- relevant_nodes$alias
  for(i in 1:nrow(relevant_nodes)){
    if(relevant_nodes$node_type[i] != "middle")
    {
      relevant_nodes$hoverinfo[i] <- paste(relevant_nodes$alias[i], "<br>Chemical probes:<br>", 
                                           (data[data$ENTRY_NAME == relevant_nodes$alias[i], 1]))
    }
  }
  
  target_connections <- network_table_temp[network_table_temp$A %in% relevant_nodes$node & network_table_temp$B %in% relevant_nodes$node, ]
  g <- graph.data.frame(as.matrix(target_connections[, c(1, 2)]), vertices = relevant_nodes, directed = FALSE)
  
  L <- layout.fruchterman.reingold(g)
  vs <- V(g)
  es <- as.data.frame(get.edgelist(g))
  Nv <- length(vs)
  Ne <- length(es[1]$V1)
  
  Xn <- L[,1]
  Yn <- L[,2]
  
  V(g)[V(g)$node_type == "probe"]$color <- "Has chemical probe"
  V(g)[V(g)$node_type == "middle"]$color <- "Intermediate node"
  V(g)[V(g)$node_type == "target"]$color <- "Target"
  
  data <- data.frame(X = Xn, Y = Yn, hovertext = V(g)$hoverinfo, alias = V(g)$alias)
  
  network <- plot_ly(data, type = "scatter", x = ~X, y = ~Y, mode = "markers", 
                     marker = list(size = 15), 
                     hoverinfo = "text", text = ~hovertext, color = V(g)$color) %>%
    add_annotations(x = data$X,
                    y = data$Y,
                    text = data$alias,
                    showarrow = FALSE)
  
  edge_shapes <- list()
  for(i in 1:Ne) {
    v0 <- es[i,]$V1
    v1 <- es[i,]$V2
    #### is not getting the right vertices
    edge_shape = list(
      type = "line",
      line = list(color = "#b2b2d8", width = 0.7),
      x0 = Xn[which(as.character(v0) == as.character(names(vs)))],
      y0 = Yn[which(as.character(v0) == as.character(names(vs)))],
      x1 = Xn[which(as.character(v1) == as.character(names(vs)))],
      y1 = Yn[which(as.character(v1) == as.character(names(vs)))]
    )
    
    edge_shapes[[i]] <- edge_shape
  }
  
  network <- layout(
    network,
    title = paste("Local pharmacology network for:", target),
    shapes = edge_shapes,
    xaxis = list(title = "", showgrid = FALSE, showticklabels = FALSE, zeroline = FALSE),
    yaxis = list(title = "", showgrid = FALSE, showticklabels = FALSE, zeroline = FALSE)
  )

  print(network) 
}