library(data.table)
library(stringr)
library(KEGGREST)
library(biomaRt)
library(reshape2)
library(igraph)
library(ggpubr)

pathways <- fread("SCOPE_PathwayOverlaps.csv")
pathways[pathways == 0] <- NA
sel_pathways <-pathways[complete.cases(pathways),]

if(!dir.exists("Plots")){
  dir.create("Plots")
}

pathway_data <- fread("SCOPE_Gene_Level.csv")

all_corrs <- c()
for(curr_file in list.files("AllCorrelations/", full.names = TRUE)){
  curr_data <- fread(curr_file)
  curr_cancer <- str_split(basename(curr_file), pattern = "_", simplify = TRUE)[1]
  curr_data$Cancer <- curr_cancer
  
  all_corrs <- rbindlist(list(all_corrs, curr_data))
}

all_data <- c()
for(curr_pathway in sel_pathways$`KEGG GeneSet`){
  print(paste0("Currently working on ", curr_pathway, "..."))
  
  query <- keggGet(curr_pathway)
  genes <- query[[1]]$GENE
  genes <- genes[seq(2, length(genes), 2)]
  genes <- str_split(genes, ";", simplify = TRUE)[,1]
  
  human <- useMart("ensembl", dataset="hsapiens_gene_ensembl")
  gene_coords <- getBM(attributes=c("ensembl_gene_id", "external_gene_name"), 
                       filters="external_gene_name", values=genes, mart=human)
  
  sel_pathway_data <- pathway_data[pathway_data$`KEGG GeneSet` == curr_pathway,]
  
  rep_info <- data.table(Cancer = unique(sel_pathway_data$Cancer),
                         times = NA)
  for(i in 1:nrow(rep_info)){
    rep_info$times[i] <- length(genes)*length(unique(sel_pathway_data$`Core Gene`[sel_pathway_data$Cancer==rep_info$Cancer[i]]))
  }
  
  unique_data <- unique(sel_pathway_data[, c("Cancer", "KEGG GeneSet", "Pathway Name", "Core Gene")])
  
  curr_data <- unique_data[rep(seq_len(nrow(unique_data)), nrow(gene_coords))]
  curr_data <- setorder(curr_data, Cancer, `KEGG GeneSet`, `Pathway Name`, `Core Gene`)
  curr_data$`Secondary Gene` <- rep(gene_coords$ensembl_gene_id, length.out = nrow(curr_data))
  curr_data$`Core Gene Name` <- sel_pathway_data$`Core Gene Name`[match(curr_data$`Core Gene`, sel_pathway_data$`Core Gene`)]
  curr_data$`Secondary Gene Name` <-gene_coords$external_gene_name[match(curr_data$`Secondary Gene`, gene_coords$ensembl_gene_id)]

  curr_data <- merge(curr_data, all_corrs[, c("Cancer", "Core Gene", "Core Gene Name",
                                              "Secondary Gene", "Secondary Gene Name",
                                              "Tumor Correlation", "Normal Correlation")],
                     by = c("Cancer", "Core Gene", "Core Gene Name",
                            "Secondary Gene", "Secondary Gene Name"), 
                     all.x = TRUE)
  
  curr_data <- curr_data[complete.cases(curr_data),]
  
  all_data <- rbindlist(list(all_data, curr_data))
}

complete_data <- melt(all_data, id.vars = c("Cancer", "Core Gene", "Core Gene Name", "Secondary Gene",
                                            "Secondary Gene Name","KEGG GeneSet", "Pathway Name"),
                      variable.name = "Tissue", value.name = "Correlation")
complete_data$Tissue <- str_remove(complete_data$Tissue, " Correlation")

core_diffs <- unique(pathway_data[, c("Cancer", "Core Gene", "adj.P.Val_Core")])
complete_data <- merge(complete_data, core_diffs, by = c("Cancer", "Core Gene"))

fwrite(complete_data, "SCOPE_Network_Plot_Data.csv")

# mm to inch
setWidth = 90*0.039370 

# font size in pt
setFontSize = 5

# 1 in R = 0.75pt, so 0.25pt is specified as 
setLwd <- 0.25/0.75

data_used <- c()
for(curr_pathway in sel_pathways$`KEGG GeneSet`){
  
  tmp <- complete_data[complete_data$`KEGG GeneSet`==curr_pathway,]
  tmp2 <- unique(tmp[, c("Cancer", "Core Gene")])
  curr_plots_data <- c()
  for(x in unique(tmp2$Cancer)){
    curr_plots_data <- rbindlist(list(curr_plots_data,
                                      tmp[tmp$Cancer == x & tmp$`Core Gene` == tmp2$`Core Gene`[tmp2$Cancer == x][1],]))
  }
  
  data_used <- rbindlist(list(data_used, curr_plots_data))
  
  curr_nodes_data <- unique(curr_plots_data[, c("Cancer", "Core Gene Name", "Secondary Gene Name")])
  curr_nodes_data <- melt(curr_nodes_data, id.vars = "Cancer", variable.name = "type", value.name = "name")
  curr_nodes_data$type <- str_split(curr_nodes_data$type, pattern = " ", simplify = TRUE)[,1]
  curr_nodes_data <- unique(curr_nodes_data)
  curr_nodes_data <- setcolorder(curr_nodes_data, c("name", "type", "Cancer"))
  
  png(file=paste0("Plots/", curr_pathway, "_full.png"), width=setWidth, height=1.5, units="in", res=300, pointsize=setFontSize)
  par(mfcol=c(2, 6),
      mai = c(0.05, 0.05, 0.05, 0))
  
  curr_col = 1
  
  for(c in unique(curr_plots_data$Cancer)){
    curr_edges <- curr_plots_data[curr_plots_data$Cancer == c & curr_plots_data$Tissue == "Tumor", c("Core Gene Name", "Secondary Gene Name", "Correlation")]
    colnames(curr_edges) <- c("from", "to", "corr")
    curr_nodes <- curr_nodes_data[curr_nodes_data$Cancer == c,]
    network <- graph_from_data_frame(d=curr_edges, vertices = curr_nodes, directed=F) 
    
    c_scale <- colorRamp(c('red','white','blue'))
    E(network)$color = apply(c_scale(((E(network)$corr)+1)/2), 1, function(x) rgb(x[1]/255,x[2]/255,x[3]/255) )
    
    plot(network, layout = layout.star,
         vertex.size = ifelse(V(network)$type == "Core", 60, 15),
         vertex.label= ifelse(V(network)$type == "Core", V(network)$name, NA),
         vertex.label.font = 2,
         vertex.label.cex = 1.5,
         vertex.label.color = "black",
         vertex.label.font = "arial",
         vertex.color = ifelse(V(network)$type == "Core", "#78d697fe", "#D3D3D3"),
         vertex.frame.color = NA,
         edge.width = abs(E(network)$corr)*2,
         lwd=setLwd,
         margin = -0.24)
    
    if(curr_col == 1) mtext(side = 2, line = -0.3, "Primary Tumor", cex = 1, font = 1, family = "sans")
    mtext(side = 3, line = -1, c, cex = 1, font = 1, family = "sans")
    
    curr_edges <- curr_plots_data[curr_plots_data$Cancer == c & curr_plots_data$Tissue == "Normal", c("Core Gene Name", "Secondary Gene Name", "Correlation")]
    colnames(curr_edges) <- c("from", "to", "corr")
    curr_nodes <- curr_nodes_data[curr_nodes_data$Cancer == c,]
    network <- graph_from_data_frame(d=curr_edges, vertices = curr_nodes, directed=F) 
    
    c_scale <- colorRamp(c('red','white','blue'))
    E(network)$color = apply(c_scale(((E(network)$corr)+1)/2), 1, function(x) rgb(x[1]/255,x[2]/255,x[3]/255) )
    
    plot(network, layout = layout.star,
         vertex.size = ifelse(V(network)$type == "Core", 60, 15),
         vertex.label= ifelse(V(network)$type == "Core", V(network)$name, NA),
         vertex.label.font = 2,
         vertex.label.cex = 1.5,
         vertex.label.color = "black",
         vertex.label.font = "arial",
         vertex.color = ifelse(V(network)$type == "Core", "#78d697fe", "#D3D3D3"),
         vertex.frame.color = NA,
         edge.width = abs(E(network)$corr)*2,
         lwd=setLwd,
         margin = -0.24)
    
    if(curr_col == 1) mtext(side = 2, line = -0.3, "Solid Tissue Normal", cex = 1, font = 1, family = "sans")
    
    curr_col = curr_col + 1
  }
  
  dev.off()
  
  tests <- c()
  for(ca in unique(curr_plots_data$Cancer)){
    tmpx <- curr_plots_data$Correlation[curr_plots_data$Cancer == ca & curr_plots_data$Tissue == "Tumor"]
    tmpy <- curr_plots_data$Correlation[curr_plots_data$Cancer == ca & curr_plots_data$Tissue == "Normal"]
    tests <- rbindlist(list(tests, data.table(Cancer = ca,
                                              group1 = "Tumor",
                                              group2 = "Normal",
                                              p = signif(ks.test(tmpx, tmpy)$p.value, 3))))
  }
  
  p <- ggboxplot(curr_plots_data, x = "Tissue", y = "Correlation", color = "Tissue", facet.by = "Cancer", 
                 xlab = FALSE, legend = "bottom", ggtheme = theme_light(),
                 legend.title = "",
                 palette = "npg",
                 order = c("Tumor", "Normal"),
                 size = 0.2) +
    stat_pvalue_manual(tests, 
                       label.size = 1.5,
                       y.position = max(curr_plots_data$Correlation) + 0.1,
                       label = "KS test, p = {p}") +
    scale_y_continuous(expand = expansion(add = c(0, 0.15))) +
    theme(text = element_text(size=5),
          axis.ticks.x = element_blank(),
          axis.text.x = element_blank(),
          strip.text.x = element_text(margin = margin(0, 0, 0.02, 0, "cm")),
          legend.margin=margin(t=-0.4, r=0, b=0, l=0, unit="cm"),
          legend.key.size = unit(3, 'mm'),
          plot.margin=unit(rep(0.05, 4),"cm"),
          axis.title.y = element_text(margin = margin(t = 0, r = 0, b = 0, l = 0)))
  

  ggsave(paste0("Plots/", curr_pathway, "_boxplot_full.png"), p,
         width = 90*0.039370, height = 2.5, units = "in", dpi = 300)
}

fwrite(data_used, "SCOPE_Network_Plots.csv")

# Only non differentially Expressed Core Genes
data_used <- c()
for(curr_pathway in sel_pathways$`KEGG GeneSet`){
  
  tmp <- complete_data[complete_data$`KEGG GeneSet`==curr_pathway,]
  tmp2 <- unique(tmp[tmp$adj.P.Val_Core > 0.05, c("Cancer", "Core Gene")])
  curr_plots_data <- c()
  for(x in unique(tmp2$Cancer)){
    curr_plots_data <- rbindlist(list(curr_plots_data,
                                      tmp[tmp$Cancer == x & tmp$`Core Gene` == tmp2$`Core Gene`[tmp2$Cancer == x][1],]))
  }
  
  data_used <- rbindlist(list(data_used, curr_plots_data))
  
  curr_nodes_data <- unique(curr_plots_data[, c("Cancer", "Core Gene Name", "Secondary Gene Name")])
  curr_nodes_data <- melt(curr_nodes_data, id.vars = "Cancer", variable.name = "type", value.name = "name")
  curr_nodes_data$type <- str_split(curr_nodes_data$type, pattern = " ", simplify = TRUE)[,1]
  curr_nodes_data <- unique(curr_nodes_data)
  curr_nodes_data <- setcolorder(curr_nodes_data, c("name", "type", "Cancer"))
  
  png(file=paste0("Plots/diff_", curr_pathway, "_full.png"), width=setWidth, height=1.5, units="in", res=300, pointsize=setFontSize)
  par(mfcol=c(2, 6),
      mai = c(0.05, 0.05, 0.05, 0))
  
  curr_col = 1
  
  for(c in unique(curr_plots_data$Cancer)){
    curr_edges <- curr_plots_data[curr_plots_data$Cancer == c & curr_plots_data$Tissue == "Tumor", c("Core Gene Name", "Secondary Gene Name", "Correlation")]
    colnames(curr_edges) <- c("from", "to", "corr")
    curr_nodes <- curr_nodes_data[curr_nodes_data$Cancer == c,]
    network <- graph_from_data_frame(d=curr_edges, vertices = curr_nodes, directed=F) 
    
    c_scale <- colorRamp(c('red','white','blue'))
    E(network)$color = apply(c_scale(((E(network)$corr)+1)/2), 1, function(x) rgb(x[1]/255,x[2]/255,x[3]/255) )
    
    plot(network, layout = layout.star,
         vertex.size = ifelse(V(network)$type == "Core", 60, 15),
         vertex.label= ifelse(V(network)$type == "Core", V(network)$name, NA),
         vertex.label.font = 2,
         vertex.label.cex = 1.5,
         vertex.label.color = "black",
         vertex.label.font = "arial",
         vertex.color = ifelse(V(network)$type == "Core", "#78d697fe", "#D3D3D3"),
         vertex.frame.color = NA,
         edge.width = abs(E(network)$corr)*2,
         lwd=setLwd,
         margin = -0.24)
    
    if(curr_col == 1) mtext(side = 2, line = -0.3, "Primary Tumor", cex = 1, font = 1, family = "sans")
    mtext(side = 3, line = -1, c, cex = 1, font = 1, family = "sans")
    
    curr_edges <- curr_plots_data[curr_plots_data$Cancer == c & curr_plots_data$Tissue == "Normal", c("Core Gene Name", "Secondary Gene Name", "Correlation")]
    colnames(curr_edges) <- c("from", "to", "corr")
    curr_nodes <- curr_nodes_data[curr_nodes_data$Cancer == c,]
    network <- graph_from_data_frame(d=curr_edges, vertices = curr_nodes, directed=F) 
    
    c_scale <- colorRamp(c('red','white','blue'))
    E(network)$color = apply(c_scale(((E(network)$corr)+1)/2), 1, function(x) rgb(x[1]/255,x[2]/255,x[3]/255) )
    
    plot(network, layout = layout.star,
         vertex.size = ifelse(V(network)$type == "Core", 60, 15),
         vertex.label= ifelse(V(network)$type == "Core", V(network)$name, NA),
         vertex.label.font = 2,
         vertex.label.cex = 1.5,
         vertex.label.color = "black",
         vertex.label.font = "arial",
         vertex.color = ifelse(V(network)$type == "Core", "#78d697fe", "#D3D3D3"),
         vertex.frame.color = NA,
         edge.width = abs(E(network)$corr)*2,
         lwd=setLwd,
         margin = -0.24)
    
    if(curr_col == 1) mtext(side = 2, line = -0.3, "Solid Tissue Normal", cex = 1, font = 1, family = "sans")
    
    curr_col = curr_col + 1
  }
  
  dev.off()
  
  tests <- c()
  for(ca in unique(curr_plots_data$Cancer)){
    tmpx <- curr_plots_data$Correlation[curr_plots_data$Cancer == ca & curr_plots_data$Tissue == "Tumor"]
    tmpy <- curr_plots_data$Correlation[curr_plots_data$Cancer == ca & curr_plots_data$Tissue == "Normal"]
    tests <- rbindlist(list(tests, data.table(Cancer = ca,
                                              group1 = "Tumor",
                                              group2 = "Normal",
                                              p = signif(ks.test(tmpx, tmpy)$p.value, 3))))
  }
  
  p <- ggboxplot(curr_plots_data, x = "Tissue", y = "Correlation", color = "Tissue", facet.by = "Cancer", 
                 xlab = FALSE, legend = "bottom", ggtheme = theme_light(),
                 legend.title = "",
                 palette = "npg",
                 order = c("Tumor", "Normal"),
                 size = 0.2) +
    stat_pvalue_manual(tests, 
                       label.size = 1.5,
                       y.position = max(curr_plots_data$Correlation) + 0.1,
                       label = "KS test, p = {p}") +
    scale_y_continuous(expand = expansion(add = c(0, 0.15))) +
    theme(text = element_text(size=5),
          axis.ticks.x = element_blank(),
          axis.text.x = element_blank(),
          strip.text.x = element_text(margin = margin(0, 0, 0.02, 0, "cm")),
          legend.margin=margin(t=-0.4, r=0, b=0, l=0, unit="cm"),
          legend.key.size = unit(3, 'mm'),
          plot.margin=unit(rep(0.05, 4),"cm"),
          axis.title.y = element_text(margin = margin(t = 0, r = 0, b = 0, l = 0)))
  
  
  ggsave(paste0("Plots/diff_", curr_pathway, "_boxplot_full.png"), p,
         width = 90*0.039370, height = 2.5, units = "in", dpi = 300)
}

fwrite(data_used, "SCOPE_Network_Plots_NonDifferential.csv")
