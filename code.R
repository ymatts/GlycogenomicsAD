
## Loading Required Libraries
library(data.table)
library(readxl)
library(ggplot2)
library(ggpie)
library(biomaRt)
library(metaviz)
library(fgsea)
library(igraph)
library(ggnet)
library(network)
library(sna)
library(ggnetwork)
library(DT)
library(Seurat)
library(CellChat)
library(nichenetr)
library(scCustomize)
library(RandomWalkRestartMH)
library(dorothea)
library(OmnipathR)
library(fedup)
library(RColorBrewer)
library(tidyverse)
library(gtools)
library(pheatmap)
library(reshape2)
library(hdWGCNA)
library(scCustomize)
library(pathwayPCA)
library(metap)

####################################
# Data Import and Initial Analysis #
####################################

gly_type <- readRDS("data/misc/glycogene.rds")
glytype_parent <- sapply(strsplit(gly_type[,1]," \\/ "),function(x)x[1])
glytype_sub <- sapply(strsplit(gly_type[,1]," \\/ "),function(x)x[2])
df <- data.frame(gene_name = gly_type[,2],parent = glytype_parent,sub = glytype_sub,check.names = F)
df$sub <- ifelse(is.na(df$sub),"NA",df$sub)
gg <- ggnestedpie(df,group_key = c("parent","sub"),
                  count_type = "full",inner_label_info = "all", inner_label_split = NULL, inner_label_size = 2,
                  outer_label_type = "circle", outer_label_pos = "out", outer_label_info = "all")
gg

# Types of glycogene
tb <- table(df$parent,df$sub)
tb

#######################################
# Meta analysis results for glycogene #
#######################################

meta_res <- readRDS("data/deg/meta_ad_ctr.rds")
glygene_deg <- meta_res[meta_res$symbol%in%gly_type[,2],]
DT::datatable(glygene_deg,glygene_deg[,c("symbol","TE.random","lower.random","upper.random","fdr.random")])

# enrichment analysis
db <- useMart("ensembl")
hg <- useDataset("hsapiens_gene_ensembl", mart = db)
stats <- dat$TE.random
gene_ids <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol"),
                  filters = "ensembl_gene_id", 
                  values = dat$ensembl_gene_id,
                  mart = hg)
names(stats) <- gene_ids$hgnc_symbol[match(dat$ensembl_gene_id,gene_ids$ensembl_gene_id)]
stats <- stats[!is.na(names(stats))]
unique_gene <- unique(names(stats))

stats2 <- rep(NA,length(unique_gene))
names(stats2) <- unique_gene
for(i in seq_along(stats2)){
  idx <- which(unique_gene[i] == names(stats))
  selidx <- idx[which.max(abs(stats[idx]))]
  stats2[i] <- stats[selidx]
}
gmt <- gmtPathways("data/misc/c2.cp.v2022.1.Hs.symbols.gmt")

gsea_result <- fgseaSimple(gmt,stats2,nperm = 1000)
gsea_sigresult <- gsea_result[gsea_result$padj <= .05,]
gsea_sigresult <- gsea_sigresult[order(gsea_sigresult$NES,decreasing = T),]
gsea_sigresult


######################################################
# core gene identification from leading edge of GSEA #
######################################################

geneset <- gsea_sigresult$leadingEdge
names(geneset) <- gsea_sigresult$pathway
edgedf <- lapply(seq_along(geneset),function(xx)data.frame(pathway = names(geneset)[xx],gene = geneset[[xx]]))
edgedf <- do.call(rbind,edgedf)
edgedf[,1] <- substr(gsub("_"," ",gsub("REACTOME_|NABA_|WP_|KEGG_","",edgedf[,1])),1,50)

edgedf <- edgedf[edgedf[,2]%in%glygene_deg$symbol,]

subdcnet <- network(edgedf,directed = F)
ggnet <- ggnetwork(subdcnet,layout="fruchtermanreingold")
ggnet$category <- ifelse(ggnet$vertex.names %in% unique(edgedf[,2]),"Gene","Pathway")
ggnet$category <- factor(ggnet$category,levels = c("Gene","Pathway"))

gg <- ggplot(ggnet, aes(x = x, y = y, xend = xend, yend = yend)) +
  geom_edges(aes(alpha=.7),
             color = "grey1",
             size = 1) + 
  scale_color_brewer(palette = "Dark2") +
  guides(alpha=FALSE) + 
  geom_nodelabel_repel(aes(label = vertex.names,color=category),
                       fontface = "bold", 
                       box.padding = unit(1, "lines")) +
  geom_nodes(aes(color = category,shape=category,alpha=.8),size=4) +
  theme_blank()
gg


# network degree
tb <- sort(table(edgedf[,2]),decreasing = T)
df <- data.frame(tb)
colnames(df) <- c("Glycogene","Degree")
gg <- ggplot(df,aes(x = Glycogene,y = Degree)) + 
  geom_bar(stat = "identity", fill=rgb(0.1,0.4,0.5,0.7)) + 
  xlab("Glycogene") + ylab("Degree") + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),text = element_text(size = 8))
gg



#####################################
# Forest plots for pathway activity #
#####################################

dat <- readRDS("data/deg/differentialExpressionSummary.rds")
subdat <- dat[dat$hgnc_symbol == "PLOD3" & dat$Comparison == "AD-CONTROL" & dat$Model=="Diagnosis",]
subdat$study <- paste(subdat$Study,subdat$Tissue,sep="_")
x <- subdat[,c("study","logFC","CI.L","CI.R")]

## Plot forest plot
gg <- ggplot(x, aes(y = study, x = logFC)) +
  geom_point(shape = 18, size = 5) +  
  geom_errorbarh(aes(xmin = CI.L, xmax = CI.R), height = 0.25) +
  geom_vline(xintercept = 0, color = "red", linetype = "dashed", cex = 1, alpha = 0.5) +
  #  scale_y_continuous(name = "", breaks=1:4, labels = dat$label, trans = "reverse") +
  xlab("Log fold change (95% CI)") + 
  ylab(" ") + 
  theme_bw() +
  theme(panel.border = element_blank(),
        panel.background = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"),
        axis.text.y = element_text(size = 12, colour = "black"),
        axis.text.x.bottom = element_text(size = 12, colour = "black"),
        axis.title.x = element_text(size = 12, colour = "black"))
gg


# colagen formaion
genes <- unlist(gsea_sigresult[grep("COLLAGEN_FORMATION",gsea_sigresult$pathway),]$leadingEdge)
subdat <- dat[dat$hgnc_symbol %in% genes & dat$Comparison == "AD-CONTROL" & dat$Model=="Diagnosis",]
subdat$study <- paste(subdat$Study,subdat$Tissue,sep="_")
x <- subdat[,c("study","logFC","CI.L","CI.R")]
x <- split(x,x$study)
x <- t(sapply(x,function(x)colMeans(x[,-1])))
x <- data.frame(study=rownames(x),x)
## Plot forest plot
gg <- ggplot(x, aes(y = study, x = logFC)) +
  geom_point(shape = 18, size = 5) +  
  geom_errorbarh(aes(xmin = CI.L, xmax = CI.R), height = 0.25) +
  geom_vline(xintercept = 0, color = "red", linetype = "dashed", cex = 1, alpha = 0.5) +
  #  scale_y_continuous(name = "", breaks=1:4, labels = dat$label, trans = "reverse") +
  xlab("Log fold change (95% CI)") + 
  ylab(" ") + 
  theme_bw() +
  theme(panel.border = element_blank(),
        panel.background = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"),
        axis.text.y = element_text(size = 12, colour = "black"),
        axis.text.x.bottom = element_text(size = 12, colour = "black"),
        axis.title.x = element_text(size = 12, colour = "black"))
gg


# ECM 
genes <- unlist(gsea_sigresult[grep("EXTRACELLULAR_MATRIX_ORGANIZATION",gsea_sigresult$pathway),]$leadingEdge)
subdat <- dat[dat$hgnc_symbol %in% genes & dat$Comparison == "AD-CONTROL" & dat$Model=="Diagnosis",]
subdat$study <- paste(subdat$Study,subdat$Tissue,sep="_")
x <- subdat[,c("study","logFC","CI.L","CI.R")]
x <- split(x,x$study)
x <- t(sapply(x,function(x)colMeans(x[,-1])))
x <- data.frame(study=rownames(x),x)
## Plot forest plot
gg <- ggplot(x, aes(y = study, x = logFC)) +
  geom_point(shape = 18, size = 5) +  
  geom_errorbarh(aes(xmin = CI.L, xmax = CI.R), height = 0.25) +
  geom_vline(xintercept = 0, color = "red", linetype = "dashed", cex = 1, alpha = 0.5) +
  #  scale_y_continuous(name = "", breaks=1:4, labels = dat$label, trans = "reverse") +
  xlab("Log fold change (95% CI)") + 
  ylab(" ") + 
  theme_bw() +
  theme(panel.border = element_blank(),
        panel.background = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"),
        axis.text.y = element_text(size = 12, colour = "black"),
        axis.text.x.bottom = element_text(size = 12, colour = "black"),
        axis.title.x = element_text(size = 12, colour = "black"))



#######################
# PathwayPCA analysis #
#######################

# ROSMAP
exprs <- readRDS("data/pathwayPCA/rosmap_exprs.rds")
ecm_gsea_sigresult <- as.data.frame(gsea_sigresult[grep("ECM|Extracellular|Collagen",gsea_sigresult$pathway,ignore.case = T),])
ecm_genes <- ecm_gsea_sigresult$leadingEdge[[2]]
common_gene <- intersect(gname,ecm_genes)
ecm_eigen <- colMeans(exprs[gname %in% ecm_genes,],na.rm = T)

gmt = read_gmt(file = "data/misc/c2.cp.reactome.v2022.1.Hs.symbols.gmt")
subexprs <- t(exprs[!rownames(exprs) %in% ecm_genes,])
pcaobj <- pcaMethods::pca(subexprs,nPcs = 5)
subexprs <- pcaobj@completeObs

df = data.frame(Sample = names(ecm_eigen),
                ecm = ecm_eigen,subexprs
)

obj = CreateOmics(
  assayData_df = df[,-2], 
  pathwayCollection_ls = gmt,
  response = df[,1:2],
  respType = "regression",
  minPathSize = 3
)

aespc_obj_ROSMAP = AESPCA_pVals(
  object = obj,
  numPCs = 1,
  parallel = TRUE,
  numCores = 20,
  numReps = 0,
  adjustment = "BH"
)

pPCA_obj_ROSMAP = getPathpVals(aespc_obj_ROSMAP, score = FALSE, numPaths = Inf)
pPCA_obj_ROSMAP


# MAYO
exprs <- readRDS("data/pathwayPCA/mayo_cbe_exprs.rds")
ecm_gsea_sigresult <- as.data.frame(gsea_sigresult[grep("ECM|Extracellular|Collagen",gsea_sigresult$pathway,ignore.case = T),])
ecm_genes <- ecm_gsea_sigresult$leadingEdge[[2]]
common_gene <- intersect(gname,ecm_genes)
ecm_eigen <- colMeans(exprs[gname %in% ecm_genes,],na.rm = T)

gmt = read_gmt(file = "data/misc/c2.cp.reactome.v2022.1.Hs.symbols.gmt")
subexprs <- t(exprs[!rownames(exprs) %in% ecm_genes,])
pcaobj <- pcaMethods::pca(subexprs,nPcs = 5)
subexprs <- pcaobj@completeObs

df = data.frame(Sample = names(ecm_eigen),
                ecm = ecm_eigen,subexprs
)


df = data.frame(Sample = colnames(case),
                ecm = ecm_eigen,
                t(case[!gname %in% ecm_genes,]))

obj = CreateOmics(
  assayData_df = df[,-2], 
  pathwayCollection_ls = gmt,
  response = df[,1:2],
  respType = "regression",
  minPathSize = 3
)

aespc_obj_MAYOCBE = AESPCA_pVals(
  object = obj,
  numPCs = 1,
  parallel = TRUE,
  numCores = 20,
  numReps = 0,
  adjustment = "BH"
)

pPCA_obj_MAYOCBE = getPathpVals(aespc_obj_MAYOCBE, score = FALSE, numPaths = Inf)
pPCA_obj_MAYOCBE

# MAYO
exprs <- readRDS("data/pathwayPCA/mayo_tcx_exprs.rds")
ecm_gsea_sigresult <- as.data.frame(gsea_sigresult[grep("ECM|Extracellular|Collagen",gsea_sigresult$pathway,ignore.case = T),])
ecm_genes <- ecm_gsea_sigresult$leadingEdge[[2]]
common_gene <- intersect(gname,ecm_genes)
ecm_eigen <- colMeans(exprs[gname %in% ecm_genes,],na.rm = T)

gmt = read_gmt(file = "data/misc/c2.cp.reactome.v2022.1.Hs.symbols.gmt")
subexprs <- t(exprs[!rownames(exprs) %in% ecm_genes,])
pcaobj <- pcaMethods::pca(subexprs,nPcs = 5)
subexprs <- pcaobj@completeObs

df = data.frame(Sample = names(ecm_eigen),
                ecm = ecm_eigen,subexprs
)


df = data.frame(Sample = colnames(case),
                ecm = ecm_eigen,
                t(case[!gname %in% ecm_genes,]))

obj = CreateOmics(
  assayData_df = df[,-2], 
  pathwayCollection_ls = gmt,
  response = df[,1:2],
  respType = "regression",
  minPathSize = 3
)

aespc_obj_MAYOTCX = AESPCA_pVals(
  object = obj,
  numPCs = 1,
  parallel = TRUE,
  numCores = 20,
  numReps = 0,
  adjustment = "BH"
)

pPCA_obj_MAYOTCX = getPathpVals(aespc_obj_MAYOTCX, score = FALSE, numPaths = Inf)
pPCA_obj_MAYOTCX

# MSSB PHG
exprs <- readRDS("data/pathwayPCA/mssb_phg_exprs.rds")
ecm_gsea_sigresult <- as.data.frame(gsea_sigresult[grep("ECM|Extracellular|Collagen",gsea_sigresult$pathway,ignore.case = T),])
ecm_genes <- ecm_gsea_sigresult$leadingEdge[[2]]
common_gene <- intersect(gname,ecm_genes)
ecm_eigen <- colMeans(exprs[gname %in% ecm_genes,],na.rm = T)

gmt = read_gmt(file = "data/misc/c2.cp.reactome.v2022.1.Hs.symbols.gmt")
subexprs <- t(exprs[!rownames(exprs) %in% ecm_genes,])
pcaobj <- pcaMethods::pca(subexprs,nPcs = 5)
subexprs <- pcaobj@completeObs

df = data.frame(Sample = names(ecm_eigen),
                ecm = ecm_eigen,subexprs
)


df = data.frame(Sample = colnames(case),
                ecm = ecm_eigen,
                t(case[!gname %in% ecm_genes,]))

obj = CreateOmics(
  assayData_df = df[,-2], 
  pathwayCollection_ls = gmt,
  response = df[,1:2],
  respType = "regression",
  minPathSize = 3
)

aespc_obj_MSSBPHG = AESPCA_pVals(
  object = obj,
  numPCs = 1,
  parallel = TRUE,
  numCores = 20,
  numReps = 0,
  adjustment = "BH"
)

pPCA_obj_MSSBPHG = getPathpVals(aespc_obj_MSSBPHG, score = FALSE, numPaths = Inf)
pPCA_obj_MSSBPHG 


# MSSB STG
exprs <- readRDS("data/pathwayPCA/mssb_stg_exprs.rds")
ecm_gsea_sigresult <- as.data.frame(gsea_sigresult[grep("ECM|Extracellular|Collagen",gsea_sigresult$pathway,ignore.case = T),])
ecm_genes <- ecm_gsea_sigresult$leadingEdge[[2]]
common_gene <- intersect(gname,ecm_genes)
ecm_eigen <- colMeans(exprs[gname %in% ecm_genes,],na.rm = T)

gmt = read_gmt(file = "data/misc/c2.cp.reactome.v2022.1.Hs.symbols.gmt")
subexprs <- t(exprs[!rownames(exprs) %in% ecm_genes,])
pcaobj <- pcaMethods::pca(subexprs,nPcs = 5)
subexprs <- pcaobj@completeObs

df = data.frame(Sample = names(ecm_eigen),
                ecm = ecm_eigen,subexprs
)


df = data.frame(Sample = colnames(case),
                ecm = ecm_eigen,
                t(case[!gname %in% ecm_genes,]))

obj = CreateOmics(
  assayData_df = df[,-2], 
  pathwayCollection_ls = gmt,
  response = df[,1:2],
  respType = "regression",
  minPathSize = 3
)

aespc_obj_MSSBSTG = AESPCA_pVals(
  object = obj,
  numPCs = 1,
  parallel = TRUE,
  numCores = 20,
  numReps = 0,
  adjustment = "BH"
)

pPCA_obj_MSSBSTG = getPathpVals(aespc_obj_MSSBSTG, score = FALSE, numPaths = Inf)
pPCA_obj_MSSBSTG

# MSSB FP
exprs <- readRDS("data/pathwayPCA/mssb_fp_exprs.rds")
ecm_gsea_sigresult <- as.data.frame(gsea_sigresult[grep("ECM|Extracellular|Collagen",gsea_sigresult$pathway,ignore.case = T),])
ecm_genes <- ecm_gsea_sigresult$leadingEdge[[2]]
common_gene <- intersect(gname,ecm_genes)
ecm_eigen <- colMeans(exprs[gname %in% ecm_genes,],na.rm = T)

gmt = read_gmt(file = "data/misc/c2.cp.reactome.v2022.1.Hs.symbols.gmt")
subexprs <- t(exprs[!rownames(exprs) %in% ecm_genes,])
pcaobj <- pcaMethods::pca(subexprs,nPcs = 5)
subexprs <- pcaobj@completeObs

df = data.frame(Sample = names(ecm_eigen),
                ecm = ecm_eigen,subexprs
)


df = data.frame(Sample = colnames(case),
                ecm = ecm_eigen,
                t(case[!gname %in% ecm_genes,]))

obj = CreateOmics(
  assayData_df = df[,-2], 
  pathwayCollection_ls = gmt,
  response = df[,1:2],
  respType = "regression",
  minPathSize = 3
)

aespc_obj_MSSBFP = AESPCA_pVals(
  object = obj,
  numPCs = 1,
  parallel = TRUE,
  numCores = 20,
  numReps = 0,
  adjustment = "BH"
)

pPCA_obj_MSSBFP = getPathpVals(aespc_obj_MSSBFP, score = FALSE, numPaths = Inf)
pPCA_obj_MSSBFP


# MSSB IFG
exprs <- readRDS("data/pathwayPCA/mssb_ifg_exprs.rds")
ecm_gsea_sigresult <- as.data.frame(gsea_sigresult[grep("ECM|Extracellular|Collagen",gsea_sigresult$pathway,ignore.case = T),])
ecm_genes <- ecm_gsea_sigresult$leadingEdge[[2]]
common_gene <- intersect(gname,ecm_genes)
ecm_eigen <- colMeans(exprs[gname %in% ecm_genes,],na.rm = T)

gmt = read_gmt(file = "data/misc/c2.cp.reactome.v2022.1.Hs.symbols.gmt")
subexprs <- t(exprs[!rownames(exprs) %in% ecm_genes,])
pcaobj <- pcaMethods::pca(subexprs,nPcs = 5)
subexprs <- pcaobj@completeObs

df = data.frame(Sample = names(ecm_eigen),
                ecm = ecm_eigen,subexprs
)


df = data.frame(Sample = colnames(case),
                ecm = ecm_eigen,
                t(case[!gname %in% ecm_genes,]))

obj = CreateOmics(
  assayData_df = df[,-2], 
  pathwayCollection_ls = gmt,
  response = df[,1:2],
  respType = "regression",
  minPathSize = 3
)

aespc_obj_MSSBIFG = AESPCA_pVals(
  object = obj,
  numPCs = 1,
  parallel = TRUE,
  numCores = 20,
  numReps = 0,
  adjustment = "BH"
)

pPCA_obj_MSSBIFG = getPathpVals(aespc_obj_MSSBIFG, score = FALSE, numPaths = Inf)
pPCA_obj_MSSBIFG

# integrating results for pathwayPCA
result <- list(
  ROSMAP = pPCA_obj_ROSMAP,
  MAYOCBE = pPCA_obj_MAYOCBE,
  MAYOTCX = pPCA_obj_MAYOTCX,
  MSSBFP = pPCA_obj_MSSBFP,
  MSSBIFG = pPCA_obj_MSSBIFG,
  MSSBSTG = pPCA_obj_MSSBSTG,
  MSSBPHG = pPCA_obj_MSSBPHG
)


gg <- result$ROSMAP[1:10,] %>%
  mutate(terms = substr(gsub("REACTOME_","",terms),1,40)) %>%
  ggplot(aes(x = reorder(terms, FDR_BH,decreasing = TRUE), y = -log10(FDR_BH))) +
  geom_bar(stat = "identity", fill = "#191970", width = 0.7) +
  labs(
    x = "Pathways",
    y = "Negative LN (FDR)",
    caption = "Yellow line: raw p-value = 0.01") +
  geom_hline(yintercept = -log(0.01), size = 1, color = "yellow") +
  coord_flip() +
  theme_bw()
gg + ggtitle("ROSMAP")

gg <- result$MAYOCBE[1:10,] %>%
  mutate(terms = substr(gsub("REACTOME_","",terms),1,40)) %>%
  ggplot(aes(x = reorder(terms, FDR_BH,decreasing = TRUE), y = -log10(FDR_BH))) +
  #geom_bar(stat = "identity", fill = "#191970", width = 0.7) +
  geom_bar(stat = "identity", fill = "#0072B2", width = 0.7) +
  labs(
    x = "Pathways",
    y = "Negative LN (FDR)",
    caption = "Yellow line: raw p-value = 0.01") +
  geom_hline(yintercept = -log(0.01), size = 1, color = "yellow") +
  coord_flip() +
  theme_bw()
gg + ggtitle("MAYOCBE")

gg <- result$MAYOTCX[1:10,] %>%
  mutate(terms = substr(gsub("REACTOME_","",terms),1,40)) %>%
  ggplot(aes(x = reorder(terms, FDR_BH,decreasing = TRUE), y = -log10(FDR_BH))) +
  #geom_bar(stat = "identity", fill = "#191970", width = 0.7) +
  geom_bar(stat = "identity", fill = "#0072B2", width = 0.7) +
  labs(
    x = "Pathways",
    y = "Negative LN (FDR)",
    caption = "Yellow line: raw p-value = 0.01") +
  geom_hline(yintercept = -log(0.01), size = 1, color = "yellow") +
  coord_flip() +
  theme_bw()
gg + ggtitle("MAYOTCX")

gg <- result$MSSBFP[1:10,] %>%
  mutate(terms = substr(gsub("REACTOME_","",terms),1,40)) %>%
  ggplot(aes(x = reorder(terms, FDR_BH,decreasing = TRUE), y = -log10(FDR_BH))) +
  #geom_bar(stat = "identity", fill = "#191970", width = 0.7) +
  geom_bar(stat = "identity", fill = "#CC79A7", width = 0.7) +
  labs(
    x = "Pathways",
    y = "Negative LN (FDR)",
    caption = "Yellow line: raw p-value = 0.01") +
  geom_hline(yintercept = -log(0.01), size = 1, color = "yellow") +
  coord_flip() +
  theme_bw()
gg + ggtitle("MSSBFP")

gg <- result$MSSBIFG[1:10,] %>%
  mutate(terms = substr(gsub("REACTOME_","",terms),1,40)) %>%
  ggplot(aes(x = reorder(terms, FDR_BH,decreasing = TRUE), y = -log10(FDR_BH))) +
  #geom_bar(stat = "identity", fill = "#191970", width = 0.7) +
  geom_bar(stat = "identity", fill = "#CC79A7", width = 0.7) +
  labs(
    x = "Pathways",
    y = "Negative LN (FDR)",
    caption = "Yellow line: raw p-value = 0.01") +
  geom_hline(yintercept = -log(0.01), size = 1, color = "yellow") +
  coord_flip() +
  theme_bw()
gg + ggtitle("MSSBIFG")


gg <- result$MSSBSTG[1:10,] %>%
  mutate(terms = substr(gsub("REACTOME_","",terms),1,40)) %>%
  ggplot(aes(x = reorder(terms, FDR_BH,decreasing = TRUE), y = -log10(FDR_BH))) +
  #geom_bar(stat = "identity", fill = "#191970", width = 0.7) +
  geom_bar(stat = "identity", fill = "#CC79A7", width = 0.7) +
  labs(
    x = "Pathways",
    y = "Negative LN (FDR)",
    caption = "Yellow line: raw p-value = 0.01") +
  geom_hline(yintercept = -log(0.01), size = 1, color = "yellow") +
  coord_flip() +
  theme_bw()
gg + ggtitle("MSSBSTG")

gg <- result$MSSBPHG[1:10,] %>%
  mutate(terms = substr(gsub("REACTOME_","",terms),1,40)) %>%
  ggplot(aes(x = reorder(terms, FDR_BH,decreasing = TRUE), y = -log10(FDR_BH))) +
  #geom_bar(stat = "identity", fill = "#191970", width = 0.7) +
  geom_bar(stat = "identity", fill = "#CC79A7", width = 0.7) +
  labs(
    x = "Pathways",
    y = "Negative LN (FDR)",
    caption = "Yellow line: raw p-value = 0.01") +
  geom_hline(yintercept = -log(0.01), size = 1, color = "yellow") +
  coord_flip() +
  theme_bw()
gg + ggtitle("MSSBPHG")



#################################################
# Single cell anlysis using Human protein atlas #
#################################################
# Please obtain data by yourself from HPA
expr <- fread("data/HPA_scRNA/rna_single_cell_type_tissue.tsv", header = TRUE,data.table=F)
expr <- expr[expr$Tissue == "brain",]
expr <- expr[grep("^COL4A[0-9]$|PLOD3",expr$`Gene name`),]
expr$id <- paste(expr$`Cell type`,expr$Cluster)
df <- spread(expr[,c("Gene name","id","nTPM")],key="Gene name",value = nTPM)
df <- df[mixedorder(df$id,decreasing = T),]
breaks <- seq(-2,2,by = .1)
cols <- colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(length(breaks))
ph <- pheatmap(t(df[,-1]),labels_col = df[,1],
               cluster_cols = T,
               cluster_rows = T,
               treeheight_col = 0,
               scale="row",
               color=cols,
               breaks=breaks)

##############################################################
# Single cell analysis ((GSE163577)) with Cellchat, Nichenet #
##############################################################

# Utility function
setIdent <- function(object, ident.use = NULL, levels = NULL){
  object@idents <- as.factor(object@meta[[ident.use]])
  if (!is.null(levels)) {
    object@idents <- factor(object@idents, levels = levels)
  }
  if (length(object@net) > 0) {
    if (all(dimnames(object@net$prob)[[1]] %in% levels(object@idents) )) {
      message("Reorder cell groups! ")
      idx <- match(dimnames(object@net$prob)[[1]], levels(object@idents))
      object@net$prob <- object@net$prob[idx, idx, ]
      object@net$pval <- object@net$pval[idx, idx, ]
      cat("The cell group order after reordering is ", dimnames(object@net$prob)[[1]],'\n')
    } else {
      message("Rename cell groups but do not change the order! ")
      cat("The cell group order before renaming is ", dimnames(object@net$prob)[[1]],'\n')
      dimnames(object@net$prob) <- list(levels(object@idents), levels(object@idents), dimnames(object@net$prob)[[3]])
      dimnames(object@net$pval) <- dimnames(object@net$prob)
      cat("The cell group order after renaming is ", dimnames(object@net$prob)[[1]],'\n')
    }
    
  }
  return(object)
}

# Please obtain data by yourself
refinfo <- read_excel("GSE163577/41586_2021_4369_MOESM3_ESM.xlsx")
refinfo <- as.data.frame(refinfo)
refinfo <- split(refinfo,f=refinfo$`Cell Type`)

stats2 <- readRDS("data/deg/metadeg_stat_randomTE.rda")
gmt <- gmtPathways("data/geneset/c2.cp.v2022.1.Hs.symbols.gmt")
gsea_result <- fgseaSimple(gmt,stats2,nperm = 1000)
gsea_sigresult <- gsea_result[gsea_result$padj <= .05,]
gsea_sigresult <- gsea_sigresult[order(gsea_sigresult$NES,decreasing = T),]

# collagen genes
col_gene <- fread("data/deg/collagen_meta.txt",data.table=F)
col_gene <- col_gene$genename[grep("COL4",col_gene$genename)]


#ECM genes
meta <- fread("data/GSE163577/atlas/meta.tsv",data.table = F)
obj <- Read10X(
  "data/GSE163577/atlas/",
  gene.column = 2,
  cell.column = 1,
  unique.features = TRUE,
  strip.suffix = FALSE
)

normObj <- CreateSeuratObject(counts = obj, project = "deconv", min.cells = 3, min.features = 200)
normObj <- NormalizeData(normObj, verbose = FALSE)
normObj <- FindVariableFeatures(normObj, selection.method = "vst",verbose = FALSE)
normObj <- ScaleData(normObj, verbose = FALSE)
normObj <- RunPCA(normObj,npcs = 30, verbose = FALSE)
normObj <- RunUMAP(normObj, reduction = "pca", dims = 1:30, verbose = FALSE)
ct <- meta$Cell_Type[match(names(Idents(normObj)),meta$Cell)]
ct <- gsub("Microglia/Mφ","Microglia",ct)
ct <- factor(ct,levels = sort(unique(ct)))
Idents(normObj) <- ct

# with condition
ct2 <- paste(meta$Cell_Type[match(colnames(normObj),meta$Cell)],meta$Treat[match(colnames(normObj),meta$Cell)],sep="_")
ct2 <- gsub("Microglia/Mφ","Microglia",ct2)
Idents(normObj) <- ct2

normObj_ad <- normObj[,grep("AD$",Idents(normObj))]
Idents(normObj_ad) <- gsub("_AD","",Idents(normObj_ad))
normObj_ctr <- normObj[,grep("Control$",Idents(normObj))]
Idents(normObj_ctr) <- gsub("_Control$","",Idents(normObj_ctr))

cellchat <- createCellChat(object = normObj_ad)
CellChatDB <- CellChatDB.human # use CellChatDB.mouse if running on mouse data
CellChatDB.use <- subsetDB(CellChatDB, search = "ECM-Receptor") # use Secreted Signaling
cellchat@DB <- CellChatDB.use
cellchat <- subsetData(cellchat)

cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)

cellchat <- computeCommunProb(cellchat)
cellchat_aggr <- aggregateNet(cellchat)

cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
groupSize <- as.numeric(table(cellchat@idents))

pathways.show <- c("COLLAGEN") 
vertex.receiver = seq(1,4) # a numeric vector. 
netVisual_aggregate(cellchat, signaling = pathways.show,  vertex.receiver = vertex.receiver)
par(mfrow=c(1,1))
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "circle")

# Neichenet
cellchat <- createCellChat(object = normObj_ad)
CellChatDB <- CellChatDB.human # use CellChatDB.mouse if running on mouse data
CellChatDB.use <- subsetDB(CellChatDB, search = "ECM-Receptor") # use Secreted Signaling

lr_network = readRDS(url("https://zenodo.org/record/7074291/files/lr_network_human_21122021.rds"))
ligand_target_matrix = readRDS(url("https://zenodo.org/record/7074291/files/ligand_target_matrix_nsga2r_final.rds"))
weighted_networks = readRDS("nichenet/weighted_networks_nsga2r_final.rds")

intCellChat <- data.frame(from = CellChatDB.use$interaction$ligand,
                          to = CellChatDB.use$interaction$receptor,
                          database = "Cellchat",
                          source = "Cellchat")
split_to <- strsplit(intCellChat$to, "_")
intCellChat <- data.frame(
  from = rep(intCellChat$from, times = sapply(split_to, length)),
  to = unlist(split_to),
  database = rep(intCellChat$source, times = sapply(split_to, length)),
  source = rep(intCellChat$source, times = sapply(split_to, length))
)

lr_network <- rbind(lr_network,intCellChat)



#####################################################
# Single cell anlysis using GSE163577 with nichenet #
#####################################################
# Please obtain data by yourself
meta <- fread("data/GSE163577/atlas/meta.tsv",data.table = F)
obj <- Read10X(
  "data/GSE163577/atlas/",
  gene.column = 2,
  cell.column = 1,
  unique.features = TRUE,
  strip.suffix = FALSE
)

normObj <- CreateSeuratObject(counts = obj, project = "deconv", min.cells = 3, min.features = 200)
normObj <- NormalizeData(normObj, verbose = FALSE)
normObj <- FindVariableFeatures(normObj, selection.method = "vst",verbose = FALSE)
normObj <- ScaleData(normObj, verbose = FALSE)
normObj <- RunPCA(normObj,npcs = 30, verbose = FALSE)
normObj <- RunUMAP(normObj, reduction = "pca", dims = 1:30, verbose = FALSE)

ct <- meta$Cell_Type[match(names(Idents(normObj)),meta$Cell)]
ct <- gsub("Microglia/Mφ","Microglia",ct)
ct <- factor(ct,levels = sort(unique(ct)))
Idents(normObj) <- ct

common_cell <- intersect(colnames(normObj),meta$Cell)
idx1 <- match(common_cell,rownames(normObj@meta.data))
idx2 <- match(common_cell,meta$Cell)

normObj <- normObj[,idx1]
meta <- meta[idx2,]
normObj@meta.data <- cbind(normObj@meta.data,meta)

seuratObj <- normObj
seuratObj@meta.data %>% head()
seuratObj %>% Idents() %>% table()

nichenet_output = nichenet_seuratobj_aggregate(
  seurat_obj = seuratObj, 
  receiver = "Astrocyte",
  condition_colname = "Treat",
  condition_oi = "AD", condition_reference = "Control", 
  sender = "all", 
  ligand_target_matrix = ligand_target_matrix,
  lr_network = lr_network,
  weighted_networks = weighted_networks)


nichenet_output$ligand_expression_dotplot
nichenet_output$ligand_differential_expression_heatmap
nichenet_output$ligand_differential_expression_heatmap
nichenet_output$ligand_target_heatmap
nichenet_output$ligand_activity_target_heatmap
nichenet_output$ligand_activities

seuratObj %>% Idents() %>% table()
nichenet_output2 = nichenet_seuratobj_aggregate(
  seurat_obj = seuratObj, 
  receiver = "Astrocyte",
  condition_colname = "Treat",
  condition_oi = "AD", condition_reference = "Control", 
  sender = "Oligo",
  ligand_target_matrix = ligand_target_matrix,
  lr_network = lr_network,
  weighted_networks = weighted_networks)

nichenet_output2$ligand_expression_dotplot
nichenet_output2$ligand_differential_expression_heatmap
nichenet_output2$top_targets
nichenet_output2$top_receptors
nichenet_output2$ligand_receptor_heatmap
nichenet_output2$ligand_target_heatmap
nichenet_output2$ligand_activity_target_heatmap
nichenet_output2$ligand_target_matrix
nichenet_output2$top_ligands

weighted_networks = readRDS("nichenet/weighted_networks_nsga2r_final.rds")
ligand_tf_matrix = readRDS(url("https://zenodo.org/record/7074291/files/ligand_tf_matrix_nsga2r_final.rds"))

lr_network = readRDS(url("https://zenodo.org/record/7074291/files/lr_network_human_21122021.rds"))
sig_network = readRDS(url("https://zenodo.org/record/7074291/files/signaling_network_human_21122021.rds"))
gr_network = readRDS(url("https://zenodo.org/record/7074291/files/gr_network_human_21122021.rds"))

ligands_all = "COL4A5" 
targets_all = c("BCL6","SGK1")

active_signaling_network = get_ligand_signaling_path(ligand_tf_matrix = ligand_tf_matrix, ligands_all = ligands_all, targets_all = targets_all, weighted_networks = weighted_networks)

active_signaling_network_min_max = active_signaling_network
active_signaling_network_min_max$sig = active_signaling_network_min_max$sig %>% mutate(weight = ((weight-min(weight))/(max(weight)-min(weight))) + 0.75)
active_signaling_network_min_max$gr = active_signaling_network_min_max$gr %>% mutate(weight = ((weight-min(weight))/(max(weight)-min(weight))) + 0.75)

library(DiagrammeRsvg)
library(DiagrammeR)
graph_min_max = diagrammer_format_signaling_graph(signaling_graph_list = active_signaling_network_min_max, ligands_all = ligands_all, targets_all = targets_all, sig_color = "indianred", gr_color = "steelblue")
gr <- render_graph(graph_min_max, layout = "tree")
export_graph(graph = gr, file_type = "png",
             file_name = "main_results/sc_bbb_nichenet_signaling.png",
             width = 8,height = 8)

normObj <- CreateSeuratObject(counts = obj, project = "deconv", min.cells = 3, min.features = 200)
normObj <- NormalizeData(normObj, verbose = FALSE)
normObj <- FindVariableFeatures(normObj, selection.method = "vst",verbose = FALSE)
normObj <- ScaleData(normObj, verbose = FALSE)
normObj <- RunPCA(normObj,npcs = 30, verbose = FALSE)
normObj <- RunUMAP(normObj, reduction = "pca", dims = 1:30, verbose = FALSE)
ct <- meta$Cell_Type[match(names(Idents(normObj)),meta$Cell)]
ct <- gsub("Microglia/Mφ","Microglia",ct)
ct <- factor(ct,levels = sort(unique(ct)))
Idents(normObj) <- ct

# receptor expression
p <- FeaturePlot_scCustom(tf_activities, features = c("CD44","SDC4","DDR2","ITGAV","ITGB8"),slot = "RNA",pt.size = .3,na_cutoff = 3,num_columns = 5)
p


###########################################################################
# Transcriptional factor activity estimation using dorothea and decoupleR #
###########################################################################
data(dorothea_hs, package = "dorothea")

regulons = dorothea_hs %>%
  filter(confidence %in% c("A", "B","C"))

sort(regulons[regulons$tf=="BCL6",]$target)

tf_activities <- dorothea::run_viper(normObj, regulons, 
                                     options =  list(method = "scale", minsize = 4, 
                                                     eset.filter = FALSE, cores = 60, 
                                                     verbose = FALSE))

p1 <- DimPlot(tf_activities, reduction = "umap", label = F) + 
  NoLegend() + ggtitle('Cell types')
p1 <- LabelClusters(p1, id = "ident",  fontface = "bold")


# BCL6
DefaultAssay(object = tf_activities) <- "dorothea"
p2 <- FeaturePlot_scCustom(tf_activities, features = "BCL6",slot = "dorothea",na_cutoff = 7.5,pt.size = .3) + ggtitle('BCL6 activity')

DefaultAssay(object = tf_activities) <- "RNA"
p3 <- FeaturePlot_scCustom(tf_activities, features = "CD44",slot = "RNA",pt.size = .3,na_cutoff = 3) + ggtitle('BCL6 expression')

p <- p2 | p3

# SGK1
p1 <- DimPlot(tf_activities, reduction = "umap", label = F) + 
  NoLegend() + ggtitle('Cell types')
p1 <- LabelClusters(p1, id = "ident",  fontface = "bold")

DefaultAssay(object = tf_activities) <- "dorothea"
p2 <- FeaturePlot_scCustom(tf_activities, features = "SGK1",slot = "dorothea",pt.size = .3,na_cutoff = 4.5) + ggtitle('SGK1 activity')

DefaultAssay(object = tf_activities) <- "RNA"
p3 <- FeaturePlot_scCustom(normObj, features = c("SGK1"),slot = "RNA",pt.size = .5) + ggtitle('SGK1 expression')

p <- p1|p3

# both
DefaultAssay(object = tf_activities) <- "RNA"
p1 <- FeaturePlot_scCustom(normObj, features = c("BCL6","SGK1"),slot = "RNA",pt.size = .3,na_cutoff = 2.5)

DefaultAssay(object = tf_activities) <- "dorothea"
p2 <- FeaturePlot_scCustom(tf_activities, features = "BCL6",slot = "dorothea",na_cutoff = 7,pt.size = .3) + ggtitle('BCL6 activity')

################################################
# Network topology identification with hdWGCNA #
################################################

meta <- fread("data/GSE163577/atlas/meta.tsv",data.table = F)
obj <- Read10X(
  "data/GSE163577/atlas/",
  gene.column = 2,
  cell.column = 1,
  unique.features = TRUE,
  strip.suffix = FALSE
)

library(Seurat)
normObj <- CreateSeuratObject(counts = obj, project = "deconv", min.cells = 3, min.features = 200)
normObj <- NormalizeData(normObj, verbose = FALSE)
normObj <- FindVariableFeatures(normObj, selection.method = "vst",verbose = FALSE)
normObj <- ScaleData(normObj, verbose = FALSE)
normObj <- RunPCA(normObj,npcs = 30, verbose = FALSE)
normObj <- RunUMAP(normObj, reduction = "pca", dims = 1:30, verbose = FALSE)
ct <- meta$Cell_Type[match(names(Idents(normObj)),meta$Cell)]
ct <- gsub("Microglia/Mφ","Microglia",ct)
ct <- factor(ct,levels = sort(unique(ct)))
Idents(normObj) <- ct

# with condition
ct2 <- paste(meta$Cell_Type[match(colnames(normObj),meta$Cell)],meta$Treat[match(colnames(normObj),meta$Cell)],sep="_")
ct2 <- gsub("Microglia/Mφ","Microglia",ct2)
#ct2 <- factor(ct2,levels = sort(unique(ct2)))
Idents(normObj) <- ct2
common_cell <- intersect(colnames(normObj),meta$Cell)
idx1 <- match(common_cell,rownames(normObj@meta.data))
idx2 <- match(common_cell,meta$Cell)

normObj <- normObj[,idx1]
meta <- meta[idx2,]
normObj@meta.data <- cbind(normObj@meta.data,meta)

normObj_ad <- normObj[,grep("AD$",Idents(normObj))]
Idents(normObj_ad) <- gsub("_AD","",Idents(normObj_ad))
normObj_ctr <- normObj[,grep("Control$",Idents(normObj))]
Idents(normObj_ctr) <- gsub("_Control$","",Idents(normObj_ctr))

seurat_obj <- SetupForWGCNA(
  normObj_ad,
  gene_select = "fraction",
  fraction = 0.05,
  wgcna_name = "scBBB"
)

# construct metacells  in each group
seurat_obj <- MetacellsByGroups(
  seurat_obj = seurat_obj,
  group.by = "Cell_Type", # specify the columns in seurat_obj@meta.data to group by
  reduction = 'umap', # select the dimensionality reduction to perform KNN on
  k = 25, # nearest-neighbors parameter
  max_shared = 10, # maximum number of shared cells between two metacells
  ident.group = 'Cell_Type' # set the Idents of the metacell seurat object
)

# normalize metacell expression matrix:
seurat_obj <- NormalizeMetacells(seurat_obj)

seurat_obj <- SetDatExpr(
  seurat_obj,
  group_name = "Astrocyte", # the name of the group of interest in the group.by column
  group.by='Cell_Type', # the metadata column containing the cell type info. This same column should have also been used in MetacellsByGroups
  assay = 'RNA', # using RNA assay
  slot = 'data' # using normalized data
)

# Test different soft powers:
seurat_obj <- TestSoftPowers(
  seurat_obj,
  networkType = 'unsigned' # you can also use "unsigned" or "signed hybrid"
)

# plot the results:
plot_list <- PlotSoftPowers(seurat_obj)

# assemble with patchwork
wrap_plots(plot_list, ncol=2)

power_table <- GetPowerTable(seurat_obj)
head(power_table)


seurat_obj <- ConstructNetwork(
  seurat_obj, soft_power=5,
  setDatExpr=FALSE,overwrite_tom = TRUE,
  tom_name = 'Astrocyte' # name of the topoligical overlap matrix written to disk
)

PlotDendrogram(seurat_obj, main='Astrocyte hdWGCNA Dendrogram')
TOM <- GetTOM(seurat_obj)

#####################################################
# Network propagation analysis (random walk restart)#
#####################################################

thresTOM <- quantile(TOM,probs = c(0.9))
subTOM <- TOM
subTOM <- ifelse(subTOM > thresTOM,1,0)
hgnc <- fread("https://ftp.ebi.ac.uk/pub/databases/genenames/hgnc/tsv/locus_types/gene_with_protein_product.txt")
idx <- rownames(subTOM) %in% hgnc$symbol[hgnc$locus_group=="protein-coding gene"]
subTOM <- subTOM[idx,idx]

saveRDS(subTOM,"main_results/hdWGCNA_subTOM.rds")
g <- graph_from_adjacency_matrix(subTOM,mode = "undirected")
MultiplexObject <- create.multiplex(list(net=g))
AdjMatrix <- compute.adjacency.matrix(MultiplexObject)
AdjMatrixNorm <- normalize.multiplex.adjacency(AdjMatrix)

data(dorothea_hs, package = "dorothea")

regulons = dorothea_hs %>%
  filter(confidence %in% c("A", "B","C"))

tf <- c("BCL6","SGK1")
SeedGene <- tf
SeedGene <- intersect(SeedGene,V(g)$name)
RWR_Results <- Random.Walk.Restart.Multiplex(AdjMatrixNorm,MultiplexObject,SeedGene)
TopResults <-
  create.multiplexNetwork.topResults(RWR_Results,MultiplexObject,k=30)

####################################################
## Enrichment analysis for gene module from hdWGCNA#
####################################################

data(pathwaysGMT)
pathwaysGMT <- gmtPathways("data/geneset/c5.go.bp.v2022.1.Hs.symbols.gmt")
geneset <- list(background = V(g)$name,topGene = setdiff(V(TopResults)$name,SeedGene))
fedupRes <- runFedup(geneset, pathwaysGMT)

fedupPlot <- fedupRes %>%
  bind_rows(.id = "set") %>%
  subset(pvalue < 0.01) %>%
  mutate(log10pvalue = -log10(pvalue)) %>%
  mutate(pathway = gsub("\\%.*", "", pathway)) %>%
  mutate(status = factor(status, levels = c("enriched", "depleted"))) %>%
  as.data.frame()

p <- plotDotPlot(
  df = fedupPlot,
  xVar = "log10pvalue",
  yVar = "pathway",
  xLab = "-log10(pvalue)",
  fillVar = "status",
  fillLab = "Enrichment\nstatus",
  sizeVar = "fold_enrichment",
  sizeLab = "Fold enrichment") +
  facet_grid("status", scales = "free", space = "free") +
  theme(strip.text.y = element_blank())
print(p)



