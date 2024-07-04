#AML data clean-up

library(maftools)

laml <- read.maf(maf = "data_mutations.txt", useAll = TRUE)
mtx <- mutCountMatrix(laml, includeSyn = FALSE, countOnly = NULL, removeNonMutated = FALSE)

#transpose mtx to have genes in columns and samples in row
mtx <- t(mtx)
#Convert counts to binary
mtx.b <- apply(mtx, 2, function(x) ifelse(x > 0, 1, x)) # So 0 = no, 1 =yes

mtx.b <- mtx.b[,-which(colSums(mtx.b) < 2)] #Remove mutations seen in less than 2 samples
mtx.b <- mtx.b[-which(rowSums(mtx.b) < 1),] #Remove patients with no mutations

###
#Run VICatMix
library(mcclust.ext)
source("VariationalMixtureModelVarSel.R")
set.seed(205)

library(foreach)
library(doParallel)
library(doRNG)
registerDoParallel(30)
amlclust <- foreach(i = 1:30) %dorng% {
  vs <- mixturemodelvarsel(mtx.b, 20, 0.01, 2, 2000, 0.00000005)
  finalresult <- list()
  #want to record cluster labels
  finalresult$VarSel <- vs$model$labels
  finalresult$SelectedVars <- vs$model$c
  finalresult
}

resultforpsm <- lapply(amlclust, "[[", 1)
p1 <- t(matrix(unlist(resultforpsm), 185, 30))
psm <- comp.psm(p1)
VIcomp <- minVI(psm, method = 'comp',max.k = 20)$cl

p2 <- t(matrix(unlist(lapply(amlclust, "[[", 2)), 151, 30))
p2[p2 <= 0.5] <- 0
p2[p2 > 0.5] <- 1

result <- vector(mode = 'numeric', length = 151)
for (i in 1:42){
  result[i] <- sum(p2[,i]) / 30
}
ThresVars <- colnames(mtx.b)[result > 0.949999]

amlclust$AvgComp <- VIcomp
amlclust$SelectedVars <- ThresVars


##
#Plots 

#Heatmaps

currentAnnotationRow <- data.frame(
  Cluster = factor(amlclust$AvgComp)
)
rownames(currentAnnotationRow) <- rownames(mtx.b) 

#Heatmap with selected variables
pheatmap(mtx.b[order(amlclust$AvgComp),amlclust$SelectedVars], cluster_rows = F,
         color = colorRampPalette(colors = c("white", "black"))(2),
         annotation_row = currentAnnotationRow, cluster_cols = F, show_rownames = F, show_colnames = T)
#Heatmap with all variables
pheatmap(mtx.b[order(amlclust$AvgComp),], cluster_rows = F,
         color = colorRampPalette(colors = c("white", "black"))(2),
         annotation_row = currentAnnotationRow, cluster_cols = F, show_rownames = F, show_colnames = T)

#ORA analysis
aml_genes6 <- c("DNMT3A", "NPM1", "FLT3", "IDH2", "RUNX1", "TP53")

organism = "org.Hs.eg.db"
BiocManager::install(organism, character.only = TRUE)
library(organism, character.only = TRUE)
library(enrichplot)
library(DOSE)

aml_entrez <- select(org.Hs.eg.db, 
                     keys = aml_genes6,
                     columns = c("ENTREZID", "SYMBOL"),
                     keytype = "SYMBOL")

all_aml_entrez <- select(org.Hs.eg.db, 
                         keys = colnames(mtx.b),
                         columns = c("ENTREZID", "SYMBOL"),
                         keytype = "SYMBOL")

library(cowplot)

AML_DO <- enrichDO(aml_entrez$ENTREZID, universe = all_aml_entrez$ENTREZID)
dotplot(AML_DO, title = "Disease Ontology ORA", font.size = 10) + xlim(0.65, 1.02) + scale_size_continuous(range  = c(3, 7), 
                                                                                                           limits = c(4, 6), 
                                                                                                           breaks = c(4, 5, 6))
AML_NCG <- enrichNCG(aml_entrez$ENTREZID, universe = all_aml_entrez$ENTREZID)
ncg <- dotplot(AML_NCG, title = "Network of Cancer Genes", font.size = 10)
AML_DGN <- enrichDGN(aml_entrez$ENTREZID, universe = all_aml_entrez$ENTREZID)
dgn <- dotplot(AML_DGN, title = "DisGeNET", font.size = 10) + xlim(0.81, 1.02) + scale_size_continuous(range  = c(3, 7), 
                                                                                                       limits = c(4, 6), 
                                                                                                       breaks = c(4, 5, 6))
cowplot::plot_grid(ncg, dgn, labels=LETTERS[1:2], rel_widths=c(.8, .8))




