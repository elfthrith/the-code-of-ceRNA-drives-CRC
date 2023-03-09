#Differential expression analysis
#input gene expression profiles:a1(cancer vs polyps)
#output Differential expression (.txt)
library(limma)
a1[,1:6]<-lapply(a1[,1:6],as.numeric)
group<-factor(c("P","C","P","C","C","P"))
library(dplyr)
design=model.matrix(~0+factor(group))
colnames(design)<-levels(group)
rownames(design)<-colnames(a1)
contrast.matrix<-makeContrasts("C-N",levels=design)
fit<-lmFit(a1,design)
fit2<-contrasts.fit(fit,contrast.matrix)
fit2<-eBayes(fit2,trend = T)
tempoutput<-topTable(fit2,coef=1,n=Inf)
tempoutput=cbind(Ensembl=row.names(tempoutput),tempoutput)
row.names(tempoutput)=NULL
output<-subset(tempoutput,P.Value<0.05&(logFC>1|logFC< -1))

#The volcano plots
#input Differential expression results. The volcano plots of the article results can be replicated through supplementary materials: DERNAs data.xlsx
#output the volcano plots
library(ggplot2)
tempoutput[which(tempoutput$logFC >= 1 & tempoutput$P.Value < 0.05),'sig'] <- 'up'
tempoutput[which(tempoutput$logFC <= -1 &tempoutput$P.Value < 0.05),'sig'] <- 'down'
tempoutput[which(abs(tempoutput$logFC ) < 1 | tempoutput$P.Value >= 0.05),'sig'] <- 'none'
rescp_select <- subset(tempoutput, sig %in% c('up', 'down'))
library(ggplot2)
p <- ggplot(data =tempoutput, aes(x =logFC, y = -log10(P.Value), color = sig)) +
  geom_point(size = 1.5) +  #»æÖÆÉ¢µãÍ¼
  scale_color_manual(values = c('red', 'black', 'green'), limits = c('up', 'none', 'down')) +  
  labs(x = 'log2 Fold Change', y = '-log10 p-value', title = 'cancer vs polyp_mRNA', color = '') + 
  theme(plot.title = element_text(hjust = 0.5, size = 14), panel.grid = element_blank(), 
        panel.background = element_rect(color = 'black', fill = 'transparent'), 
        legend.key = element_rect(fill = 'transparent')) +
  geom_vline(xintercept = c(-1, 1), lty = 3, color = 'black') +
  geom_hline(yintercept = -log10(0.05), lty = 3, color = 'black') +
  xlim(-7,7) + ylim(0,6) 

# the pheatmap plots
#input gene expression profiles of DERNAs
#output the pheatmap plots
library(gplots)
fix(a1)
cg=rescp_select[,1]
group<-factor(c("P","C","P","C","C","P"))
names(group)=colnames(a1)
group=as.data.frame(group)
fix(group)
DEG_expr_matr<-a1[cg,]
DEG_expr_matr<-DEG_expr_matr[,c(1,3,6,4,2,5)]
z_score_matrix<-t(scale(t(DEG_expr_matr)))
annotation_color=list(group=c(C="orange",P="green"))
library(pheatmap)
p<-pheatmap(z_score_matrix,cluter_rows=T,show_rownames=F,show_colnames=T,fontsize_row=12,cluster_cols=T,angle_col=0,
            annotation_col=group,annotation_colors=annotation_color,color=colorRampPalette(c("blue", "black", "red"))(50),
            legend=T,main = "cancer vs polyp_mRNA")

#2.3	GSEA Enrichment Analysis
#input geneList: DERNAs
library(ReactomePA)
library(tidyverse)
library(data.table)
library(org.Hs.eg.db)
library(clusterProfiler)
library(biomaRt)
library(enrichplot)
Go_gseresult <- gseGO(geneList, 'org.Hs.eg.db', keyType = "ENTREZID", ont="all", nPerm = 1000, minGSSize = 10, maxGSSize = 1000, pvalueCutoff=1)

#2.4	WGCNA Analysis and Identification of Clinically Significant Modules
#input gene expression profiles of validated lncRNA+mRNA of TCGA samples
library(data.table)
datExpr0<-fread("datExpr0_log2(tpm+1).txt")
datExpr0<-as.data.frame(datExpr0)
rownames(datExpr0)<-datExpr0[,1]
datExpr0<-datExpr0[,-1]
library(WGCNA) 
options(stringsAsFactors = FALSE) 
data.mat<-as.data.frame(t(datExpr0))
gsg = goodSamplesGenes(data.mat, verbose = 3);
gsg$allOK
sampleTree <- hclust(dist(data.mat), method = "average")
sizeGrWindow(12,9) 
par(cex = 0.6);
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
     cex.axis = 1.5, cex.main =2)
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)
library(stringr)
traitDate<-as.data.frame(rownames(data.mat))
traitDate$group<-ifelse(as.numeric(str_sub(traitDate$`rownames(data.mat)`,14,15)) < 10,'tumor','normal')
femaleSamples = rownames(data.mat);
traitRows = match(femaleSamples, traitDate$`rownames(data.mat)`);
collectGarbage()
allowWGCNAThreads()
type <- "unsigned"
powers <- c(1:10, seq(from = 12, to=30, by=2))
sft <- pickSoftThreshold(
  data.mat, powerVector=powers, 
  networkType=type, verbose=3
)
sizeGrWindow(9, 5)
par(mfrow = c(1,2))
cex1 = 0.9 
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
abline(h=0.85,col="red") 
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
sft$powerEstimate
cor = WGCNA::cor
net = blockwiseModules(data.mat, power = 12,
                       TOMType = "unsigned", minModuleSize = 30,
                       reassignThreshold = 0, mergeCutHeight = 0.1,
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       saveTOMs = TRUE,
                       saveTOMFileBase = "TOM",
                       verbose = 3)
cor = stats::cor 
table(net$colors)
labels2colors(net$colors)
moduleColors <- labels2colors(net$colors)
sizeGrWindow(12,9)
library(ggplot2)
plotDendroAndColors(
  net$dendrograms[[1]],
  moduleColors[net$blockGenes[[1]]],
  "Module colors",
  dendroLabels = FALSE,
  hang = 0.03,
  addGuide = TRUE,
  guideHang = 0.05
)
nGenes = ncol(data.mat);
nSamples = nrow(data.mat);
MEs0 = moduleEigengenes(data.mat, moduleColors)$eigengenes
MEs = orderMEs(MEs0)
datTraits<-traitDate
design <- model.matrix(~0 + datTraits$group)
dimnames(design) <- list(datTraits$`rownames(data.mat)`, sort(unique(datTraits$group)))
design <- design[rownames(MEs),]
modTraitCor <- cor(MEs, design, use = "p")
modTraitP <- corPvalueStudent(modTraitCor, dim(datTraits)[1])
textMatrix <- paste0(signif(modTraitCor, 2), "\n(", signif(modTraitP, 1), ")")
dim(textMatrix) <- dim(modTraitCor)
sizeGrWindow(12,9)
par(mar = c(4,7,4,4))
labeledHeatmap(
  Matrix = modTraitCor,
  xLabels = colnames(design),
  yLabels = colnames(MEs),
  cex.lab = 0.9,
  ySymbols = colnames(MEs),
  colorLabels = FALSE,
  colors = blueWhiteRed(50),
  textMatrix = textMatrix,
  setStdMargins = FALSE,
  cex.text = 1,
  zlim = c(-1, 1),
  main = paste("Module-trait relationships")
)

geneModuleMembership <- cor(data.mat, MEs, use = "p")
MMPvalue <- corPvalueStudent(geneModuleMembership, nSamples)
geneSignificanceCor <- cor(data.mat, design, use = "p")
geneSignificanceP <- corPvalueStudent(geneSignificanceCor, nSamples)
module <- "brown"
column <- paste0("ME", module)
moduleGenes <- names(net$colors)[which(moduleColors == module)]
MM <- abs(geneModuleMembership[moduleGenes, column])
GS <- abs(geneSignificanceCor[moduleGenes, 1])
verboseScatterplot(
  MM, GS,
  xlab = paste("Module Membership in", module, "module"),
  ylab = "Gene significance for proliferating",
  main = paste("Module membership vs. gene significance\n"),
  abline = TRUE,
  pch = 21,
  cex.main = 1.2,
  cex.lab = 1.2,
  cex.axis = 1.2,
  col = "black",
  bg = module
)
moduleGenes[(GS > 0.7 & MM > 0.7)]
