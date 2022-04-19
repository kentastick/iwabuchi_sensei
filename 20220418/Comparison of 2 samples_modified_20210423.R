# #install foreach packages
# 
# install.packages("foreach")
# install.packages("doParallel")
# install.packages("doMPI")
# install.packages("doSNOW")
# install.packages('Seurat')
# install.packages('robustbase')
# install.packages("ddalpha")
# install.packages("jackstraw")

#Nx1-seq analysis with Seurat package

library(foreach)
library(doParallel)
cl = detectCores()
cl = makeCluster(cl)
registerDoParallel(cl)

library(Seurat)
library(dplyr)
library(Matrix)
library(tidyverse)

#cellcycle <- readRDS(file = "Cellcycle/all.rds")
#Session -> Set working directory -> Choose directory#
# Load and label the time-course dataset
S1 = read.table("table_th1000andOver_list_kohara-NV2-tdt-bd_S1.txt", header=TRUE, row.names=1, sep="\t", quote="")
#S1 <- S1[!(rownames(S1) %in% cellcycle),]
#knee plot
tmp = colSums(S1)
tmp = sort(tmp, decreasing = TRUE)
plot(tmp, log="xy")
abline(h=0)
S1=S1[,colSums(S1)>0]  # 20000 reads

# detect cells of which min.genes over 200 
# set threshold as in the Seurat analysis (value of min.genes)

temp_func = function(i,j){sum(i[,j]>0)}
tmp = NULL
tmp = foreach(i=1:ncol(S1), .combine='c') %dopar% temp_func(S1, i)
tmp1 = tmp >= 0 
S1 = S1[,tmp1]

#check survivied cells
ncol(S1)
oridinalS1 = colnames(S1)
colnames(S1) = paste('S1_', c(1:ncol(S1)), sep='')
S1 = as.data.frame(S1)

#Comparing data???
S2 = read.table("table_th1000andOver_list_kohara-NP2-tdt-bd_S2.txt", header=TRUE, row.names=1, sep="\t", quote="")
#S2 <- S2[!(rownames(S2) %in% cellcycle),]
tmp = colSums(S2)
tmp = sort(tmp, decreasing = TRUE)
plot(tmp, log="xy")
abline(h=0)
S2 = S2[,colSums(S2)>0]

# detect cells of which min.genes over 200 
# set threshold as in the Seurat analysis (value of min.genes)
tmp = NULL
tmp = foreach(i=1:ncol(S2), .combine='c') %dopar% temp_func(S2, i)
tmp1 = tmp >= 0 
S2 = S2[,tmp1]

#check survivied cells
ncol(S2)
oridinalS2 = colnames(S2)
colnames(S2) = paste('S2_', c(1:ncol(S2)), sep='')
S2 = as.data.frame(S2)

bar_code = c(oridinalS1, oridinalS2)


#ここまでで、全てのデータをS1やS2としてmergeできるプログラムがあると良いのですが、スカスカなデータをたす意味を考えてました#
# combine tables  based on rownames(gene symbol)
merged_data = merge(S1, S2, by = "row.names", all = T)
merged_data[is.na(merged_data)] = 0
merged_data_1 = merged_data[,2:ncol(merged_data)]
rownames(merged_data_1) = merged_data[,1]

SAM_metadata = cbind(bar_code, colnames(merged_data_1))
SAM = CreateSeuratObject(raw.data = merged_data_1, min.cells = 0, min.genes = 0)
mito.genes = grep(pattern = "^MT.", x = rownames(x = SAM@data), value = TRUE)
percent.mito = Matrix::colSums(SAM@raw.data[mito.genes, ])/Matrix::colSums(SAM@raw.data)
SAM = AddMetaData(object = SAM, metadata = percent.mito, col.name = "percent.mito")

par(mfrow = c(1, 2))
GenePlot(object = SAM, gene1 = "nGene", gene2 = "percent.mito")

SAM = NormalizeData(object = SAM)
SAM = ScaleData(object = SAM, vars.to.regress = c("percent.mito"))
SAM = FindVariableGenes(object = SAM, mean.function = ExpMean, dispersion.function = LogVMR,
                        
                        x.low.cutoff = 0.03, x.high.cutoff = 8, y.cutoff = 0.8)

length(x = SAM@var.genes)

# perform PCA
SAM = RunPCA(object = SAM, pc.genes = SAM@var.genes, do.print = FALSE)
SAM = ProjectPCA(object = SAM, do.print = FALSE)
SAM = JackStraw(object = SAM, num.replicate = 100)
JackStrawPlot(object = SAM, PCs = 1:18)

#plot each principal component-associated genes
PCHeatmap(object = SAM, pc.use = 1:9, cells.use = 500, do.balanced = TRUE,
          
          label.columns = FALSE, use.full = FALSE)

PCHeatmap(object = SAM, pc.use = 10:18, cells.use = 500, do.balanced = TRUE,
          
          label.columns = FALSE, use.full = FALSE)

PCElbowPlot(object = SAM)

#clustering and detection of cellular subsets
#dims.use = 1: X, X could be read by Standard Deviation of PC
SAM = FindClusters(object = SAM, reduction.type = "pca", dims.use = 1:20,
                   resolution = 1.0, print.output = 0, save.SNN = TRUE, force.recalc = TRUE)
SAM = RunTSNE(object = SAM, dims.use = 1:20, do.fast = TRUE)
TSNEPlot(object = SAM, pt.size = 1.5)
TSNEPlot(object = SAM, cells.use = colnames(S1), pt.size = 1.5)
TSNEPlot(object = SAM, cells.use = colnames(S2), pt.size = 1.5)

#create oridinal cell barcode-condition-cell cluster table
cluster = SAM@meta.data
colnames(cluster) = c("nGene", "nReads", "Condition", "percent.mito", "Cell_Type_Seurat")
tmp = is.element(SAM_metadata[,2],rownames(cluster))
SAM_metadata = SAM_metadata[tmp,]
colnames(SAM_metadata)=c("cellSAM_oridinal", "CellIndex")
SAM_metadata = as.data.frame(SAM_metadata)
Cell_cluster_info = cbind(SAM_metadata, cluster)
rawData = t(SAM@raw.data)
Cell_cluster_info = cbind(Cell_cluster_info, rawData)

#Sort by Cell_Type_Seurat (Cell cluster)
tmp = Cell_cluster_info[order(Cell_cluster_info[,7]),] 
tmp = t(tmp)
tmp = as.data.frame(tmp)

#Export cell metadata x gene expression table

write.table(tmp,"Cellmetadata.txt",
            row.names=T, col.names=F, sep="\t", quote=F)

p1 <- TSNEPlot(object = SAM, do.return = TRUE, pt.size = 0.5)
plot_grid(p1)

tsne_coords = FetchData(SAM, c("tSNE_1", "tSNE_2"))
write.table(tsne_coords, "seurat_tSNE.txt", quote=F, sep="\t")


SAM.markers = FindAllMarkers(object = SAM, test.use = "bimod", only.pos = TRUE, min.pct = 0.25)
SAM.markers = SAM.markers [!grepl("^MT-|^RPL|^RPS|^RNA|^MTRNR", SAM.markers$gene),]
top10 = SAM.markers %>% group_by(cluster) %>% top_n(10, avg_logFC)
DoHeatmap(object = SAM, genes.use = top10$gene, slim.col.label = TRUE, remove.key = FALSE)
current.cluster.ids = c(0:17)

new.cluster.ids = c("0","1", "2", "3", "4", "5","6","7","8","9","10","11","12","13","14","15", "16", "17", "18", "19")
SAM@ident = plyr::mapvalues(x = SAM@ident, from = current.cluster.ids, to = new.cluster.ids)
TSNEPlot(object = SAM, do.label = TRUE, pt.size = 1.5, label.size = 6)

top10 = SAM.markers %>% group_by(cluster) %>% top_n(10, avg_logFC)

top10 = as.data.frame(top10)
gene = rownames(merged_data_1)
tmp = cbind (gene, merged_data_1)
marker_gene_table =  merge(top10, tmp, by="gene")    #
write.table(marker_gene_table,"~seurat_marker_genes.txt",
            row.names=F, sep="\t", quote=F)


#plot changes of conditions on t-SNE plot
cell_info = as.vector(SAM@cell.names)
tmp = cell_info
tmp_length=length(tmp)
tmp=strsplit(tmp, "\\_")
tmp=unlist(tmp)
tmp1 = seq(1,2*tmp_length,by=2)
experiment.info=data.frame(tmp[tmp1],row.names = SAM@cell.names)

#experiment.info = data.frame(c(rep(S1, ncol(S1)), rep(S2, ncol(S2))),row.names = SAM@cell.names)

colnames(experiment.info) = "Condition"
SAM = AddMetaData(object = SAM, metadata = experiment.info)
SAM = StashIdent(object = SAM, save.name = "CellType")
SAM = SetAllIdent(object = SAM, id = "Condition")
TSNEPlot(object = SAM, pt.size = 1.5, colors.use = c("black", "red"))
TSNEPlot(object = SAM, cells.use = colnames(S1), pt.size = 1.0)
TSNEPlot(object = SAM, cells.use = colnames(S2), pt.size = 1.0)

SAM = SetAllIdent(object = SAM, id = "CellType")     # Switch back to cell type labels

#export the number of each subsets
tmp = table(SAM@ident, SAM@meta.data$Condition)
write.table(tmp,"PASS.txt", row.names=T, col.names=F, sep="\t", quote=F)

#plot expression pattern of each marker genes on t-SNE plot
tmp = top10$gene


for (i in c(1:length(tmp))) {
  file.name = sprintf("%s.tiff", tmp[i])
  tiff(file.name, width = 512, height = 512)
  FeaturePlot(object = SAM, features.plot = tmp[i], cols.use = c("grey", "red"),
              reduction.use = "tsne", no.legend = FALSE)
  dev.off()
}

# extract RAW expression matrix for all subsets (perhaps, to load into another package)

for (cellType in new.cluster.ids) {         
  file.name = sprintf("%s.txt", cellType)
  raw.data = as.matrix(x = SAM@raw.data[, WhichCells(object = SAM, ident = cellType)])
  gene = rownames(SAM@raw.data)
  raw.data = cbind(gene, raw.data)
  write.table(raw.data, file.name, row.names=F, sep="\t", quote=F)      
}  

saveRDS(object = SAM, file = "SAM.rds")
readRDS(file="SAM.rds")

FeaturePlot(object = SAM, features.plot = c("Adgre1", "Cd68", "Cd9", "Cd163", "Trem2", "Hmox1", "Mmp9", "Ly6c1"), 
            cols.use = c("grey", "red"), reduction.use = "tsne")
VlnPlot(object = SAM, features.plot = c("Adgre1", "Cd68", "Cd9", "Cd163", "Trem2", "Hmox1", "Mmp9", "Ly6c1"), use.raw = TRUE, y.log = TRUE)   


FeaturePlot(object = SAM, features.plot = c("Cd3e", "Cd4", "Cd8a", "Lyz1", "Cd74", "Cd14", "Igj", "Fn1", "Apoe"), 
            cols.use = c("grey", "red"), reduction.use = "tsne")

VlnPlot(object = SAM, features.plot = c("Cd3e", "Cd4", "Cd8a", "Lyz1", "Cd74", "Cd14", "Igj", "Fn1", "Apoe"), use.raw = TRUE, y.log = TRUE)   

SAM = SetAllIdent(object = SAM, id = "Cluster")
FindAllMarkers(object = SAM,genes.use = c("Cd74"), logfc.threshold = 0, min.diff.pct = 0, min.pct = 0)

SAM = SetAllIdent(object = SAM, id = "CellType")
FindAllMarkers(object = SAM,genes.use = c("Cd74"), logfc.threshold = 0, min.diff.pct = 0, min.pct = 0)


FeaturePlot(object = SAM, features.plot = c("Adgre1"), 
            cols.use = c("grey", "red"), reduction.use = "tsne", pt.size = 1.5)

VlnPlot(object = SAM, features.plot = c("Adgre1"), use.raw = TRUE, y.log = TRUE)   

FeaturePlot(object = SAM, features.plot = c("Cd68"), 
            cols.use = c("grey", "red"), reduction.use = "tsne", pt.size = 1.5)

VlnPlot(object = SAM, features.plot = c("Cd68"), use.raw = TRUE, y.log = TRUE)   

FeaturePlot(object = SAM, features.plot = c("Cd9"), 
            cols.use = c("grey", "red"), reduction.use = "tsne", pt.size = 1.5)

VlnPlot(object = SAM, features.plot = c("Cd9"), use.raw = TRUE, y.log = TRUE)   

FeaturePlot(object = SAM, features.plot = c("Cd163"), 
            cols.use = c("grey", "red"), reduction.use = "tsne", pt.size = 1.5)

VlnPlot(object = SAM, features.plot = c("Cd163"), use.raw = TRUE, y.log = TRUE)   

FeaturePlot(object = SAM, features.plot = c("Trem2"), 
            cols.use = c("grey", "red"), reduction.use = "tsne", pt.size = 1.5)

VlnPlot(object = SAM, features.plot = c("Trem2"), use.raw = TRUE, y.log = TRUE)   

FeaturePlot(object = SAM, features.plot = c("Hmox1"), 
            cols.use = c("grey", "red"), reduction.use = "tsne", pt.size = 1.5)

VlnPlot(object = SAM, features.plot = c("Mmp9"), use.raw = TRUE, y.log = TRUE)   

FeaturePlot(object = SAM, features.plot = c("Marco"), 
            cols.use = c("grey", "red"), reduction.use = "tsne", pt.size = 1.5)

VlnPlot(object = SAM, features.plot = c("Marco"), use.raw = TRUE, y.log = TRUE)   

FeaturePlot(object = SAM, features.plot = c("Mmp9"), 
            cols.use = c("grey", "red"), reduction.use = "tsne", pt.size = 3.0)

SAM = SetAllIdent(object = SAM, id = "Condition")
FeaturePlot(object = SAM, features.plot = c("Mmp9"), 
            cols.use = c("grey", "red"), reduction.use = "tsne", pt.size = 3.0)
FeaturePlot(object = SAM, features.plot = c("Mmp9"), 
            cols.use = c("grey", "red"), reduction.use = "tsne", cells.use = colnames(S1), pt.size = 3.0)
FeaturePlot(object = subcluster, features.plot = c("Mmp9"), 
            cols.use = c("grey", "red"), reduction.use = "tsne", cells.use = colnames(S2), pt.size = 3.0)


SAM = SetAllIdent(object = SAM, id = "CellType")
subcluster <- SubsetData(SAM, ident.use = "6")
subcluster = SetAllIdent(object = subcluster, id = "Condition")

subcluster@meta.data$Mmp9 = FetchData(subcluster, "Mmp9")

subcluster@meta.data$Mmp9_group = subcluster@meta.data %>% mutate(Mmp9_group = case_when(Mmp9>0&Condition == "S1"~"Mmp9_posi_S1",
                                                                                         Mmp9==0&Condition == "S1"~"Mmp9_nega_S1",
                                                                                         Mmp9>0&Condition == "S2"~"Mmp9_posi_S2",
                                                                                         Mmp9==0&Condition == "S2"~"Mmp9_nega_S2"
))  %>%    pull(Mmp9_group)  %>% 
  fct_relevel(c("Mmp9_posi_S1", "Mmp9_nega_S1", "Mmp9_posi_S2", "Mmp9_nega_S2"))

subcluster = SetAllIdent(object = subcluster, id = "Mmp9_group")

TSNEPlot(object = subcluster, colors.use = c("grey", "grey", "black", "red"),  pt.size = 2.0) 

TSNEPlot(object = subcluster, plot.order = c("Mmp9_posi_S1", "Mmp9_posi_S2", "Mmp9_nega_S1", "Mmp9_nega_S2"),
         colors.use = c("grey", "grey", "red", "black"), pt.size = 2.0) 



FeaturePlot(object = SAM, features.plot = c("Marco"),no.axes =TRUE,overlay = FALSE, dark.theme = TRUE,
            cols.use = c("grey", "red"), reduction.use = "tsne", pt.size = 1.5)

FeaturePlot(object = SAM, features.plot = c("Marco"),　cols.use = c("grey", "red"), reduction.use = "tsne", pt.size = 1.5)+theme_classic()

FeaturePlot(object = SAM, features.plot = c("Marco"),　cols.use = c("grey", "red"), reduction.use = "tsne", pt.size = 1.5)+theme(panel.background = element_rect("white"))


subcluster@meta.data$Mmp9 = FetchData(subcluster, "Mmp9")

subcluster@meta.data$Mmp9_group = subcluster@meta.data %>% mutate(Mmp9_group = case_when(Mmp9>0&Condition == "S1"~"Mmp9_posi_S1",
                                                                                         Mmp9==0&Condition == "S1"~"Mmp9_nega_S1",
                                                                                         Mmp9>0&Condition == "S2"~"Mmp9_posi_S2",
                                                                                         Mmp9==0&Condition == "S2"~"Mmp9_nega_S2"
)) > pull(Mmp9_group)> 
  fct_relevel(c("Mmp9_posi_S1", "Mmp9_nega_S1", "Mmp9_posi_S2", "Mmp9_nega_S2"))

subcluster = SetAllIdent(object = subcluster, id = "Mmp9_group")

TSNEPlot(object = subcluster, colors.use = c("#FF4040", "#FFB6C1", "#00CD00", "#9BCD9B"),  pt.size = 2.0) 



subcluster = SetAllIdent(object = subcluster, id = "Condition")
VlnPlot(object = subcluster, features.plot = c("Mmp9"),use.raw = TRUE, y.log = TRUE)  
FeaturePlot(object = subcluster, features.plot = c("Mmp9"), 
            cols.use = c("grey", "red"), reduction.use = "tsne", pt.size = 3.0)


VlnPlot(object = SAM, features.plot = c("Mmp9"), use.raw = TRUE, y.log = TRUE)   

FeaturePlot(object = SAM, features.plot = c("Ly6c1"), 
            cols.use = c("grey", "red"), reduction.use = "tsne", pt.size = 1.5)

VlnPlot(object = SAM, features.plot = c("Ly6c1"), use.raw = TRUE, y.log = TRUE)   




##########################################################################
#pan macrophagel#
FeaturePlot(object = SAM, features.plot = c("Lyz1", "Lyz2", "Adgre1"), 
            cols.use = c("grey", "red"), reduction.use = "tsne", overlay = FALSE) 

#Macrophage M1: Cd80, Cd86, Cd68, H2-Ab1, Il6, Il15, Nos2, Tlr2, Tlr4, Ccl4, Ccl8, Cxcl11#
FeaturePlot(object = SAM, features.plot = c("Cd80", "Cd86", "Cd68", "H2-Ab1", "Il6", "Il15", "Nos2", "Tlr2", "Tlr4", "Ccl4", "Ccl8", "Cxcl11"), 
            cols.use = c("grey", "red"), reduction.use = "tsne", overlay = FALSE) 

#Marcophage M2a: Cd163, Cxcr1, Cxcr2, Il10, Tgm2   Il12 and MHC class low#
FeaturePlot(object = SAM, features.plot = c("Cd163", "Cxcr2", "Il10", "Tgm2"), 
            cols.use = c("grey", "red"), reduction.use = "tsne", overlay = FALSE) 

#Marcophage M2b: high classII, Cd86, Ccl1, Il4, Irf4, SOCS3#
FeaturePlot(object = SAM, features.plot = c("Cd86", "Ccl1", "Il4", "Irf4", "Socs3"), 
            cols.use = c("grey", "red"), reduction.use = "tsne", overlay = FALSE) 

#Marcophage M2c: Cd163, Tlr1, Trl8, Ccr2, Ccl16#
FeaturePlot(object = SAM, features.plot = c("Cd163", "Tlr1", "Tlr8", "Ccr2"), 
            cols.use = c("grey", "red"), reduction.use = "tsne", overlay = FALSE) 

#Marcophage TAM: Ccl2, Cd31(pcam1), Cd68, PCNA, VEGF8(vegfa)#
FeaturePlot(object = SAM, features.plot = c("Ccl2", "Pecam1", "Cd68", "Pcna", "Vegfa"), 
            cols.use = c("grey", "red"), reduction.use = "tsne", overlay = FALSE) 

#DC marker: Cd11b(Itgam), Cd11c(Itgax), Cd24, Cd115(Csf1r), Cd117(Kit), Tlr7, Tlr9, Ly-6c
FeaturePlot(object = SAM, features.plot = c("Itgam", "Itgax", "Csf1r", "Kit", "Tlr7", "Ly6c1"), 
            cols.use = c("grey", "red"), reduction.use = "tsne", overlay = FALSE) 

#Monocyte marker: Cd11b, Cd16(Fcgr3), Cd31, Cd43(Spn), Cd44, Cd45(Ptprc)
FeaturePlot(object = SAM, features.plot = c("Fcgr3", "Pecam1", "Spn", "Cd44", "Ptprc"), 
            cols.use = c("grey", "red"), reduction.use = "tsne", overlay = FALSE) 

#T-cell#
FeaturePlot(object = SAM, features.plot = c("Cd3e", "Cd4", "Cd8a", "Nkg7"), 
            cols.use = c("grey", "red"), reduction.use = "tsne", overlay = FALSE) 

#Th1: Ifng, Stat4, T-bet(Tbx21)#
#Th2: Il4, Stat6, Gata3#
#Treg: tgfb1(TGFbeta), Foxp3#
#Th17: Il17a, RORgT(Rorc)#
FeaturePlot(object = SAM, features.plot = c("Stat4", "Tbx21", "Il4", "Stat6", "Gata3", "Tgfb1", "Foxp3", "Rorc"), 
            cols.use = c("grey", "red"), reduction.use = "tsne", overlay = FALSE) 

#B-cell: Cd19, Cd20(mature)(Ms4a1), Cd22(mature), Cd23(activated)(Fcer2a), Cd40(B cell), Cd138(plasma cell)(Sdc1)#
#Memory B-cell: Cd21(Cr2), Cd24(Cd24a),Cd1c#
FeaturePlot(object = SAM, features.plot = c("Cd22", "Fcer2a", "Cd40", "Sdc1"), 
            cols.use = c("grey", "red"), reduction.use = "tsne", overlay = FALSE) 

#Fibroblast: Vim, aSMA(Acta2), Hsp47(Serpinh1), S100a4, Cd26(Dpp4), Dlk1, Pdgfra#
FeaturePlot(object = SAM, features.plot = c("Acta2", "Vim", "Serpinh1", "S100a4", "Dpp4", "Pdgfra"), 
            cols.use = c("grey", "red"), reduction.use = "tsne", overlay = FALSE)

#Endothelial: Cd31(Pecam1), Cd34, Icam1, Cd45(Ptprc), Cd36, Cd39(Entpd1)#
FeaturePlot(object = SAM, features.plot = c("Pecam1", "Cd34", "Icam1", "Ptprc", "Cd36", "Entpd1"), 
            cols.use = c("grey", "red"), reduction.use = "tsne", overlay = FALSE)

#Epithelial cell: Krt19, Epcam#
FeaturePlot(object = SAM, features.plot = c("Krt19", "Epcam"), 
            cols.use = c("grey", "red"), reduction.use = "tsne", overlay = FALSE)

#Squamouse cell: Krt14, Epithelai carcinoma: Survivin(Birc5) #
#Proliferating malignant epithelial cells: Krt18#
FeaturePlot(object = SAM, features.plot = c("Krt14", "Birc5", "Krt18"), 
            cols.use = c("grey", "red"), reduction.use = "tsne", overlay = FALSE)

#Smooth muscle cells: aSMA(Acta1), VE-Cadherin(Cdh5), Cald1, Calponin1(Cnn1), Hexim1, Tagln#
FeaturePlot(object = SAM, features.plot = c("Acta1", "Cdh5", "Cald1", "Cnn1", "Hexim1", "Tagln"), 
            cols.use = c("grey", "red"), reduction.use = "tsne", overlay = FALSE)




SAM = SetAllIdent(object = SAM, id = "Condition")
VlnPlot(object = SAM, features.plot = c("Ndrg1"), use.raw = TRUE, y.log = TRUE, cols.use = c("black", "red"))   
VlnPlot(object = SAM, features.plot = c("Fabp5"), use.raw = TRUE, y.log = TRUE, cols.use = c("black", "red"))   

SAM = SetAllIdent(object = SAM, id = "CellType")
VlnPlot(object = SAM, features.plot = c("Fabp5"), use.raw = TRUE, y.log = TRUE)  
VlnPlot(object = SAM, features.plot = c("Ndrg1"), use.raw = TRUE, y.log = TRUE)  

library(tidyverse)
SAM = SetAllIdent(object = SAM, id = "Condition")
TSNEPlot(object = SAM, colors.use = c ("black","red"))
comparison_avg = FindMarkers(object = SAM, ident.1 = "S1", ident.2 = "S2")
comparison_avg = comparison_avg %>% rownames_to_column(var = "gene")
comparison_avg = comparison_avg %>% mutate(cluster = if_else(avg_logFC>0, "id1_posi", "id2_posi")) 
comparison_avg = comparison_avg [!grepl("^MT-|^Rpl|^Rps|^Rna|^Mtrnr|^LOC", comparison_avg$gene),]
comparison_avg = comparison_avg %>% group_by(cluster) %>% do(head(., 20)) %>% arrange(cluster)
DoHeatmap(object = SAM, cells.use = WhichCells(object = SAM, ident = c("S1", "S2")), genes.use = comparison_avg$gene,
          slim.col.label = TRUE, remove.key = T)

SAM = SetAllIdent(object = SAM, id = "CellType")
subcluster <- SubsetData(SAM, ident.use = "6")
subcluster = SetAllIdent(object = subcluster, id = "Condition")
TSNEPlot(object = subcluster, colors.use = c ("black","red"), pt.size = 2.0)
FeaturePlot(object = subcluster, features.plot = c("Mmp9"), 
            cols.use = c("grey", "red"), reduction.use = "tsne",pt.size = 3.0)
FeaturePlot(object = subcluster, features.plot = c("Mmp9"), 
            cells.use = subcluster@cell.names[subcluster@meta.data$Condition == "S1"],
            cols.use = c("grey", "red"),
            reduction.use = "tsne", pt.size = 3.0)


comparison_avg = FindMarkers(object = subcluster, ident.1 = "S1", ident.2 = "S2")
comparison_avg = comparison_avg %>% rownames_to_column(var = "gene")
comparison_avg = comparison_avg %>% mutate(cluster = if_else(avg_logFC>0, "id1_posi", "id2_posi")) 
comparison_avg = comparison_avg [!grepl("^MT-|^Rpl|^Rps|^Rna|^Mtrnr|^LOC", comparison_avg$gene),]
comparison_avg = comparison_avg %>% group_by(cluster) %>% do(head(., 10)) %>% arrange(cluster)
DoHeatmap(object = subcluster, cells.use = WhichCells(object = subcluster, ident = c("S1", "S2")), genes.use = comparison_avg$gene,
          slim.col.label = TRUE, remove.key = T)



VlnPlot(object = SAM, features.plot = c("Krt19", "Cd74", "H2-Aa", "H2-Eb1", "Gbp2b", "B2m", "Cebpb", "Tpt1", "Cxcl9"), 
            cols.use = c("black", "red"), use.raw = TRUE, y.log = TRUE)

label = TSNEPlot(SAM, do.identify = TRUE)
SUBCLUSTER = SetIdent(SAM, cells.use = label, ident.use = "SUB1")
TSNEPlot(SUBCLUSTER)
label = TSNEPlot(SUBCLUSTER, do.identify = TRUE)
SUBCLUSTER = SetIdent(SUBCLUSTER, cells.use = label, ident.use = "SUB2")
TSNEPlot(SUBCLUSTER)
comparison_avg = FindMarkers(object =SUBCLUSTER, ident.1 = "SUB1", ident.2 = "SUB2")
comparison_avg = comparison_avg %>% rownames_to_column(var = "gene")
comparison_avg = comparison_avg [!grepl("^MT-|^RPL|^RPS|^RNA|^MTRNR", comparison_avg$gene),]
comparison_avg = comparison_avg %>% mutate(cluster = if_else(avg_logFC > 0, "SUB1", "SUB2")) 
comparison_avg = comparison_avg %>% group_by(cluster) %>% do(head(., 30)) %>% arrange(cluster)
DoHeatmap(object = SUBCLUSTER, cells.use = WhichCells(object = SUBCLUSTER, ident = c("SUB1", "SUB2")), genes.use = comparison_avg$gene,
          slim.col.label = TRUE, remove.key = T)

FeaturePlot(object = SAM, features.plot = c("Ndrg1"), 
            cols.use = c("grey", "red"), reduction.use = "tsne", pt.size = 1.0)



# count gene percent ------------------------------------------------------

FetchData(SAM, c("Adgre1", "orig.ident", "res.1")) %>% group_by(orig.ident, res.1) %>% 
  summarise(count = sum(Adgre1>0), total = n(), ratio = count/total)

library(clipr)
FetchData(SAM, vars.all = c("Adgre1", "orig.ident","res.1")) %>% group_by(orig.ident, res.1) %>%
  summarise(count = sum(Adgre1>0), total = n(), ratio = round(count/total, 3)) %>% print(as_tibble(iris), n = 100) %>% write_clip　

# S1 vs S2 comparison in cluster 5 -----------------------------------------


SAM = SetAllIdent(object = SAM, id = "CellType")
subcluster <- SubsetData(SAM, ident.use = "5")

#これが追加したコードです
subcluster = SetAllIdent(object = subcluster, id = "Condition")

TSNEPlot(object = subcluster, colors.use = c ("black","red"))
comparison_avg = FindMarkers(object = subcluster, ident.1 = "S1", ident.2 = "S2")
comparison_avg = comparison_avg %>% rownames_to_column(var = "gene")
comparison_avg = comparison_avg %>% mutate(cluster = if_else(avg_logFC>0, "id1_posi", "id2_posi"))
comparison_avg = comparison_avg [!grepl("^MT-|^Rpl|^Rps|^Rna|^Mtrnr|^LOC", comparison_avg$gene),]
comparison_avg = comparison_avg %>% group_by(cluster) %>% do(head(., 20)) %>% arrange(cluster)
DoHeatmap(object = subcluster, cells.use = WhichCells(object = subcluster, ident = c("S1", "S2")), genes.use = comparison_avg$gene,
          slim.col.label = TRUE, remove.key = T)






