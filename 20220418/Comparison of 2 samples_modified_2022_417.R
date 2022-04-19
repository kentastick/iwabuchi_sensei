#install foreach packages

# install.packages("foreach")
# install.packages("doParallel")
# install.packages("doMPI")
# install.packages("doSNOW")
# # install.packages('Seurat')
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

#S1 = read.table("table_th1000andOver_list_kohara-NV2-tdt-bd_S1.txt", header=TRUE, row.names=1, sep="\t", quote="")


# data という名のフォルダ入っている処理するデータのリストを作る
# マトリックスデータはdata というフォルダに入れてください
setwd("20220418/")

file_list = list.files("data", full.names = T)
file_list


# 画像フォルダの作成 ここにimage file が保存されます
dir.create("image")


# それぞれのファイルに名前をつける、
names(file_list) = c("S1", "S10", "S2", "S3", "S4", "S5", "S6", "S7", "S8", "S9")


# 処理を関数にしてまとめました、これでいくつものデータを使用したときに処理が楽になるのではと思います
do_process <- function(path_to_data, name) {
  df = read.table(path_to_data, header=TRUE, row.names=1, sep="\t", quote="")
  tmp = colSums(df)
  tmp = sort(tmp, decreasing = TRUE)
  
  # count plot の保存、image フォルダに保存されます
  jpeg(file.path("image", paste0(name, "_count.jpg")))
  plot(tmp, log="xy")
  abline(h=0)
  dev.off()
  
  df=df[,colSums(df)>0]  # 20000 reads
  
  # detect cells of which min.genes over 200 
  # set threshold as in the Seurat analysis (value of min.genes)
  
  temp_func = function(i,j){sum(i[,j]>0)}
  tmp = NULL
  tmp = foreach(i=1:ncol(df), .combine='c') %dopar% temp_func(df, i)
  tmp1 = tmp >= 0 
  df = df[,tmp1]
  
  
  #check survivied cells
  ncol(df)
  oridinaldf = colnames(df)
  colnames(df) = paste(name,'_', c(1:ncol(df)), sep='')
  df = as.data.frame(df)
  seurat_obj = CreateSeuratObject(raw.data = df, project = name)
  return(seurat_obj)
}


# 上記の処理を行ってseurat object を作る
seurat_obj_list = map2(file_list, names(file_list), do_process)

# seurat object のリストをマージする
SAM <- reduce(seurat_obj_list,MergeSeurat,do.normalize = FALSE)

# マージされたデータの内訳
table(SAM@meta.data$orig.ident)

# ラベリング　NAFLD (S1-S6) vs NASH(S7-S9) vs HCC(S10)
SAM@meta.data$label = fct_collapse(SAM@meta.data$orig.ident, NAFLD = c("S1","S2","S3", "S4","S5","S6"), NASH = c("S7","S8","S9"), HCC = c("S10"))


mito.genes = grep(pattern = "^MT.", x = rownames(x = SAM@data), value = TRUE)
percent.mito = Matrix::colSums(SAM@raw.data[mito.genes, ])/Matrix::colSums(SAM@raw.data)
SAM = AddMetaData(object = SAM, metadata = percent.mito, col.name = "percent.mito")


jpeg(filename = file.path("image", "geneplot.jpg"))
GenePlot(object = SAM, gene1 = "nGene", gene2 = "percent.mito")
dev.off()

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

PCElbowPlot(object = SAM, )

#clustering and detection of cellular subsets
#dims.use = 1: X, X could be read by Standard Deviation of PC
SAM = FindClusters(object = SAM, reduction.type = "pca", dims.use = 1:20,
                   resolution = 1.0, print.output = 0, save.SNN = TRUE, force.recalc = TRUE)
SAM = RunTSNE(object = SAM, dims.use = 1:20, do.fast = TRUE)

TSNEPlot(object = SAM, pt.size = 1.5)
TSNEPlot(object = SAM, pt.size = 1.5, group.by = "orig.ident")
TSNEPlot(object = SAM, pt.size = 1.5, group.by = "label")





DotPlot(SAM, genes.plot = "")


# NAFLD vs NASH vs HCC comparison -----------------------------------------


SAM = SetAllIdent(object = SAM, id = "label")

TSNEPlot(object = SAM, colors.use = c ("black","red"))
comparison_avg = FindMarkers(object = SAM, ident.1 = "new_S1", ident.2 = "new_S2")
comparison_avg = comparison_avg %>% rownames_to_column(var = "gene")
comparison_avg = comparison_avg %>% mutate(cluster = if_else(avg_logFC>0, "id1_posi", "id2_posi"))
comparison_avg = comparison_avg [!grepl("^MT-|^Rpl|^Rps|^Rna|^Mtrnr|^LOC", comparison_avg$gene),]
comparison_avg = comparison_avg %>% group_by(cluster) %>% do(head(., 20)) %>% arrange(cluster)
DoHeatmap(object = SAM, cells.use = WhichCells(object = SAM, ident = c("new_S1", "new_S2")), genes.use = comparison_avg$gene,
          slim.col.label = TRUE, remove.key = T)


SAM = SetAllIdent(object = SAM, id = "orig.ident")

# S1 ,S2, S3 vs S4 
vs_S4_list = map(c("S1","S2","S3"),
    ~{
    comparison_avg = FindMarkers(object = SAM, ident.1 = .x, ident.2 = "S4")
    comparison_avg = comparison_avg %>% rownames_to_column(var = "gene")
    comparison_avg = comparison_avg %>% mutate(cluster = if_else(avg_logFC>0, "id1_posi", "id2_posi"))
    comparison_avg = comparison_avg [!grepl("^MT-|^Rpl|^Rps|^Rna|^Mtrnr|^LOC", comparison_avg$gene),]
    comparison_avg = comparison_avg %>% group_by(cluster) %>% do(head(., 20)) %>% arrange(cluster)
    DoHeatmap(object = SAM, cells.use = WhichCells(object = SAM, ident = c(.x, "S4")), genes.use = comparison_avg$gene,
              slim.col.label = TRUE, remove.key = T)
    
    }
    )
.Last.value




