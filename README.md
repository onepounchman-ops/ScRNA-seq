E12.5 <- Read10X_h5(file.path(dataset_loc, ids[1], "GSM4504959_E12.5.h5"), use.names = T)
E15.5 <- Read10X_h5(file.path(dataset_loc, ids[1], "GSM4504960_E15.5.h5"), use.names = T)
E17.5 <- Read10X_h5(file.path(dataset_loc, ids[1], "GSM4504961_E17.5.h5"), use.names = T)
colnames(E12.5) <- paste(colnames(E12.5),"E12.5",sep = "_")
colnames(E15.5) <- paste(colnames(E15.5),"E15.5",sep = "_")
colnames(E17.5) <- paste(colnames(E17.5),"E17.5",sep = "_")
sam.name <- "multi"
if(!dir.exists(sam.name)){
  dir.create(sam.name)
}
experiment.aggregate <- CreateSeuratObject(
  experiment.data,
  project = "multi", 
  min.cells = 10,
  min.features = 200,
  names.field = 2,
  names.delim = "_")
save(experiment.aggregate,file=paste0("./",sam.name,"/",sam.name,"_raw_SeuratObject.RData"))
slotNames(experiment.aggregate)
experiment.aggregate@assays
dim(experiment.aggregate@meta.data)
View(experiment.aggregate@meta.data)
experiment.aggregate.matrix <- as.matrix(experiment.aggregate@assays$RNA@counts)
experiment.aggregate[["percent.mt"]] <- PercentageFeatureSet(experiment.aggregate, 
                                                             pattern = "^mt-")
pdf(paste0("./",sam.name,"/QC-VlnPlot.pdf"),width = 8,height = 4.5)
VlnPlot(experiment.aggregate, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
        ncol = 3)
dev.off()
gene.freq <- do.call("cbind", tapply(experiment.aggregate@meta.data$nFeature_RNA,experiment.aggregate@meta.data$orig.ident,quantile,probs=seq(0,1,0.05)))
rna.freq <- do.call("cbind", tapply(experiment.aggregate@meta.data$nCount_RNA,experiment.aggregate@meta.data$orig.ident,quantile,probs=seq(0,1,0.05)))
mt.freq <- do.call("cbind", tapply(experiment.aggregate@meta.data$percent.mt,experiment.aggregate@meta.data$orig.ident,quantile,probs=seq(0,1,0.05)))
freq.combine <- as.data.frame(cbind(gene.freq,rna.freq,mt.freq))
colnames(freq.combine) <- c(paste(colnames(gene.freq),"Gene",sep = "_"),
                            paste(colnames(rna.freq),"RNA",sep = "_"),
                            paste(colnames(mt.freq),"mt",sep = "_"))
write.table(freq.combine,file = paste0(sam.name,"/QC-gene_frequency.txt"),quote = F,sep = "\t")
rm(gene.freq,rna.freq,mt.freq)
View(freq.combine)
plot1 <- FeatureScatter(experiment.aggregate, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(experiment.aggregate, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
pdf(paste0("./",sam.name,"/QC-FeatureScatter.pdf"),width = 8,height = 4.5)
CombinePlots(plots = list(plot1, plot2),legend = "none")
dev.off()
rm(plot1,plot2)
cat("Before filter :",nrow(experiment.aggregate@meta.data),"cells\n")
experiment.aggregate <- subset(experiment.aggregate, 
                               subset = 
                                 nFeature_RNA > 500 & 
                                 nCount_RNA > 1000 & 
                                 nCount_RNA < 20000 &
                                 percent.mt < 5)
cat("After filter :",nrow(experiment.aggregate@meta.data),"cells\n")
experiment.aggregate <- NormalizeData(experiment.aggregate, 
                                      normalization.method = "LogNormalize",
                                      scale.factor = 10000)
experiment.aggregate <- FindVariableFeatures(experiment.aggregate, 
                                             selection.method = "vst",
                                             nfeatures = 1000)
top10 <- head(x = VariableFeatures(experiment.aggregate), 10)
plot1 <- VariableFeaturePlot(experiment.aggregate)
plot2 <- LabelPoints(plot = plot1, points = top10)
pdf(file = paste0(sam.name,"/Norm-feature_variable_plot.pdf"),width = 8,height = 5)
CombinePlots(plots = list(plot1, plot2),legend = "none")
dev.off()

experiment.aggregate <- ScaleData(
  object = experiment.aggregate,
  do.scale = FALSE,
  do.center = FALSE,
  vars.to.regress = c("orig.ident","percent.mt"))

experiment.aggregate <- RunPCA(object = experiment.aggregate, 
                               features = VariableFeatures(experiment.aggregate),
                               verbose = F,npcs = 50)

pdf(paste0("./",sam.name,"/PCA-VizDimLoadings.pdf"),width = 7,height = 5)
VizDimLoadings(experiment.aggregate, dims = 1:2, reduction = "pca")
dev.off()

pdf(paste0("./",sam.name,"/PCA-DimPlot.pdf"),width = 5,height = 4)
DimPlot(experiment.aggregate, reduction = "pca")
dev.off()

pdf(paste0("./",sam.name,"/PCA-DimHeatmap.pdf"))
DimHeatmap(experiment.aggregate, dims = 1:6, cells = 500, balanced = TRUE)
dev.off()

experiment.aggregate <- JackStraw(experiment.aggregate, num.replicate = 100,dims = 40)
experiment.aggregate <- ScoreJackStraw(experiment.aggregate, dims = 1:40)
pdf(paste0("./",sam.name,"/PCA-JackStrawPlot_40.pdf"),width = 6,height = 5)
JackStrawPlot(object = experiment.aggregate, dims = 1:40)
dev.off()

pdf(paste0("./",sam.name,"/PCA-ElbowPlot.pdf"),width = 8,height = 5)
ElbowPlot(experiment.aggregate,ndims = 40)
dev.off()

dim.use <- 1:14

experiment.aggregate <- FindNeighbors(experiment.aggregate, dims = dim.use)
experiment.aggregate <- FindClusters(experiment.aggregate, resolution = 0.5)

experiment.aggregate <- RunTSNE(experiment.aggregate, dims = dim.use, 
                                do.fast = TRUE)
pdf(paste0("./",sam.name,"/CellCluster-TSNEPlot_res0.5_",max(dim.use),"PC.pdf"),width = 8.5,height = 7)
DimPlot(object = experiment.aggregate, pt.size=0.5,label = T) 
dev.off()

pdf(paste0("./",sam.name,"/CellCluster-TSNEPlot_SamGroup_",max(dim.use),"PC.pdf"),width = 8,height = 7)
DimPlot(object = experiment.aggregate, 
        group.by="orig.ident", 
        pt.size=0.5,reduction = "tsne")
dev.off()

pdf(paste0("./",sam.name,"/CellCluster-TSNEPlot_SamGroup_slipt_",max(dim.use),"PC.pdf"),width = 10,height = 4)
DimPlot(object = experiment.aggregate, 
        split.by ="orig.ident", 
        pt.size=0.5,reduction = "tsne")
dev.off()

table(experiment.aggregate@meta.data$orig.ident)

all.markers <- FindAllMarkers(experiment.aggregate, only.pos = TRUE, 
                              min.pct = 0.3, logfc.threshold = 0.25)
write.table(all.markers,
            file=paste0("./",sam.name,"/",sam.name,"_total_marker_genes_tsne_",max(dim.use),"PC.txt"),
            sep="\t",quote = F,row.names = F)

Selective.marker.genes<-c("Epcam","Ager","Hopx","Sftpb","Sftpc","S100a8",
                          "S100a9","Hba-a2","Hbb-bt","Tgfbi","Actg2",
                          "Acta2","Myh11","Tagln","Col1a1","Col13a1",
                          "Col14a1","Scgb3a2","Scgb1a1","Cyp2f2",
                          "Upk3b","Wt1","Ube2c","Hmmr","Cldn5",
                          "Emcn","Pecam1","Cd68","C1qa","Pdgfrb",
                          "Higd1b","Gucy1a3","Cox4i2","Postn",
                          "Pdzd2","Cspg4","Ptn","Lum","Prrx1")
Epithelium.marker.genes<-c("Epcam","Ager","Hopx","Sftpb","Sftpc")  
  #VlnPlot
  pvn <- VlnPlot(experiment.aggregate, features = Selective.marker.genes,ncol = 3)
  pdf(paste0("./",sam.name,"/MarkerGene-VlnPlot_cluster","_tsne_",max(dim.use),"PC.pdf"),width = 7,height = 6)
  print(pvn)
  dev.off()
  
  #Feather plot 
  pvn <- FeaturePlot(experiment.merged,features=Epithelium.marker.genes,ncol = 3,cols = brewer.pal(8,"YlOrBr"))
  pdf(paste0("./",sam.name,"/MarkerGene-FeaturePlot_cluster","_tsne_",max(dim.use),"PC.pdf"),width = 10,height = 6)
  print(pvn)
  dev.off()
  
  #RidgePlot
  pvn<-RidgePlot(experiment.aggregate, features = cl4.genes, ncol = 2)
  pdf(paste0("./",sam.name,"/MarkerGene-RidgePlot_cluster",cluster_id,"_tsne_",max(dim.use),"PC.pdf"),width = 7,height = 6)
  print(pvn)
  dev.off()
}
