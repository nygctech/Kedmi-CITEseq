library(Seurat)
library(dplyr)
library(ggplot2)
library(RColorBrewer)
grd.col <- brewer.pal(11,"RdYlGn")[11:1]

#set the relevant folder

RNA_matrix = Read10X(data.dir = "./cDNA1_counts_matrix/") 
ADT_dir='./ADT1_counts_matrix'
ADT_umi_matrix <- readMM(paste0(ADT_dir, "/genecounts.mtx"))
ADT_features <- readLines(paste0(ADT_dir, '/genecounts.genes.txt'))
ADT_bcs <-  readLines(paste0(ADT_dir, '/genecounts.barcodes.txt'))
rownames(ADT_umi_matrix) <- paste0(ADT_bcs, "-1")
colnames(ADT_umi_matrix) <- ADT_features
ADT_umi_matrix <- t(ADT_umi_matrix)

HTO_dir='./HTO1_counts_matrix_2'
HTO_umi_matrix <- readMM(paste0(HTO_dir, '/genecounts.mtx'))
HTO_features <- readLines(paste0(HTO_dir, '/genecounts.genes.txt'))
HTO_bcs <-  readLines(paste0(HTO_dir, '/genecounts.barcodes.txt'))
rownames(HTO_umi_matrix) <-  paste0(HTO_bcs, "-1") 
colnames(HTO_umi_matrix) <- HTO_features
HTO_umi_matrix <- t(HTO_umi_matrix)


cells.use = Reduce(intersect, list(colnames(RNA_matrix) , colnames(ADT_umi_matrix) , colnames(HTO_umi_matrix)))
obj.1 <- CreateSeuratObject(counts = RNA_matrix[, cells.use])
obj.1[["ADT"]] <- CreateAssayObject(counts = ADT_umi_matrix[, cells.use])
obj.1[["HTO"]] <- CreateAssayObject(counts = HTO_umi_matrix[, cells.use])

obj.1[['percent.mito']] <- PercentageFeatureSet(obj.1, pattern = "^mt-")
obj.1$Lane <-"Lane1"

VlnPlot(obj.1, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), pt.size = 0)
obj.1 <- subset(x = obj.1, subset = nFeature_RNA > 700 & nFeature_RNA < 5000 & percent.mito < 6)
 

# remove adt and hto too high cells
VlnPlot(obj.1, features = c("nCount_ADT", "nCount_HTO"), pt.size = 0, log = T)
VlnPlot(obj.1, features = c("nCount_ADT", "nCount_HTO"), pt.size = 0.1 )

obj.1 <- subset(obj.1, subset = nCount_ADT < 5000 & nCount_HTO < 4000)


RNA_matrix = Read10X(data.dir = "./cDNA2_counts_matrix/") 
ADT_dir='./ADT2_counts_matrix'
ADT_umi_matrix <- readMM(paste0(ADT_dir, "/genecounts.mtx"))
ADT_features <- readLines(paste0(ADT_dir, '/genecounts.genes.txt'))
ADT_bcs <-  readLines(paste0(ADT_dir, '/genecounts.barcodes.txt'))
rownames(ADT_umi_matrix) <- paste0(ADT_bcs, "-1")
colnames(ADT_umi_matrix) <- ADT_features
ADT_umi_matrix <- t(ADT_umi_matrix)

HTO_dir='./HTO2_counts_matrix_2'
HTO_umi_matrix <- readMM(paste0(HTO_dir, '/genecounts.mtx'))
HTO_features <- readLines(paste0(HTO_dir, '/genecounts.genes.txt'))
HTO_bcs <-  readLines(paste0(HTO_dir, '/genecounts.barcodes.txt'))
rownames(HTO_umi_matrix) <-  paste0(HTO_bcs, "-1") 
colnames(HTO_umi_matrix) <- HTO_features
HTO_umi_matrix <- t(HTO_umi_matrix)


cells.use = Reduce(intersect, list(colnames(RNA_matrix) , colnames(ADT_umi_matrix) , colnames(HTO_umi_matrix)))
 
obj.2 <- CreateSeuratObject(counts = RNA_matrix[, cells.use])
obj.2[["ADT"]] <- CreateAssayObject(counts = ADT_umi_matrix[, cells.use])
obj.2[["HTO"]] <- CreateAssayObject(counts = HTO_umi_matrix[, cells.use])

obj.2[['percent.mito']] <- PercentageFeatureSet(obj.2, pattern = "^mt-")
obj.2$Lane <-"Lane2"

obj.2 <- subset(x = obj.2, subset = nFeature_RNA > 700 & nFeature_RNA < 5000 & percent.mito < 6)
 
VlnPlot(obj.2, features = c("nCount_ADT", "nCount_HTO"), pt.size = 0, log = T)
VlnPlot(obj.2, features = c("nCount_ADT", "nCount_HTO"), pt.size = 0.1 )
obj.2 <- subset(obj.2, subset = nCount_ADT < 5000 & nCount_HTO < 4000)



obj.combined <- merge(x = obj.1, y = obj.2, add.cell.ids = c("Lane1", "Lane2"))
 

#ADT and HTO normalization and scaling

obj.combined <- NormalizeData(obj.combined, assay = "HTO", normalization.method = "CLR")
obj.combined <- ScaleData(obj.combined, assay = "HTO")
obj.combined <- HTODemux(obj.combined, positive.quantile = 0.99)

VlnPlot(obj.combined, features = "M1" )
print (table(obj.combined$HTO_classification.global))
obj.combined <- subset(obj.combined, subset = HTO_classification.global != "Negative" )

obj.single <- subset(obj.combined, subset = HTO_classification.global == "Singlet" )

DefaultAssay(obj.single) <- "RNA"
obj.single <- SCTransform(obj.single)%>%RunPCA()%>%RunUMAP(dims = 1:20)
DimPlot(obj.single, reduction = "umap")


## when normalize ADT by CLR, the margin need to be set 2
DefaultAssay(obj.single) <- "ADT"
obj.single <- NormalizeData(obj.single, normalization.method = "CLR", margin = 2)
VariableFeatures(obj.single) <- rownames(obj.single[["ADT"]])
# when you have a big protein panel, usually set do.scale to FALSE 
obj.single <- ScaleData(obj.single, do.scale = FALSE)


obj.single <- RunPCA(obj.single, reduction.name = "pca.adt", reduction.key = "PCadt_")
obj.single <- RunUMAP(obj.single, reduction = "pca.adt", dims = 1:20, reduction.name = "umap.adt", reduction.key = "Uadt_")
DimPlot(obj.single, reduction = "umap.adt")


## WNN analysis
obj.single <- FindMultiModalNeighbors(object = obj.single, reduction.list = list("pca", "pca.adt"), dims.list = list(1:20, 1:20))
obj.single <- RunUMAP(obj.single, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "Uwnn_" )
obj.single <- FindClusters(obj.single, graph.name = "wsnn", algorithm = 3)
DimPlot(obj.single, reduction = "wnn.umap", label = T)


## looking at cell fractions, based on genes express mainly in cytoplasm or nucleus
mygenescyto <- c("Cox8a","Fth1","Ftl1","Ndufa4","Rpl37a","Rpl41","Rps29","Serf2","Ubb","Cox6c","Cox7c","Eif1","Itm2b","Prdx1","Rpl35","Rpl37","Rps19","Rps27","Atp5l","Cox4i1","Cox7a2","Ndufa13","Uqcr10","Uqcr11")
obj.single <- AddModuleScore(object = obj.single, features = list(mygenescyto), ctrl = 5, assay = "SCT", name = "cytoplam")

mygenesnuc <- c("Neat1","Zbtb20","Malat1","Pde4d","Rbm39","Esrrg","Ghr","Rsrp1","Pkhd1","Ank3","Magi2","Mecom","Spag5","Arl15","Dach1","Gm21092","Kap","Lrba","Magi1","Plcb1","Ptprd","Rbms1","Slc34a1","Slit2")
obj.single <- AddModuleScore(object = obj.single, features = list(mygenesnuc), ctrl = 5, assay = "SCT", name = "nucleus")


DefaultAssay(obj.single) <- "SCT"
grd.col <- brewer.pal(11,"RdYlGn")[11:1]
FeaturePlot(obj.single, features = c("cytoplam1", "nucleus1"), cols = grd.col, pt.size = 0.03, reduction = "wnn.umap")
 
Idents(obj.single) <- obj.single$wsnn_res.0.8
obj.single <- subset(obj.single, invert = TRUE, idents = c(19,20,12,16,17)) # remove cell fractions and not relevant subsets (dont include MHCII)

Idents(obj.single) <- obj.single$HTO_maxID
obj.single <- RenameIdents(object = obj.single,'M1' = 'zDC GFP', 'M2' = 'zDC GFP', 'M3' = 'zDC GFP;RORgt tdt', 'M4' = 'zDC GFP;RORgt tdt', 'M5' = 'RORgt tdt','M6' = 'RORgt tdt','M7' = 'CD11c tdt','M8' = 'CD11c tdt')
obj.single$genotype <- Idents(obj.single)

Idents(obj.single) <- obj.single$wsnn_res.0.8

DimPlot(obj.single, reduction = "wnn.umap", label = T, repel = T)#, split.by = 'genotype') 
DimPlot(obj.single, reduction = "wnn.umap", label = T, repel = T, split.by = "genotype" ) 

obj.single <- RenameIdents(object = obj.single,'0' = 'pDC','1' = 'migratory cDC2', '9' = 'migratory cDC2','2' = 'migratory cDC1', '8' = 'ILC3', '13' = 'ILC1','14' = 'Aire+ ILC-like', '3' = 'CD11b+ DC','4' = 'CD11b+ DC','5' = 'CD11b+ DC','10' = 'proliferating CD11b+ DC','6' = 'cDC1','7' = 'cDC1', '18' = 'monocytes', '11' = 'B cells', '15' = 'plasmablasts')
#VlnPlot(CleanCD11c_tdt,  features = c("IA/IE"), pt.size = 0.01)#, max.cutoff = "q95", min.cutoff = "q05", cols = grd.col)#, split.by = "HTO_maxID")
FeaturePlot(CleanCD11c_tdt, reduction = "wnn.umap",  features = c("Igha"), pt.size = 0.2 , cols = grd.col, max.cutoff = "q95", min.cutoff = "q05")#, split.by = "genotype")#,)
grd.col <- brewer.pal(11,"RdYlBu")[11:1]

# change cell type color
obj.single$celltype <- Idents(obj.single)

set.seed(123)
random.order <- sample(unique(obj.single$celltype ))
obj.single$celltype  <- factor(obj.single$celltype , levels = random.order)
DimPlot(obj.single, group.by = "celltype", reduction = "wnn.umap") + xlab("UMAP_1") + ylab("UMAP_2")
DimPlot(obj.single, reduction = "wnn.umap") + xlab("UMAP_1") + ylab("UMAP_2")
VlnPlot(obj.single,  features = c("Itgav"), pt.size = 0.01, grouped.by = 'celltype')



Idents(obj.single) <- obj.single$genotype

CD11c_tdt <- subset(x = obj.single, idents = c("CD11c tdt"))
Idents(CD11c_tdt) <- CD11c_tdt$celltype
DimPlot(CD11c_tdt, reduction = "wnn.umap", label = F, repel = T, pt.size = 0.05) + xlab("UMAP_1") + ylab("UMAP_2")
VlnPlot(CD11c_tdt, features = "adt_IA/IE", pt.size = 0.01) + NoLegend()
FeaturePlot(CD11c_tdt, reduction = "wnn.umap", features = "adt_IA/IE", pt.size = 0.01) + NoLegend()


saveRDS(CD11c_tdt, "./CD11c_tdt-11721.rds")



###########flowjo
#get UMAP coordinates


umapcoord <- as.data.frame(CD11c_tdt@reductions$wnn.umap@cell.embeddings)
#get ADT data
adt.dataflo <- GetAssayData(object = CD11c_tdt, assay = "ADT", slot = "scale.data")
adt.distflo2 <- as.data.frame(t(adt.dataflo)) #flip row and column to merge
#merge UMAP and ADT
merge1 <- merge(adt.distflo2, umapcoord, by = "row.names")
rownames(merge1) <- merge1$Row.names

#add HTO info
#htodata <- GetAssayData(object = citeseq, assay = "HTO", slot = "scale.data")
#htodata2 <- as.data.frame(t(htodata))

#merge again
#merge3 <- merge(merge2, htodata2, by = "row.names")
#rownames(merge3) <- merge3$Row.names
#add hto maxid
#htomaxid <- as.data.frame(citeseq$HTO_maxID)
#merge4 <- merge(merge3, htomaxid, by = "row.names")
#rownames(merge4) <- merge4$Row.names

#add genotype
genotype <- as.data.frame(CD11c_tdt$genotype)
merge2 <- merge(merge1, genotype, by = "row.names")
rownames(merge2) <- merge2$Row.names

#add cluster identity
clustercoord <- as.data.frame(CD11c_tdt$wsnn_res.0.8)
merge3 <- merge(merge2, clustercoord, by = "row.names")
rownames(merge3) <- merge3$Row.names

#add ADT umap identity
DefaultAssay(CD11c_tdt) <- 'ADT'
adtumapcoord <- as.data.frame(CD11c_tdt@reductions$umap@cell.embeddings)
merge4 <- merge(merge3, adtumapcoord, by = "row.names")

head(merge4)
#delete first rowname columns
merge5 <- merge4[, -c(1:4)]
#write to csv
write.table(merge5, file="CD11c_tdt.csv", quote = F, sep = ',', row.names = F)

#fo zbtbgfp RORgt Fatemap
RORgt_fatemap <- subset(x = obj.single, idents = c("CD11c tdt"), invert = TRUE)
Idents(RORgt_fatemap) <- RORgt_fatemap$celltype
DimPlot(RORgt_fatemap, reduction = "wnn.umap", label = F, repel = T, pt.size = 0.05, split.by = 'genotype') + xlab("UMAP_1") + ylab("UMAP_2") + NoLegend()
VlnPlot(CD11c_tdt, features = "sct_Itgav", pt.size = 0.01) + NoLegend()
FeaturePlot(obj.single, reduction = "wnn.umap", features = "adt_CD51", pt.size = 0.05, cols = brewer.pal(11,"RdYlBu")[11:1], min.cutoff = 'q03', max.cutoff = 'q99') 

saveRDS(RORgt_fatemap, "./RORgt_fatemap-11721.rds")

