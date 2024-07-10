#load required libraries
library(Seurat)
library(cowplot)
library(dplyr)
library(scran)
library(scCATCH)
library(SingleR)
library(ImmClassifier)
library(ggplot2)
library(viridis)

#Analyze all patients except patient 8 as shown below using patient 5 as example:
#read patient 5 pre data file from its directory
patient5pre <- Read10X(data.dir = "C:/Users/sszymura/Desktop/NewData/Patient5/Pre/filtered_feature_bc_matrix")
#read patient 5 post data file from its directory
patient5post <- Read10X(data.dir = "C:/Users/sszymura/Desktop/NewData/Patient5/Post/filtered_feature_bc_matrix")
#create Seurat object for pre file
patient5preSeurat <- CreateSeuratObject(counts = patient5pre, project = "pre")
#create Seurat object for post file
patient5postSeurat <- CreateSeuratObject(counts = patient5post, project = "post")
#filter based on number on mRNA molecules
patient5preSeurat <- subset(patient5preSeurat, subset = nFeature_RNA > 200 & nFeature_RNA < 2500)
#filter based on number of mRNA molecules
patient5postSeurat <- subset(patient5postSeurat, subset = nFeature_RNA > 200 & nFeature_RNA < 2500)
#determine percentage of mitochondrial RNA
patient5preSeurat[["percent.mt"]] <- PercentageFeatureSet(patient5preSeurat, pattern = "^MT-")
#determine percentage of mitochondrial DNA
patient5postSeurat[["percent.mt"]] <- PercentageFeatureSet(patient5postSeurat, pattern = "^MT-")
#check visually
VlnPlot(patient5preSeurat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
#check visually
VlnPlot(patient5postSeurat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
#sort out based on mitochondrial RNA
patient5preSeurat <- subset(patient5preSeurat, subset = percent.mt < 20)
#sort out based on mitochondrial RNA)
patient5postSeurat <- subset(patient5postSeurat, subset = percent.mt <20)
#add metadata to Seurat object pre: these steps will label this Seurat object as pre treatment
patient5preSeurat <- AddMetaData(object = patient5preSeurat, metadata = rep("PrePatient5", length(Cells(patient5preSeurat))), col.name = "Treatment")
#add metadata to Seurat object post: these steps will label this Seurat object as post treatment  
patient5postSeurat <- AddMetaData(object = patient5postSeurat, metadata = rep("PostPatient5", length(Cells(patient5postSeurat))), col.name = "Treatment")
#normalize and find variable features
patient5preSeurat <- NormalizeData(patient5preSeurat)
patient5postSeurat <- NormalizeData(patient5postSeurat)
#optional:find doublets
#pre5 <- doubletCells(GetAssayData(patient5preSeurat))
#post5 <- doubletCells(GetAssayData(patient5postSeurat))
#Cells with low scores are singlets
#pre5filtered <- pre5 < 100
#post5filtered <- post5 < 100
#cells5pre <- Cells(patient5preSeurat)
#cells5post <- Cells(patient5postSeurat)
#cells5prefiltered <- cells5pre[pre5filtered]
#cells5postfiltered <- cells5post[post5filtered]
#patient5prefiltered <- subset(patient5preSeurat, cells = cells5prefiltered)
#patient5postfiltered <- subset(patient5postSeurat, cells = cells5postfiltered)
#save Seurat objects
patient5pre <- FindVariableFeatures(patient5preSeurat, selection.method = "vst", nfeatures = 2000)
patient5post <- FindVariableFeatures(patient5postSeurat, selection.method = "vst", nfeatures = 2000)
save(patient5pre, file = "C:/Users/sszymura/Desktop/patient5pre.Robj")
save(patient5post, file = "C:/Users/sszymura/Desktop/patient5post.Robj")

#read patient 8 pre data file from its directory
patient8 <- Read10X(data.dir = "C:/Users/sszymura/Desktop/AllPatientData/Patient8/filtered_feature_bc_matrix")
patient8 <- patient8[[1]]
#create Seurat object for pre file
patient8 <- CreateSeuratObject(counts = patient8, project = "Patient8")
#filter based on number on mRNA molecules
patient8cells <- Cells(patient8)
write.csv(patient8cells, file = "C:/Users/sszymura/Desktop/patient8cells.csv")
patient8 <- RenameCells(patient8, new.names = cells$x)
patient8 <- subset(patient8, subset = nFeature_RNA > 200 & nFeature_RNA < 2500)
#determine percentage of mitochondrial RNA
patient8[["percent.mt"]] <- PercentageFeatureSet(patient8, pattern = "^MT-")
#check visually
VlnPlot(patient8, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
#sort out based on mitochondrial RNA
patient8 <- subset(patient8, subset = percent.mt < 20)
#split object into pre and post:
cells8pre <- read.csv(file = "C:/Users/sszymura/Desktop/NewData/Patient8/treatment_pre_patient8.csv")
cells8post <- read.csv(file = "C:/Users/sszymura/Desktop/NewData/Patient8/treatment_post_patient8.csv")
patient8pre <- subset(patient8, cells = as.character(cells8pre$Barcode))
patient8post <- subset(patient8, cells = as.character(cells8post$Barcode))
#add metadata to Seurat object pre: these steps will label this Seurat object as pre treatment
patient8pre <- AddMetaData(object = patient8pre, metadata = rep("PrePatient8", length(Cells(patient8pre))), col.name = "Treatment")
#add metadata to Seurat object post: these steps will label this Seurat object as post treatment  
patient8post <- AddMetaData(object = patient8post, metadata = rep("PostPatient8", length(Cells(patient8post))), col.name = "Treatment")
#normalize and find variable features
patient8pre <- NormalizeData(patient8pre)
patient8post <- NormalizeData(patient8post)
#optional:find doublets
#pre8 <- doubletCells(GetAssayData(patient8pre))
#post8 <- doubletCells(GetAssayData(patient8post))
#Cells with low scores are singlets
#pre8filtered <- pre8 < 100
#post8filtered <- post8 < 100
#cells8pre <- Cells(patient8pre)
#cells8post <- Cells(patient8post)
#cells8prefiltered <- cells8pre[pre8filtered]
#cells8postfiltered <- cells8post[post8filtered]
#patient8prefiltered <- subset(patient8pre, cells = cells8prefiltered)
#patient8postfiltered <- subset(patient8post, cells = cells8postfiltered)
#save Seurat objects
patient8pre <- FindVariableFeatures(patient8prefiltered, selection.method = "vst", nfeatures = 2000)
patient8post <- FindVariableFeatures(patient8postfiltered, selection.method = "vst", nfeatures = 2000)
save(patient8pre, file = "C:/Users/sszymura/Desktop/patient8pre.Robj")
save(patient8post, file = "C:/Users/sszymura/Desktop/patient8post.Robj")

#filter BCR clonotypes from raw data
setwd("C:/Users/sszymura/Desktop/ProcessedData/RawVDJ/BCR/RawData")
files <- list.files()
for (file in files) {
  x <- read.csv(file, stringsAsFactors = FALSE);
  x <- split(x, x$chain);
  x <- full_join(rbind(x$IGK, x$IGL), x$IGH, by = "barcode");
  y <- filter(x, x$productive.x == TRUE);
  y <- filter(y, y$productive.y == TRUE);
  z <- subset(x, !(x$barcode %in% y$barcode))
  x <- rbind(y, z)
  x <- filter(x, !(x$productive.x == FALSE & x$productive.y == FALSE))
  x <- filter(x, !(x$productive.x == FALSE & x$productive.y == 'None'))
  x <- filter(x, !(x$productive.x == 'None' & x$productive.y == 'None'))
  x <- filter(x, !(x$productive.x == 'None' & x$productive.y == FALSE))
  x <- filter(x, !(duplicated(x$barcode)))
  write.csv(x, file = paste0("C:/Users/sszymura/Desktop/ProcessedData/RawVDJ/BCR/", file))
}

#filter TCR clonotypes from raw data
setwd("C:/Users/sszymura/Desktop/ProcessedData/RawVDJ/TCR/RawData")
files <- list.files()
for (file in files) {
  x <- read.csv(file, stringsAsFactors = FALSE);
  x <- split(x, x$chain);
  x <- full_join(x$TRA, x$TRB, by = "barcode");
  y <- filter(x, x$productive.x == TRUE);
  y <- filter(y, y$productive.y == TRUE);
  z <- subset(x, !(x$barcode %in% y$barcode))
  x <- rbind(y, z)
  x <- filter(x, !(x$productive.x == FALSE & x$productive.y == FALSE))
  x <- filter(x, !(x$productive.x == FALSE & x$productive.y == 'None'))
  x <- filter(x, !(x$productive.x == 'None' & x$productive.y == 'None'))
  x <- filter(x, !(x$productive.x == 'None' & x$productive.y == FALSE))
  x <- filter(x, !(duplicated(x$barcode)))
  write.csv(x, file = paste0("C:/Users/sszymura/Desktop/ProcessedData/RawVDJ/TCR/", file))
}

#split patient8 sample
setwd("C:/Users/sszymura/Desktop/ProcessedData")
patient8cells <- read.csv("C:/Users/sszymura/Desktop/ProcessedData/patient8totalcellsnewanalysis.csv")
tcr <- read.csv(file = "C:/Users/sszymura/Desktop/ProcessedData/RawVDJ/TCR/8_tcr.csv")
bcr <- read.csv(file = "C:/Users/sszymura/Desktop/ProcessedData/RawVDJ/BCR/8_bcr.csv")
cells <- split(patient8cells, patient8cells$Treatment)
cells8post <- as.data.frame(cells$Post)
cells8pre <- as.data.frame(cells$Pre)
names(cells8post) <- c("barcode")
names(cells8pre) <- c("barcode")
tcr8pre <- semi_join(tcr, cells8pre, by = "barcode")
tcr8post <- semi_join(tcr, cells8post, by = "barcode")
bcr8pre <- semi_join(bcr, cells8pre, by = "barcode")
bcr8post <- semi_join(bcr, cells8post, by = "barcode")
write.csv(bcr8pre, file = "C:/Users/sszymura/Desktop/ProcessedData/RawVDJ/BCR/8_pre_bcr.csv")
write.csv(bcr8post, file = "C:/Users/sszymura/Desktop/ProcessedData/RawVDJ/BCR/8_post_bcr.csv")
write.csv(tcr8pre, file = "C:/Users/sszymura/Desktop/ProcessedData/RawVDJ/TCR/8_pre_tcr.csv")
write.csv(tcr8post, file = "C:/Users/sszymura/Desktop/ProcessedData/RawVDJ/TCR/8_post_tcr.csv")

############# Process all TCR and BCR pre and post files for each patient in similar manner as described below #############################################
setwd("C:/Users/sszymura/Desktop/NewTCRanalysis")
post <- read.csv("1post_tcr.csv", stringsAsFactors = FALSE)
pre <- read.csv("1pre_tcr.csv", stringsAsFactors = FALSE)
library(dplyr)

for (i in 1:length(pre$barcode)) {
  if (pre$productive.x[i] == "FALSE") {
    pre$cdr3.x[i] <- "None"
  }
  if (pre$productive.y[i] == "FALSE") {
    pre$cdr3.y[i] <- "None"
  }
}

for (i in 1:length(post$barcode)) {
  if (post$productive.x[i] == "FALSE") {
    post$cdr3.x[i] <- "None"
  }
  if (post$productive.y[i] == "FALSE") {
    post$cdr3.y[i] <- "None"
  }
}

patient1 <- rbind(pre, post)
patient1 <- as.data.frame(cbind(patient1$cdr3.x, patient1$cdr3.y))
patient1 <- cbind(patient1, rep(0, length(patient1$V1)))
names(patient1)[3] <- "clonotype"
patient1clonotypes <- distinct(patient1)

for (i in 1:length(patient1clonotypes$clonotype)) {patient1clonotypes$clonotype[i] <- paste0("clonotype", i)}

post <- cbind(post, rep(0, length(post$barcode)))
pre <- cbind(pre, rep(0, length(pre$barcode)))
names(post)[36] <- "NewClonotype"
names(pre)[36] <- "NewClonotype"

for (i in 1:length(post$barcode)) {
  for (j in 1:length(patient1clonotypes$clonotype))
    if ((post$cdr3.x[i] == patient1clonotypes$V1[j]) & (post$cdr3.y[i] == patient1clonotypes$V2[j])){
      post$NewClonotype[i] <- patient1clonotypes$clonotype[j]
    }
}

for (i in 1:length(pre$barcode)) {
  for (j in 1:length(patient1clonotypes$clonotype))
    if ((pre$cdr3.x[i] == patient1clonotypes$V1[j]) & (pre$cdr3.y[i] == patient1clonotypes$V2[j])){
      pre$NewClonotype[i] <- patient1clonotypes$clonotype[j]
    }
}

write.csv(pre, file = "C:/Users/sszymura/Desktop/NewestTCRanalysis/patient1pre.csv")
write.csv(post, file = "C:/Users/sszymura/Desktop/NewestTCRanalysis/patient1post.csv")


#Add BCR and TCR data to Seurat object as metadata. As an example:
load("patient1post.Robj")
load("patient1pre.Robj")
load("patient2post.Robj")  
load("patient2pre.Robj")
load("patient3post.Robj")
load("patient3pre.Robj")
load("patient4post.Robj")
load("patient4pre.Robj")
load("patient5post.Robj")
load("patient5pre.Robj")
load("patient6pre.Robj")
load("patient7post.Robj")
load("patient7pre.Robj")
load("patient8post.Robj")
load("patient8pre.Robj")
load("patient9post.Robj")
load("patient9pre.Robj")

meta = as.data.frame(cbind(Cells(patient1post)))
names(meta) = "barcode"
meta = left_join(meta, post, by = "barcode")
patient1post = AddMetaData(patient1post, metadata = meta$NewClonotype, rownames = "TCR") #can also add CDR, VDJ information

#Creating files for analysis of changes in clonotype abundancies pre and post vaccination. Identical for BCR and TCR analysis
Post1bcr <- read.csv("post1.csv", stringsAsFactors = FALSE)
Pre1bcr <- read.csv("pre1.csv", stringsAsFactors = FALSE)
Post2bcr <- read.csv("post2.csv", stringsAsFactors = FALSE)
Pre2bcr <- read.csv("pre2.csv", stringsAsFactors = FALSE)
Post3bcr <- read.csv("post3.csv", stringsAsFactors = FALSE) 
Pre3bcr <- read.csv("pre3.csv", stringsAsFactors = FALSE)
Post4bcr <- read.csv("post4.csv", stringsAsFactors = FALSE) 
Pre4bcr <- read.csv("pre4.csv", stringsAsFactors = FALSE)
Post5bcr <- read.csv("post5.csv", stringsAsFactors = FALSE)
Pre5bcr <- read.csv("pre5.csv", stringsAsFactors = FALSE)
Pre6bcr <- read.csv("pre6.csv", stringsAsFactors = FALSE)
Post7bcr <- read.csv("post7.csv", stringsAsFactors = FALSE)
Pre7bcr <- read.csv("pre7.csv", stringsAsFactors = FALSE)
Post8bcr <- read.csv("post8.csv", stringsAsFactors = FALSE)
Pre8bcr <- read.csv("pre8.csv", stringsAsFactors = FALSE)
Post9bcr <- read.csv("post9.csv", stringsAsFactors = FALSE)
Pre9bcr <- read.csv("pre9.csv", stringsAsFactors = FALSE)

Patient1 <- full_join(Post1bcr, Pre1bcr, by = "Var1")
Patient2 <- full_join(Post2bcr, Pre2bcr, by = "Var1")
Patient3 <- full_join(Post3bcr, Pre3bcr, by = "Var1")
Patient4 <- full_join(Post4bcr, Pre4bcr, by = "Var1")
Patient5 <- full_join(Post5bcr, Pre5bcr, by = "Var1")
Patient7 <- full_join(Post7bcr, Pre7bcr, by = "Var1")
Patient8 <- full_join(Post8bcr, Pre8bcr, by = "Var1")
Patient9 <- full_join(Post9bcr, Pre9bcr, by = "Var1")

write.csv(Patient1, file = "patient1.csv")
write.csv(Patient2, file = "patient2.csv")
write.csv(Patient3, file = "patient3.csv")
write.csv(Patient4, file = "patient4.csv")
write.csv(Patient5, file = "patient5.csv")
write.csv(Patient7, file = "patient7.csv")
write.csv(Patient8, file = "patient8.csv")
write.csv(Patient9, file = "patient9.csv")

#merging patient data files
setwd("C:/Users/sszymura/Desktop/Currentanalysis/ProcessedData")
load("patient1post.Robj")
load("patient1pre.Robj")
load("patient2post.Robj")  
load("patient2pre.Robj")
load("patient3post.Robj")
load("patient3pre.Robj")
load("patient4post.Robj")
load("patient4pre.Robj")
load("patient5post.Robj")
load("patient5pre.Robj")
load("patient6pre.Robj")
load("patient7post.Robj")
load("patient7pre.Robj")
load("patient8post.Robj")
load("patient8pre.Robj")
load("patient9post.Robj")
load("patient9pre.Robj")

immune.anchors <- FindIntegrationAnchors(object.list = list(patient1post, patient1pre, patient2post, patient2pre, patient3post, patient3pre, patient4post, patient4pre, patient5post, patient5pre, patient6pre, patient7post, patient7pre, patient8post, patient8pre, patient9post, patient9pre), dims = 1:30)
Patientscombined <- IntegrateData(anchorset = immune.anchors, dims = 1:30)
save(Patientscombined, file = "C:/Users/sszymura/Desktop/Patientscombined.Robj")
DefaultAssay(PatientscoSmbined) <- "integrated"
Patientscombined <- ScaleData(Patientscombined, verbose = FALSE, features = rownames(data))
Patientscombined <- RunPCA(Patientscombined, npcs = 30, verbose = FALSE)
ElbowPlot(Patientscombined)
Patientscombined <- FindNeighbors(Patientscombined, reduction = "pca", dims = 1:20)
Patientscombined <- FindClusters(Patientscombined, resolution = 0.5)
Patientscombined <- RunUMAP(Patientscombined, reduction = "pca", dims = 1:20)
DimPlot(Patientscombined, reduction = "umap", label = TRUE)
DimPlot(Patientscombined, reduction = "umap", split.by = "Treatment")
#save merged seurat object
save(Patientscombined, file = "C:/Users/sszymura/Desktop/Patientscombined2.Robj")

#Analysis using singleR package
ExpressionIntegrated <- GetAssayData(object = Patientscombined, assay = "integrated")
mid <- MonacoImmuneData()
pred55 <- SingleR(test = matrix, ref = mid, labels = mid$label.fine)
Patientscombined <- AddMetaData(Patientscombined, monaco2$labels, col.name = "Monaco")

#######Analyze bone marrow samples from human cell atlas################

atlas <- read.table("C:/Users/sszymura/Desktop/1MimmuneCellsMeta.tsv", header = TRUE, sep = "\t")
males50 <- filter(atlas, atlas$specimen_from_organism.organ_parts.ontology_label == "bone marrow" & atlas$donor_organism.organism_age == 50)
females52 <- filter(atlas, atlas$specimen_from_organism.organ_parts.ontology_label == "bone marrow" & atlas$donor_organism.organism_age == 52)
males50 <- as.character(males50$file_uuid)
females52 <- as.character(females52$file_uuid)

setwd("C:/Users/sszymura/Desktop/EntireAtlas")
files <- list.files()
for (file in files) {
  x <- readRDS(file);
  if (as.character(x$file_uuid)[1] %in% females52) {saveRDS(x, file = paste0("C:/Users/sszymura/Desktop/", file))
  }}

for (file in files) {
  x <- readRDS(file);
  if (as.character(x$file_uuid)[1] %in% males50) {saveRDS(x, file = paste0("C:/Users/sszymura/Desktop/", file))
  }}

#QC filter in similar manner as patient data
bmlist <- c(bm1, bm107, bm109, bm117, bm119, bm120, bm15, bm3, bm40, bm49, bm63, bm84, bm92, bm94, bm96)
for (i in 1:length(bmlist)) {
  x <- bmlist[[i]]
  x <- subset(x, subset = nFeature_RNA > 200 & nFeature_RNA < 2500);
  x[["percent.mt"]] <- PercentageFeatureSet(x, pattern = "^MT-");
  x <- subset(x, subset = percent.mt < 20);
  x <- AddMetaData(object = x, metadata = rep(paste0("bm", 10), length(Cells(x))), col.name = "Treatment")
  x <- NormalizeData(x)
  #z <- doubletCells(GetAssayData(x))
  #y <- z < 100
  #cellsx <- Cells(x)
  #cellsxfiltered <- cellsx[y]
  #x <- subset(x, cells = cellsxfiltered)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
  saveRDS(x, file = paste0("C:/Users/sszymura/Desktop/bm_new_", 10, ".Robj"))
}

#Merge bone marrow Human Cell Atlas data to LPL data:
bm1 <- readRDS("C:/Users/sszymura/Desktop/Currentanalysis/Data/bm1.Robj")
bm2 <- readRDS("C:/Users/sszymura/Desktop/Currentanalysis/Data/bm2.Robj")
bm3 <- readRDS("C:/Users/sszymura/Desktop/Currentanalysis/Data/bm3.Robj")
bm4 <- readRDS("C:/Users/sszymura/Desktop/Currentanalysis/Data/bm4.Robj")
bm5 <- readRDS("C:/Users/sszymura/Desktop/Currentanalysis/Data/bm5.Robj")
bm6 <- readRDS("C:/Users/sszymura/Desktop/Currentanalysis/Data/bm6.Robj")
bm7 <- readRDS("C:/Users/sszymura/Desktop/Currentanalysis/Data/bm7.Robj")
bm8 <- readRDS("C:/Users/sszymura/Desktop/Currentanalysis/Data/bm8.Robj")
bm9 <- readRDS("C:/Users/sszymura/Desktop/Currentanalysis/Data/bm9.Robj")
bm10 <- readRDS("C:/Users/sszymura/Desktop/Currentanalysis/Data/bm10.Robj")
bm11 <- readRDS("C:/Users/sszymura/Desktop/Currentanalysis/Data/bm11.Robj")
bm12 <- readRDS("C:/Users/sszymura/Desktop/Currentanalysis/Data/bm12.Robj")
bm13 <- readRDS("C:/Users/sszymura/Desktop/Currentanalysis/Data/bm13.Robj")
bm14 <- readRDS("C:/Users/sszymura/Desktop/Currentanalysis/Data/bm14.Robj")
bm15 <- readRDS("C:/Users/sszymura/Desktop/Currentanalysis/Data/bm15.Robj")

#Merge sequentially with LPL data
immune.anchors <- FindIntegrationAnchors(object.list = list(Patientscombinedbonemarrow, bm1), dims = 1:30) #repeat for all bm datasets
Patientscombinedbonemarrow <- IntegrateData(anchorset = immune.anchors, dims = 1:30)
saveRDS(Patientscombinedbonemarrow, file = "C:/Users/sszymura/Desktop/Patientscombinedbonemarrow.Robj")
DefaultAssay(Patientscombinedbonemarrow) <- "integrated"
Patientscombinedbonemarrow <- ScaleData(Patientscombinedbonemarrow, verbose = FALSE)
Patientscombinedbonemarrow <- RunPCA(Patientscombinedbonemarrow, npcs = 30, verbose = FALSE)
ElbowPlot(Patientscombinedbonemarrow)
Patientscombinedbonemarrow <- FindNeighbors(Patientscombinedbonemarrow, reduction = "pca", dims = 1:20)
Patientscombinedbonemarrow <- FindClusters(Patientscombinedbonemarrow, resolution = 0.5)
Patientscombinedbonemarrow <- RunUMAP(Patientscombinedbonemarrow, reduction = "pca", dims = 1:20)

#subset into separate main populations
Bcells <- Cells(Patientscombinedbonemarrow)[Patientscombinedbonemarrow$seurat_clusters == "1" | Patientscombinedbonemarrow$seurat_clusters == "9" | Patientscombinedbonemarrow$seurat_clusters == "10" | Patientscombinedbonemarrow$seurat_clusters == "14" ]
Tcells <- Cells(Patientscombinedbonemarrow)[Patientscombinedbonemarrow$seurat_clusters == "3" | Patientscombinedbonemarrow$seurat_clusters == "5" | Patientscombinedbonemarrow$seurat_clusters == "6" | Patientscombinedbonemarrow$seurat_clusters == "19" | Patientscombinedbonemarrow$seurat_clusters == "2" | Patientscombinedbonemarrow$seurat_clusters == "11" | Patientscombinedbonemarrow$seurat_clusters == "8" | Patientscombinedbonemarrow$seurat_clusters == "0" | Patientscombinedbonemarrow$seurat_clusters == "12"]
Monocytes <- Cells(Patientscombinedbonemarrow)[Patientscombinedbonemarrow$seurat_clusters == "17" | Patientscombinedbonemarrow$seurat_clusters == "21" | Patientscombinedbonemarrow$seurat_clusters == "16" | Patientscombinedbonemarrow$seurat_clusters == "7" | Patientscombinedbonemarrow$seurat_clusters == "4" | Patientscombinedbonemarrow$seurat_clusters == "18"]

#Create separate Seurat objects for B_cells, T_cells and Monocytes
Patientscombinedbonemarrow_bcells <- subset(Patientscombinedbonemarrow, cells = Bcells)
Patientscombinedbonemarrow_tcells <- subset(Patientscombinedbonemarrow, cells = Tcells)
Patientscombinedbonemarrow_monocytes <- subset(Patientscombinedbonemarrow, cells = Monocytes)

#Analyze T cells separately
Patientscombinedbonemarrow_tcells <- RunPCA(Patientscombinedbonemarrow_tcells, npcs = 30, verbose = FALSE)
Patientscombinedbonemarrow_tcells <- FindNeighbors(Patientscombinedbonemarrow_tcells, reduction = "pca", dims = 1:20)
Patientscombinedbonemarrow_tcells <- FindClusters(Patientscombinedbonemarrow_tcells, resolution = 0.5)
Patientscombinedbonemarrow_tcells <- RunUMAP(Patientscombinedbonemarrowtcells, reduction = "pca", dims = 1:20)
DimPlot(Patientscombinedbonemarrowtcells, reduction = 'umap', label = TRUE, repel = TRUE)
#Remove cluster that contains CD4/CD8 double negative cells.
Patientscombinedbonemarrow_tcells = subset(Patientscombinedbonemarrow_tcells, cells = Cells(Patientscombinedbonemarrow_tcells)[Patientscombinedbonemarrow_tcells$seurat_clusters != "11"])

#Analyze B cells separately
Patientscombinedbonemarrow_bcells <- RunPCA(Patientscombinedbonemarrow_bcells, npcs = 30, verbose = FALSE)
Patientscombinedbonemarrow_bcells <- FindNeighbors(Patientscombinedbonemarrow_bcells, reduction = "pca", dims = 1:20)
Patientscombinedbonemarrow_bcells <- FindClusters(Patientscombinedbonemarrow_bcells, resolution = 0.5)
Patientscombinedbonemarrow_bcells <- RunUMAP(Patientscombinedbonemarrowbcells, reduction = "pca", dims = 1:20)
DimPlot(Patientscombinedbonemarrowbcells, reduction = 'umap', label = TRUE, repel = TRUE)
Patientscombinedbonemarrow_bcells_data = GetAssatData(Patientscombinedbonemarrow_bcells)
Patientscombinedbonemarrow_bcells = AddMetaData(Patientscombinedbonemarrow_bcells, metadata = Patientscombinedbonemarrow_bcells_data['ENSG00000160654'], col.name = "CD3G")
Patientscombinedbonemarrow_bcells = AddMetaData(Patientscombinedbonemarrow_bcells, metadata = Patientscombinedbonemarrow_bcells_data['ENSG00000167286'], col.name = "CD3D")
Patientscombinedbonemarrow_bcells = AddMetaData(Patientscombinedbonemarrow_bcells, metadata = Patientscombinedbonemarrow_bcells_data['ENSG00000198851'], col.name = "CD3E")
#Filter out CD3 positive cells
Patientscombinedbonemarrow_bcells = subset(Patientscombinedbonemarrow_bcells, cells = Cells(Patientscombinedbonemarrow_bcells)[Patientscombinedbonemarrow_bcells$CD3G > 0])
Patientscombinedbonemarrow_bcells = subset(Patientscombinedbonemarrow_bcells, cells = Cells(Patientscombinedbonemarrow_bcells)[Patientscombinedbonemarrow_bcells$CD3D > 0])
Patientscombinedbonemarrow_bcells = subset(Patientscombinedbonemarrow_bcells, cells = Cells(Patientscombinedbonemarrow_bcells)[Patientscombinedbonemarrow_bcells$CD3E > 0])

#Adding information about tumor BCR clone based on the idiotype sequence used for vaccination. Create table with list of clonotypes that are tumor called "tumorclones"
bcellstumor = as.data.frame(cbind(Cells(bcells), bcells$BCR, rep(0, length(bcells$BCR))))
names(bcellstumor)[2] = "clonotypes"
names(bcellstumor)[3] = "tumorclone"
for (i in 1:length(bcellstumor$clonotypes)) {if (bcellstumor$clonotypes[i] %in% tumorclones$tumor) {bcellstumor$tumorclone[i] <- "Tumor"} else {bcellstumor$tumorclone[i] <- "None"}}
bcells <- AddMetaData(bcells, metadata = bcellstumor$tumorclone, col.name = "Tumor")

#Analyze monocyte cells separately
Patientscombinedbonemarrow_monocytes <- RunPCA(Patientscombinedbonemarrow_monocytes, npcs = 30, verbose = FALSE)
Patientscombinedbonemarrow_monocytes <- FindNeighbors(Patientscombinedbonemarrow_monocytes, reduction = "pca", dims = 1:20)
Patientscombinedbonemarrow_monocytes <- FindClusters(Patientscombinedbonemarrow_monocytes, resolution = 0.5)
Patientscombinedbonemarrow_monocytes <- RunUMAP(Patientscombinedbonemarrowmonocytes, reduction = "pca", dims = 1:20)
DimPlot(Patientscombinedbonemarrowmonocytes, reduction = 'umap', label = TRUE, repel = TRUE)

#Perform all downstream analysis on patient samples only
Patientscombinedmonocytes <- subset(Patientscombinedbonemarrow_monocytes, idents = c("PostPatient1", "PostPatient2", "PostPatient3", "PostPatient4", "PostPatient5", "PostPatient7", "PostPatient8", "PostPatient9", "PrePatient1", "PrePatient2", "PrePatient3", "PrePatient4", "PrePatient5", "PrePatient6", "PrePatient7", "PrePatient8", "PrePatient9"))
Patientscombinedbcells <- subset(Patientscombinedbonemarrow_bcells, idents = c("PostPatient1", "PostPatient2", "PostPatient3", "PostPatient4", "PostPatient5", "PostPatient7", "PostPatient8", "PostPatient9", "PrePatient1", "PrePatient2", "PrePatient3", "PrePatient4", "PrePatient5", "PrePatient6", "PrePatient7", "PrePatient8", "PrePatient9"))
Patientscombinedtcells <- subset(Patientscombinedbonemarrow_tcells, idents = c("PostPatient1", "PostPatient2", "PostPatient3", "PostPatient4", "PostPatient5", "PostPatient7", "PostPatient8", "PostPatient9", "PrePatient1", "PrePatient2", "PrePatient3", "PrePatient4", "PrePatient5", "PrePatient6", "PrePatient7", "PrePatient8", "PrePatient9"))

#Example of generation of heatmap for gene markers
bcells_markers <- FindAllMarkers(bcells, only.pos = TRUE, min.pct = 0.25, logfc.thershold = 0.25)
top30_bcells <- bcells_markers %>% group_by(cluster) %>% top_n(n = 30, wt = avg_logFC)
DoHeatmap(bcells, features = top30_bcells$gene) + NoLegend()

#Example of the DEG analysis and volcano plots
#Create cell annotations based on seurat cluster and treatment status (pre or post)
meta = as.data.frame(cbind(Cells(tcells), tcells$Treatment, tcells$seurat_cluster))
#Combine treatment and seurat cluster annotations into single annotation called "SeuratCondition" that will be used to analyze differential expression between different clusters
DefaultAssay(tcells) = "RNA"
Idents(tcells) = "SeuratCondition"
for (i in 0:11){clust1 = paste0("Pre_", i); clust2 = paste0("Post_", i); diff = FindMarkers(myeloid_vaccine, ident.1 = clust2, ident.2 = clust1, verbose = FALSE, min.pct = 0.05, logfc.threshold = 0, test.use = "wilcox"); name = paste0("/Users/sszymura/Desktop/vaccine_paper/tcells_wilcox/tcells_clust", i, ".csv"); write.csv(diff, file = name)}
#This last line will create differential expression files for each cluster pre vs post in a folder that you specify. Then need to read this file back into R. This list is used to load into IPA software to identify differential pathways and biological processess
#Example of creation of volcano plots
b.interferon.response = read.csv("/Users/sszymura/Desktop/vaccine_paper/tcells_wilcox/tcells_clust0.csv")

library("EnhancedVolcano")
library("AnnotationDbi")
library("org.Hs.eg.db")
library("EnsDb.Hsapiens.v75")
library("rtracklayer")

b.interferon.response$ens75 = mapIds(EnsDb.Hsapiens.v75,
                                     keys=b.interferon.response$X, 
                                     column="SYMBOL",
                                     keytype="GENEID",
                                     multiVals="first")
EnhancedVolcano(b.interferon.response, lab=b.interferon.response$ens75, x="avg_log2FC", y="p_val", FCcutoff = log2(1.25), pCutoff = 0.05, pointSize = 1.2, labSize = 4, colAlpha = 1, title=NULL, subtitle=NULL, caption=NULL, xlim = c(-1.5,1.5)) #legendPosition="right", legendPosition = "none") 

#Create heatmaps of pathways and biological processess identified by IPA using file created in IPA software.
library(pheatmap)
data <- read.csv("C:/Users/sszymura/Desktop/bcells_post_vs_pre_responders_040521.csv", header = TRUE, row.names = 1)
pheatmap(as.matrix(data)[1:80, 1:4], cluster_rows = FALSE, cluster_cols = FALSE)

#Example of creating Heatmaps and DotPlots of HLA gene marker
tumor_bcells = subset(bcells, cells = Cells(bcells)[bcells$Tumor == "Tumor"])
Idents(tumor_bcells) = "seurat_clusters"
#Read "features" file from any of the original patient samples to obtain a list of gene names with ENSEMBL symbols
features = read.table("/Users/sszymura/Desktop/features.tsv", sep = '\t')

genes <- as.data.frame(c("HLA-DMA", "HLA-DMB", "HLA-DOA", "HLA-DOB", "HLA-DPA1", "HLA-DPB1", "HLA-DQA1", "HLA-DQA2", "HLA-DQB1", "HLA-DQB2", "HLA-DRA", "HLA-DRB1", "HLA-DRB3", "HLA-DRB4", "HLA-DRB5", "HLA-A", "HLA-B", "HLA-C", "HLA-E", "HLA-F", "HLA-G"))

library(dplyr)
names(genes) <- "V2"
genes = left_join(genes, features, by="V2")
DefaultAssay(bcells_vaccine_pre) = "RNA"
DoHeatmap(bcells_vaccine_pre, group.by = "Treatment", features = genes$V1)
DotPlot(bcells_cluster2, features= genes$V1, cols = c("darkred", "white", "darkblue"))

DoHeatmap(tumor_bcells, features = c("HLA-DMA", "HLA-DMB", "HLA-DOA", "HLA-DOB", "HLA-DPA1", "HLA-DPB1", "HLA-DQA1", "HLA-DQA2", "HLA-DQB1", "HLA-DRA", "HLA-DRB1", "HLA-DRB5", "CD74"), assay = "RNA", cells = WhichCells(tumor_bcells, idents = c("0", "1", "2", "5", "10")))

#UMAP plot of T cells with TCR clonotypes
DimPlot(tcells, cells.highlight = Cells(tcells)[tcells$TCR != "None"])

#Plotting TCR correlation plots
plot(jitter(log10(Patient1$Freq.x + 1), 2) ~ jitter(log10(Patient1$Freq.y + 1), 2), xlim = c(0, 0.4), ylim = c(0, 0.4), pch = 1)
abline(0,1)

#Plotting top 20 clonotypes pre and post vaccine
all_tcrs <- repLoad("/Users/sszymura/Desktop/TCR/untitledfolder/", .mode = "paired")
library(dplyr)
pre = as.data.frame(all_tcrs$data[[1]]$CDR3.aa)
names(pre) = "clonotypes"
post = as.data.frame(all_tcrs$data[[2]]$CDR3.aa)
names(post) = "clonotypes"
clonotypes = full_join(pre, post, by = "clonotypes")
tc2 <- trackClonotypes(all_tcrs$data[c(7,8)], .which = list(1, 20), .col = "aa+v")
p2 <- vis(tc2)
p2

#Example of the analysis of shannon entropy
devtools::install_github("kylebittinger/abdiv")
setwd("/Users/sszymura/Desktop/Final_TCR_three_methods2/")
lpl1 = read.csv("lpl1.csv", stringsAsFactors = FALSE)
shannon(lpl1$Pre_Immunarch, base = exp(1))
shannon(lpl1$Post_Immunarch, base = exp(1))
##

#UMAP of top 20 clonotypes
separated_patients <- SplitObject(Patientscombinedbonemarrow_tcells, split.by = "Patient")
setwd("C:/Users/sszymura/Desktop/NewestTCRanalysis/")
files <- list.files("C:/Users/sszymura/Desktop/NewestTCRanalysis/adjusted/")
patient1 <- read.csv("Patient1.csv", stringsAsFactors = FALSE)
patient2 <- read.csv("Patient2.csv", stringsAsFactors = FALSE)
patient3 <- read.csv("Patient3.csv", stringsAsFactors = FALSE)
patient4 <- read.csv("Patient4.csv", stringsAsFactors = FALSE)
patient5 <- read.csv("Patient5.csv", stringsAsFactors = FALSE)
patient7 <- read.csv("Patient7.csv", stringsAsFactors = FALSE)
patient8 <- read.csv("Patient8.csv", stringsAsFactors = FALSE)
patient9 <- read.csv("Patient9.csv", stringsAsFactors = FALSE)
patient1_top20 <- patient1[1:20,]
patient2_top20 <- patient2[1:20,]
patient3_top20 <- patient3[1:20,]
patient4_top20 <- patient4[1:20,]
patient5_top20 <- patient5[1:20,]
patient7_top20 <- patient7[1:20,]
patient8_top20 <- patient8[1:20,]
patient9_top20 <- patient9[1:20,]


#Violin plots of top20 clonotypes gene expression pre and post vaccination
expression_test_bcells = function(seurat_object, cluster_pre, cluster_post, gene){
  gene = features[features$V2 == gene, ]$V1
  myeloid_8_pre = subset(seurat_object, cells = Cells(seurat_object)[seurat_object$Condition == cluster_pre])
  myeloid_8_post = subset(seurat_object, cells = Cells(seurat_object)[seurat_object$Condition == cluster_post])
  myeloid_8_pre_expression = GetAssayData(myeloid_8_pre, assay = "RNA")
  myeloid_8_post_expression = GetAssayData(myeloid_8_post, assay = "RNA")
  test = wilcox.test(myeloid_8_pre_expression[(rownames(myeloid_8_pre_expression) == gene), ], myeloid_8_post_expression[(rownames(myeloid_8_post_expression) == gene), ], alternative = "two.sided")
  print(test)}

difference_tcells = function(seurat_object, cluster_pre, cluster_post, gene){
  gene = features[features$V2 == gene, ]$V1
  myeloid_8_pre = subset(seurat_object, cells = Cells(seurat_object)[seurat_object$Condition == cluster_pre])
  myeloid_8_post = subset(seurat_object, cells = Cells(seurat_object)[seurat_object$Condition == cluster_post])
  myeloid_8_pre_expression = GetAssayData(myeloid_8_pre, assay = "RNA")
  myeloid_8_post_expression = GetAssayData(myeloid_8_post, assay = "RNA")
  post = mean(myeloid_8_post_expression[(rownames(myeloid_8_post_expression) == gene), ])
  pre = mean(myeloid_8_pre_expression[(rownames(myeloid_8_pre_expression) == gene), ])
  print(pre)
  print(post)}

Idents(tcells_top20) = "Condition"
Idents(tcells_not_top20) = "Condition"
expression_test_bcells(tcells_top20, "Pre", "Post", "ENTPD1")
difference_tcells(tcells_top20, "Pre", "Post", "ENTPD1")


#Only single Tm unique expanded, shared expanded
tcr_processing <- function(lpl001){
  for (i in 1:length(lpl001$Var1)) {
    if (lpl001$Freq.x[i] <= 1 && lpl001$Freq.y[i] <= 1) {lpl001$analysis[i] = "Background"}
    else if (lpl001$Freq.x[i] == 0 && lpl001$Freq.y[i] > 1) {lpl001$analysis[i] = "Emerged"}
    else if (lpl001$Freq.x[i] != 0 && (lpl001$Freq2.y[i]/lpl001$Freq2.x[i] > 1)) {lpl001$analysis[i] = "Expanded"}
    else if (lpl001$Freq.y[i] == 0 && lpl001$Freq.x[i] > 1) {lpl001$analysis[i] = "Dissapeared"}
    else if (lpl001$Freq.y[i] != 0 && (lpl001$Freq2.x[i]/lpl001$Freq2.y[i] > 1)) {lpl001$analysis[i] = "Shrunk"}
    else (lpl001$analysis[i] == "Unknown")
  }; return(lpl001)}

#Non top 20 clonotypes
load("/Users/sszymura/Desktop/tcells_vaccine_subset_020124.Robj")
VlnPlot(tcells_top20, features = "ENSG00000138185", group.by = "Condition")
top20_barcodes = as.data.frame(Cells(tcells_vaccine_subset_scaled_top20))
tcells_not_top20 = subset(tcells, cells = Cells(tcells)[!Cells(tcells) %in% top20_barcodes$`Cells(tcells_top20)`])
VlnPlot(tcells_not_top20, features = "ENSG00000138185", group.by = "Condition")

