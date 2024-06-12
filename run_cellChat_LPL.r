library(CellChat)
library(ggplot2)
library(ggalluvial)
options(stringsAsFactors = FALSE)
library(extrafont)
loadfonts(device = "pdf", quiet = F)


data.input <- GetAssayData(LPL_scRNA_data, assay = "RNA", slot = "data")
labels <- Idents(LPL_scRNA_data)
identity <- data.frame(group = labels, row.names = names(labels)) # create a dataframe of the cell labels
cellchatLPL <- createCellChat(data = data.input)
cellchatLPL <- addMeta(cellchatLPL, meta = identity, meta.name = "labels")
cellchatLPL <- setIdent(cellchatLPL, ident.use = "labels") # set "labels" as default cell identity
levels(cellchatLPL@idents) # show factor levels of the cell labels
CellChatDB <- CellChatDB.human # use CellChatDB.human if running on human data 
CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling") # use Secreted Signaling for cell-cell communication analysis
cellchatLPL@DB <- CellChatDB.use # set the used database in the objec
cellchatLPL <- subsetData(cellchatLPL) # subset the expression data of signaling genes for saving computation cost
options(future.globals.maxSize= 1691289600)
future::plan("multiprocess", workers = 6) # do parallel
cellchatLPL <- identifyOverExpressedGenes(cellchatLPL)
cellchatLPL <- identifyOverExpressedInteractions(cellchatLPL)
cellchatLPL <- projectData(cellchatLPL, LPLI.human)
cellchatLPL <- computeCommunProb(cellchatLPL)
cellchatLPL <- computeCommunProbPathway(cellchatLPL)
cellchatLPL <- aggregateNet(cellchatLPL)
groupSize <- as.numeric(table(cellchatLPL@idents)) # number of cells in each cell group
netVisual_circle(cellchatLPL@net$count, vertex.size = groupSize, weight.scale = T, label.edge= T, edge.label.cex = 0.8, vertex.label.cex = 1)
netVisual_circle(cellchatLPL@net$count, vertex.size = groupSize, weight.scale = T, label.edge= F, edge.label.cex = 0.8, vertex.label.cex = 1)
vertex.receiver = seq(1,3)
##output all the significant pathway

cellchatLPL@netP$pathways
pathways.show <- "TGFb"
netVisual_aggregate(cellchatLPL, signaling = pathways.show,  vertex.receiver = vertex.receiver, vertex.size = groupSize)
netVisual_aggregate(cellchatLPL, signaling = c("TGFb"), layout = "circle", vertex.size = groupSize)
netAnalysis_contribution(cellchatLPL, signaling = pathways.show)


cellchatLPL_Pre<-readRDS("cellchatLPL_Pre.rds")
cellchatPL_Post<-readRDS("cellchatPL_Post.rds")
cellchat <- mergeCellChat(list(cellchatLPL_Pre, cellchatPL_Post), add.names = c("Pre","Post"))

cellchat <- computeNetSimilarityPairwise(cellchat, type = "structural")
cellchat <- netEmbedding(cellchat, type = "structural")
cellchat <- netClustering(cellchat, type = "structural")
# Visualization
netVisual_embeddingPairwise(cellchat, type = "structural")
netVisual_embeddingPairwiseZoomIn(cellchat, type = "structural")
rankSimilarity(cellchat, type = "structural")
rankNet(cellchat, mode = "comparison")




