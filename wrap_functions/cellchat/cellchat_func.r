
###### load zebrafish cellchatdb database ######

zebrafishcellchatdb <- function(){
	print("Set the ligand-receptor interaction database")
	CellChatDB <- CellChatDB.zebrafish # use CellChatDB.mouse if running on mouse data
	showDatabaseCategory(CellChatDB)
	dplyr::glimpse(CellChatDB$interaction)
	
	# use a subset of CellChatDB for cell-cell communication analysis
	CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling") # use Secreted Signaling
	
	## 去除斑马鱼中重复的ligand-receptors
	new_interaction <- NULL
	all_pathway_name <- unique(CellChatDB.use$interaction$pathway_name)
	
	lrnum <- matrix(0, length(all_pathway_name), 4)
	for (i in seq_along(all_pathway_name)){
		print(all_pathway_name[i])
		temp <- dplyr::filter(CellChatDB.use$interaction, pathway_name == all_pathway_name[i])
		lrnum[i, c(1, 2)] <- dim(temp)
		temp <- temp %>% distinct(interaction_name_2, .keep_all = TRUE)
		lrnum[i, c(3, 4)] <- dim(temp)
		new_interaction <- rbind(new_interaction, temp)
	}
	
	rownames(lrnum) <- all_pathway_name
	colnames(lrnum) <- c("pre_nrow", "pre_ncol", "rm_duplicate_nrow", "rm_duplicate_ncol")
	
	CellChatDB.use$interaction <- new_interaction
	cellchat@DB <- CellChatDB.use
	
	saveRDS(CellChatDB.use, "/home/zhengjh/scripts/cellchat/cellchatdb/rm_duplicate_Secreted_Signaling_lr_zebrafish_cellchatdb.rds")
}

ratcellchatdb <- function(){
	
	library(CellChat)
	print("Set the ligand-receptor interaction database")
	CellChatDB <- CellChatDB.mouse # use CellChatDB.mouse if running on mouse data
	showDatabaseCategory(CellChatDB)
	dplyr::glimpse(CellChatDB$interaction)
	
	# use a subset of CellChatDB for cell-cell communication analysis
	CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling") # use Secreted Signaling
	
	## 去除斑马鱼中重复的ligand-receptors
	new_interaction <- NULL
	all_pathway_name <- unique(CellChatDB.use$interaction$pathway_name)
	
	lrnum <- matrix(0, length(all_pathway_name), 4)
	for (i in seq_along(all_pathway_name)){
		print(all_pathway_name[i])
		temp <- dplyr::filter(CellChatDB.use$interaction, pathway_name == all_pathway_name[i])
		lrnum[i, c(1, 2)] <- dim(temp)
		temp <- temp %>% distinct(interaction_name_2, .keep_all = TRUE)
		lrnum[i, c(3, 4)] <- dim(temp)
		new_interaction <- rbind(new_interaction, temp)
	}
	
	rownames(lrnum) <- all_pathway_name
	colnames(lrnum) <- c("pre_nrow", "pre_ncol", "rm_duplicate_nrow", "rm_duplicate_ncol")
	
	CellChatDB.use$interaction <- new_interaction

	saveRDS(CellChatDB.use, "/home/zhengjh/scripts/cellchat/cellchatdb/rm_duplicate_Secreted_Signaling_lr_rat_cellchatdb.rds")
}

monkeycellchatdb <- function(){
	
	library(CellChat)
	print("Set the ligand-receptor interaction database")
	CellChatDB <- CellChatDB.human # use CellChatDB.mouse if running on mouse data
	showDatabaseCategory(CellChatDB)
	dplyr::glimpse(CellChatDB$interaction)
	
	# use a subset of CellChatDB for cell-cell communication analysis
	CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling") # use Secreted Signaling
	
	## 去除斑马鱼中重复的ligand-receptors
	new_interaction <- NULL
	all_pathway_name <- unique(CellChatDB.use$interaction$pathway_name)
	
	lrnum <- matrix(0, length(all_pathway_name), 4)
	for (i in seq_along(all_pathway_name)){
		print(all_pathway_name[i])
		temp <- dplyr::filter(CellChatDB.use$interaction, pathway_name == all_pathway_name[i])
		lrnum[i, c(1, 2)] <- dim(temp)
		temp <- temp %>% distinct(interaction_name_2, .keep_all = TRUE)
		lrnum[i, c(3, 4)] <- dim(temp)
		new_interaction <- rbind(new_interaction, temp)
	}
	
	rownames(lrnum) <- all_pathway_name
	colnames(lrnum) <- c("pre_nrow", "pre_ncol", "rm_duplicate_nrow", "rm_duplicate_ncol")
	
	CellChatDB.use$interaction <- new_interaction
	
	saveRDS(CellChatDB.use, "/home/zhengjh/scripts/cellchat/cellchatdb/rm_duplicate_Secreted_Signaling_lr_monkey_cellchatdb.rds")
}



#### CreateCellChatObject #######

CreateCellChatObject <- function(seurat_rds, CellChatDB, thresh_p){
	
	data.input = seurat_rds@assays$RNA@data # normalized data matrix
	meta <- seurat_rds@meta.data # a dataframe with rownames containing cell mata data
	meta$labels <- meta$celltype_assign
	unique(meta$labels)
	#### Create a CellChat object ####
	print("Create a CellChat object")
	cellchat <- createCellChat(object = data.input, meta = meta, group.by = "labels")
	cellchat <- addMeta(cellchat, meta = meta)
	cellchat <- setIdent(cellchat, ident.use = "labels") # set "labels" as default cell identity
	print(levels(cellchat@idents)) # show factor levels of the cell labels

	print("load the ligand-receptor interaction database")
	
	cellchat@DB <- CellChatDB
	# duplicated(CellChatDB.use$interaction$interaction_name_2)
	### Preprocessing the expression data for cell-cell communication analysis #####
	print("Preprocessing the expression data for cell-cell communication analysis")
	cellchat <- subsetData(cellchat) # This step is necessary even if using the whole database
	cellchat <- identifyOverExpressedGenes(cellchat, thresh.p = thresh_p)
	cellchat <- identifyOverExpressedInteractions(cellchat)
	
	return(cellchat)
}


AggregateNetWrap <- function(cellchat, thresh_p){
	
	folder_name="1.InferenceCellCellCommunicationNetwork"
	if (file.exists(folder_name)){
		print("1.InferenceCellCellCommunicationNetwork existed.")
	}else{
		dir.create(folder_name)
	}
	
	print("Compute the communication probability and infer cellular communication network")
	cellchat <- computeCommunProb(cellchat)
	cellchat@options$parameter
	cellchat <- filterCommunication(cellchat, min.cells = 10)
	print("Infer the cell-cell communication at a signaling pathway level")
	cellchat <- computeCommunProbPathway(cellchat, thresh = thresh_p)
	print("Calculate the aggregated cell-cell communication network")
	cellchat <- aggregateNet(cellchat, thresh = thresh_p)
	
	file_name1 <- "Number_of_interactions_strength.pdf"
	pdf(file_name1, 10, 7)
	par(mfrow = c(1,2), xpd=TRUE)
	print(netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions"))
	print(netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength"))
	dev.off()
	
	mat <- cellchat@net$weight
	file_name2 <- "netVisual_circle_each_celltype.pdf"
	pdf(file_name2, 12, 6)
	par(mfrow = c(2,4), xpd=TRUE)
	for (i in 1:nrow(mat)) {
		mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
		mat2[i, ] <- mat[i, ]
		netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
	}
	dev.off()
	
	result_files <- c(file_name1, file_name2)
	file.copy(result_files, folder_name, overwrite = T)
	file.remove(result_files) 
	
	return(cellchat)
}


VisualizeSignalPathway <- function(cellchat, pathways.show){
	
	folder_name="2.VisualizationCellCellCommunicationNetwork"
	if (file.exists(folder_name)){
		print("2.VisualizationCellCellCommunicationNetwork existed.")
	}else{
		dir.create(folder_name)
	}
	
	# Circle plot
	file_name1 <- "Using_CirclePlot_visualizing_each_signaling_pathway.pdf"
	pdf(file_name1, 20, 40)
	par(mfrow=c(ceiling(length(pathways.show)/4), 4))
	for(i in 1:length(pathways.show)){
		netVisual_aggregate(cellchat, signaling = pathways.show[i], layout = "circle")
	}
	dev.off()
	
	# Chord diagram
	file_name2 <- "Using_ChordPlot_visualizing_each_signaling_pathway.pdf"
	pdf(file_name2, 20, 40)
	par(mfrow=c(ceiling(length(pathways.show)/4), 4))
	for(i in 1:length(pathways.show)){
		netVisual_aggregate(cellchat, signaling = pathways.show[i], layout = "chord")
	}
	dev.off()
	
	# Heatmap
	file_name3 <- "Using_HeatmapPlot_visualizing_each_signaling_pathway.pdf"
	pdf(file_name3, 10, 10)
	par(mfrow=c(ceiling(length(pathways.show)/4), 4))
	for(i in 1:length(pathways.show)){
		print(netVisual_heatmap(cellchat, signaling = pathways.show[i], color.heatmap = "Reds"))
	}
	dev.off()
	
	for (i in 1:length(pathways.show)) {
		print(pathways.show[i])
		# Visualize communication network associated with both signaling pathway and individual L-R pairs
		netVisual(cellchat, signaling = pathways.show[i], layout = "circle", out.format = "pdf")
		# Compute and visualize the contribution of each ligand-receptor pair to the overall signaling pathway
		gg <- netAnalysis_contribution(cellchat, signaling = pathways.show[i])
		ggsave(filename=paste0(folder_name, "/", pathways.show[i], "_L-R_contribution.pdf"), plot=gg, width = 6, height = 4, units = 'in', dpi = 300)
	}
	
	result_files <- c(file_name1, file_name2, file_name3)
	file.copy(result_files, folder_name, overwrite = T)
	file.remove(result_files) 
	
}


VisualizeLRs <- function(cellchat, folder_name, pathways.show){
	
	### Bubble plot ###
	file_name1 <- "BubblePlot_multiple_ligand_receptor.pdf"
	pdf(file_name1, 12, 10)
	par(mfrow=c(ceiling(length(levels(cellchat@idents))/2), 2))
	for (i in 1:length(levels(cellchat@idents))){
		print(netVisual_bubble(cellchat, sources.use = i, targets.use = 1:length(levels(cellchat@idents)), remove.isolate = FALSE))
	}
	dev.off()
	
	##### Plot the signaling gene expression distribution using violin/dot plot ######
	file_name2 <- "plotGeneExpression.pdf"
	pdf(file_name2, 12, 10)
	for ( i in 1:length(pathways.show)){
		print(plotGeneExpression(cellchat, signaling = pathways.show[i]))
	}
	dev.off()
	
	file_name3 <- "netVisual_chord_gene_LR.pdf"
	pdf(file_name3, 16, 16)
	netVisual_chord_gene(cellchat, lab.cex = 0.5,legend.pos.y = 30)
	dev.off()
	
	file_name4 <- "cellchat_lr.csv"
	df.net <- subsetCommunication(cellchat)
	write.csv(df.net, file_name4, quote = F,sep = ',')

	result_files <- c(file_name1, file_name2, file_name3, file_name4)
	file.copy(result_files, folder_name, overwrite = T)
	file.remove(result_files)
	
}

SystemsAnalysisCellCellNetwork <- function(cellchat, nPatterns, pathways.show){
	
	print("3.SystemsAnalysisCellCellCommunicationNetwork")
	folder_name="3.SystemsAnalysisCellCellCommunicationNetwork"
	if (file.exists(folder_name)){
		print("3.SystemsAnalysisCellCellCommunicationNetwork existed.")
	}else{
		dir.create(folder_name)
	}
	
	# Compute the network centrality scores
	cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP") # the slot 'netP' means the inferred intercellular communication network of signaling pathways
	
	# Visualize the computed centrality scores using heatmap, allowing ready identification of major signaling roles of cell groups
	file_name1 <- "netAnalysis_signalingRole_network.pdf"
	pdf(file_name1, 8, 6)
	for ( i in 1:length(pathways.show)){
		print(netAnalysis_signalingRole_network(cellchat, signaling = pathways.show[i], width = 8, height = 2.5, font.size = 10))
	}
	dev.off()
	
	#### Identify signals contributing most to outgoing or incoming signaling of certain cell groups ####
	# Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
	# ht1 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "outgoing")
	# ht2 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "incoming")
	
	file_name2 <- "netAnalysis_signalingRole_heatmap.pdf"
	pdf(file_name2, 12, 8)
	print(ht1 + ht2)
	dev.off()
	
	#### Identify global communication patterns to explore how multiple cell types and signaling pathways coordinate together ####
	
	#selectK(cellchat, pattern = "outgoing")
	#nPatterns=3
	cellchat <- identifyCommunicationPatterns(cellchat, pattern = "outgoing", k = nPatterns)
	
	file_name3 <- "Outgoing_communication_patterns.pdf"
	pdf(file_name3, 10, 8)
	# river plot
	print(netAnalysis_river(cellchat, pattern = "outgoing"))
	# dot plot
	print(netAnalysis_dot(cellchat, pattern = "outgoing"))
	dev.off()
	
	#### Identify and visualize incoming communication pattern of target cells ######
	# selectK(cellchat, pattern = "incoming")
	# nPatterns = 3
	cellchat <- identifyCommunicationPatterns(cellchat, pattern = "incoming", k = nPatterns)
	
	file_name4 <- "Incoming_communication_patterns.pdf"
	pdf(file_name4, 10, 8)
	# river plot
	print(netAnalysis_river(cellchat, pattern = "incoming"))
	# dot plot
	print(netAnalysis_dot(cellchat, pattern = "incoming"))
	dev.off()
	
	#### Identify signaling groups based on their functional similarity #####
	
	# cellchat <- computeNetSimilarity(cellchat, type = "functional")
	# cellchat <- netEmbedding(cellchat, type = "functional")
	#> Manifold learning of the signaling networks for a single dataset
	# cellchat <- netClustering(cellchat, type = "functional")
	#> Classification learning of the signaling networks for a single dataset
	# Visualization in 2D-space
	# file_name5 <- "FunctionalSimilarity.pdf"
	# pdf(file_name5, 10, 8)
	# print(netVisual_embedding(cellchat, type = "functional", label.size = 3.5))
	# dev.off()
	
	#### Identify signaling groups based on structure similarity ####
	
	# cellchat <- computeNetSimilarity(cellchat, type = "structural")
	# cellchat <- netEmbedding(cellchat, type = "structural")
	#> Manifold learning of the signaling networks for a single dataset
	# cellchat <- netClustering(cellchat, type = "structural")
	#> Classification learning of the signaling networks for a single dataset
	# Visualization in 2D-space
	
	# file_name6 <- "StructuralSimilarity.pdf"
	# pdf(file_name6, 10, 8)
	# print(netVisual_embedding(cellchat, type = "structural", label.size = 3.5))
	# dev.off()
	
	result_files <- c(file_name1, file_name2, file_name3, file_name4)
	file.copy(result_files, folder_name, overwrite = T)
	file.remove(result_files)
	
	return(cellchat)
}


