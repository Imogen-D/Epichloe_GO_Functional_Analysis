#Read in file from http://ekhidna2.biocenter.helsinki.fi/sanspanz/
#This code based on https://github.com/dwinter/genome_factory/wiki/Taking-a-gene-list-and-PANNZER-output-through-to-GO-term-enrichment

library(topGO)

AmarillansAnnoout <- read.delim("~/CoxExtension/GO output/AmarillansAnnoout.txt", stringsAsFactors=FALSE, na.strings = "n.d.")

go_filt <- AmarillansAnnoout[which(AmarillansAnnoout$PPV >= 0.5),]   #only with high PPV                          
go_filt$id <- paste0('GO:', go_filt$id)
panzer_to_golist <- function(panzer_df){
  go_df <- aggregate( id ~ qpid, data=panzer_df, FUN=c)
  structure(go_df$id, .Names=go_df$qpid)
}
all_golist <- panzer_to_golist(go_filt)
str(head(all_golist))

make_topGO_DO <- function(gene_list, ontology, gene2GO_list){
  topGO_data <- new("topGOdata", ontology = ontology, allGenes = gene_list,
                    annot = annFUN.gene2GO, gene2GO = gene2GO_list)
  fishers_result <- runTest(topGO_data, algorithm = "elim", statistic = "fisher")
  fishers_table <- GenTable(topGO_data, Fishers = fishers_result, useLevels = TRUE)
  fishers_table$Ontology <- ontology
  fishers_table
}

amarillansuniquegeneid <- read.table("~/CoxExtension/LineageSpecificity/amarillansuniquegeneid.tsv", quote="\"", comment.char="", stringsAsFactors=FALSE)
ortho_long <- read.delim("~/CoxExtension/ortho_long.tsv", header=FALSE, stringsAsFactors=FALSE)
Eam702 <- ortho_long[which(ortho_long$V3 == "Eam702"),]
uniquegenes <- amarillansuniquegeneid$V1
geneList <- factor(as.integer(Eam702$V2 %in% uniquegenes))
names(geneList) <- Eam702$V2
str(geneList)

topGO_BP_table <- make_topGO_DO(geneList, "BP", all_golist)
topGO_MF_table <- make_topGO_DO(geneList, "MF", all_golist)
topGO_CC_table <- make_topGO_DO(geneList, "CC", all_golist)

topGO_all_table <- (rbind(topGO_BP_table, topGO_MF_table, topGO_CC_table))
topGO_all_table <- topGO_all_table[order(topGO_all_table$Fishers),]
View(topGO_all_table)

write.csv(topGO_all_table, file = "topGOfull.csv", quote = FALSE)


#for bromicola


BromicolaAnnoout <- read.delim("~/CoxExtension/GO output/BromicolaAnnoout.txt", stringsAsFactors=FALSE, na.strings =  "n.d.")
subsetbromGO <- BromicolaAnnoout[which(BromicolaAnnoout$PPV >= 0.5), c(1, 5)]

go_filt_brom <- BromicolaAnnoout[which(BromicolaAnnoout$PPV >= 0.5),]                             
go_filt_brom$id <- paste0('GO:', go_filt_brom$id)
all_golist_brom <- panzer_to_golist(go_filt_brom)
str(head(all_golist_brom))

bromicolauniquegeneid <- read.table("~/CoxExtension/LineageSpecificity/bromicolauniquegeneid.tsv", quote="\"", comment.char="", stringsAsFactors=FALSE)
ortho_long <- read.delim("~/CoxExtension/ortho_long.tsv", header=FALSE, stringsAsFactors=FALSE)
Ebrom <- ortho_long[which(ortho_long$V3 == "EbroNfe1"),]
uniquegenes_brom <- bromicolauniquegeneid$V1
geneList_brom <- factor(as.integer(Ebrom$V2 %in% uniquegenes_brom))
names(geneList_brom) <- Ebrom$V2
str(geneList_brom)

topGO_BP_table_brom <- make_topGO_DO(geneList_brom, "BP", all_golist_brom)
topGO_MF_table_brom <- make_topGO_DO(geneList_brom, "MF", all_golist_brom)
topGO_CC_table_brom <- make_topGO_DO(geneList_brom, "CC", all_golist_brom)

topGO_all_table_brom <- (rbind(topGO_BP_table_brom, topGO_MF_table_brom, topGO_CC_table_brom))
topGO_all_table_brom <- topGO_all_table_brom[order(topGO_all_table_brom$Fishers),]
View(topGO_all_table_brom)

write.csv(topGO_all_table_brom, file = "topGOfull_brom.csv", quote = FALSE)
