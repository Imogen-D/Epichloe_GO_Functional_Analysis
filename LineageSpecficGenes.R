#read in long orthos

orthos <- read.table("ortho_long.tsv", stringsAsFactors=FALSE, col.names=c("ortho_id", "gene_id", "strain"))
# convert to named list, with strains represente by each ortho
by_ortho <- tapply(orthos$strain, orthos$ortho_id, FUN=unique)


strains <- unique(orthos$strain)

#fxn to find orthos with only one set of IDs
lineage_specific <- function(ortho_set, lineage_ids){
  all(ortho_set %in% lineage_ids) & length(ortho_set) == length(lineage_ids)                                                                                      
  }
# which indices in by_ortho are amarillans speific.
all_amarillans <- sapply(by_ortho, lineage_specific, c("EamaE4668", "Eam702", "EamaE57"))
sum(all_amarillans) #142
long_amarillans <- sapply(by_ortho, lineage_specific, c("Eam702"))
sum(long_amarillans) #302

trueamaog <- data.frame(which(long_amarillans == TRUE))
colnames(trueamaog) <- "num"

shared_null <- function(i = sample(strains, 3)){
  random <- sapply(by_ortho, lineage_specific, c(i))
  shared <- sum(random)
  return(shared)
}

random_shared <- replicate(1000, shared_null())
mean(random_shared) #0.394

mean(random_shared >= 142) #0.001

all_bromicola <- sapply(by_ortho, lineage_specific, c("EbroAL0426", "EbroAL0434", "EbroNfe1"))
sum(all_bromicola) #32
long_bromicola <- sapply(by_ortho, lineage_specific, c("EbroNfe1"))
sum(long_bromicola) #394
truebromog <- data.frame(which(long_bromicola == TRUE))
colnames(truebromog) <- "num"

mean(random_shared >= 32) #0.003
hist(random_shared)
abline(v=c(32, 142))

amaorthonum <- by_ortho[trueamaog$num]
amaorthos <- data.frame(names(amaorthonum))

bromorthonum <- by_ortho[truebromog$num]
bromorthos <- data.frame(names(bromorthonum))

amageneID <- orthos[c(which(orthos$ortho_id %in% amaorthos$names.amaorthonum.)),]
write.table(amageneID$gene_id, file = "amarillansuniquegeneid.tsv", row.names = FALSE, col.names = FALSE, quote = FALSE)
bromgeneID <- orthos[c(which(orthos$ortho_id %in% bromorthos$names.bromorthonum.)),]
write.table(bromgeneID$gene_id, file = "bromicolauniquegeneid.tsv", row.names = FALSE, col.names = FALSE, quote = FALSE)
