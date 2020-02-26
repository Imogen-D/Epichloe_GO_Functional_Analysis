ortho_long <- read.delim("ortho_long.tsv", header=FALSE, stringsAsFactors=FALSE)
#quant files from full quantification in https://github.com/Imogen-D/Cyclic-Peptide-Expression
Eamquant <- read.delim("Eamquant.sf", stringsAsFactors=FALSE)

Eam702 <- ortho_long[which(ortho_long$V3 == "Eam702"),]
colnames(Eam702) <- c("ortholog", "Name", "Strain")
mergedTPMandid <- merge(Eam702, Eamquant)

amarillansuniquegeneid <- read.table("amarillansuniquegeneid.tsv", quote="\"", comment.char="", stringsAsFactors=FALSE)
uniquegenes <- amarillansuniquegeneid$V1
geneList <- factor(as.integer(Eam702$Name %in% uniquegenes))
names(geneList) <- Eam702$Name
TFgenes <- data.frame(TF = geneList, Name = Eam702$Name, row.names = NULL)
mergedTPMandTF <- merge(mergedTPMandid, TFgenes)

Amaclosest_AT <- read.delim("Amaclosest_AT.tsv", header=FALSE, stringsAsFactors=FALSE)
library(stringr)
Amaclosest_AT$V1 <- str_remove(Amaclosest_AT$V1, "ID=")
Amaclosest_AT$V1 <- str_replace(Amaclosest_AT$V1, ";", "-T1")
colnames(Amaclosest_AT) <- c("Name", "DistanceAT")

Amarillansfinalframe <- merge(Amaclosest_AT, mergedTPMandTF)

write.csv(Amarillansfinalframe, file = "AmarillansData.csv")



##for Bromicola

EbroNfe1 <- ortho_long[which(ortho_long$V3 == "EbroNfe1"),]
colnames(EbroNfe1) <- c("ortholog", "Name", "Strain")
bromicolauniquegeneid <- read.table("bromicolauniquegeneid.tsv", quote="\"", comment.char="", stringsAsFactors=FALSE)
uniquegenes_brom <- bromicolauniquegeneid$V1
geneList_brom <- factor(as.integer(EbroNfe1$Name %in% uniquegenes_brom))
names(geneList_brom) <- EbroNfe1$Name
TFgenes_brom <- data.frame(TF = geneList_brom, Name = EbroNfe1$Name, row.names = NULL)
merged_brom <- merge(EbroNfe1, TFgenes_brom)
#no length

Ebroclosest_AT <- read.table("BromicolaATclosest.tsv", quote="\"", stringsAsFactors=FALSE)
colnames(Ebroclosest_AT) <- c("Name", "DistanceAT")
BromicolaLengths <- read.delim("BromicolaLengths.csv", header=FALSE, stringsAsFactors=FALSE)
colnames(BromicolaLengths) <- c("Name", "Length")

BromicolaAT <- merge(Ebroclosest_AT, merged_brom)
BromicolaFinalFrame <- merge(BromicolaAT, BromicolaLengths)
write.csv(BromicolaFinalFrame, file = "BromicolaData.csv")

Amarillansfinalframe$TF <- str_replace(Amarillansfinalframe$TF, "0", "FALSE")
Amarillansfinalframe$TF <- str_replace(Amarillansfinalframe$TF, "1", "TRUE")

BromicolaFinalFrame$TF <- str_replace(BromicolaFinalFrame$TF, "0", "FALSE")
BromicolaFinalFrame$TF <- str_replace(BromicolaFinalFrame$TF, "1", "TRUE")

#graphs
boxplot(Amarillansfinalframe$Length~Amarillansfinalframe$TF,
        data=Amarillansfinalframe,
        main="Boxplot of Lineage Specificity and Protein Length",
        ylab="Lineage Specific",
        xlab="Length",
        col= c("orange","red"),
        border="Black",
        horizontal = TRUE
)

boxplot(log(Amarillansfinalframe$TPM)~Amarillansfinalframe$TF,
        data=Amarillansfinalframe,
        main="Boxplot of Lineage Specificity and log(TPM)",
        ylab="Lineage Specific",
        xlab="log(TPM)",
        col= c("orange","red"),
        border="Black",
        horizontal = TRUE
)

boxplot(Amarillansfinalframe$DistanceAT~Amarillansfinalframe$TF,
        data=Amarillansfinalframe,
        main="Boxplot of Lineage Specificity and Distance to AT region",
        ylab="Lineage Specific",
        xlab="Distance to AT",
        col= c("orange","red"),
        border="Black",
        horizontal = TRUE
)

##bromicola
boxplot(BromicolaFinalFrame$DistanceAT~BromicolaFinalFrame$TF,
        data=Amarillansfinalframe,
        main="Boxplot of Bromicola Lineage Specificity and Distance to AT region",
        ylab="Lineage Specific",
        xlab="Distance to AT",
        col= c("orange","red"),
        border="Black",
        horizontal = TRUE
)

boxplot(BromicolaFinalFrame$Length~BromicolaFinalFrame$TF,
        data=Amarillansfinalframe,
        main="Boxplot of Bromicola Lineage Specificity and Protein Length",
        ylab="Lineage Specific",
        xlab="Length",
        col= c("orange","red"),
        border="Black",
        horizontal = TRUE
)

#stats
mod <- glm(TPM < 2 ~ TF, data=Amarillansfinalframe, family = "binomial")
summary(mod) #SIGNFICANT
mod <- glm(DistanceAT < 2000 ~ TF, data=Amarillansfinalframe, family = "binomial")
summary(mod) #SIGNIFICANT
mod <- glm(Length < 500 ~ TF, data=Amarillansfinalframe, family = "binomial")
summary(mod) #SIGNIFICANT

mod <- glm(Length < 500 ~ TF, data=BromicolaFinalFrame, family = "binomial")
summary(mod) #SIGNIFICANT
mod <- glm(DistanceAT < 2000 ~ TF, data=BromicolaFinalFrame, family = "binomial")
summary(mod) #SIGNIFICANT
