library(limma)
library(Biobase)

#read in necessary files:
df <- read.delim("C:/Users/david/workspace/Ovarian_cancer/TCGA_489_UE.txt")
pf<- read.delim("C:/Users/david/workspace/Ovarian_cancer/phenotype_matrix.txt")

#set row names to genes
row.names(df) <- df$X
row.names(pf) <- pf$Sample_ID

#removes first column of characters (gene names)
df<- df[-1]
pf<- pf[-1]

#subsets phenotype_frame (pf) to get data with platinum status of resistant or sensitive and 
#of known ovarian cancer subtype
fallopian_pm <- subset(pf, pf$Platinum.Status %in% c("Resistant", "Sensitive") & pf$Subtype %in% c("Fallopian"))
proliferative_pm <- subset(pf, pf$Platinum.Status %in% c("Resistant", "Sensitive") & pf$Subtype %in% c("Proliferative"))
mesenchymal_pm <- subset(pf, pf$Platinum.Status %in% c("Resistant", "Sensitive") & pf$Subtype %in% c("Mesenchymal"))
immunoreactive_pm <- subset(pf, pf$Platinum.Status %in% c("Resistant", "Sensitive") & pf$Subtype %in% c("Immunoreactive"))

#converts phenotype_matrix to data frame
#phenotype_matrix<-data.frame(phenotype_matrix)
fallopian_pm <- data.frame(fallopian_pm)
proliferative_pm <- data.frame(proliferative_pm)
mesenchymal_pm <- data.frame(mesenchymal_pm)
immunoreactive_pm <- data.frame(immunoreactive_pm)


#converting platinum status to factor
fallopian_pm$Platinum.Status <- factor(fallopian_pm$Platinum.Status)
proliferative_pm$Platinum.Status <- factor(proliferative_pm$Platinum.Status)
mesenchymal_pm$Platinum.Status <- factor(mesenchymal_pm$Platinum.Status)
immunoreactive_pm$Platinum.Status <- factor(immunoreactive_pm$Platinum.Status)

#subsetting expression data frame(df) to get only columns with relevant 
#platinum status 
fallopian_list <- names(df) %in% row.names(fallopian_pm)
proliferative_list <-names(df) %in% row.names(proliferative_pm)
mesenchymal_list <- names(df) %in% row.names(mesenchymal_pm)
immunoreactive_list <- names(df) %in% row.names(immunoreactive_pm)

data_fallopian <- df[fallopian_list]
data_proliferative <- df[proliferative_list]
data_mesenchymal <- df[mesenchymal_list]
data_immunoreactive<- df[immunoreactive_list]

#creating design matrix
design_fallopian <- model.matrix(~ fallopian_pm$Platinum.Status - 1)
colnames(design_fallopian) <- c("Resistant", "Sensitive")

design_proliferative <- model.matrix(~ proliferative_pm$Platinum.Status - 1)
colnames(design_proliferative) <- c("Resistant", "Sensitive")

design_mesenchymal <- model.matrix(~ mesenchymal_pm$Platinum.Status - 1)
colnames(design_mesenchymal) <- c("Resistant", "Sensitive")

design_immunoreactive <- model.matrix(~ immunoreactive_pm$Platinum.Status - 1)
colnames(design_immunoreactive) <- c("Resistant", "Sensitive")


#generating linear models
fallopian_fit <- lmFit(data_fallopian, design_fallopian)
proliferative_fit <- lmFit(data_proliferative, design_proliferative)
mesenchymal_fit <- lmFit(data_mesenchymal, design_mesenchymal)
immunoreactive_fit <- lmFit(data_immunoreactive, design_immunoreactive)

fallopian_contrast <- makeContrasts("Resistant-Sensitive", levels = design_fallopian)
proliferative_contrast <- makeContrasts("Resistant-Sensitive", levels = design_proliferative)
mesenchymal_contrast <- makeContrasts("Resistant-Sensitive", levels = design_mesenchymal)
immunoreactive_contrast <- makeContrasts("Resistant-Sensitive", levels = design_immunoreactive)

final_fallopian_fit <- contrasts.fit(fallopian_fit, fallopian_contrast)
final_proliferative_fit <- contrasts.fit(proliferative_fit, proliferative_contrast)
final_mesenchymal_fit <- contrasts.fit(mesenchymal_fit, mesenchymal_contrast)
final_immunoreactive_fit <- contrasts.fit(immunoreactive_fit, immunoreactive_contrast)

final_fallopian_fit <- eBayes(final_fallopian_fit)
final_proliferative_fit <- eBayes(final_proliferative_fit)
final_mesenchymal_fit <- eBayes(final_mesenchymal_fit)
final_immunoreactive_fit <- eBayes(final_immunoreactive_fit)

#Generating outputs
fallopian_output <- data.frame(topTable(final_fallopian_fit, number=11865))
proliferative_output <- data.frame(topTable(final_proliferative_fit, number=11865))
mesenchymal_output <- data.frame(topTable(final_mesenchymal_fit, number=11865))
immunoreactive_output <- data.frame(topTable(final_immunoreactive_fit, number=11865))

#appending Gene Names
gene_fallopian <- data.frame(row.names(fallopian_output))
gene_proliferative <- data.frame(row.names(proliferative_output))
gene_mesenchymal <- data.frame(row.names(mesenchymal_output))
gene_immunoreactive <- data.frame(row.names(immunoreactive_output))

fallopian_output[,"gene name"] <- gene_fallopian
proliferative_output[,"gene name"] <- gene_proliferative
mesenchymal_output[,"gene name"] <- gene_mesenchymal
immunoreactive_output[,"gene name"] <- gene_immunoreactive

#removing adj-p values from output data
fallopian_output$adj.P.Val <- NULL
proliferative_output$adj.P.Val <- NULL
mesenchymal_output$adj.P.Val <- NULL
immunoreactive_output$adj.P.Val <- NULL

#cleaning workspace
rm(gene_fallopian)
rm(gene_proliferative)
rm(gene_mesenchymal)
rm(gene_immunoreactive)

topTable(final_fallopian_fit)
