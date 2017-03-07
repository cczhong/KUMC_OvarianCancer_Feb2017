library(limma)
library(Biobase)

#read in necessary files:
df <- read.delim("<insert file path>")
pf<- read.delim("<insert file path>")

#set row names to genes
row.names(df) <- df$X
row.names(pf) <- pf$Sample_ID

#removes first column of characters (gene names)
df<- df[-1]
pf<- pf[-1]

#subsets phenotype_frame (pf) to get data with platinum status of resistant or sensitive and 
#of known ovarian cancer subtype
phenotype_matrix<- subset(pf, pf$Platinum.Status %in% c("Resistant", "Sensitive") & pf$Subtype %in% c("Fallopian", "Proliferative", "Mesenchymal", "Immunoreactive"))

#converts phenotype_matrix to data frame
phenotype_matrix<-data.frame(phenotype_matrix)

#converting platinum status to factor
phenotype_matrix$Platinum.Status <- factor(phenotype_matrix$Platinum.Status)

#subsetting expression data frame(df) to get only columns with relevant 
#platinum status 
samples <- names(df) %in% row.names(phenotype_matrix)
final_data_set <- df[samples]

#creating design matrix
design1 <- model.matrix(~ phenotype_matrix$Platinum.Status)


design2 <- model.matrix(~ phenotype_matrix$Platinum.Status - 1)
colnames(design2) <- c("Resistant", "Sensitive")

#generating linear models======================================================
#-----------------------
# Model 1: Simple Parametrization
#-----------------------

fit1 <- lmFit(final_data_set, design1)
fit1 <- eBayes(fit1)


#-----------------------
# Model 2: Two Parameter
#-----------------------

fit2 <- lmFit(final_data_set, design2)

#generating contrast matrix
contrast.matrix <- makeContrasts("Resistant-Sensitive", levels = design2)

fit2_with_cont <- contrasts.fit(fit2, contrast.matrix)
fit2_with_cont <- eBayes(fit2_with_cont) 


#==============================================================================
#generating output tables
model1_output <- data.frame(topTable(fit1, number=50))
model2_output <- data.frame(topTable(fit2_with_cont, number=50))

#writing output tables to files
write.table(model1_output, "<insertfilepath>/<insert desired filename>")
write.table(model2_output, "<insertfilepath>/<insert desired filename>")
