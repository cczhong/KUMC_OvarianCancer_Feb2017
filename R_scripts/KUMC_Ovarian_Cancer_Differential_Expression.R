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

#subsetting pf to sensitive only
s_fallopian_pm <- subset(pf, pf$Platinum.Status %in% c("Sensitive") & pf$Subtype %in% c("Fallopian"))
s_proliferative_pm <- subset(pf, pf$Platinum.Status %in% c("Sensitive") & pf$Subtype %in% c("Proliferative"))
s_mesenchymal_pm <- subset(pf, pf$Platinum.Status %in% c("Sensitive") & pf$Subtype %in% c("Mesenchymal"))
s_immunoreactive_pm <- subset(pf, pf$Platinum.Status %in% c("Sensitive") & pf$Subtype %in% c("Immunoreactive"))

#subsetting pf to resistant only
r_fallopian_pm <- subset(pf, pf$Platinum.Status %in% c("Resistant") & pf$Subtype %in% c("Fallopian"))
r_proliferative_pm <- subset(pf, pf$Platinum.Status %in% c("Resistant") & pf$Subtype %in% c("Proliferative"))
r_mesenchymal_pm <- subset(pf, pf$Platinum.Status %in% c("Resistant") & pf$Subtype %in% c("Mesenchymal"))
r_immunoreactive_pm <- subset(pf, pf$Platinum.Status %in% c("Resistant") & pf$Subtype %in% c("Immunoreactive"))

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

s_fallopian_list <- names(df) %in% row.names(s_fallopian_pm)
s_proliferative_list <- names(df) %in% row.names(s_proliferative_pm)
s_mesenchymal_list <- names(df) %in% row.names(s_mesenchymal_pm)
s_immunoreactive_list <- names(df) %in% row.names(s_immunoreactive_pm)

r_fallopian_list <- names(df) %in% row.names(r_fallopian_pm)
r_proliferative_list <- names(df) %in% row.names(r_proliferative_pm)
r_mesenchymal_list <- names(df) %in% row.names(r_mesenchymal_pm)
r_immunoreactive_list <- names(df) %in% row.names(r_immunoreactive_pm)

#subsetting df for sensitive platinum status
s_data_fallopian <- df[s_fallopian_list]
s_data_proliferative <- df[s_proliferative_list]
s_data_mesenchymal <- df[s_mesenchymal_list]
s_data_immunoreactive <- df[s_immunoreactive_list]

#subsetting df for resistant platinum status
r_data_fallopian <- df[r_fallopian_list]
r_data_proliferative <- df[r_proliferative_list]
r_data_mesenchymal <- df[r_mesenchymal_list]
r_data_immunoreactive <- df[r_immunoreactive_list]

#subsetting df based on cancper subtype
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


#removing adj-p values from output data
#fallopian_output$adj.P.Val <- NULL
#proliferative_output$adj.P.Val <- NULL
#mesenchymal_output$adj.P.Val <- NULL
#immunoreactive_output$adj.P.Val <- NULL

#removing B values from output data
fallopian_output$B <- NULL
proliferative_output$B <- NULL
mesenchymal_output$B <- NULL
immunoreactive_output$B<- NULL

#generating output sets with significant p values.
f_fallopian_output <- subset(fallopian_output, fallopian_output$P.Value <= 0.05)
f_immunoreactive_output <- subset(immunoreactive_output, immunoreactive_output$P.Value <= 0.05)
f_proliferative_output <- subset(proliferative_output, proliferative_output$P.Value <= 0.05)
f_mesenchymal_output <- subset(mesenchymal_output, mesenchymal_output$P.Value <= 0.05)

#generate Gene Sets with significant p values
fallopian_gene_set <- row.names(f_fallopian_output)
proliferative_gene_set <- row.names(f_proliferative_output)
mesenchymal_gene_set <- row.names(f_mesenchymal_output)
immunoreactive_gene_set <- row.names(f_immunoreactive_output)

#writing full output to files
write.table(fallopian_output, "C:/Users/david/Documents/GitHub/KUMC_OvarianCancer_Feb2017/Outputs/full_fallopian_output.txt", sep = '\t')
write.table(proliferative_output, "C:/Users/david/Documents/GitHub/KUMC_OvarianCancer_Feb2017/Outputs/full_proliferative_output.txt", sep = '\t')
write.table(mesenchymal_output, "C:/Users/david/Documents/GitHub/KUMC_OvarianCancer_Feb2017/Outputs/full_mesenchymal_output.txt", sep = '\t')
write.table(immunoreactive_output, "C:/Users/david/Documents/GitHub/KUMC_OvarianCancer_Feb2017/Outputs/full_immunoreactive_output.txt", sep = '\t')

#writing output with significant p values
write.table(f_fallopian_output, "C:/Users/david/Documents/GitHub/KUMC_OvarianCancer_Feb2017/Outputs/final_fallopian_output.txt", sep = '\t')
write.table(f_proliferative_output, "C:/Users/david/Documents/GitHub/KUMC_OvarianCancer_Feb2017/Outputs/final_proliferative_output.txt", sep='\t')
write.table(f_mesenchymal_output, "C:/Users/david/Documents/GitHub/KUMC_OvarianCancer_Feb2017/Outputs/final_mesenchymal_output.txt", sep='\t')
write.table(f_immunoreactive_output, "C:/Users/david/Documents/GitHub/KUMC_OvarianCancer_Feb2017/Outputs/final_immunoreactive_output.txt", sep='\t')

#writing significant gene set output
write.table(fallopian_gene_set, "C:/Users/david/Documents/GitHub/KUMC_OvarianCancer_Feb2017/Outputs/fallopian_gene_set.txt", row.names = FALSE)
write.table(proliferative_gene_set, "C:/Users/david/Documents/GitHub/KUMC_OvarianCancer_Feb2017/Outputs/proliferative_gene_set.txt", row.names = FALSE)
write.table(mesenchymal_gene_set, "C:/Users/david/Documents/GitHub/KUMC_OvarianCancer_Feb2017/Outputs/mesenchymal_gene_set.txt", row.names = FALSE)
write.table(immunoreactive_gene_set, "C:/Users/david/Documents/GitHub/KUMC_OvarianCancer_Feb2017/Outputs/immunoreactive_gene_set.txt", row.names = FALSE)


#creating data-frames for boxplots
#f = fallopian, imm = immuno, prolif = proliferative, mesen=mesenchymal

names <- data.frame(row.names(df))
cols <- c("Sample.ID", "Platinum.Status")
names(names) <- Gene.Name

for(x in names$Gene.Name)
{
  cols<- append(cols,x)
}

f <-read.delim("C:/Users/david/Documents/Github/KUMC_OvarianCancer_Feb2017/R_data_inputs/fallopian.txt", header=FALSE)
imm <- read.delim("C:/Users/david/Documents/GitHub/KUMC_OvarianCancer_Feb2017/R_data_inputs/immunoreactive.txt", header=FALSE)
mesen <- read.delim("C:/Users/david/Documents/GitHub/KUMC_OvarianCancer_Feb2017/R_data_inputs/mesenchymal.txt", header=FALSE)
prolif <- read.delim("C:/Users/david/Documents/GitHub/KUMC_OvarianCancer_Feb2017/R_data_inputs/proliferative.txt", header=FALSE)

colnames(f) <- cols
colnames(imm) <- cols
colnames(mesen) <- cols
colnames(prolif) <- cols

#Boxplots to verify results -- top 5 genes for each subtype
boxplot(f$C10ORF95~f$Platinum.Status, col=c("red","blue"), horizontal=TRUE, xlab="Expression Level", ylab="Platinum Status", main="C10RF95 Expression in Fallopian Subtypes")
boxplot(f$EPYC~f$Platinum.Status, col=c("red","blue"), horizontal=TRUE, xlab="Expression Level", ylab="Platinum Status", main= "EPYC Expression in Fallopian Subtype")
boxplot(f$CX3CL1~f$Platinum.Status, col=c("red","blue"), horizontal=TRUE, xlab="Expression Level", ylab="Platinum Status", main= "CX3CL1 Expression in Fallopian Subtype")
boxplot(f$XPNPEP1~f$Platinum.Status, col=c("red","blue"), horizontal=TRUE, xlab="Expression Level", ylab="Platinum Status", main= "XPNPEP1 Expression in Fallopian Subtype")
boxplot(f$SCEL~f$Platinum.Status, col=c("red","blue"), horizontal=TRUE, xlab="Expression Level", ylab="Platinum Status", main= "SCEL Expression in Fallopian Subtype")

boxplot(imm$TRAM1~imm$Platinum.Status, col=c("red","blue"), horizontal=TRUE, xlab="Expression Level", ylab= "Platinum Status", main= "TRAM1 Expression in Immunoreactive Subtype")
boxplot(imm$GBL~imm$Platinum.Status, col=c("red","blue"), horizontal=TRUE, xlab="Expression Level", ylab= "Platinum Status", main= "GBL Expression in Immunoreactive Subtype")
boxplot(imm$PQBP1~imm$Platinum.Status, col=c("red","blue"), horizontal=TRUE, xlab="Expression Level", ylab= "Platinum Status", main= "PQBP1 Expression in Immunoreactive Subtype")
boxplot(imm$OGFRL1~imm$Platinum.Status, col=c("red","blue"), horizontal=TRUE, xlab="Expression Level", ylab= "Platinum Status", main= "OGFRL1 Expression in Immunoreactive Subtype")
boxplot(imm$HDAC6~imm$Platinum.Status, col=c("red","blue"), horizontal=TRUE, xlab="Expression Level", ylab= "Platinum Status", main= "HDAC6 Expression in Immunoreactive Subtype")

boxplot(mesen$C16ORF5~mesen$Platinum.Status, col=c("red", "blue"), horizontal=TRUE, xlab="Expression Level", ylab="Platinum Status", main= "C16ORF5 Expression in Mesenchymal Subtype")
boxplot(mesen$MAD1L1~mesen$Platinum.Status, col=c("red", "blue"), horizontal=TRUE, xlab="Expression Level", ylab="Platinum Status", main= "MAD1L1 Expression in Mesenchymal Subtype")
boxplot(mesen$NLGN1~mesen$Platinum.Status, col=c("red", "blue"), horizontal=TRUE, xlab="Expression Level", ylab="Platinum Status", main= "NLGN1 Expression in Mesenchymal Subtype")
boxplot(mesen$SERGEF~mesen$Platinum.Status, col=c("red", "blue"), horizontal=TRUE, xlab="Expression Level", ylab="Platinum Status", main= "SERGEF Expression in Mesenchymal Subtype")
boxplot(mesen$ARMC9~mesen$Platinum.Status, col=c("red", "blue"), horizontal=TRUE, xlab="Expression Level", ylab="Platinum Status", main= "ARMC9 Expression in Mesenchymal Subtype")

boxplot(prolif$MXRA5~prolif$Platinum.Status, col=c("red", "blue"), horizontal=TRUE, xlab="Expression Level", ylab="Platinum Status", main="MXRA5 Expression Level in Proliferative Subtype")
boxplot(prolif$DKFZP566H0824~prolif$Platinum.Status, col=c("red", "blue"), horizontal=TRUE, xlab="Expression Level", ylab="Platinum Status", main="DKFZP566H0824 Expression Level in Proliferative Subtype")
boxplot(prolif$PRKCDBP~prolif$Platinum.Status, col=c("red", "blue"), horizontal=TRUE, xlab="Expression Level", ylab="Platinum Status", main="PRKCDBP Expression Level in Proliferative Subtype")
boxplot(prolif$MMP1~prolif$Platinum.Status, col=c("red", "blue"), horizontal=TRUE, xlab="Expression Level", ylab="Platinum Status", main="MMP1 Expression Level in Proliferative Subtype")
boxplot(prolif$ADAMTS1~prolif$Platinum.Status, col=c("red", "blue"), horizontal=TRUE, xlab="Expression Level", ylab="Platinum Status", main="ADAMTS1 Expression Level in Proliferative Subtype")
