#First thing you need to do is set the working directory using the setwd() function. Put the filepath in quotations!!
#it should be in the same directory as youre in at the explorer window
setwd("/home/pierce/Desktop/Lab/Biostat/BLRNASEQ") #change this to your own working directory
source("https://bioconductor.org/biocLite.R") #load the bioconductor repository
biocLite("edgeR") #make sure you have the edgeR package installed
library(edgeR)
library("edgeR") #load the edgeR library
#Make sure you add an extra column at the end of the spreadsheet that contains your count values in .csv format
counts <- read.csv(file="PL_CK_vs_PL_flg22_30m.csv", header=TRUE, sep=",") #Then create the count matrix using read.csv()
cmat <- counts[ , -c(1,ncol(counts)) ] #format the data
rownames(cmat) <- counts[ , 1 ] # add gene names to new object
libsize <- colSums(cmat) # calculate library size in each sample
libmr <- libsize/1e06 #library size in millions of reads
cmat <- cmat[rowSums(cmat > 1) >= 3,] #keep only those rows where there are at least 10 counts in at least 3 samples
sample.description <- data.frame(colnames(cmat))
genotype=regmatches(colnames(cmat),regexpr("PL",colnames(cmat))) #define genotype vectors
treatment=regmatches(colnames(cmat),regexpr("CK|flg22",colnames(cmat))) #define treatment vectors
sample.description$genotype <- genotype #assign genotype vectors to sample.description dataframe
sample.description$treatment <- treatment #assign treatment vectors to sample.description dataframe
sample.description$group <- paste(sample.description$genotype,sample.description$treatment,sep="_")
sample.description$treatment <- as.factor(sample.description$treatment) #convert treatment data classes from characters to factors
sample.description$treatment <- relevel(sample.description$treatment,ref="CK") #set the reference treatment to "Ctrl"


#Perform full-model DE analysis for gene *AND* treatment 
dge <- DGEList(cmat, group=sample.description$group ) # Create object
dge <- calcNormFactors(dge, method='TMM') # Normalize library sizes using TMM
dge <- dge[rowSums(1e+06 * dge$counts/expandAsMatrix(dge$samples$lib.size, dim(dge)) > 1) >= 3, ] #keep only those genes that have at least 1 read per million in at least 3 samples.
design <- model.matrix(~treatment,data = sample.description)# Create design matrix for glm
rownames(design) <- sample.description$colnames.cmat.
dge <- estimateGLMCommonDisp(dge, design, verbose=TRUE) # Estimate common dispersion
dge <- estimateGLMTrendedDisp(dge,design) #Estimate Trended dispersion
dge <- estimateGLMTagwiseDisp(dge, design) # Estimate tagwise dispersion
fit <- glmFit(dge,design) # Fit general linear model
lrt <- glmLRT(fit, coef = "treatmentflg22") # Likelihood ratio test
out <- topTags(lrt, n=Inf) #create object which is organizing by decscending log ratio
keep <- out$table$FDR<= 0.05 & abs(out$table$logFC) >= 2 #eliminate all genes with a false discovery rate of .05 or lower, and whose |FC| is greater than or equal to 2
PL_CK_vs_PL_flg22_30m <- out[keep,]  #create matrix filled with results
plotMDS(dge, method = "bcv")  #make multidimensional scaling plot
plotBCV(dge) #plot biological coefficient of variation
decide <- decideTestsDGE(lrt, adjust.method = "fdr",p.value=0.05, lfc=2) #create summary of numbers of upregulated and downregulated genes
plotSmear(lrt, de.tags=decide) #create smear plot 
abline(h=c(-2,2), col="blue") #create blue cutofflines on the smear plot
write.csv(PL_CK_vs_PL_flg22_30m, file="PL_CK_vs_PL_flg22_30m_DE_lfc2.csv") #write the results to a .csv file
summary(decide)
  
