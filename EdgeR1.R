#First thing you need to do is set the working directory using the setwd() function. Put the filepath in quotations!!
#it should be in the same directory as youre in at the explorer window
setwd("/home/pierce/Desktop/Biostat/RNA Seq 8") #change this to your own working directory
source("https://bioconductor.org/biocLite.R") #load the bioconductor repository
biocLite("edgeR") #make sure you have the edgeR package installed
library(edgeR) #load the edgeR library
#Make sure you add an extra column at the end of the spreadsheet that contains your count values in .csv format
counts <- read.csv(file="aggie5.csv", header=TRUE, sep=",") #Then create the count matrix using read.csv()
cmat <- counts[ , -c(1,ncol(counts)) ] #format the data
rownames(cmat) <- counts[ , 1 ] # add gene names to new object
libsize <- colSums(cmat) # calculate library size in each sample
libmr <- libsize/1e06 #library size in millions of reads
cmat <- cmat[rowSums(cmat > 10) >= 3,] #keep only those rows where there are at least 10 counts in at least 3 samples
sample.description <- data.frame(sample=colnames(cmat))
genotype=regmatches(colnames(cmat),regexpr("WT|Aggie5",colnames(cmat))) #define genotypes for each column
treatment=regmatches(colnames(cmat),regexpr("Ctrl|Flg22",colnames(cmat))) #define treatment variable for each column

#time to build the object for edgeR
batches = c("WT", "WT", "WT", "WT", "WT", "WT", "OEaggie5", "OEaggie5", "OEaggie5", "OEaggie5", "OEaggie5", "OEaggie5") #Make batch vector
conditions = c("ctrl", "ctrl", "ctrl", "flg22", "flg22", "flg22", "ctrl", "ctrl", "ctrl", "flg22", "flg22", "flg22") #Make condition vector

#The pipeline could look like this:
dge <- DGEList(cmat, group = batches ) # Create object
dge <- calcNormFactors(dge, method='TMM') # Normalize library sizes using TMM
dge <- dge[rowSums(1e+06 * dge$counts/expandAsMatrix(dge$samples$lib.size, dim(dge)) > 1) >= 6, ] #keep only those genes that have at least 1 read per million in at least 6 samples.
design <- model.matrix(~0+conditions+batches) # Create design matrix for glm
dge <- estimateGLMCommonDisp(dge, design) # Estimate common dispersion
dge <- estimateGLMTrendedDisp(dge,design) #Estimate Trended dispersion
dge <- estimateGLMTagwiseDisp(dge, design) # Estimate tagwise dispersion
fit <- glmFit(dge,design) # Fit general linear model
pair_vector <- sprintf("%s-%s", "conditionsctrl", "conditionsflg22") # Conditions to be compared
pair_contrast <- makeContrasts(contrasts=pair_vector, levels=design) # Make contrast
lrt <- glmLRT(fit, contrast=pair_contrast) # Likelihood ratio test
out <- topTags(lrt, n=Inf) #create object which is organizing by decscending log ratio
keep <- out$table$FDR<= 0.05 & abs(out$table$logFC) >= 2 #eliminate all genes with a false discovery rate of .05 or lower, and whose |FC| is greater than or equal to 2
aggie5_flg22vsctrl <- out[keep,]  #create matrix filled with results
plotMDS(dge) #make multidimensional scaling plot
plotBCV(dge) #plot biological coefficient of variation
decide <- decideTestsDGE(lrt) #create summary of numbers of upregulated and downregulated genes
plotSmear(lrt, de.tags=decide) #create smear plot 
abline(h=c(-1,1), col="blue") #create blue cutofflines on the smear plot
write.csv(aggie5_flg22vsctrl, file="CL_aggie5_flg22vsctrl.csv") #write the results to a .csv file

#Create gene by treatment interaction model matrix
design.gbtinteraction <- model.matrix(~genotype*treatment, data = sample.description)
rownames(design.gbtinteraction) <- sample.description$sample
sample.description$group <- paste(sample.description$genotype,sample.description$treatment,sep="_") #paste the treatment and genotype columns together to get a group identifier



summary(decide)
