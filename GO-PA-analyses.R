library(tidyverse)
library(tximport)
library(edgeR)
source("/Users/sunshine2.0/Desktop/WORK/ColumnCombFunction-1.R")
library(pheatmap)

### RUN DE SCRIPT FIRST ***

# Assign unnormalized quantification
unnorm<-salmon_quant$counts #<-from tximport

## NORMALIZE WITH RELATIVE LOG EXPRESSSION METHOD

#Create a new DGEList object and call it RLE
#Calculate the normalization factors. There should now be 2 total: unnorm, RLE

RLE<-unnorm

# Use April Female with the slimmest emus as the reference because it corresponds to the null condition. 
RLE<-DGEList(counts=RLE,samples=samples)
RLE<-calcNormFactors(RLE,method="RLE",refColumn=5)

# Take a look at the normalization factors of RLE They should be different from 1
print(RLE$samples$norm.factors)

# Convert unnormalized list to DGElist
unnorm<-DGEList(counts=unnorm,samples=samples)

#import gene length------------------------------------------------
library("EDASeq")
head(geneslist$gene)
geneslist<-data.frame("gene"=rownames(salmon_quant$counts), stringsAsFactors = FALSE)
geneslist$gene = substr(geneslist$gene,1,nchar(geneslist$gene)-2) # remove .n from gene ID
genelength<-getGeneLengthAndGCContent(geneslist$gene, "dno")  # ****choose mart 1 for emu*******
#genelength<-genelength # export genes copy and keep copy because this took forever to download
# write.csv(genelength, "/Users/sunshine2.0/Desktop/WORK/output/genelength.csv")
genelength<-as.data.frame(genelength)
#------------------------------------------------------------------ OR below code if already done this chunk

# Import saved gene length file from above code 
genelength<-read.table("/Users/sunshine2.0/Desktop/WORK/output-without pvalues on GO PA/genelength.csv", header=TRUE, sep=",")

## NORMALIZE AS TRANSCRIPTS PER MILLION FOR GENE LENGTH
#TPM calculation (https://gist.github.com/slowkow/6e34ccb4d1311b8fe62e)
tpm <- function(counts, lengths) {
  rate <- counts / lengths
  rate / sum(rate) * 1e6
}
tpms <- apply((RLE$counts+0.25), 2, function(x) tpm(x, genelength$length)) #prior count of 0.25
log2tpm <- log2(tpms)
log2tpm<-as.data.frame(log2tpm)
head(log2tpm)

# CREATE PCA PLOTS

#Before performing the PCA, we must determine what 500 genes have the highest
#between sample variation. Note that in plotMDS(), the prior.count added to 
#each read count before computing the log2 CPM is 2, so we also add that 
#prior.count here. I do not know why they chose 2, but they did!
# Log2cpmRLE<-cpm(RLE,normalized.lib.sizes=TRUE,log=TRUE,prior.count=2)

#Then, for each gene (row), compute the variance in the log2 CPM across all
#samples.
geneVariance<-rowMeans((log2tpm-rowMeans(log2tpm))^2)

#Then, order the genes in decreasing order of their variance and take the 
#first 500 genes after ordering.
top500VarGenes<-order(geneVariance,decreasing=TRUE)[1:500]

#The function prcomp() performs a PCA on a matrix. Note that because of 
#the way the function is written in R, the principle components are computed 
#for each row of the matrix rather than each column. Therefore, to perform a 
#PCA on the samples, rather than the genes, the log 2 CPM matrix needs to be 
#transposed using the t() function. The arguments center=TRUE and scale=FALSE
#mean that we want the log2 CPM values centered so that the mean log 2 CPM 
#for each gene across all samples is zero, but not scaled.
pcaRLE<-prcomp(t(log2tpm[top500VarGenes,]),center=TRUE,scale=FALSE)

#Note that the sum of the squares of each PC across all the samples gives the 
#variance explained by that PC, and so the proportion of variance 
#explained by each PC can be found by taking that sum and dividing by the total
#variance. We do that below:
totalVariance<-sum(colSums(pcaRLE$x**2))
pcVariance<-colSums(pcaRLE$x**2)
propVarPC<-pcVariance/totalVariance

round(propVarPC,digits=3)

#The output of prcomp() is a list. The element of the list called "x" gives
#the principle components as magnitudes. Each row corresponds to a sample and
#each column to a principle component. The first two principle components (PCs),
#namely those that explain the most variance, are the first two columns. 
#We plot PC 1 on the x axis and PC 2 on the y axis below, labeling the points
#with the samples they correspond to. coloured by month.

# PCA coloured by month

png("/Users/sunshine2.0/Desktop/WORK/output/emu_PCA_month.png", width=9,height=9,units="in",res=1200)

plot(pcaRLE$x[,1],pcaRLE$x[,2],type="n",
     xlab=paste0("PC 1 (",round(propVarPC,digits=3)[[1]]*100,"%)"), ylab=paste0("PC 2 (",round(propVarPC,digits=3)[[2]]*100,"%)"),
     xlim=c(-20,51),ylim=c(-32,20))

text(pcaRLE$x[,1],pcaRLE$x[,2],labels=rownames(pcaRLE$x),col=c(rep('red',8),rep('#009300',8),rep('blue',8)))

legend("bottomright", inset=.02,
       c("April","June","Nov"), fill=c("red","#009300","blue"), horiz=TRUE, cex=0.8)

dev.off()

# PCA coloured by gender

png("/Users/sunshine2.0/Desktop/WORK/output/emu_PCA_gender.png",width=9,height=9,units="in",res=1200)

plot(pcaRLE$x[,1],pcaRLE$x[,2],type="n",
     xlab=paste0("PC 1 (",round(propVarPC,digits=3)[[1]]*100,"%)"), ylab=paste0("PC 2 (",round(propVarPC,digits=3)[[2]]*100,"%)"),
     xlim=c(-20,51),ylim=c(-32,20))

text(pcaRLE$x[,1],pcaRLE$x[,2],labels=rownames(pcaRLE$x),col=c(rep('#660033',4),rep('purple',4)))

legend("bottomright", inset=.02,
       c("Male","Female"), fill=c("#660033","purple"), horiz=TRUE, cex=0.8)

dev.off()


# Edit metadata-----------------------------------
samples<-samples %>%
  mutate("grp"= paste0(month, ".", gender)) 
row.names(samples) <- samples$sample
samples[1] <- NULL
rownames(samples)==colnames(topDat)


## HEATMAP HIERARCHICAL GENE CLUSTERING 

# Take list of DE genes less than adjusted p-value cut off of 0.05 and 1 < LFC < -1 for clustering
length(unique(c(DEgenelist1, DEgenelist2, DEgenelist3, DEgenelist4, DEgenelist5, DEgenelist6, DEgenelist7, DEgenelist8, DEgenelist9)))
DElist0.05<-unique(c(DEgenelist1, DEgenelist2, DEgenelist3, DEgenelist4, DEgenelist5, DEgenelist6, DEgenelist7, DEgenelist8, DEgenelist9))

# Scaled data of topGenes
topDat <- log2tpm[rownames(log2tpm)%in% DElist0.05, ]
nrow(topDat )

geneC.dis <- dist(topDat, method = "euclidean")

geneC.hc.a <- hclust(geneC.dis, method = "average")

plot(geneC.hc.a, labels = FALSE, main = "Hierarchical with Average Linkage", xlab = "")

# Heatmap padj < 0.05

save_pheatmap_pdf <- function(x, filename, width=7, height=7) {
  stopifnot(!missing(x))
  stopifnot(!missing(filename))
  pdf(filename, width=width, height=height)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}

pheatmap0.05<-pheatmap(topDat, cluster_rows = TRUE, scale = "none", clustering_method = "average", 
         clustering_distance_cols = "euclidean", clustering_distance_rows = "euclidean", 
         annotation = samples[, c("month", "gender","grp")], show_rownames = FALSE)

save_pheatmap_pdf(pheatmap0.05, "/Users/sunshine2.0/Desktop/WORK/output/heatmap0.05.pdf")


# Take list of DE genes less than adjusted p-value cut off of 0.01 and 1 < LFC < -1  for clustering

length(unique(c(DEgenelist1.01, DEgenelist2.01, DEgenelist3.01, DEgenelist4.01, DEgenelist5.01, DEgenelist6.01, DEgenelist7.01, DEgenelist8.01, DEgenelist9.01)))
DElist0.01<-unique(c(DEgenelist1.01, DEgenelist2.01, DEgenelist3.01, DEgenelist4.01, DEgenelist5.01, DEgenelist6.01, DEgenelist7.01, DEgenelist8.01, DEgenelist9.01))

# Scaled data of topGenes 
topDat <- log2tpm[rownames(log2tpm)%in% DElist0.01, ]
nrow(topDat )

geneC.dis <- dist(topDat, method = "euclidean")

geneC.hc.a <- hclust(geneC.dis, method = "average")

png("/Users/sunshine2.0/Desktop/WORK/output/hierarchical_clustering0.01.png",width=10,height=7,units="in",res=1200)

plot(geneC.hc.a, labels = FALSE, main = "Hierarchical with Average Linkage", xlab = "")

dev.off()

# Heatmap padj < 0.01

pheatmap0.01<-pheatmap(topDat, cluster_rows = TRUE, scale = "none", clustering_method = "average", 
         clustering_distance_cols = "euclidean", clustering_distance_rows = "euclidean", 
         annotation = samples[, c("month", "gender","grp")], show_rownames = FALSE)

save_pheatmap_pdf(pheatmap0.01, "/Users/sunshine2.0/Desktop/WORK/output/heatmap0.01.pdf")


# Take list of DE genes less than adjusted p-value cut off of 0.01 and 1 < LFC  for clustering (upregulated)

length(unique(c(DEgenelist1.up01, DEgenelist2.up01, DEgenelist3.up01, DEgenelist4.up01, DEgenelist5.up01, DEgenelist6.up01, DEgenelist7.up01, DEgenelist8.up01, DEgenelist9.up01)))
DElist0.up01<-unique(c(DEgenelist1.up01, DEgenelist2.up01, DEgenelist3.up01, DEgenelist4.up01, DEgenelist5.up01, DEgenelist6.up01, DEgenelist7.up01, DEgenelist8.up01, DEgenelist9.up01))

# Scaled data of topGenes 
topDat <- log2tpm[rownames(log2tpm)%in% DElist0.up01, ]
nrow(topDat )

# Heatmap padj < 0.01, upregulated

pheatmap0.up01<-pheatmap(topDat, cluster_rows = TRUE, scale = "none", clustering_method = "average", 
                       clustering_distance_cols = "euclidean", clustering_distance_rows = "euclidean", 
                       annotation = samples[, c("month", "gender","grp")], show_rownames = FALSE)

save_pheatmap_pdf(pheatmap0.up01, "/Users/sunshine2.0/Desktop/WORK/output/heatmap0.01upregulated.pdf")


# Take list of DE genes less than adjusted p-value cut off of 0.01 and LFC < -1  for clustering (downregulated)

length(unique(c(DEgenelist1.down01, DEgenelist2.down01, DEgenelist3.down01, DEgenelist4.down01, DEgenelist5.down01, DEgenelist6.down01, DEgenelist7.down01, DEgenelist8.down01, DEgenelist9.down01)))
DElist0.down01<-unique(c(DEgenelist1.down01, DEgenelist2.down01, DEgenelist3.down01, DEgenelist4.down01, DEgenelist5.down01, DEgenelist6.down01, DEgenelist7.down01, DEgenelist8.down01, DEgenelist9.down01))

# Scaled data of topGenes 
topDat <- log2tpm[rownames(log2tpm)%in% DElist0.down01, ]
nrow(topDat )

# Heatmap padj < 0.01, downregulated

pheatmap0.down01<-pheatmap(topDat, cluster_rows = TRUE, scale = "none", clustering_method = "average", 
                       clustering_distance_cols = "euclidean", clustering_distance_rows = "euclidean", 
                       annotation = samples[, c("month", "gender","grp")], show_rownames = FALSE)

save_pheatmap_pdf(pheatmap0.down01, "/Users/sunshine2.0/Desktop/WORK/output/heatmap0.01downregulated.pdf")


# Take list of overlap of both upregulated and downregulated genes padj < 0.01

length(DElist0.up01) +length(DElist0.down01)-length(DElist0.01)
up0.01<-as.data.frame(DElist0.up01, stringsAsFactors=FALSE)  %>%
  rename( DElist0.up01= "x")
down0.01<-as.data.frame("x"=DElist0.down01, stringsAsFactors=FALSE) %>%
  rename( DElist0.down01= "x")
overlap0.01<-up0.01 %>%
  inner_join(down0.01, by= "x")

# Scaled data of topGenes 
topDat <- log2tpm[rownames(log2tpm)%in% overlap0.01$x, ]
nrow(topDat )

# Heatmap overlap of both upregulated and downregulated genes padj < 0.01

pheatmap0.overlap01<-pheatmap(topDat, cluster_rows = TRUE, scale = "none", clustering_method = "average", 
                           clustering_distance_cols = "euclidean", clustering_distance_rows = "euclidean", 
                           annotation = samples[, c("month", "gender","grp")], show_rownames = FALSE)

save_pheatmap_pdf(pheatmap0.overlap01, "/Users/sunshine2.0/Desktop/WORK/output/heatmap0.01overlap.pdf")

# Male and female both down and up regulated genes list for GO bar plots---------------------------------------********
DElist0.down01F<-unique(c( DEgenelist4.down01, DEgenelist5.down01, DEgenelist6.down01))
DElist0.down01M<-unique(c( DEgenelist7.down01, DEgenelist8.down01, DEgenelist9.down01))
DElist0.up01F<-unique(c( DEgenelist4.up01, DEgenelist5.up01, DEgenelist6.up01))
DElist0.up01M<-unique(c( DEgenelist7.up01, DEgenelist8.up01, DEgenelist9.up01))

# Female overlap gene list
up0.01F<-as.data.frame(DElist0.up01F, stringsAsFactors=FALSE)  %>%
  rename( DElist0.up01F = "x")
down0.01F<-as.data.frame(DElist0.down01F, stringsAsFactors=FALSE) %>%
  rename( DElist0.down01F = "x")
overlap0.01F<-up0.01F %>%
  inner_join(down0.01F, by= "x")
nrow(overlap0.01F)

# Male overlap gene list
up0.01M<-as.data.frame(DElist0.up01M, stringsAsFactors=FALSE)  %>%
  rename( DElist0.up01M= "x")
down0.01M<-as.data.frame("x"=DElist0.down01M, stringsAsFactors=FALSE) %>%
  rename( DElist0.down01M = "x")
overlap0.01M<-up0.01M %>%
  inner_join(down0.01M, by= "x")
nrow(overlap0.01M)

# Overlap by month
#J vs A female
JvsA_female_overlap_genes<-as.data.frame(DEgenelist4.01, stringsAsFactors=FALSE)  %>%
  rename( DEgenelist4.01 = "x") %>%
  filter(x %in% overlap0.01F$x) #overlap subset
#N vs J female
NvsJ_female_overlap_genes<-as.data.frame(DEgenelist5.01, stringsAsFactors=FALSE)  %>%
  rename(DEgenelist5.01 = "x") %>%
  filter(x %in% overlap0.01F$x) #overlap subset
#A vs N female
AvsN_female_overlap_genes<-as.data.frame(DEgenelist6.01, stringsAsFactors=FALSE)  %>%
  rename(DEgenelist6.01 = "x") %>%
  filter(x %in% overlap0.01F$x) #overlap subset

#J vs A male
JvsA_male_overlap_genes<-as.data.frame(DEgenelist7.01, stringsAsFactors=FALSE)  %>%
  rename(DEgenelist7.01  = "x") %>%
  filter(x %in% overlap0.01M$x) #overlap subset
#N vs J male
NvsJ_male_overlap_genes<-as.data.frame(DEgenelist8.01, stringsAsFactors=FALSE)  %>%
  rename(DEgenelist8.01  = "x") %>%
  filter(x %in% overlap0.01M$x) #overlap subset
#A vs N male
AvsN_male_overlap_genes<-as.data.frame(DEgenelist9.01, stringsAsFactors=FALSE)  %>%
  rename( DEgenelist9.01  = "x") %>%
  filter(x %in% overlap0.01M$x) #overlap subset

# Save gene lists
write.csv(overlap0.01F, "/Users/sunshine2.0/Desktop/WORK/output/overlap_genelist_0.01_F.csv")
write.csv(overlap0.01M, "/Users/sunshine2.0/Desktop/WORK/output/overlap_genelist_0.01_M.csv")

write.csv(JvsA_female_overlap_genes, "/Users/sunshine2.0/Desktop/WORK/output/JvsA_female_overlap_genes_0.01.csv")
write.csv(NvsJ_female_overlap_genes, "/Users/sunshine2.0/Desktop/WORK/output/NvsJ_female_overlap_genes_0.01.csv")
write.csv(AvsN_female_overlap_genes, "/Users/sunshine2.0/Desktop/WORK/output/AvsN_female_overlap_genes_0.01.csv")
write.csv(JvsA_male_overlap_genes, "/Users/sunshine2.0/Desktop/WORK/output/JvsA_male_overlap_genes_0.01.csv")
write.csv(NvsJ_male_overlap_genes, "/Users/sunshine2.0/Desktop/WORK/output/NvsJ_male_overlap_genes_0.01.csv")
write.csv(AvsN_male_overlap_genes, "/Users/sunshine2.0/Desktop/WORK/output/AvsN_male_overlap_genes_0.01.csv")

write.csv(DEgenelist4.01, "/Users/sunshine2.0/Desktop/WORK/output/JvsA_female_DE_genes_0.01.csv")
write.csv(DEgenelist5.01, "/Users/sunshine2.0/Desktop/WORK/output/NvsJ_female_DE_genes_0.01.csv")
write.csv(DEgenelist6.01, "/Users/sunshine2.0/Desktop/WORK/output/AvsN_female_DE_genes_0.01.csv")
write.csv(DEgenelist7.01, "/Users/sunshine2.0/Desktop/WORK/output/JvsA_male_DE_genes_0.01.csv")
write.csv(DEgenelist8.01, "/Users/sunshine2.0/Desktop/WORK/output/NvsJ_male_DE_genes_0.01.csv")
write.csv(DEgenelist9.01, "/Users/sunshine2.0/Desktop/WORK/output/AvsN_male_DE_genes_0.01.csv")

# Remove .n from gene ID
overlap0.01F$x = substr(overlap0.01F$x,1,nchar(overlap0.01F$x)-2) 
overlap0.01M$x = substr(overlap0.01M$x,1,nchar(overlap0.01M$x)-2)

# Obtain external gene names and other ID types
require('biomaRt')

# mart <- useMart('ENSEMBL_MART_ENSEMBL')
# mart <- useDataset('dnovaehollandiae_gene_ensembl', mart)
mart <- useEnsembl(biomart = "ensembl", dataset = "dnovaehollandiae_gene_ensembl", version = "102")

annotLookupF <-getBM(
  mart = mart,
  attributes = c(
    'ensembl_gene_id',
    'gene_biotype',
    'external_gene_name',
    'uniprot_gn_symbol',
    'uniprot_gn_id',
    'entrezgene_id'),
  filter = 'ensembl_gene_id',
  values = overlap0.01F$x,
  uniqueRows = TRUE)
write.csv(annotLookupF$external_gene_name, "/Users/sunshine2.0/Desktop/WORK/output/annotLookup$external_gene_name_overlap0.01F.csv")

annotLookupM <-getBM(
  mart = mart,
  attributes = c(
    'ensembl_gene_id',
    'gene_biotype',
    'external_gene_name',
    'uniprot_gn_symbol',
    'uniprot_gn_id',
    'entrezgene_id'),
  filter = 'ensembl_gene_id',
  values = overlap0.01M$x,
  uniqueRows = TRUE)
write.csv(annotLookupM$external_gene_name, "/Users/sunshine2.0/Desktop/WORK/output/annotLookup$external_gene_name_overlap0.01M.csv")

# Create dataframe of up and down regulated genes for the 6 pair-wise comparisons for each sex
ensembl_gene_id<-c(DEgenelist4.up01, DEgenelist4.down01, DEgenelist5.up01, DEgenelist5.down01, DEgenelist6.up01, DEgenelist6.down01)
regulated<-c(rep("up", length(DEgenelist4.up01)), rep("down", length(DEgenelist4.down01)), rep("up", length(DEgenelist5.up01)), rep("down", length(DEgenelist5.down01)), rep("up", length(DEgenelist6.up01)), rep("down", length(DEgenelist6.down01)))
month_comparison<-c(rep("June vs. April", length(DEgenelist4.up01)+ length(DEgenelist4.down01)), rep("November vs. June", length(DEgenelist5.up01)+ length(DEgenelist5.down01)), rep("April vs. November", length(DEgenelist6.up01)+ length(DEgenelist6.down01)))
Female_DE_genes.1<-data.frame(ensembl_gene_id, regulated, month_comparison, stringsAsFactors = FALSE)

ensembl_gene_id<-c(DEgenelist7.up01, DEgenelist7.down01, DEgenelist8.up01, DEgenelist8.down01, DEgenelist9.up01, DEgenelist9.down01)
regulated<-c(rep("up", length(DEgenelist7.up01)), rep("down", length(DEgenelist7.down01)), rep("up", length(DEgenelist8.up01)), rep("down", length(DEgenelist8.down01)), rep("up", length(DEgenelist9.up01)), rep("down", length(DEgenelist9.down01)))
month_comparison<-c(rep("June vs. April", length(DEgenelist7.up01)+ length(DEgenelist7.down01)), rep("November vs. June", length(DEgenelist8.up01)+ length(DEgenelist8.down01)), rep("April vs. November", length(DEgenelist9.up01)+ length(DEgenelist9.down01)))
Male_DE_genes.1<-data.frame(ensembl_gene_id, regulated, month_comparison, stringsAsFactors = FALSE)

# Remove .n from gene ID
Female_DE_genes.1$ensembl_gene_id = substr(Female_DE_genes.1$ensembl_gene_id,1,nchar(Female_DE_genes.1$ensembl_gene_id)-2) 
Male_DE_genes.1$ensembl_gene_id = substr(Male_DE_genes.1$ensembl_gene_id,1,nchar(Male_DE_genes.1$ensembl_gene_id)-2)

# Obtain external gene names for overlap DE genes and add to dataframe
annotLookupF_DE <-getBM(
  mart = mart,
  attributes = c(
    'ensembl_gene_id',
    'external_gene_name'),
  filter = 'ensembl_gene_id',
  values = Female_DE_genes.1$ensembl_gene_id,
  uniqueRows = TRUE)

annotLookupM_DE <-getBM(
  mart = mart,
  attributes = c(
    'ensembl_gene_id',
    'external_gene_name'),
  filter = 'ensembl_gene_id',
  values = Male_DE_genes.1$ensembl_gene_id,
  uniqueRows = TRUE)

Female_DE_genes.1<- full_join(Female_DE_genes.1, annotLookupF_DE, key="ensembl_gene_id")
Female_DE_genes<-Female_DE_genes.1 %>%
  mutate("regulated.month_comparison" =paste(Female_DE_genes.1$regulated, "-", Female_DE_genes.1$month_comparison)) %>%
  dplyr::select(-regulated, -month_comparison) %>%
  arrange(ensembl_gene_id) %>%
  group_by(ensembl_gene_id) %>%
  mutate("regulated.month_comparison" = paste0( regulated.month_comparison, collapse = " & ")) %>%
  unique() %>%
  mutate("regulated.month_comparison" = paste0("(", regulated.month_comparison, ")"))

Male_DE_genes.1<-full_join(Male_DE_genes.1, annotLookupM_DE, key="ensembl_gene_id")
Male_DE_genes<-Male_DE_genes.1 %>%
  mutate("regulated.month_comparison" =paste(Male_DE_genes.1$regulated, "-", Male_DE_genes.1$month_comparison)) %>%
  dplyr::select(-regulated, -month_comparison) %>%
  arrange(ensembl_gene_id) %>%
  group_by(ensembl_gene_id) %>%
  mutate("regulated.month_comparison" = paste0( regulated.month_comparison, collapse = " & ")) %>%
  unique() %>%
  mutate("regulated.month_comparison" = paste0("(", regulated.month_comparison, ")")) 

# export overlap DE genes with monthly up and down regulation
Female_DE_overlap_genes<-Female_DE_genes %>%
  filter(ensembl_gene_id %in% annotLookupF$ensembl_gene_id) %>%
  as.data.frame()

Male_DE_overlap_genes<-Male_DE_genes %>%
  filter(ensembl_gene_id %in% annotLookupM$ensembl_gene_id) %>%
  as.data.frame()

require(openxlsx)
sheets <- list("Female_DE_genes" = Female_DE_overlap_genes, "Male_DE_genes" = Male_DE_overlap_genes)
write.xlsx(sheets, "/Users/sunshine2.0/Desktop/WORK/output/DE_overlap_genes.xlsx")


## GENE ONTOLOGY ANALYSIS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~____________________---------------------------------

## BUILD plot of distribution of DE gene p-values  ***************  move below all exressed genes after this chunk<----
# Females
names(annotLookupF)[names(annotLookupF) == "ensembl_gene_id"] <- "gene"

# remove .n from gene ID
DEgenelist4.01p$gene = substr(DEgenelist4.01p$gene,1,nchar(DEgenelist4.01p$gene)-2) 
DEgenelist5.01p$gene = substr(DEgenelist5.01p$gene,1,nchar(DEgenelist5.01p$gene)-2) 
DEgenelist6.01p$gene = substr(DEgenelist6.01p$gene,1,nchar(DEgenelist6.01p$gene)-2) 

# add external gene name to list of DE 0.01 genes with pvalues
DEgenelist4.01pe<-left_join(DEgenelist4.01p, annotLookupF, by = "gene") %>%
  dplyr::select(c(gene, padj, month.comparison, external_gene_name))
DEgenelist5.01pe<-left_join(DEgenelist5.01p, annotLookupF, by = "gene") %>%
  dplyr::select(c(gene, padj, month.comparison, external_gene_name))
DEgenelist6.01pe<-left_join(DEgenelist6.01p, annotLookupF, by = "gene") %>%
  dplyr::select(c(gene, padj, month.comparison, external_gene_name))

# list of cherry picked genes with relevant functions to analysis
CherryF <- c("ANGPTL4", "LPL")

# combine list of DE 0.01 genes with pvalues and filter for overlap genes
Female_DE_genes_p<-rbind(DEgenelist4.01pe, DEgenelist5.01pe, DEgenelist6.01pe) %>%
  filter(gene %in% overlap0.01F$x) %>%
  mutate("cherry" = ifelse(external_gene_name %in% CherryF, paste0("yes"), paste0("no")))
Female_DE_genes_p$month.comparison <-factor(Female_DE_genes_p$month.comparison, levels= unique(Female_DE_genes_p$month.comparison))

# make plots of p-value ditribution

pvalueFDEplot<-ggplot(Female_DE_genes_p, aes(month.comparison, padj, label = external_gene_name, colour= cherry))+
  geom_text()+
  ggtitle("DE Female p-value distribution")

outfile<-"/Users/sunshine2.0/Desktop/WORK/output/pvalue_distribution_DE_F.png"

ggsave(outfile, 
       dpi=300, 
       device="png",
       width=10, 
       height=25,
       units="cm")

# Males
names(annotLookupM)[names(annotLookupM) == "ensembl_gene_id"] <- "gene"

# remove .n from gene ID
DEgenelist7.01p$gene = substr(DEgenelist7.01p$gene,1,nchar(DEgenelist7.01p$gene)-2) 
DEgenelist8.01p$gene = substr(DEgenelist8.01p$gene,1,nchar(DEgenelist8.01p$gene)-2) 
DEgenelist9.01p$gene = substr(DEgenelist9.01p$gene,1,nchar(DEgenelist9.01p$gene)-2) 

# add external gene name to list of DE 0.01 genes with pvalues
DEgenelist7.01pe<-left_join(DEgenelist7.01p, annotLookupM, by = "gene") %>%
  dplyr::select(c(gene, padj, month.comparison, external_gene_name))
DEgenelist8.01pe<-left_join(DEgenelist8.01p, annotLookupM, by = "gene") %>%
  dplyr::select(c(gene, padj, month.comparison, external_gene_name))
DEgenelist9.01pe<-left_join(DEgenelist9.01p, annotLookupM, by = "gene") %>%
  dplyr::select(c(gene, padj, month.comparison, external_gene_name))

# list of cherry picked genes with relevant functions to analysis
CherryM <- c("ITGB2", "LUM", "OGN", "ALDOB", "SLC37A2", "MOGAT2", "LPL", "DGAT2")

# combine list of DE 0.01 genes with pvalues and filter for overlap genes
Male_DE_genes_p<-rbind(DEgenelist7.01pe, DEgenelist8.01pe, DEgenelist9.01pe) %>%
  filter(gene %in% overlap0.01M$x) %>%
  mutate("cherry" = ifelse(external_gene_name %in% CherryM, paste0("yes"), paste0("no")))
Male_DE_genes_p$month.comparison <-factor(Male_DE_genes_p$month.comparison, levels= unique(Male_DE_genes_p$month.comparison))

# make plots of p-value ditribution

pvalueMDEplot<-ggplot(Male_DE_genes_p, aes(month.comparison, padj, label = external_gene_name, colour= cherry))+
  geom_text()+
  ggtitle("DE Male p-value distribution")

outfile<-"/Users/sunshine2.0/Desktop/WORK/output/pvalue_distribution_DE_M.png"

ggsave(outfile, 
       dpi=300, 
       device="png",
       width=10, 
       height=25,
       units="cm")

# export dataframe
sheets <- list("Female_DE_genes" = Female_DE_genes_p, "Male_DE_genes" = Male_DE_genes_p)
write.xlsx(sheets, "/Users/sunshine2.0/Desktop/WORK/output/DE_overlap_genes_pvalue.xlsx")

## Barplot of expression levels of cherry picked genes

CherryF <- c("ANGPTL4", "LPL", "MOGAT2", "DGAT2")
# obtain median log2tpm of each monthly comparison for each gene

grp<-unique(samples$grp)
samples$samples <- rownames(samples)

# A.M
A.M_meta<-samples%>%
  filter(grp == "A.M")
A.M_meta$samples
A.M<-log2tpm%>%
  dplyr::select(A.M_meta$samples)
A.M_median<-as.data.frame(apply(A.M, 1, median))
A.M_median$gene<- rownames(A.M_median)

# A.F
A.F_meta<-samples%>%
  filter(grp == "A.F")
A.F_meta$samples
A.F<-log2tpm%>%
  dplyr::select(A.F_meta$samples)
A.F_median<-as.data.frame(apply(A.F, 1, median))
A.F_median$gene<- rownames(A.F_median)

# J.M
J.M_meta<-samples%>%
  filter(grp == "J.M")
J.M_meta$samples
J.M<-log2tpm%>%
  dplyr::select(J.M_meta$samples)
J.M_median<-as.data.frame(apply(J.M, 1, median))
J.M_median$gene<- rownames(J.M_median)

# J.F
J.F_meta<-samples%>%
  filter(grp == "J.F")
J.F_meta$samples
J.F<-log2tpm%>%
  dplyr::select(J.F_meta$samples)
J.F_median<-as.data.frame(apply(J.F, 1, median))
J.F_median$gene<- rownames(J.F_median)

# N.M
N.M_meta<-samples%>%
  filter(grp == "N.M")
N.M_meta$samples
N.M<-log2tpm%>%
  dplyr::select(N.M_meta$samples)
N.M_median<-as.data.frame(apply(N.M, 1, median))
N.M_median$gene<- rownames(N.M_median)

# N.F
N.F_meta<-samples%>%
  filter(grp == "N.F")
N.F_meta$samples
N.F<-log2tpm%>%
  dplyr::select(N.F_meta$samples)
N.F_median<-as.data.frame(apply(N.F, 1, median))  
N.F_median$gene<- rownames(N.F_median)

# Combine each group medians
full_median<-full_join(A.M_median, A.F_median, by= "gene") %>%
  full_join(J.M_median, by= "gene" ) %>%
  full_join(J.F_median, by= "gene" ) %>%
  full_join(N.M_median, by= "gene" ) %>%
  full_join(N.F_median, by= "gene" ) 

# remove .n from gene ID
full_median$gene = substr(full_median$gene,1,nchar(full_median$gene)-2) 

rownames(full_median) <- full_median$gene
full_median$gene<-NULL
colnames(full_median) <- grp
full_median$gene<-rownames(full_median)

# split median log2tpm dataframe by sex, add external gene name, and filter for DE genes only
Female_median<-full_median %>%
  dplyr::select("gene", "A.F", "J.F", "N.F") %>%
  right_join( annotLookupF, by= "gene") %>%
  dplyr::select("gene", "A.F", "J.F", "N.F", "external_gene_name")
Male_median<-full_median %>%
  dplyr::select("gene","A.M", "J.M", "N.M") %>%
  right_join( annotLookupM, by= "gene") %>%
  dplyr::select("gene", "A.M", "J.M", "N.M", "external_gene_name")

# filter for cherry picked genes only
Female_median_cherry <- Female_median %>%
  filter(external_gene_name %in% CherryF)
colnames(Female_median_cherry)[2:4]<- c("April", "June", "November")
Male_median_cherry <- Male_median %>%
  filter(external_gene_name %in% CherryM)
colnames(Male_median_cherry)[2:4]<- c("April", "June", "November")

# The median absolute deviation is used to draw the error bars on the graph
A.F_mad<-as.data.frame(apply(A.F, 1, mad))
J.F_mad<-as.data.frame(apply(J.F, 1, mad))
N.F_mad<-as.data.frame(apply(N.F, 1, mad))
A.M_mad<-as.data.frame(apply(A.M, 1, mad))
J.M_mad<-as.data.frame(apply(J.M, 1, mad))
N.M_mad<-as.data.frame(apply(N.M, 1, mad))
A.F_mad$gene<- rownames(A.F_mad)
J.F_mad$gene<- rownames(J.F_mad)
N.F_mad$gene<- rownames(N.F_mad)
A.M_mad$gene<- rownames(A.M_mad)
J.M_mad$gene<- rownames(J.M_mad)
N.M_mad$gene<- rownames(N.M_mad)
# Combine each group mads
full_mad<-full_join(A.M_mad, A.F_mad, by= "gene") %>%
  full_join(J.M_mad, by= "gene" ) %>%
  full_join(J.F_mad, by= "gene" ) %>%
  full_join(N.M_mad, by= "gene" ) %>%
  full_join(N.F_mad, by= "gene" ) 

rownames(full_mad) <- full_mad$gene
full_mad$gene<-NULL
colnames(full_mad) <- grp
full_mad$gene<-rownames(full_mad)
# remove .n from gene ID
full_mad$gene = substr(full_mad$gene,1,nchar(full_mad$gene)-2) 

# split mad log2tpm dataframe by sex, add external gene name
Female_mad<-full_mad %>%
  dplyr::select("gene", "A.F", "J.F", "N.F")
colnames(Female_mad)[2:4]<- c("April.mad", "June.mad", "November.mad")
Male_mad<-full_mad %>%
  dplyr::select("gene","A.M", "J.M", "N.M") 
colnames(Male_mad)[2:4]<- c("April.mad", "June.mad", "November.mad")

#dataframe median log2tpm, mad, padj, gene name
full_med_mad<-left_join(full_median, full_mad, by= "gene")
colnames(full_med_mad)[1:6]<- c("Male_April_median_log2(TPM)", "Female_April_median_log2(TPM)", "Male_June_median_log2(TPM)", "Female_June_median_log2(TPM)","Male_November_median_log2(TPM)", "Female_November_median_log2(TPM)")
colnames(full_med_mad)[8:13]<- c("Male_April_MAD", "Female_April_MAD", "Male_June_MAD", "Female_June_MAD","Male_November_MAD", "Female_November_MAD")

DE_padj4<-DE_padj4%>% dplyr::select("gene", "padj")
colnames(DE_padj4)[2]<-"Female_June_vs_April_pvalue"
DE_padj4$gene = substr(DE_padj4$gene,1,nchar(DE_padj4$gene)-2) 
DE_padj5<-DE_padj5%>% dplyr::select("gene", "padj")
colnames(DE_padj5)[2]<-"Female_November_vs_June_pvalue"
DE_padj5$gene = substr(DE_padj5$gene,1,nchar(DE_padj5$gene)-2) 
DE_padj6<-DE_padj6%>% dplyr::select("gene", "padj")
colnames(DE_padj6)[2]<-"Female_April_vs_November_pvalue"
DE_padj6$gene = substr(DE_padj6$gene,1,nchar(DE_padj6$gene)-2) 
DE_padj7<-DE_padj7%>% dplyr::select("gene", "padj")
colnames(DE_padj7)[2]<-"Male_June_vs_April_pvalue"
DE_padj7$gene = substr(DE_padj7$gene,1,nchar(DE_padj7$gene)-2) 
DE_padj8<-DE_padj8%>% dplyr::select("gene", "padj")
colnames(DE_padj8)[2]<-"Male_November_vs_June_pvalue"
DE_padj8$gene = substr(DE_padj8$gene,1,nchar(DE_padj8$gene)-2) 
DE_padj9<-DE_padj9%>% dplyr::select("gene", "padj")
colnames(DE_padj9)[2]<-"Male_April_vs_November_pvalue"
DE_padj9$gene = substr(DE_padj9$gene,1,nchar(DE_padj9$gene)-2) 

full_padj<-full_join(DE_padj4, DE_padj5, by= "gene") %>%
  full_join(DE_padj6, by= "gene" ) %>%
  full_join(DE_padj7, by= "gene" ) %>%
  full_join(DE_padj8, by= "gene" ) %>%
  full_join(DE_padj9, by= "gene" ) 
colnames(full_padj)[2:7]<-c("Female_June_vs_April_pvalue", "Female_November_vs_June_pvalue", "Female_April_vs_November_pvalue", "Male_June_vs_April_pvalue", "Male_November_vs_June_pvalue", "Male_April_vs_November_pvalue")

full_med_mad_padj<-left_join(full_med_mad, full_padj, by= "gene")

mart <- useEnsembl(biomart = "ensembl", dataset = "dnovaehollandiae_gene_ensembl", version = "102")
full_gene_names <-getBM(
  mart = mart,
  attributes = c(
    'ensembl_gene_id',
    'external_gene_name'),
  filter = 'ensembl_gene_id',
  values = full_med_mad_padj$gene,
  uniqueRows = TRUE)
colnames(full_gene_names)[1]<-"gene"

full_med_mad_padj_names<-left_join(full_gene_names, full_med_mad_padj, by= "gene")
F_overlap_full_med_mad_padj_names<-full_med_mad_padj_names%>%
  filter(gene %in% overlap0.01F$x)
M_overlap_full_med_mad_padj_names<-full_med_mad_padj_names%>%
  filter(gene %in% overlap0.01M$x)

colnames(DE_padj1)[2]<-"April_Male_vs_Female_pvalue"
DE_padj1$gene = substr(DE_padj1$gene,1,nchar(DE_padj1$gene)-2) 
colnames(DE_padj2)[2]<-"June_Male_vs_Female_pvalue"
DE_padj2$gene = substr(DE_padj2$gene,1,nchar(DE_padj2$gene)-2) 
colnames(DE_padj3)[2]<-"November_Male_vs_Female_pvalue"
DE_padj3$gene = substr(DE_padj3$gene,1,nchar(DE_padj3$gene)-2) 

full_padj_sex<-full_join(DE_padj1, DE_padj2, by= "gene") %>%
  full_join(DE_padj3, by= "gene" )

full_sex_padj_name<-left_join(full_gene_names, full_padj_sex, by= "gene")
DE_full_sex_padj_name<-full_sex_padj_name%>%
  filter(gene %in% DE_genes$ensembl_gene_id)
  
sheet_final<-list("All genes"=full_med_mad_padj_names, "Female DE genes"=F_overlap_full_med_mad_padj_names, "Male DE genes"= M_overlap_full_med_mad_padj_names, "All genes Sex"=full_sex_padj_name, "DE genes Sex"=DE_full_sex_padj_name)
write.xlsx(sheet_final, "/Users/sunshine2.0/Desktop/WORK/output/Emu_expression_and_pvalues.xlsx")




# combine mad and median dataframes
Female_median_cherry_mad<-left_join(Female_median_cherry,Female_mad, by= "gene")
Male_median_cherry_mad<-left_join(Male_median_cherry,Male_mad, by= "gene")

# tidy data for plotting
library(reshape2)
dm1F <- melt(Female_median_cherry_mad[,c("gene","April","June","November","external_gene_name")], id=c("gene","external_gene_name"))
dm2F <- melt(Female_median_cherry_mad[,c("gene","April.mad","June.mad","November.mad","external_gene_name")], id=c("gene","external_gene_name"))
colnames(dm2F) <- c("gene", "external_gene_name", "variable2", "value2")
dm1F$gene.month <-paste0(dm1F$external_gene_name, dm1F$variable)
dm2F$gene.month <-paste0(dm1F$external_gene_name, dm1F$variable)
dmF <- full_join(dm1F, dm2F, by="gene.month") %>%
  dplyr::select("gene.x", "external_gene_name.x", "variable", "value", "value2")
colnames(dmF)<-c("gene", "external_gene_name", "month", "median_log2tpm", "mad")

dm1M <- melt(Male_median_cherry_mad[,c("gene","April","June","November","external_gene_name")], id=c("gene","external_gene_name"))
dm2M <- melt(Male_median_cherry_mad[,c("gene","April.mad","June.mad","November.mad","external_gene_name")], id=c("gene","external_gene_name"))
colnames(dm2M) <- c("gene", "external_gene_name", "variable2", "value2")
dm1M$gene.month <-paste0(dm1M$external_gene_name, dm1M$variable)
dm2M$gene.month <-paste0(dm1M$external_gene_name, dm1M$variable)
dmM <- full_join(dm1M, dm2M, by="gene.month") %>%
  dplyr::select("gene.x", "external_gene_name.x", "variable", "value", "value2")
colnames(dmM)<-c("gene", "external_gene_name", "month", "median_log2tpm", "mad")

# create barplot of cherry picked gene expression levels
ggplot(dmF, aes(month, median_log2tpm, fill=external_gene_name)) +
  geom_col(position="dodge")+
  ggtitle("Female: Expression Level by Month") +
  geom_errorbar(aes(ymin=median_log2tpm-mad, ymax=median_log2tpm+mad), width=.2, position=position_dodge(.9)) 

outfile<-"/Users/sunshine2.0/Desktop/WORK/output/cherry_picked_genes_expression_F.png"

ggsave(outfile, 
       dpi=300, 
       device="png",
       width=15, 
       height=15,
       units="cm")

ggplot(dmM, aes(month, median_log2tpm, fill=external_gene_name)) +
  geom_col(position="dodge")+
  ggtitle("Male: Expression Level by Month") +
  geom_errorbar(aes(ymin=median_log2tpm-mad, ymax=median_log2tpm+mad), width=.2, position=position_dodge(.9)) 

outfile<-"/Users/sunshine2.0/Desktop/WORK/output/cherry_picked_genes_expression_M.png"

ggsave(outfile, 
       dpi=300, 
       device="png",
       width=15, 
       height=15,
       units="cm")

#INDIVIDUAL EXPRESSION PLOTS FOR PATHWAY MECHANISMS

## obtain p-values for plot

#OBTAIN DATA ON ALL EXPRESSED GENES

#female
RLE$samples$samples<-rownames(RLE$samples) #filter for female only in metadata
RLE_F_meta<-RLE$samples %>%
  filter(gender=="F") 

RLEcounts <-as.data.frame(RLE$counts)
RLEcounts$gene<-rownames(RLE$counts) #filter for female only in counts
RLE_F_counts<-RLEcounts%>%
  dplyr::select(c(RLE_F_meta$samples, gene)) %>%
  filter_all(all_vars(. > 1)) #filter for counts >1
nrow(RLE$counts)
nrow(RLE_F_counts)

# remove .n from gene ID
RLE_F_counts$gene = substr(RLE_F_counts$gene,1,nchar(RLE_F_counts$gene)-2) 
require('biomaRt')

#get gene symbol of emsembl ID's
# mart <- useMart('ENSEMBL_MART_ENSEMBL')
# mart <- useDataset('dnovaehollandiae_gene_ensembl', mart)
mart <- useEnsembl(biomart = "ensembl", dataset = "dnovaehollandiae_gene_ensembl", version = "102")

annotLookupallF <-getBM( 
  #tried to get uniprot annotation for emu ***emu had no or few annotation for uniprot*--> used chicken uniprot annotation on panther
  mart = mart,
  attributes = c(
    'ensembl_gene_id',
    'gene_biotype',
    'external_gene_name',
    'uniprot_gn_symbol',
    'uniprot_gn_id'),
  filter = 'ensembl_gene_id',
  values = RLE_F_counts$gene,
  uniqueRows = TRUE)

names(annotLookupallF)[names(annotLookupallF) == "ensembl_gene_id"] <- "gene"

#add month comparisons
df_log_4<-df_log_4%>%
  mutate("month.comparison"="June vs. April")
df_log_5<-df_log_5%>%
  mutate("month.comparison"="November vs. June")
df_log_6<-df_log_6%>%
  mutate("month.comparison"="April vs. November")

# remove .n from gene ID
df_log_4$gene = substr(df_log_4$gene,1,nchar(df_log_4$gene)-2)   #i think already did this above somewhere
df_log_5$gene = substr(df_log_5$gene,1,nchar(df_log_5$gene)-2) 
df_log_6$gene = substr(df_log_6$gene,1,nchar(df_log_6$gene)-2) 

# add external gene name to list of genes with pvalues
df_log_4pe<-left_join(df_log_4, annotLookupallF, by = "gene") %>%
  dplyr::select(c(gene, padj, month.comparison, external_gene_name))
df_log_5pe<-left_join(df_log_5, annotLookupallF, by = "gene") %>%
  dplyr::select(c(gene, padj, month.comparison, external_gene_name))
df_log_6pe<-left_join(df_log_6, annotLookupallF, by = "gene") %>%
  dplyr::select(c(gene, padj, month.comparison, external_gene_name))

#list of female DE genes and month comparisons
Female_DE_genes_pe<-Female_DE_genes.1 %>%
  mutate("gene.month_comparison" =paste(Female_DE_genes.1$external_gene_name, "-", Female_DE_genes.1$month_comparison))
#list of female genes and month comparisons
Female_genes_pe_all<- rbind(df_log_4pe, df_log_5pe, df_log_6pe) %>%
  mutate("gene.month_comparison" =paste(external_gene_name, "-", month.comparison)) 
Female_genes_pe_all<-Female_genes_pe_all %>%
  mutate("overlap" = ifelse(((Female_genes_pe_all$gene.month_comparison %in% Female_DE_genes_pe$gene.month_comparison) & (gene %in% overlap0.01F$x)), paste0("*"), paste0(""))) %>%
  mutate("p_value"= paste0(signif(padj, digits=3), overlap)) %>%
  mutate("p_values"=paste(month.comparison, "=", p_value)) %>%
  mutate("month1"= str_split_fixed(month.comparison, " vs. ", 2)[,1])%>%
  mutate("month2"= str_split_fixed(month.comparison, " vs. ", 2)[,2])

#CREATE PATHVISIO p-value datframe
#filter genes with p-values for fatty acid synthesis pathway genes from wikipathways https://www.wikipathways.org/index.php/Pathway:WP782
FA_genes<-c("PPAP2C", "PPAP2B", "PPAP2A", "PLPP1", "PNPLA2", "MOGAT3", "MOGAT2", "MOGAT1", "LPL", "LIPF", "LIPE", "LIPC", "GPD1", "GPAM", "GNPAT", "GK2", "GK", "DGAT2", "DGAT1", "AYR1", "AGPS", "AGPAT5", "AGPAT4", "AGPAT3", "AGPAT2", "AGPAT1")
FA_sythesis_pathway<-Female_genes_pe_all %>%
  filter(external_gene_name %in% FA_genes)
FA_sythesis_pathway%>%
  dplyr::select(gene, external_gene_name)%>%
  unique()

#add ensembl chicken ID, entrez gene, and HGNC id to dataset for Pathvisio (entrz gene id only imported into Pathvisio)
##did not use only needed external gene name which I already had
library("biomartr")
library("Biostrings")

# mart <- useMart('ENSEMBL_MART_ENSEMBL')
# mart <- useDataset('dnovaehollandiae_gene_ensembl', mart)
mart <- useEnsembl(biomart = "ensembl", dataset = "dnovaehollandiae_gene_ensembl", version = "102")
#a<-getAttributes(dataset="dnovaehollandiae_gene_ensembl", mart = 'ENSEMBL_MART_ENSEMBL')
ggallus_allF1 <-getBM(
  mart = mart,
  attributes = c(
    'ensembl_gene_id',
    'external_gene_name',
    'ggallus_homolog_ensembl_gene'),
  filter = 'ensembl_gene_id',
  values = Female_genes_pe_all$gene,
  uniqueRows = TRUE)

# mart <- useMart('ENSEMBL_MART_ENSEMBL')
# mart <- useDataset('ggallus_gene_ensembl', mart)
a<-getAttributes(dataset="ggallus_gene_ensembl", mart = 'ENSEMBL_MART_ENSEMBL')
mart <- useEnsembl(biomart = "ensembl", dataset = "ggallus_gene_ensembl", version = "102")

ggallus_allF <-getBM(
  mart = mart,
  attributes = c(
    'ensembl_gene_id',
    'entrezgene_id'),
  filter = 'ensembl_gene_id',
  values = ggallus_allF1$ggallus_homolog_ensembl_gene,
    uniqueRows = TRUE)
#rearrange dataframe with p-values for each month split up in new columns (ex. June vs. April p-values)
Female_all_pvalues_monthly<-Female_genes_pe_all%>%
  dplyr::select(gene, external_gene_name, padj, month.comparison)%>%
  pivot_wider(names_from = month.comparison, 
              values_from = padj)

names(Female_all_pvalues_monthly)[3:5]<-c("June vs. April p-value", "November vs. June p-value", "April vs. November p-value")
names(ggallus_allF)[1]<-"ggallus_homolog_ensembl_gene"

ggallus_allF_join<-full_join(ggallus_allF1, ggallus_allF, by="ggallus_homolog_ensembl_gene") %>%
  unique()
names(ggallus_allF_join)[1]<-"gene"
#add ggallus ensembl id to dataframe 
chicken_pvalues_monthly_F<-full_join(Female_all_pvalues_monthly, ggallus_allF_join, by="gene") %>%
  unique()
#remove ENSDNVG00000018159	since is a duplicate
chicken_pvalues_monthly_F<-filter(chicken_pvalues_monthly_F, ggallus_homolog_ensembl_gene !="ENSGALG00000050616")

#chicken_pvalues_monthly_F$entrezgene_id[chicken_pvalues_monthly_F$external_gene_name=='GPD1']<-424263 #add missing entrez id
#chicken_pvalues_monthly_F$entrezgene_id[chicken_pvalues_monthly_F$external_gene_name=='GNPAT']<-421548

#export dataset of all female p-values for Pathvisio https://projects.bigcat.unimaas.nl/ebi-roadshow/pathway-analysis/instructions/
write.csv(chicken_pvalues_monthly_F, "/Users/sunshine2.0/Desktop/WORK/output/emu_entrez_pvalue_pathvisio.csv")
#https://bridgedb.github.io/data/gene_database/
#110448 rows of data were imported succesfully.
#45540 identifiers were not recognized
#new 20433 rows of data were imported succesfully.9981 identifiers were not recognized


#Male pathvisio

#add month comparisons
df_log_7<-df_log_7%>%
  mutate("month.comparison"="June vs. April")
df_log_8<-df_log_8%>%
  mutate("month.comparison"="November vs. June")
df_log_9<-df_log_9%>%
  mutate("month.comparison"="April vs. November")

# remove .n from gene ID
df_log_7$gene = substr(df_log_7$gene,1,nchar(df_log_7$gene)-2)  
df_log_8$gene = substr(df_log_8$gene,1,nchar(df_log_8$gene)-2) 
df_log_9$gene = substr(df_log_9$gene,1,nchar(df_log_9$gene)-2) 

# add external gene name to list of genes with pvalues
df_log_7pe<-left_join(df_log_7, annotLookupallF, by = "gene") %>%
  dplyr::select(c(gene, padj, month.comparison, external_gene_name))
df_log_8pe<-left_join(df_log_8, annotLookupallF, by = "gene") %>%
  dplyr::select(c(gene, padj, month.comparison, external_gene_name))
df_log_9pe<-left_join(df_log_9, annotLookupallF, by = "gene") %>%
  dplyr::select(c(gene, padj, month.comparison, external_gene_name))

#list of female DE genes and month comparisons
Male_DE_genes_pe<-Male_DE_genes.1 %>%
  mutate("gene.month_comparison" =paste(Male_DE_genes.1$external_gene_name, "-", Male_DE_genes.1$month_comparison))
#list of female genes and month comparisons
Male_genes_pe_all<- rbind(df_log_7pe, df_log_8pe, df_log_9pe) %>%
  mutate("gene.month_comparison" =paste(external_gene_name, "-", month.comparison)) 
Male_genes_pe_all<-Male_genes_pe_all %>%
  mutate("overlap" = ifelse(((Male_genes_pe_all$gene.month_comparison %in% Male_DE_genes_pe$gene.month_comparison) & (gene %in% overlap0.01M$x)), paste0("*"), paste0(""))) %>%
  mutate("p_value"= paste0(signif(padj, digits=3), overlap)) %>%
  mutate("p_values"=paste(month.comparison, "=", p_value)) %>%
  mutate("month1"= str_split_fixed(month.comparison, " vs. ", 2)[,1])%>%
  mutate("month2"= str_split_fixed(month.comparison, " vs. ", 2)[,2])

#CREATE PATHVISIO p-value datframe
#filter genes with p-values for fatty acid synthesis pathway genes from wikipathways https://www.wikipathways.org/index.php/Pathway:WP782
FA_genes<-c("PPAP2C", "PPAP2B", "PPAP2A", "PLPP1", "PNPLA2", "MOGAT3", "MOGAT2", "MOGAT1", "LPL", "LIPF", "LIPE", "LIPC", "GPD1", "GPAM", "GNPAT", "GK2", "GK", "DGAT2", "DGAT1", "AYR1", "AGPS", "AGPAT5", "AGPAT4", "AGPAT3", "AGPAT2", "AGPAT1")
FA_sythesis_pathway_M<-Male_genes_pe_all %>%
  filter(external_gene_name %in% FA_genes)
FA_sythesis_pathway_M%>%
  dplyr::select(gene, external_gene_name)%>%
  unique()


#add ensembl chicken ID, entrez gene, and HGNC id to dataset for Pathvisio (entrz gene id only imported into Pathvisio)
##did not use only needed external gene name which I already had
library("biomartr")
library("Biostrings")

# mart <- useMart('ENSEMBL_MART_ENSEMBL')
# mart <- useDataset('dnovaehollandiae_gene_ensembl', mart)
mart <- useEnsembl(biomart = "ensembl", dataset = "dnovaehollandiae_gene_ensembl", version = "102")
#a<-getAttributes(dataset="dnovaehollandiae_gene_ensembl", mart = 'ENSEMBL_MART_ENSEMBL')
ggallus_allM1 <-getBM(
  mart = mart,
  attributes = c(
    'ensembl_gene_id',
    'external_gene_name',
    'ggallus_homolog_ensembl_gene'),
  filter = 'ensembl_gene_id',
  values = Male_genes_pe_all$gene,
  uniqueRows = TRUE)

# mart <- useMart('ENSEMBL_MART_ENSEMBL')
# mart <- useDataset('ggallus_gene_ensembl', mart)
a<-getAttributes(dataset="ggallus_gene_ensembl", mart = 'ENSEMBL_MART_ENSEMBL')
mart <- useEnsembl(biomart = "ensembl", dataset = "ggallus_gene_ensembl", version = "102")

ggallus_allM <-getBM(
  mart = mart,
  attributes = c(
    'ensembl_gene_id',
    'entrezgene_id'),
  filter = 'ensembl_gene_id',
  values = ggallus_allM1$ggallus_homolog_ensembl_gene,
  uniqueRows = TRUE)
#rearrange dataframe with p-values for each month split up in new columns (ex. June vs. April p-values)
Male_all_pvalues_monthly<-Male_genes_pe_all%>%
  dplyr::select(gene, external_gene_name, padj, month.comparison)%>%
  pivot_wider(names_from = month.comparison, 
              values_from = padj)

names(Male_all_pvalues_monthly)[3:5]<-c("June vs. April p-value", "November vs. June p-value", "April vs. November p-value")
names(ggallus_allM)[1]<-"ggallus_homolog_ensembl_gene"

ggallus_allM_join<-full_join(ggallus_allM1, ggallus_allM, by="ggallus_homolog_ensembl_gene") %>%
  unique()
names(ggallus_allM_join)[1]<-"gene"
#add ggallus ensembl id to dataframe 
chicken_pvalues_monthly_M<-full_join(Male_all_pvalues_monthly, ggallus_allM_join, by="gene") %>%
  unique()
#remove ENSDNVG00000018159	since is a duplicate
chicken_pvalues_monthly_M<-filter(chicken_pvalues_monthly_M, ggallus_homolog_ensembl_gene !="ENSGALG00000050616")

#chicken_pvalues_monthly_F$entrezgene_id[chicken_pvalues_monthly_F$external_gene_name=='GPD1']<-424263 #add missing entrez id
#chicken_pvalues_monthly_F$entrezgene_id[chicken_pvalues_monthly_F$external_gene_name=='GNPAT']<-421548

#export dataset of all female p-values for Pathvisio https://projects.bigcat.unimaas.nl/ebi-roadshow/pathway-analysis/instructions/
write.csv(chicken_pvalues_monthly_M, "/Users/sunshine2.0/Desktop/WORK/output/emu_entrez_pvalue_pathvisio_males.csv")
#https://bridgedb.github.io/data/gene_database/
#110448 rows of data were imported succesfully.
#45540 identifiers were not recognized
#new 20433 rows of data were imported succesfully.9981 identifiers were not recognized
Male_genes_pe_all%>%
  filter(external_gene_name %in% FA_genes)%>% unique() # only MOGAT2
library(ggpubr)

#MOGAT2 MALE
MOGAT2<-dmM%>% filter(external_gene_name == "MOGAT2")
# plot items
MOGAT2_pvalues_M<-Male_genes_pe_all%>% #pvalues
  filter(external_gene_name == "MOGAT2") %>%
  dplyr::select(month1, month2, p_value) %>%
  as_tibble()
names(MOGAT2_pvalues_M)<-c("group1", "group2", "p.adj")
MOGAT2_pvalues_M$group1[MOGAT2_pvalues_M$group1 == "Nov"] <- "November"#change Nov to November
MOGAT2_pvalues_M$group2[MOGAT2_pvalues_M$group2 == "Nov"] <- "November"

library('ggpubr')
# plt <- 
ggplot(MOGAT2, aes(month, median_log2tpm, fill="red")) +
  geom_col()+
  labs(title = "MOGAT2 Expression") + theme(plot.title = element_text(hjust = 0.5))+
  geom_errorbar(aes(ymin=median_log2tpm-mad, ymax=median_log2tpm+mad), width=.2, size =1) +
  guides(fill=FALSE)+
  theme_bw()+
  theme(plot.title = element_text(size =26), axis.title=element_text(size = 27), axis.text=element_text(size = 25))+
  ylim(0, 15)+
  ylab("Median log2(TPM)")+
  xlab("Month")+
  stat_pvalue_manual(
    MOGAT2_pvalues_M, 
    y.position = 7, step.increase = 0.4,
    label = "p.adj",
    size=9,
    bracket.size=1
  )+
  scale_fill_grey()

outfile<-"/Users/sunshine2.0/Desktop/WORK/output/expression_M_MOGAT2_brackets.png"
ggsave(outfile, 
       dpi=300, 
       device="png",
       width=13, 
       height=15,
       units="cm")

#LPL MALE 
dmMM_med_mad<-full_med_mad_padj_names %>%
  filter(external_gene_name %in% c("LPL", "DGAT2"))%>%
  dplyr::select(gene, external_gene_name, 
                "Male_April_median_log2(TPM)", "Male_June_median_log2(TPM)", "Male_November_median_log2(TPM)",
                Male_April_MAD, Male_June_MAD, Male_November_MAD)
colnames(dmMM_med_mad)[3:8]<-c("April", "June", "November", "April.mad", "June.mad", "November.mad")

dm1MM <- melt(dmMM_med_mad[,c("gene","April","June","November","external_gene_name")], id=c("gene","external_gene_name"))
dm2MM <- melt(dmMM_med_mad[,c("gene","April.mad","June.mad","November.mad","external_gene_name")], id=c("gene","external_gene_name"))
colnames(dm2MM) <- c("gene", "external_gene_name", "variable2", "value2")
dm1MM$gene.month <-paste0(dm1MM$external_gene_name, dm1MM$variable)
dm2MM$gene.month <-paste0(dm1MM$external_gene_name, dm1MM$variable)
dmMM <- full_join(dm1MM, dm2MM, by="gene.month") %>%
  dplyr::select("gene.x", "external_gene_name.x", "variable", "value", "value2")
colnames(dmMM)<-c("gene", "external_gene_name", "month", "median_log2tpm", "mad")




dmMM_pvalues<-full_med_mad_padj_names %>%
  filter(external_gene_name %in% c("LPL", "DGAT2"))%>%
  dplyr::select(gene, external_gene_name, 
                Male_June_vs_April_pvalue, Male_November_vs_June_pvalue, Male_April_vs_November_pvalue )
colnames(dmMM_pvalues)[3:5]<-c("June.April", "November.June", "April.November")
dmMMM1_pvalues <- melt(dmMM_pvalues[,c("gene","external_gene_name","June.April","November.June","April.November")], id=c("gene","external_gene_name"))%>%
  separate(variable, c("month1", "month2"))
colnames(dmMMM1_pvalues)[5]<-"p_value"
# dmMMM_pvalues$p_value<-signif(dmMMM_pvalues$p_value, digits=3)
#dmMMM_pvalues$p_value<-formatC(dmMMM_pvalues$p_value, format = "e", digits=3)
# dmMMM_pvalues<-dmMMM_pvalues%>%
#   mutate(p_value= ifelse(p_value<0.01, paste(signif(dmMMM_pvalues$p_value, digits=3),"*"), paste(signif(dmMMM_pvalues$p_value, digits=3))))
# dmMMM_pvalues<-dmMMM_pvalues%>%
  # mutate(p_value= ifelse(p_value<0.01, paste(formatC(p_value, format = "e", digits =2),"*"), paste(formatC(p_value, format = "e", digits =2))))
dmMMM_pvalues<-dmMMM1_pvalues%>%
  mutate("ast" = ifelse((p_value<0.01), paste0("*"), paste0(""))) %>%
  mutate("p_value"= paste0(signif(p_value, digits=3), ast)) 
dmMMM_pvalues$p_value[4]<-paste0(formatC(dmMMM1_pvalues$p_value[4], format = "e", digits=2), dmMMM_pvalues$ast[4])
dmMMM_pvalues$p_value[5]<-formatC(dmMMM1_pvalues$p_value[5], format = "e", digits=2)
dmMMM_pvalues$p_value[6]<-formatC(dmMMM1_pvalues$p_value[6], format = "e", digits=2)

#filter for LPL
LPL<-dmMM%>% filter(external_gene_name == "LPL")
# plot items
LPL_pvalues<-dmMMM_pvalues%>% #pvalues
  filter(external_gene_name == "LPL") %>%
  dplyr::select(month1, month2, p_value) %>%
  as_tibble()
names(LPL_pvalues)<-c("group1", "group2", "p.adj")
LPL_pvalues$group1[LPL_pvalues$group1 == "Nov"] <- "November"#change Nov to November
LPL_pvalues$group2[LPL_pvalues$group2 == "Nov"] <- "November"

# plt <-
ggplot(LPL, aes(month, median_log2tpm, fill="red")) +
  geom_col()+
  labs(title = "LPL Expression") + theme(plot.title = element_text(hjust = 0.5))+
  geom_errorbar(aes(ymin=median_log2tpm-mad, ymax=median_log2tpm+mad), width=.2, size =1) +
  guides(fill=FALSE)+
  theme_bw()+
  theme(plot.title = element_text(size =30), axis.title=element_text(size = 27), axis.text=element_text(size = 25))+
  ylim(0, 15)+
  ylab("Median log2(TPM)")+
  xlab("Month")+
  stat_pvalue_manual(
    LPL_pvalues,
    y.position = 11.5, step.increase = 0.35,
    label = "p.adj",
    size=9,
    bracket.size=1
  )+
  scale_fill_grey()

outfile<-"/Users/sunshine2.0/Desktop/WORK/output/expression_M_LPL_brackets.png"
ggsave(outfile,
       dpi=300,
       device="png",
       width=13,
       height=15,
       units="cm")

#DGAT2 MALE
DGAT2<-dmMM%>% filter(external_gene_name == "DGAT2")
# plot items
DGAT2_pvalues<-dmMMM_pvalues%>% #pvalues
  filter(external_gene_name == "DGAT2") %>%
  dplyr::select(month1, month2, p_value) %>%
  as_tibble()
names(DGAT2_pvalues)<-c("group1", "group2", "p.adj")
DGAT2_pvalues$group1[DGAT2_pvalues$group1 == "Nov"] <- "November"#change Nov to November
DGAT2_pvalues$group2[DGAT2_pvalues$group2 == "Nov"] <- "November"

# plt <-
ggplot(DGAT2, aes(month, median_log2tpm, fill="red")) +
  geom_col()+
  labs(title = "DGAT2 Expression") + theme(plot.title = element_text(hjust = 0.5))+
  geom_errorbar(aes(ymin=median_log2tpm-mad, ymax=median_log2tpm+mad), width=.2, size =1) +
  guides(fill=FALSE)+
  theme_bw()+
  theme(plot.title = element_text(size =26), axis.title=element_text(size = 27), axis.text=element_text(size = 25))+
  ylim(0, 15)+
  ylab("Median log2(TPM)")+
  xlab("Month")+
  stat_pvalue_manual(
    DGAT2_pvalues,
    y.position = 11.5, step.increase = 0.26,
    label = "p.adj",
    size=9,
    bracket.size=1
  )+
  scale_fill_grey()

outfile<-"/Users/sunshine2.0/Desktop/WORK/output/expression_M_DGAT2_brackets.png"
ggsave(outfile,
       dpi=300,
       device="png",
       width=13,
       height=15,
       units="cm")

#LPL FEMALE
#filter for LPL
LPL<-dmF%>% filter(external_gene_name == "LPL")
# plot items
LPL_pvalues<-Female_genes_pe_all%>% #pvalues
  filter(external_gene_name == "LPL") %>%
  dplyr::select(month1, month2, p_value) %>%
  as_tibble()
names(LPL_pvalues)<-c("group1", "group2", "p.adj")
LPL_pvalues$group1[LPL_pvalues$group1 == "Nov"] <- "November"#change Nov to November
LPL_pvalues$group2[LPL_pvalues$group2 == "Nov"] <- "November"

# plt <- 
ggplot(LPL, aes(month, median_log2tpm, fill="red")) +
  geom_col()+
  labs(title = "LPL Expression") + theme(plot.title = element_text(hjust = 0.5))+
  geom_errorbar(aes(ymin=median_log2tpm-mad, ymax=median_log2tpm+mad), width=.2, size =1) +
  guides(fill=FALSE)+
  theme_bw()+
  theme(plot.title = element_text(size =30), axis.title=element_text(size = 27), axis.text=element_text(size = 25))+
  ylim(0, 15)+
  ylab("Median log2(TPM)")+
  xlab("Month")+
  stat_pvalue_manual(
    LPL_pvalues, 
    y.position = 11.5, step.increase = 0.28,
    label = "p.adj",
    size=9, 
    bracket.size=1
  )+
  scale_fill_grey()

outfile<-"/Users/sunshine2.0/Desktop/WORK/output/expression_F_LPL_brackets.png"
ggsave(outfile, 
       dpi=300, 
       device="png",
       width=13, 
       height=15,
       units="cm")

#MOGAT2 Female
MOGAT2<-dmF%>% filter(external_gene_name == "MOGAT2")
# plot items
MOGAT2_pvalues<-Female_genes_pe_all%>% #pvalues
  filter(external_gene_name == "MOGAT2") %>%
  dplyr::select(month1, month2, p_value) %>%
  as_tibble()
names(MOGAT2_pvalues)<-c("group1", "group2", "p.adj")
MOGAT2_pvalues$group1[MOGAT2_pvalues$group1 == "Nov"] <- "November"#change Nov to November
MOGAT2_pvalues$group2[MOGAT2_pvalues$group2 == "Nov"] <- "November"

# plt <- 
ggplot(MOGAT2, aes(month, median_log2tpm, fill="red")) +
  geom_col()+
  labs(title = "MOGAT2 Expression") + theme(plot.title = element_text(hjust = 0.5))+
  geom_errorbar(aes(ymin=median_log2tpm-mad, ymax=median_log2tpm+mad), width=.2, size =1) +
  guides(fill=FALSE)+
  theme_bw()+
  theme(plot.title = element_text(size =26), axis.title=element_text(size = 27), axis.text=element_text(size = 25))+
  ylim(0, 15)+
  ylab("Median log2(TPM)")+
  xlab("Month")+
  stat_pvalue_manual(
    MOGAT2_pvalues, 
    y.position = 8, step.increase = 0.28,
    label = "p.adj",
    size=9,
    bracket.size=1
  )+
  scale_fill_grey()

outfile<-"/Users/sunshine2.0/Desktop/WORK/output/expression_F_MOGAT2_brackets.png"
ggsave(outfile, 
       dpi=300, 
       device="png",
       width=13, 
       height=15,
       units="cm")


#DGAT2 Female
DGAT2<-dmF%>% filter(external_gene_name == "DGAT2")
# plot items
DGAT2_pvalues<-Female_genes_pe_all%>% #pvalues
  filter(external_gene_name == "DGAT2") %>%
  dplyr::select(month1, month2, p_value) %>%
  as_tibble()
names(DGAT2_pvalues)<-c("group1", "group2", "p.adj")
DGAT2_pvalues$group1[DGAT2_pvalues$group1 == "Nov"] <- "November"#change Nov to November
DGAT2_pvalues$group2[DGAT2_pvalues$group2 == "Nov"] <- "November"

# plt <- 
ggplot(DGAT2, aes(month, median_log2tpm, fill="red")) +
  geom_col()+
  labs(title = "DGAT2 Expression") + theme(plot.title = element_text(hjust = 0.5))+
  geom_errorbar(aes(ymin=median_log2tpm-mad, ymax=median_log2tpm+mad), width=.2, size =1) +
  guides(fill=FALSE)+
  theme_bw()+
  theme(plot.title = element_text(size =26), axis.title=element_text(size = 27), axis.text=element_text(size = 25))+
  ylim(0, 15)+
  ylab("Median log2(TPM)")+
  xlab("Month")+
  stat_pvalue_manual(
    DGAT2_pvalues, 
    y.position = 12, step.increase = 0.19,
    label = "p.adj",
    size=9,
    bracket.size=1
  )+
  scale_fill_grey()

outfile<-"/Users/sunshine2.0/Desktop/WORK/output/expression_F_DGAT2_brackets.png"
ggsave(outfile, 
       dpi=300, 
       device="png",
       width=13, 
       height=15,
       units="cm")

# #INDIVIDUAL EXPRESSION PLOTS FOR PATHWAY MECHANISMS with tables
# 
# ## obtain p-values for plot
# 
# #OBTAIN DATA ON ALL EXPRESSED GENES
# 
# #female
# RLE$samples$samples<-rownames(RLE$samples) #filter for female only in metadata
# RLE_F_meta<-RLE$samples %>%
#   filter(gender=="F") 
# 
# RLEcounts <-as.data.frame(RLE$counts)
# RLEcounts$gene<-rownames(RLE$counts) #filter for female only in counts
# RLE_F_counts<-RLEcounts%>%
#   dplyr::select(c(RLE_F_meta$samples, gene)) %>%
#   filter_all(all_vars(. > 1)) #filter for counts >1
# nrow(RLE$counts)
# nrow(RLE_F_counts)
# 
# # remove .n from gene ID
# RLE_F_counts$gene = substr(RLE_F_counts$gene,1,nchar(RLE_F_counts$gene)-2) 
# require('biomaRt')
# 
# #get gene symbol of emsembl ID's
# mart <- useMart('ENSEMBL_MART_ENSEMBL')
# mart <- useDataset('dnovaehollandiae_gene_ensembl', mart)
# 
# annotLookupallF <-getBM( 
#   #tried to get uniprot annotation for emu ***emu had no or few annotation for uniprot*--> used chicken uniprot annotation on panther
#   mart = mart,
#   attributes = c(
#     'ensembl_gene_id',
#     'gene_biotype',
#     'external_gene_name',
#     'uniprot_gn_symbol',
#     'uniprot_gn_id'),
#   filter = 'ensembl_gene_id',
#   values = RLE_F_counts$gene,
#   uniqueRows = TRUE)
# 
# names(annotLookupallF)[names(annotLookupallF) == "ensembl_gene_id"] <- "gene"
# 
# df_log_4<-df_log_4%>%
#   mutate("month.comparison"="June vs. April")
# df_log_5<-df_log_5%>%
#   mutate("month.comparison"="Nov vs. June")
# df_log_6<-df_log_6%>%
#   mutate("month.comparison"="April vs. Nov")
# 
# # remove .n from gene ID
# df_log_4$gene = substr(df_log_4$gene,1,nchar(df_log_4$gene)-2) 
# df_log_5$gene = substr(df_log_5$gene,1,nchar(df_log_5$gene)-2) 
# df_log_6$gene = substr(df_log_6$gene,1,nchar(df_log_6$gene)-2) 
# 
# # add external gene name to list of DE 0.01 genes with pvalues
# df_log_4pe<-left_join(df_log_4, annotLookupallF, by = "gene") %>%
#   dplyr::select(c(gene, padj, month.comparison, external_gene_name))
# df_log_5pe<-left_join(df_log_5, annotLookupallF, by = "gene") %>%
#   dplyr::select(c(gene, padj, month.comparison, external_gene_name))
# df_log_6pe<-left_join(df_log_6, annotLookupallF, by = "gene") %>%
#   dplyr::select(c(gene, padj, month.comparison, external_gene_name))
# 
# #list of female DE genes and month comparisons
# Female_DE_genes_pe<-Female_DE_genes.1 %>%
#   mutate("gene.month_comparison" =paste(Female_DE_genes.1$external_gene_name, "-", Female_DE_genes.1$month_comparison))
# #list of female genes and month comparisons
# Female_genes_pe_all<- rbind(df_log_4pe, df_log_5pe, df_log_6pe) %>%
#   mutate("gene.month_comparison" =paste(external_gene_name, "-", month.comparison)) 
# Female_genes_pe_all<-Female_genes_pe_all %>%
#   mutate("overlap" = ifelse(((Female_genes_pe_all$gene.month_comparison %in% Female_DE_genes_pe$gene.month_comparison) & (gene %in% overlap0.01F$x)), paste0("*"), paste0(""))) %>%
#   mutate("p_value"= paste0(signif(padj, digits=3), overlap)) %>%
#   mutate("p_values"=paste(month.comparison, "=", p_value))
# 
# library(ggpubr)
# # p-values
# stat.test <- tibble::tribble(
#   ~group1, ~group2,   ~p.adj,
#   "0.5",     "1", 2.54e-07,
#   "0.5",     "2", 1.32e-13,
#   "1",     "2", 1.91e-05
# )
# stat.test
# #diff dfile names
# 
# #LPL Female
# LPL<-dmF%>% filter(external_gene_name == "LPL")
# # plot items
# library(gridExtra)
# library(grid)
# LPL_pvalues<-Female_genes_pe_all%>% #pvalues
#   filter(external_gene_name == "LPL") %>%
#   dplyr::select(month.comparison, p_value)%>%
#   as.matrix()
# my.tGrob <- tableGrob(LPL_pvalues)
# plt <- ggplot(LPL, aes(month, median_log2tpm, fill="red")) +
#   geom_col()+
#   labs(title = "Female: LPL Expression") + theme(plot.title = element_text(hjust = 0.5))+
#   geom_errorbar(aes(ymin=median_log2tpm-mad, ymax=median_log2tpm+mad), width=.2) +
#   guides(fill=FALSE)+
#   theme_bw(15)+
#   ylim(0, 15)
# 
# # layout
# vp.layout <- grid.layout(nrow=2, ncol=1, heights=unit(1, "null"),
#                          widths=unit(c(1,1), c("null","line")) )
# 
# png("/Users/sunshine2.0/Desktop/WORK/output/expression_F_LPL.png",width=4,height=5,units="in",res=1200)
# # start drawing
# grid.newpage()
# pushViewport(viewport(layout=vp.layout, name="layout"))
# # plot
# pushViewport(viewport(layout.pos.row=1, layout.pos.col=1, name="plot"))
# print(plt, newpage=FALSE)
# upViewport()
# # table
# pushViewport(viewport(layout.pos.row=2, layout.pos.col=1, name="table"))
# grid.draw(my.tGrob)
# upViewport()
# dev.off()
# 
# #MOGAT2 Female
# MOGAT2<-dmF%>% filter(external_gene_name == "MOGAT2")
# # plot items
# MOGAT2_pvalues<-Female_genes_pe_all%>% #pvalues
#   filter(external_gene_name == "MOGAT2") %>%
#   dplyr::select(month.comparison, p_value)%>%
#   as.matrix()
# my.tGrob <- tableGrob(MOGAT2_pvalues)
# plt <- ggplot(MOGAT2, aes(month, median_log2tpm, fill="red")) +
#   geom_col()+
#   labs(title = "Female: MOGAT2 Expression" ) + theme(plot.title = element_text(hjust = 0.5))+
#   geom_errorbar(aes(ymin=median_log2tpm-mad, ymax=median_log2tpm+mad), width=.2) +
#   guides(fill=FALSE)+
#   theme_bw(15)+
#   ylim(0, 15)
# 
# # layout
# vp.layout <- grid.layout(nrow=2, ncol=1, heights=unit(1, "null"),
#                          widths=unit(c(1,1), c("null","line")) )
# 
# png("/Users/sunshine2.0/Desktop/WORK/output/expression_F_MOGAT2.png",width=4,height=5,units="in",res=1200)
# # start drawing
# grid.newpage()
# pushViewport(viewport(layout=vp.layout, name="layout"))
# # plot
# pushViewport(viewport(layout.pos.row=1, layout.pos.col=1, name="plot"))
# print(plt, newpage=FALSE)
# upViewport()
# # table
# pushViewport(viewport(layout.pos.row=2, layout.pos.col=1, name="table"))
# grid.draw(my.tGrob)
# upViewport()
# dev.off()
# 
# 
# #DGAT2 Female
# DGAT2<-dmF%>% filter(external_gene_name == "DGAT2")
# # plot items
# DGAT2_pvalues<-Female_genes_pe_all%>% #pvalues
#   filter(external_gene_name == "DGAT2") %>%
#   dplyr::select(month.comparison, p_value)%>%
#   as.matrix()
# my.tGrob <- tableGrob(DGAT2_pvalues)
# plt <- ggplot(DGAT2, aes(month, median_log2tpm, fill="red")) +
#   geom_col()+
#   labs(title = "Female: DGAT2 Expression" ) + theme(plot.title = element_text(hjust = 0.5))+
#   geom_errorbar(aes(ymin=median_log2tpm-mad, ymax=median_log2tpm+mad), width=.2) +
#   guides(fill=FALSE)+
#   theme_bw(15)+
#   ylim(0, 15)
# 
# # layout
# vp.layout <- grid.layout(nrow=2, ncol=1, heights=unit(1, "null"),
#                          widths=unit(c(1,1), c("null","line")) )
# 
# png("/Users/sunshine2.0/Desktop/WORK/output/expression_F_DGAT2.png",width=4,height=5,units="in",res=1200)
# # start drawing
# grid.newpage()
# pushViewport(viewport(layout=vp.layout, name="layout"))
# # plot
# pushViewport(viewport(layout.pos.row=1, layout.pos.col=1, name="plot"))
# print(plt, newpage=FALSE)
# upViewport()
# # table
# pushViewport(viewport(layout.pos.row=2, layout.pos.col=1, name="table"))
# grid.draw(my.tGrob)
# upViewport()
# dev.off()



##Barplots KIM paper genes only

# change column name
names(annotLookupallF)[names(annotLookupallF) == "ensembl_gene_id"] <- "gene"
names(annotLookupallM)[names(annotLookupallM) == "ensembl_gene_id"] <- "gene"

# split median log2tpm dataframe by sex, add external gene name, and filter for all genes 
Female_median_all<-full_median %>%
  dplyr::select("gene", "A.F", "J.F", "N.F") %>%
  right_join( annotLookupallF, by= "gene") %>%
  dplyr::select("gene", "A.F", "J.F", "N.F", "external_gene_name")
Male_median_all<-full_median %>%
  dplyr::select("gene","A.M", "J.M", "N.M") %>%
  right_join( annotLookupallM, by= "gene") %>%
  dplyr::select("gene", "A.M", "J.M", "N.M", "external_gene_name")

# filter for Kim genes only
Kim<-c("FABP4", "SCD", "LEPR", "LEP", "ACTB")  #either no annotation or expression for ACTB, LEP
Female_median_kim <- Female_median_all %>%
  filter(external_gene_name %in% Kim)
colnames(Female_median_kim)[2:4]<- c("April", "June", "November")
Male_median_kim <- Male_median_all %>%
  filter(external_gene_name %in% Kim)
colnames(Male_median_kim)[2:4]<- c("April", "June", "November")

# The median absolute deviation is used to draw the error bars on the graph
# combine mad and kim median dataframes
Female_median_kim_mad<-left_join(Female_median_kim,Female_mad, by= "gene")
Male_median_kim_mad<-left_join(Male_median_kim,Male_mad, by= "gene")

# tidy data for plotting
library(reshape2)
dm1F <- melt(Female_median_kim_mad[,c("gene","April","June","November","external_gene_name")], id=c("gene","external_gene_name"))
dm2F <- melt(Female_median_kim_mad[,c("gene","April.mad","June.mad","November.mad","external_gene_name")], id=c("gene","external_gene_name"))
colnames(dm2F) <- c("gene", "external_gene_name", "variable2", "value2")
dm1F$gene.month <-paste0(dm1F$external_gene_name, dm1F$variable)
dm2F$gene.month <-paste0(dm1F$external_gene_name, dm1F$variable)
dmF <- full_join(dm1F, dm2F, by="gene.month") %>%
  dplyr::select("gene.x", "external_gene_name.x", "variable", "value", "value2")
colnames(dmF)<-c("gene", "external_gene_name", "month", "median_log2tpm", "mad")

dm1M <- melt(Male_median_kim_mad[,c("gene","April","June","November","external_gene_name")], id=c("gene","external_gene_name"))
dm2M <- melt(Male_median_kim_mad[,c("gene","April.mad","June.mad","November.mad","external_gene_name")], id=c("gene","external_gene_name"))
colnames(dm2M) <- c("gene", "external_gene_name", "variable2", "value2")
dm1M$gene.month <-paste0(dm1M$external_gene_name, dm1M$variable)
dm2M$gene.month <-paste0(dm1M$external_gene_name, dm1M$variable)
dmM <- full_join(dm1M, dm2M, by="gene.month") %>%
  dplyr::select("gene.x", "external_gene_name.x", "variable", "value", "value2")
colnames(dmM)<-c("gene", "external_gene_name", "month", "median_log2tpm", "mad")


###check if annotation for those non genes bc will be asked****************
annotLookupF
# create barplot of cherry picked gene expression levels
ggplot(dmF, aes(month, median_log2tpm, fill=external_gene_name)) +
  geom_col(position="dodge")+
  ggtitle("Female: Expression Level by Month") +
  geom_errorbar(aes(ymin=median_log2tpm-mad, ymax=median_log2tpm+mad), width=.2, position=position_dodge(.9)) 

outfile<-"/Users/sunshine2.0/Desktop/WORK/output/kim_picked_genes_expression_F.png"

ggsave(outfile, 
       dpi=300, 
       device="png",
       width=15, 
       height=15,
       units="cm")

ggplot(dmM, aes(month, median_log2tpm, fill=external_gene_name)) +
  geom_col(position="dodge")+
  ggtitle("Male: Expression Level by Month") +
  geom_errorbar(aes(ymin=median_log2tpm-mad, ymax=median_log2tpm+mad), width=.2, position=position_dodge(.9)) 

outfile<-"/Users/sunshine2.0/Desktop/WORK/output/kim_picked_genes_expression_M.png"

ggsave(outfile, 
       dpi=300, 
       device="png",
       width=15, 
       height=15,
       units="cm")









names(annotLookupF)[names(annotLookupF) == "gene"] <- "ensembl_gene_id"
names(annotLookupM)[names(annotLookupM) == "gene"] <- "ensembl_gene_id"
names(annotLookupallF)[names(annotLookupallF) == "gene"] <- "ensembl_gene_id"
names(annotLookupallM)[names(annotLookupallM) == "gene"] <- "ensembl_gene_id"

#********************************************************
library("gprofiler2")
# BP = biological process, M F= molecular function, CC = cellular components
# GO ANALYSIS FOR ALL EXPRESSED GENES


#OBTAIN DATA ON ALL EXPRESSED GENES

#female
RLE$samples$samples<-rownames(RLE$samples) #filter for female only in metadata
RLE_F_meta<-RLE$samples %>%
  filter(gender=="F") 

RLEcounts <-as.data.frame(RLE$counts)
RLEcounts$gene<-rownames(RLE$counts) #filter for female only in counts
RLE_F_counts<-RLEcounts%>%
  dplyr::select(c(RLE_F_meta$samples, gene)) %>%
  filter_all(all_vars(. > 1)) #filter for counts >1
nrow(RLE$counts)
nrow(RLE_F_counts)

#male
RLE$samples$samples<-rownames(RLE$samples) #filter for male only in metadata
RLE_M_meta<-RLE$samples %>%
  filter(gender=="M") 

RLE_M_counts<-RLEcounts%>% #filter for male only in counts
  dplyr::select(c(RLE_M_meta$samples, gene)) %>%
  filter_all(all_vars(. > 1)) #filter for counts >1
nrow(RLE$counts)
nrow(RLE_M_counts)

# remove .n from gene ID
RLE_F_counts$gene = substr(RLE_F_counts$gene,1,nchar(RLE_F_counts$gene)-2) 
RLE_M_counts$gene = substr(RLE_M_counts$gene,1,nchar(RLE_M_counts$gene)-2) 
require('biomaRt')

#get gene symbol of emsembl ID's
# mart <- useMart('ENSEMBL_MART_ENSEMBL')
# mart <- useDataset('dnovaehollandiae_gene_ensembl', mart, version = "102")
mart <- useEnsembl(biomart = "ensembl", dataset = "dnovaehollandiae_gene_ensembl", version = "102")


annotLookupallF <-getBM( 
  #tried to get uniprot annotation for emu ***emu had no or few annotation for uniprot*--> used chicken uniprot annotation on panther
  mart = mart,
  attributes = c(
    'ensembl_gene_id',
    'gene_biotype',
    'external_gene_name',
    'uniprot_gn_symbol',
    'uniprot_gn_id'),
  filter = 'ensembl_gene_id',
  values = RLE_F_counts$gene,
  uniqueRows = TRUE)
write.csv(annotLookupallF$external_gene_name, "/Users/sunshine2.0/Desktop/WORK/output/annotLookup$external_gene_name_all0.01F.csv")

annotLookupallM <-getBM(
  mart = mart,
  attributes = c(
    'ensembl_gene_id',
    'gene_biotype',
    'external_gene_name',
    'uniprot_gn_symbol',
    'uniprot_gn_id'),
  filter = 'ensembl_gene_id',
  values = RLE_M_counts$gene,
  uniqueRows = TRUE)
write.csv(annotLookupallM$external_gene_name, "/Users/sunshine2.0/Desktop/WORK/output/annotLookup$external_gene_name_all0.01M.csv")


#BP FEMALE all
gprofiler_Female_GO_BP_all <-gost(
  annotLookupallF$external_gene_name,
  organism = "ggallus",
  ordered_query = FALSE,
  multi_query = FALSE,
  significant = FALSE,
  exclude_iea = FALSE,
  measure_underrepresentation = FALSE,
  evcodes = FALSE,
  user_threshold = 0.05,
  correction_method = "fdr",
  domain_scope = "annotated",
  custom_bg = NULL,
  sources = "GO:BP",
  as_short_link = FALSE)

# Slim down GO categories using parent nodes
length(gprofiler_Female_GO_BP_all$result$term_id)
library("GOfuncR") #biocmanager::install
parent_Female_BP_all<-get_parent_nodes(gprofiler_Female_GO_BP_all$result$term_id, term_df = NULL, graph_path_df = NULL, godir = NULL)
slim_Female_BP_all<-parent_Female_BP_all %>% filter(parent_go_id %in% c("GO:0032502", "GO:0032501", "GO:0009987", "GO:0000003", "GO:0008283", "GO:0051179", "GO:0022414", "GO:0051704", "GO:0022610", "GO:0002376", "GO:0071840", "GO:0065007", "GO:0040007", "GO:0023052", "GO:0008152", "GO:0050896", "GO:0043473", "GO:0110148", "GO:0044848", "GO:0007610", "GO:0048511", "GO:0040011" ) )
#parent_Female_BP%>% filter(child_go_id %in% slim_Female_BP$child_go_id ==FALSE)
length(unique(slim_Female_BP_all$parent_name))
length(unique(slim_Female_BP_all$child_go_id))


#MF FEMALE all
gprofiler_Female_GO_MF_all <-gost(
  annotLookupallF$external_gene_name,
  organism = "ggallus",
  ordered_query = FALSE,
  multi_query = FALSE,
  significant = FALSE,
  exclude_iea = FALSE,
  measure_underrepresentation = FALSE,
  evcodes = FALSE,
  user_threshold = 0.05,
  correction_method = "fdr",
  domain_scope = "annotated",
  custom_bg = NULL,
  sources = "GO:MF",
  as_short_link = FALSE)

# Slim down GO categories using parent nodes
length(gprofiler_Female_GO_MF_all$result$term_id)
parent_Female_MF_all<-get_parent_nodes(gprofiler_Female_GO_MF_all$result$term_id, term_df = NULL, graph_path_df = NULL, godir = NULL)
slim_Female_MF_all<-parent_Female_MF_all %>% filter(parent_go_id %in% c("GO:0045182", "GO:0140110", "GO:0060089", "GO:0005488", "GO:0005198", "GO:0098772", "GO:0003824", "GO:0005215", "GO:0038024" , "GO:0016209","GO:0140104") ) #slim terms
length(unique(slim_Female_MF_all$parent_name))
length(unique(slim_Female_MF_all$child_go_id))



#CC FEMALE all
gprofiler_Female_GO_CC_all <-gost(
  annotLookupallF$external_gene_name,
  organism = "ggallus",
  ordered_query = FALSE,
  multi_query = FALSE,
  significant = FALSE,
  exclude_iea = FALSE,
  measure_underrepresentation = FALSE,
  evcodes = FALSE,
  user_threshold = 0.05,
  correction_method = "fdr",
  domain_scope = "annotated",
  custom_bg = NULL,
  sources = "GO:CC",
  as_short_link = FALSE)

# Slim down GO categories using parent nodes
length(unique(gprofiler_Female_GO_CC_all$result$term_id))
parent_Female_CC_all<-get_parent_nodes(gprofiler_Female_GO_CC_all$result$term_id, term_df = NULL, graph_path_df = NULL, godir = NULL)
slim_Female_CC_all<-parent_Female_CC_all %>% filter(parent_go_id %in% c("GO:0044456", "GO:0044425", "GO:0016020", "GO:0045202", "GO:0044422", "GO:0044421", "GO:0030054", "GO:0031974", "GO:0032991", "GO:0099080", "GO:0005576", "GO:0005623", "GO:0044464", "GO:0043226" ,"GO:0005622", "GO:0005737", "GO:0110165") ) #slim terms
length(unique(slim_Female_CC_all$parent_name))
length(unique(slim_Female_CC_all$child_go_id))
parent_Female_CC_all%>% filter(child_go_id %in% slim_Female_CC_all$child_go_id ==FALSE)

#BP MALE all
gprofiler_Male_GO_BP_all <-gost(
  annotLookupallM$external_gene_name,
  organism = "ggallus",
  ordered_query = FALSE,
  multi_query = FALSE,
  significant = FALSE,
  exclude_iea = FALSE,
  measure_underrepresentation = FALSE,
  evcodes = FALSE,
  user_threshold = 0.05,
  correction_method = "fdr",
  domain_scope = "annotated",
  custom_bg = NULL,
  sources = "GO:BP",
  as_short_link = FALSE)

# Slim down GO categories using parent nodes
length(gprofiler_Male_GO_BP_all$result$term_id)
library("GOfuncR") #biocmanager::install
parent_Male_BP_all<-get_parent_nodes(gprofiler_Male_GO_BP_all$result$term_id, term_df = NULL, graph_path_df = NULL, godir = NULL)
slim_Male_BP_all<-parent_Male_BP_all %>% filter(parent_go_id %in% c("GO:0032502", "GO:0032501", "GO:0009987", "GO:0000003", "GO:0008283", "GO:0051179", "GO:0022414", "GO:0051704", "GO:0022610", "GO:0002376", "GO:0071840", "GO:0065007", "GO:0040007", "GO:0023052", "GO:0008152", "GO:0050896", "GO:0043473", "GO:0110148", "GO:0044848", "GO:0007610", "GO:0048511", "GO:0040011" ) )
#parent_Female_BP%>% filter(child_go_id %in% slim_Female_BP$child_go_id ==FALSE)
length(unique(slim_Male_BP_all$parent_name))
length(unique(slim_Male_BP_all$child_go_id))

#MF MALE all
gprofiler_Male_GO_MF_all <-gost(
  annotLookupallM$external_gene_name,
  organism = "ggallus",
  ordered_query = FALSE,
  multi_query = FALSE,
  significant = FALSE,
  exclude_iea = FALSE,
  measure_underrepresentation = FALSE,
  evcodes = FALSE,
  user_threshold = 0.05,
  correction_method = "fdr",
  domain_scope = "annotated",
  custom_bg = NULL,
  sources = "GO:MF",
  as_short_link = FALSE)
 
# Slim down GO categories using parent nodes
length(gprofiler_Male_GO_MF_all$result$term_id)
parent_Male_MF_all<-get_parent_nodes(gprofiler_Male_GO_MF_all$result$term_id, term_df = NULL, graph_path_df = NULL, godir = NULL)
slim_Male_MF_all<-parent_Male_MF_all %>% filter(parent_go_id %in% c("GO:0045182", "GO:0140110", "GO:0060089", "GO:0005488", "GO:0005198", "GO:0098772", "GO:0003824", "GO:0005215", "GO:0038024" ,"GO:0016209","GO:0140104") ) #slim terms
#parent_Male_MF_all%>% filter(child_go_id %in% slim_Male_MF_all$child_go_id ==FALSE)
length(unique(slim_Male_MF_all$parent_name))
length(unique(slim_Male_MF_all$child_go_id))

#CC MALE all
gprofiler_Male_GO_CC_all <-gost(
  annotLookupallM$external_gene_name,
  organism = "ggallus",
  ordered_query = FALSE,
  multi_query = FALSE,
  significant = FALSE,
  exclude_iea = FALSE,
  measure_underrepresentation = FALSE,
  evcodes = FALSE,
  user_threshold = 0.05,
  correction_method = "fdr",
  domain_scope = "annotated",
  custom_bg = NULL,
  sources = "GO:CC",
  as_short_link = FALSE)

# Slim down GO categories using parent nodes
length(gprofiler_Male_GO_CC_all$result$term_id)
parent_Male_CC_all<-get_parent_nodes(gprofiler_Male_GO_CC_all$result$term_id, term_df = NULL, graph_path_df = NULL, godir = NULL)
slim_Male_CC_all<-parent_Male_CC_all %>% filter(parent_go_id %in% c("GO:0044456", "GO:0044425", "GO:0016020", "GO:0045202", "GO:0044422", "GO:0044421", "GO:0030054", "GO:0031974", "GO:0032991", "GO:0099080", "GO:0005576", "GO:0005623", "GO:0044464", "GO:0043226" ,"GO:0005622", "GO:0005737", "GO:0110165") ) #slim terms
#parent_Male_CC_all%>% filter(child_go_id %in% slim_Male_CC_all$child_go_id ==FALSE)
length(unique(slim_Male_CC_all$parent_name))
length(unique(slim_Male_CC_all$child_go_id))

# REFORMAT DATA FOR GO PLOT ALL
# Add male and female column
slim_Female_BP_all$sex<-"Female" 
slim_Male_BP_all$sex<-"Male" 
slim_Female_MF_all$sex<-"Female" 
slim_Male_MF_all$sex<-"Male" 
slim_Female_CC_all$sex<-"Female" 
slim_Male_CC_all$sex<-"Male" 

# Combine classification tables together and summarize counts and add percentages
BP_all<-as.data.frame(rbind(slim_Female_BP_all, slim_Male_BP_all)%>%
                    dplyr::select(child_go_id, parent_go_id, parent_name, sex) %>%
                    group_by(parent_go_id, parent_name, sex) %>%
                    summarise("n"=n_distinct(child_go_id))) 
BP_all<- BP_all%>%
  mutate("percentage_of_total"=ifelse(sex=="Female",(n/sum(BP_all%>%filter(sex=="Female") %>% dplyr::select(n))*100), (n/sum(BP_all%>%filter(sex=="Male") %>% dplyr::select(n))*100)))
MF_all<-as.data.frame(rbind(slim_Female_MF_all,slim_Male_MF_all) %>%
                    dplyr::select(child_go_id, parent_go_id, parent_name, sex) %>%
                    group_by(parent_go_id, parent_name, sex) %>%
                    summarise("n"=n_distinct(child_go_id))) 
MF_all<-MF_all%>%
  mutate("percentage_of_total"=ifelse(sex=="Female",(n/sum(MF_all%>%filter(sex=="Female") %>% dplyr::select(n))*100), (n/sum(MF_all%>%filter(sex=="Male") %>% dplyr::select(n))*100))) 
CC_all<-as.data.frame(rbind(slim_Female_CC_all, slim_Male_CC_all) %>%
                    dplyr::select(child_go_id, parent_go_id, parent_name, sex) %>%
                    group_by(parent_go_id, parent_name, sex) %>%
                    summarise("n"=n_distinct(child_go_id))) 
CC_all<-CC_all%>%
  mutate("percentage_of_total"=ifelse(sex=="Female",(n/sum(CC_all%>%filter(sex=="Female") %>% dplyr::select(n))*100), (n/sum(CC_all%>%filter(sex=="Male") %>% dplyr::select(n))*100)))

names(BP_all)[names(BP_all) == "parent_name"] <- "category"
names(MF_all)[names(MF_all) == "parent_name"] <- "category"
names(CC_all)[names(CC_all) == "parent_name"] <- "category"

## X^2 TESTS
# Statistical similarity with X^2 tests and Bonferroni correction factor of 6 and 85 
library(MASS)

BP_Pivot_all<- BP_all %>% pivot_wider(id_cols=category, #reformat dataframe
                                      names_from= sex,
                                      values_from= n) 
BP_Pivot_all$Male[is.na(BP_Pivot_all$Male) ] <- 0
BP_Pivot_all$Female[is.na(BP_Pivot_all$Female) ] <- 0

MF_Pivot_all<- MF_all %>% pivot_wider(id_cols=category, #reformat dataframe
                                      names_from= sex,
                                      values_from= n) 
MF_Pivot_all$Female[is.na(MF_Pivot_all$Female) ] <- 0

CC_Pivot_all<- CC_all %>% pivot_wider(id_cols=category, #reformat dataframe
                                      names_from= sex,
                                      values_from= n) 
CC_Pivot_all$Female[is.na(CC_Pivot_all$Female) ] <- 0

# Create dataframe of p-values from statistical tests
prop_test_all<-data.frame("type"= c("BP", "MF", "CC", rep("BP", 22),rep("MF", 11),rep("CC", 11)), 
                          "category"=c(rep("all", 3),BP_Pivot_all$category[1:22],MF_Pivot_all$category[1:11],CC_Pivot_all$category[1:11]),
                          "pvalue"=c(prop.test(x = as.matrix(BP_Pivot_all[2:3]),correct=FALSE)[3]$p.value,
                                                                                 prop.test(x = as.matrix(MF_Pivot_all[,2:3]),correct=FALSE)[3]$p.value, 
                                                                                 prop.test(x = as.matrix(CC_Pivot_all[,2:3]),correct=FALSE)[3]$p.value,
                                                                                 prop.test(x = as.matrix(BP_Pivot_all[1,2:3]),correct=FALSE)[3]$p.value,
                                                                                 prop.test(x = as.matrix(BP_Pivot_all[2,2:3]),correct=FALSE)[3]$p.value,
                                                                                 prop.test(x = as.matrix(BP_Pivot_all[3,2:3]),correct=FALSE)[3]$p.value,
                                                                                 prop.test(x = as.matrix(BP_Pivot_all[4,2:3]),correct=FALSE)[3]$p.value,
                                                                                 prop.test(x = as.matrix(BP_Pivot_all[5,2:3]),correct=FALSE)[3]$p.value,
                                                                                 prop.test(x = as.matrix(BP_Pivot_all[6,2:3]),correct=FALSE)[3]$p.value,
                                                                                 prop.test(x = as.matrix(BP_Pivot_all[7,2:3]),correct=FALSE)[3]$p.value,
                                                                                 prop.test(x = as.matrix(BP_Pivot_all[8,2:3]),correct=FALSE)[3]$p.value,
                                                                                 prop.test(x = as.matrix(BP_Pivot_all[9,2:3]),correct=FALSE)[3]$p.value,
                                                                                 prop.test(x = as.matrix(BP_Pivot_all[10,2:3]),correct=FALSE)[3]$p.value,
                                                                                 prop.test(x = as.matrix(BP_Pivot_all[11,2:3]),correct=FALSE)[3]$p.value,
                                                                                 prop.test(x = as.matrix(BP_Pivot_all[12,2:3]),correct=FALSE)[3]$p.value,
                                                                                 prop.test(x = as.matrix(BP_Pivot_all[13,2:3]),correct=FALSE)[3]$p.value,
                                                                                 prop.test(x = as.matrix(BP_Pivot_all[14,2:3]),correct=FALSE)[3]$p.value,
                                                                                 prop.test(x = as.matrix(BP_Pivot_all[15,2:3]),correct=FALSE)[3]$p.value,
                                                                                 prop.test(x = as.matrix(BP_Pivot_all[16,2:3]),correct=FALSE)[3]$p.value,
                                                                                 prop.test(x = as.matrix(BP_Pivot_all[17,2:3]),correct=FALSE)[3]$p.value,
                                                                                 prop.test(x = as.matrix(BP_Pivot_all[18,2:3]),correct=FALSE)[3]$p.value,
                                                                                 prop.test(x = as.matrix(BP_Pivot_all[19,2:3]),correct=FALSE)[3]$p.value,
                                                                                 prop.test(x = as.matrix(BP_Pivot_all[20,2:3]),correct=FALSE)[3]$p.value,
                                                                                 prop.test(x = as.matrix(BP_Pivot_all[21,2:3]),correct=FALSE)[3]$p.value,
                                                                                 prop.test(x = as.matrix(BP_Pivot_all[22,2:3]),correct=FALSE)[3]$p.value,
                                                                                 prop.test(x = as.matrix(MF_Pivot_all[1,2:3]),correct=FALSE)[3]$p.value,
                                                                                 prop.test(x = as.matrix(MF_Pivot_all[2,2:3]),correct=FALSE)[3]$p.value,
                                                                                 prop.test(x = as.matrix(MF_Pivot_all[3,2:3]),correct=FALSE)[3]$p.value,
                                                                                 prop.test(x = as.matrix(MF_Pivot_all[4,2:3]),correct=FALSE)[3]$p.value,
                                                                                 prop.test(x = as.matrix(MF_Pivot_all[5,2:3]),correct=FALSE)[3]$p.value,
                                                                                 prop.test(x = as.matrix(MF_Pivot_all[6,2:3]),correct=FALSE)[3]$p.value,
                                                                                 prop.test(x = as.matrix(MF_Pivot_all[7,2:3]),correct=FALSE)[3]$p.value,
                                                                                 prop.test(x = as.matrix(MF_Pivot_all[8,2:3]),correct=FALSE)[3]$p.value,
                                                                                 prop.test(x = as.matrix(MF_Pivot_all[9,2:3]),correct=FALSE)[3]$p.value,
                                                                                 prop.test(x = as.matrix(MF_Pivot_all[10,2:3]),correct=FALSE)[3]$p.value,
                                                                                 prop.test(x = as.matrix(MF_Pivot_all[11,2:3]),correct=FALSE)[3]$p.value,
                                                                                 prop.test(x = as.matrix(CC_Pivot_all[1,2:3]),correct=FALSE)[3]$p.value,
                                                                                 prop.test(x = as.matrix(CC_Pivot_all[2,2:3]),correct=FALSE)[3]$p.value,
                                                                                 prop.test(x = as.matrix(CC_Pivot_all[3,2:3]),correct=FALSE)[3]$p.value,
                                                                                 prop.test(x = as.matrix(CC_Pivot_all[4,2:3]),correct=FALSE)[3]$p.value,
                                                                                 prop.test(x = as.matrix(CC_Pivot_all[5,2:3]),correct=FALSE)[3]$p.value,
                                                                                 prop.test(x = as.matrix(CC_Pivot_all[6,2:3]),correct=FALSE)[3]$p.value,
                                                                                 prop.test(x = as.matrix(CC_Pivot_all[7,2:3]),correct=FALSE)[3]$p.value,
                                                                                 prop.test(x = as.matrix(CC_Pivot_all[8,2:3]),correct=FALSE)[3]$p.value,
                                                                                 prop.test(x = as.matrix(CC_Pivot_all[9,2:3]),correct=FALSE)[3]$p.value,
                                                                                 prop.test(x = as.matrix(CC_Pivot_all[10,2:3]),correct=FALSE)[3]$p.value,
                                                                                 prop.test(x = as.matrix(CC_Pivot_all[11,2:3]),correct=FALSE)[3]$p.value)) %>%
 
  #Sidak padj=1-(1-p)^n 
  mutate( "Sidak"=c(1-(1-pvalue[1:3])^6,1-(1-pvalue[4:47])^83)) %>%
  mutate("reject"=ifelse(Sidak<0.05, "reject", "fail to reject"))


BP_all_pvalue<-prop_test_all%>% filter(type== "BP", category!= "all") %>% dplyr::select(category, Sidak, reject)%>% full_join(BP_all, by= "category")%>%
  mutate("categoryPvalue"= ifelse(Sidak<0.9995, paste0(category, " (", formatC(Sidak, format = "e", digits =2), ")", ifelse(reject=="reject",paste0("*"), paste0(""))), paste0(category)) )
MF_all_pvalue<-prop_test_all%>% filter(type== "MF", category!= "all") %>% dplyr::select(category, Sidak, reject)%>% full_join(MF_all, by= "category")%>%
  mutate("categoryPvalue"= ifelse(Sidak<0.9995, paste0(category, " (", formatC(Sidak, format = "e",digits = 2), ")", ifelse(reject=="reject",paste0("*"), paste0(""))), paste0(category)) )
CC_all_pvalue<-prop_test_all%>% filter(type== "CC", category!= "all") %>% dplyr::select(category, Sidak, reject)%>% full_join(CC_all, by= "category")%>%
  mutate("categoryPvalue"= ifelse(Sidak<0.9995, paste0(category, " (", formatC(Sidak, format = "e",digits = 2), ")", ifelse(reject=="reject",paste0("*"), paste0(""))), paste0(category)) )


# GO PLOT FOR ALL EXPRESSED GENES
# Plot bar plot of percent of gene hit against total number of process hits
BP_plot_all<-BP_all_pvalue%>%
  ggplot(aes(sex, percentage_of_total, fill=category))+
  geom_bar(position="stack", stat="identity") +
  ylab("Percentage of Total") +
  xlab("Sex") +
  ggtitle(paste0("Biological Processes")) +
  theme_bw()+
  theme(axis.text.x = element_text(size=31),
           axis.text.y = element_text(size=31),
           axis.title.x = element_text(size=33),
           axis.title.y = element_text(size=33),
           plot.title = element_text(size=38),
           legend.text = element_text(size=33),
           legend.title = element_text(size=33))  +
  labs(fill="Category")+
  guides(color = guide_legend(override.aes = list(size = 15)))

MF_plot_all<-MF_all_pvalue%>%
  ggplot(aes(sex, percentage_of_total, fill=category))+
  geom_bar(position="stack", stat="identity") +
  ylab("Percentage of Total") +
  xlab("Sex") +
  ggtitle(paste0("Molecular Functions")) +
  theme_bw()+
  theme(axis.text.x = element_text(size=31),
        axis.text.y = element_text(size=31),
        axis.title.x = element_text(size=33),
        axis.title.y = element_text(size=33),
        plot.title = element_text(size=38),
        legend.text = element_text(size=33),
        legend.title = element_text(size=33))  +
  labs(fill="Category                                                                                                                            ")+
  guides(color = guide_legend(override.aes = list(size = 15)))

CC_plot_all<-CC_all_pvalue%>%
  ggplot(aes(sex, percentage_of_total, fill=category))+
  geom_bar(position="stack", stat="identity") +
  ylab("Percentage of Total") +
  xlab("Sex") +
  ggtitle(paste0("Cellular Components")) +
  theme_bw()+
  theme(axis.text.x = element_text(size=31),
        axis.text.y = element_text(size=31),
        axis.title.x = element_text(size=33),
        axis.title.y = element_text(size=33),
        plot.title = element_text(size=38),
        legend.text = element_text(size=33),
        legend.title = element_text(size=33))  +
  labs(fill="Category                                                                                                                           ")+
  guides(color = guide_legend(override.aes = list(size = 15)))




# Arrange GO bar plots together
GO_figure_all <-ggarrange(BP_plot_all, MF_plot_all, CC_plot_all,
                          widths = rep(5), 
                          hjust = -0.5,
                          vjust = 1,
                          ncol = 1, nrow = 3,
                          align = "hv") 

outfile.GO<-"/Users/sunshine2.0/Desktop/WORK/output/GO_bar_gProfiler_all.png"
ggsave(outfile.GO, 
       dpi=300, 
       device="png",
       width=68, 
       height=58,
       units="cm")




# GO ANALYSIS FOR OVERLAP OF BOTH UP AND DOWN REGULATED GENES
set_base_url("http://biit.cs.ut.ee/gprofiler_archive3/e102_eg49_p15")
#BP FEMALE overlap
gprofiler_Female_GO_BP <-gost(
  annotLookupF$external_gene_name,
  organism = "ggallus",
  ordered_query = FALSE,
  multi_query = FALSE,
  significant = FALSE,
  exclude_iea = FALSE,
  measure_underrepresentation = FALSE,
  evcodes = FALSE,
  user_threshold = 0.05,
  correction_method = "fdr",
  domain_scope = "annotated",
  custom_bg = NULL,
  sources = "GO:BP",
  as_short_link = FALSE)

library("GOfuncR") #biocmanager::install

# Slim down GO categories using parent nodes
length(gprofiler_Female_GO_BP$result$term_id)
parent_Female_BP<-get_parent_nodes(gprofiler_Female_GO_BP$result$term_id, term_df = NULL, graph_path_df = NULL, godir = NULL)
slim_Female_BP<-parent_Female_BP %>% filter(parent_go_id %in% c("GO:0032502", "GO:0032501", "GO:0009987", "GO:0000003", "GO:0008283", "GO:0051179", "GO:0022414", "GO:0051704", "GO:0022610", "GO:0002376", "GO:0071840", "GO:0065007", "GO:0040007", "GO:0023052", "GO:0008152", "GO:0050896", "GO:0043473", "GO:0110148", "GO:0044848", "GO:0007610", "GO:0048511", "GO:0040011" ) )
#parent_Female_BP%>% filter(child_go_id %in% slim_Female_BP$child_go_id ==FALSE)
length(unique(slim_Female_BP$parent_name))
length(unique(slim_Female_BP$child_go_id))


#MF FEMALE overlap
gprofiler_Female_GO_MF <-gost(
  annotLookupF$external_gene_name,
  organism = "ggallus",
  ordered_query = FALSE,
  multi_query = FALSE,
  significant = FALSE,
  exclude_iea = FALSE,
  measure_underrepresentation = FALSE,
  evcodes = FALSE,
  user_threshold = 0.05,
  correction_method = "fdr",
  domain_scope = "annotated",
  custom_bg = NULL,
  sources = "GO:MF",
  as_short_link = FALSE)

# Slim down GO categories using parent nodes
length(gprofiler_Female_GO_MF$result$term_id)
parent_Female_MF<-get_parent_nodes(gprofiler_Female_GO_MF$result$term_id, term_df = NULL, graph_path_df = NULL, godir = NULL)
slim_Female_MF<-parent_Female_MF %>% filter(parent_go_id %in% c("GO:0045182", "GO:0140110", "GO:0060089", "GO:0005488", "GO:0005198", "GO:0098772", "GO:0003824", "GO:0005215", "GO:0038024","GO:0016209" ,"GO:0140104") ) #slim terms
# parent_Female_MF%>% filter(child_go_id %in% slim_Female_MF$child_go_id ==FALSE)
length(unique(slim_Female_MF$parent_name))
length(unique(slim_Female_MF$child_go_id))

#CC FEMALE overlap
gprofiler_Female_GO_CC <-gost(
  annotLookupF$external_gene_name,
  organism = "ggallus",
  ordered_query = FALSE,
  multi_query = FALSE,
  significant = FALSE,
  exclude_iea = FALSE,
  measure_underrepresentation = FALSE,
  evcodes = FALSE,
  user_threshold = 0.05,
  correction_method = "fdr",
  domain_scope = "annotated",
  custom_bg = NULL,
  sources = "GO:CC",
  as_short_link = FALSE)

# Slim down GO categories using parent nodes
length(gprofiler_Female_GO_CC$result$term_id)
parent_Female_CC<-get_parent_nodes(gprofiler_Female_GO_CC$result$term_id, term_df = NULL, graph_path_df = NULL, godir = NULL)
slim_Female_CC<-parent_Female_CC %>% filter(parent_go_id %in% c("GO:0044456", "GO:0044425", "GO:0016020", "GO:0045202", "GO:0044422", "GO:0044421", "GO:0030054", "GO:0031974", "GO:0032991", "GO:0099080", "GO:0005576", "GO:0005623", "GO:0044464", "GO:0043226", "GO:0005622", "GO:0005737", "GO:0110165" ) ) #slim terms
#parent_Female_CC%>% filter(child_go_id %in% slim_Female_CC$child_go_id ==FALSE)
length(unique(slim_Female_CC$parent_name))
length(unique(slim_Female_CC$child_go_id))

#BP MALE overlap
gprofiler_Male_GO_BP <-gost(
  annotLookupM$external_gene_name,
  organism = "ggallus",
  ordered_query = FALSE,
  multi_query = FALSE,
  significant = FALSE,
  exclude_iea = FALSE,
  measure_underrepresentation = FALSE,
  evcodes = FALSE,
  user_threshold = 0.05,
  correction_method = "fdr",
  domain_scope = "annotated",
  custom_bg = NULL,
  sources = "GO:BP",
  as_short_link = FALSE)

# Slim down GO categories using parent nodes
length(gprofiler_Male_GO_BP$result$term_id)
parent_Male_BP<-get_parent_nodes(gprofiler_Male_GO_BP$result$term_id, term_df = NULL, graph_path_df = NULL, godir = NULL)
slim_Male_BP<-parent_Male_BP %>% filter(parent_go_id %in% c("GO:0032502", "GO:0032501", "GO:0009987", "GO:0000003", "GO:0008283", "GO:0051179", "GO:0022414", "GO:0051704", "GO:0022610", "GO:0002376", "GO:0071840", "GO:0065007", "GO:0040007", "GO:0023052", "GO:0008152", "GO:0050896", "GO:0043473", "GO:0110148", "GO:0044848", "GO:0007610", "GO:0048511", "GO:0040011" ) )
#parent_Male_BP%>% filter(child_go_id %in% slim_Male_BP$child_go_id ==FALSE) 
length(unique(slim_Male_BP$parent_name))
length(unique(slim_Male_BP$child_go_id))

#*****************************************************************
# Find genes that make up immune system process for males

# Create function to obtain all immune system process GO insignificant results for a given gene

genes_male_immune <-function(gene){I1<-gost( 
  gene,
  organism = "ggallus",
  ordered_query = FALSE,
  multi_query = FALSE,
  significant = FALSE,
  exclude_iea = FALSE,
  measure_underrepresentation = FALSE,
  evcodes = FALSE,
  user_threshold = 0.05,
  correction_method = "fdr",
  domain_scope = "annotated",
  custom_bg = NULL,
  sources = "GO:BP",
  as_short_link = FALSE)
parent_I1<-get_parent_nodes(I1$result$term_id, term_df = NULL, graph_path_df = NULL, godir = NULL)
slim_I1<-parent_I1 %>% filter(parent_go_id %in% "GO:0002376" )
slim_I1}

genes_males <- annotLookupM$external_gene_name[-which(annotLookupM$external_gene_name == "") ] #remove NA blanks
# Obtain all immune system process GO insignificant results for each gene inputted 
# genes_male_immune(annotLookupM$external_gene_name[1])
genes_male_immune(genes_males[1])
genes_male_immune(genes_males[2])
genes_male_immune(genes_males[3])
genes_male_immune(genes_males[4])#*
genes_male_immune(genes_males[5])
genes_male_immune(genes_males[6])
genes_male_immune(genes_males[7])
genes_male_immune(genes_males[8])
genes_male_immune(genes_males[9])
genes_male_immune(genes_males[10])
genes_male_immune(genes_males[11])
genes_male_immune(genes_males[12])
genes_male_immune(genes_males[13])
genes_male_immune(genes_males[14])
genes_male_immune(genes_males[15])
genes_male_immune(genes_males[16])
genes_male_immune(genes_males[17])
genes_male_immune(genes_males[18])
genes_male_immune(genes_males[19])
genes_male_immune(genes_males[20])#*
genes_male_immune(genes_males[21])
genes_male_immune(genes_males[22])
genes_male_immune(genes_males[23])
genes_male_immune(genes_males[24])
genes_male_immune(genes_males[25])
genes_male_immune(genes_males[26])
genes_male_immune(genes_males[27])
genes_male_immune(genes_males[28])
genes_male_immune(genes_males[29])
genes_male_immune(genes_males[30])
genes_male_immune(genes_males[31])
genes_male_immune(genes_males[32])
genes_male_immune(genes_males[33])
genes_male_immune(genes_males[34])
genes_male_immune(genes_males[35])
genes_male_immune(genes_males[36])
#for(i in seq_along(1:length(genes_males))){genes_male_immune(genes_males[i]) slim_[i]<-slim_I1 }**
 # if(nrow(slim_I1)>0){immune.genes.M<-rbind(immune.genes.M,slim_I1)} } **
#Genes that are in immune system process
genes_males[4]
genes_males[20]

# 
# slim_males<-data.frame(child_go_id=character(0),parent_go_id=character(0), parent_name=character(0),distance=integer(0), gene=character(0))
# #slim_males<-matrix("double", 1, 5)
# for (i in 1:36){slim_I1<-data.frame(child_go_id=character(0),parent_go_id=character(0), parent_name=character(0),distance=integer(0), gene=character(0))
#   slim_I1<-genes_male_immune(genes_males[i])
#   if (nrow(slim_I1)!=0 ) {  #is not null
#   slim_I1<-slim_I1 %>% mutate("gene"= rep(genes_males[i]))
#   slim_males<-rbind(slim_males,slim_I1)
#   slim_I1<-0
#   }}
# slim_males
# 
# slim_I1<-data.frame(child_go_id=character(0),parent_go_id=character(0), parent_name=character(0),distance=integer(0), gene=character(0))
# slim_I1<-genes_male_immune(genes_males[10])
# if (nrow(slim_I1)!=0 ) {  #is not null
#   slim_I1<-slim_I1 %>% mutate("gene"= rep(genes_males[4]))
#   slim_males<-rbind(slim_males,slim_I1)
#   slim_I1<-0}
# 
# slim_I1<-genes_male_immune(genes_males[10])
# nrow(slim_I1)!=0 
# 
# slim_I2<-data.frame(child_go_id=character(0),parent_go_id=character(0), parent_name=character(0),distance=integer(0), gene=character(0))
# nrow(slim_I2)!=0 


# do the same for cellular GO term
genes_male_cellular <-function(gene){I1<-gost(
  gene,
  organism = "ggallus",
  ordered_query = FALSE,
  multi_query = FALSE,
  significant = FALSE,
  exclude_iea = FALSE,
  measure_underrepresentation = FALSE,
  evcodes = FALSE,
  user_threshold = 0.05,
  correction_method = "fdr",
  domain_scope = "annotated",
  custom_bg = NULL,
  sources = "GO:BP",
  as_short_link = FALSE)
parent_I1<-get_parent_nodes(I1$result$term_id, term_df = NULL, graph_path_df = NULL, godir = NULL)
slim_I1<-parent_I1 %>% filter(parent_go_id %in% "GO:0009987" )
slim_I1}
genes_male_cellular(genes_males[1])#*
genes_male_cellular(genes_males[2])
genes_male_cellular(genes_males[3])#*
genes_male_cellular(genes_males[4])#*
genes_male_cellular(genes_males[5])#*
genes_male_cellular(genes_males[6])#*
genes_male_cellular(genes_males[7])#*
genes_male_cellular(genes_males[8])#*
genes_male_cellular(genes_males[9])#*
genes_male_cellular(genes_males[10])#*
genes_male_cellular(genes_males[11])#*
genes_male_cellular(genes_males[12])#*
genes_male_cellular(genes_males[13])#*
genes_male_cellular(genes_males[14])#*
genes_male_cellular(genes_males[15])#*
genes_male_cellular(genes_males[16])#*
genes_male_cellular(genes_males[17])
genes_male_cellular(genes_males[18])#*
genes_male_cellular(genes_males[19])#*
genes_male_cellular(genes_males[20])#*
genes_male_cellular(genes_males[21])
genes_male_cellular(genes_males[22])#*
genes_male_cellular(genes_males[23])
genes_male_cellular(genes_males[24])#*
genes_male_cellular(genes_males[25])#*
genes_male_cellular(genes_males[26])#*
genes_male_cellular(genes_males[27])#*
genes_male_cellular(genes_males[28])#*
genes_male_cellular(genes_males[29])#*
genes_male_cellular(genes_males[30])#*
genes_male_cellular(genes_males[31])#*
genes_male_cellular(genes_males[32])
genes_male_cellular(genes_males[33])#*
genes_male_cellular(genes_males[34])#*
genes_male_cellular(genes_males[35])
genes_male_cellular(genes_males[36])#*
genes_males[c(1, 3:16, 18:20, 22, 24:31, 33, 34, 36) ]

# do the same for localization GO term
genes_male_localization <-function(gene){I1<-gost(
  gene,
  organism = "ggallus",
  ordered_query = FALSE,
  multi_query = FALSE,
  significant = FALSE,
  exclude_iea = FALSE,
  measure_underrepresentation = FALSE,
  evcodes = FALSE,
  user_threshold = 0.05,
  correction_method = "fdr",
  domain_scope = "annotated",
  custom_bg = NULL,
  sources = "GO:BP",
  as_short_link = FALSE)
parent_I1<-get_parent_nodes(I1$result$term_id, term_df = NULL, graph_path_df = NULL, godir = NULL)
slim_I1<-parent_I1 %>% filter(parent_go_id %in% "GO:0051179" )
slim_I1}
genes_male_localization(genes_males[1])
genes_male_localization(genes_males[2])
genes_male_localization(genes_males[3])
genes_male_localization(genes_males[4])
genes_male_localization(genes_males[5])
genes_male_localization(genes_males[6])
genes_male_localization(genes_males[7])
genes_male_localization(genes_males[8])#*
genes_male_localization(genes_males[9])
genes_male_localization(genes_males[10])
genes_male_localization(genes_males[11])#*
genes_male_localization(genes_males[12])
genes_male_localization(genes_males[13])
genes_male_localization(genes_males[14])
genes_male_localization(genes_males[15])#*
genes_male_localization(genes_males[16])
genes_male_localization(genes_males[17])
genes_male_localization(genes_males[18])
genes_male_localization(genes_males[19])
genes_male_localization(genes_males[20])#*
genes_male_localization(genes_males[21])
genes_male_localization(genes_males[22])
genes_male_localization(genes_males[23])
genes_male_localization(genes_males[24])
genes_male_localization(genes_males[25])#*
genes_male_localization(genes_males[26])
genes_male_localization(genes_males[27])
genes_male_localization(genes_males[28])
genes_male_localization(genes_males[29])#*
genes_male_localization(genes_males[30])
genes_male_localization(genes_males[31])#*
genes_male_localization(genes_males[32])#*
genes_male_localization(genes_males[33])
genes_male_localization(genes_males[34])
genes_male_localization(genes_males[35])
genes_male_localization(genes_males[36])
genes_males[c(8, 11, 15, 20, 25, 29, 31, 32)]


# do the same for females immune GO term
genes_females <- annotLookupF$external_gene_name[-which(annotLookupF$external_gene_name == "") ]
#genes_male_immune(annotLookupM$external_gene_name[1])
genes_male_immune(genes_females[1])
genes_male_immune(genes_females[2])
genes_male_immune(genes_females[3])
genes_male_immune(genes_females[4])
genes_male_immune(genes_females[5])
genes_male_immune(genes_females[6])#*
genes_male_immune(genes_females[7])
genes_male_immune(genes_females[8])
genes_male_immune(genes_females[9])
genes_male_immune(genes_females[10])
genes_male_immune(genes_females[11])
genes_male_immune(genes_females[12])
genes_male_immune(genes_females[13])
genes_male_immune(genes_females[14])
genes_male_immune(genes_females[15])
genes_male_immune(genes_females[16])
genes_male_immune(genes_females[17])
genes_male_immune(genes_females[18])
genes_male_immune(genes_females[19])
genes_male_immune(genes_females[20])
genes_male_immune(genes_females[21])
genes_male_immune(genes_females[22])
genes_male_immune(genes_females[23])
genes_male_immune(genes_females[24])
genes_male_immune(genes_females[25])
genes_male_immune(genes_females[26])
genes_male_immune(genes_females[27])
genes_male_immune(genes_females[28])
genes_females[6]

# do the same for cellular GO term
genes_male_cellular(genes_females[1])#*
genes_male_cellular(genes_females[2])#*
genes_male_cellular(genes_females[3])#*
genes_male_cellular(genes_females[4])
genes_male_cellular(genes_females[5])#*
genes_male_cellular(genes_females[6])#*
genes_male_cellular(genes_females[7])#*
genes_male_cellular(genes_females[8])#*
genes_male_cellular(genes_females[9])#*
genes_male_cellular(genes_females[10])#*
genes_male_cellular(genes_females[11])#*
genes_male_cellular(genes_females[12])#*
genes_male_cellular(genes_females[13])#*
genes_male_cellular(genes_females[14])#*
genes_male_cellular(genes_females[15])
genes_male_cellular(genes_females[16])#*
genes_male_cellular(genes_females[17])#*
genes_male_cellular(genes_females[18])
genes_male_cellular(genes_females[19])
genes_male_cellular(genes_females[20])#*
genes_male_cellular(genes_females[21])
genes_male_cellular(genes_females[22])#*
genes_male_cellular(genes_females[23])#*
genes_male_cellular(genes_females[24])#*
genes_male_cellular(genes_females[25])#*
genes_male_cellular(genes_females[26])#*
genes_male_cellular(genes_females[27])#*
genes_male_cellular(genes_females[28])#*
genes_females[c(1:3, 5:14, 16, 17, 20, 22:28)]

# do the same for localization GO term
genes_male_localization(genes_females[1])
genes_male_localization(genes_females[2])#*
genes_male_localization(genes_females[3])#*
genes_male_localization(genes_females[4])
genes_male_localization(genes_females[5])#*
genes_male_localization(genes_females[6])
genes_male_localization(genes_females[7])
genes_male_localization(genes_females[8])
genes_male_localization(genes_females[9])#*
genes_male_localization(genes_females[10])
genes_male_localization(genes_females[11])
genes_male_localization(genes_females[12])#*
genes_male_localization(genes_females[13])#*
genes_male_localization(genes_females[14])
genes_male_localization(genes_females[15])
genes_male_localization(genes_females[16])
genes_male_localization(genes_females[17])
genes_male_localization(genes_females[18])
genes_male_localization(genes_females[19])
genes_male_localization(genes_females[20])
genes_male_localization(genes_females[21])
genes_male_localization(genes_females[22])
genes_male_localization(genes_females[23])
genes_male_localization(genes_females[24])#*
genes_male_immune(genes_females[25])
genes_male_immune(genes_females[26])
genes_male_immune(genes_females[27])
genes_male_immune(genes_females[28])


genes_males_all <- annotLookupallM$external_gene_name[-which(annotLookupallM$external_gene_name == "") ] #remove NA blanks
# Obtain all immune system process GO insignificant results for each gene inputted 
# genes_male_immune(annotLookupM$external_gene_name[1])
genes_male_immune(genes_males[1])
genes_male_immune(genes_males[2])
genes_male_immune(genes_males[3])
genes_male_immune(genes_males[4])#*
genes_male_immune(genes_males[5])
genes_male_immune(genes_males[6])
genes_male_immune(genes_males[7])
genes_male_immune(genes_males[8])
genes_male_immune(genes_males[9])
genes_male_immune(genes_males[10])
genes_male_immune(genes_males[11])
genes_male_immune(genes_males[12])
genes_male_immune(genes_males[13])
genes_male_immune(genes_males[14])
genes_male_immune(genes_males[15])
genes_male_immune(genes_males[16])
genes_male_immune(genes_males[17])
genes_male_immune(genes_males[18])
genes_male_immune(genes_males[19])
genes_male_immune(genes_males[20])#*
genes_male_immune(genes_males[21])
genes_male_immune(genes_males[22])
genes_male_immune(genes_males[23])
genes_male_immune(genes_males[24])
genes_male_immune(genes_males[25])
genes_male_immune(genes_males[26])
genes_male_immune(genes_males[27])
genes_male_immune(genes_males[28])
genes_male_immune(genes_males[29])
genes_male_immune(genes_males[30])
genes_male_immune(genes_males[31])
genes_male_immune(genes_males[32])
genes_male_immune(genes_males[33])
genes_male_immune(genes_males[34])
genes_male_immune(genes_males[35])
genes_male_immune(genes_males[36])

# did not do other significantly different GO terms

# genes_GO_category <-function(gene){I1<-gost( #function creation
#   gene,
#   organism = "ggallus",
#   ordered_query = FALSE,
#   multi_query = FALSE,
#   significant = FALSE,
#   exclude_iea = FALSE,
#   measure_underrepresentation = FALSE,
#   evcodes = FALSE,
#   user_threshold = 0.05,
#   correction_method = "fdr",
#   domain_scope = "annotated",
#   custom_bg = NULL,
#   sources = "GO:BP",
#   as_short_link = FALSE)
# parent_I1<-get_parent_nodes(I1$result$term_id, term_df = NULL, graph_path_df = NULL, godir = NULL)
# slim_I1<-parent_I1 %>% filter(parent_go_id %in% "GO:0002376" )
# slim_I1
# }
#
# slim_males<-data.frame(child_go_id=character(0),parent_go_id=character(0), parent_name=character(0),distance=integer(0), gene=character(0))
# #slim_males<-matrix("double", 1, 5)
# for (i in 1:20){slim_I1<-genes_GO_category(genes_males[i])
#   slim_I1<-slim_I1 %>% mutate("gene"= rep(genes_males[i]))
#   slim_males<-rbind(slim_males,slim_I1)
#   }
# slim_males
# genes_GO_category(genes_males[20])
# genes_male_immune(genes_males[4])#*


#******************************************************************



#MF MALE overlap
gprofiler_Male_GO_MF <-gost(
  annotLookupM$external_gene_name,
  organism = "ggallus",
  ordered_query = FALSE,
  multi_query = FALSE,
  significant = FALSE,
  exclude_iea = FALSE,
  measure_underrepresentation = FALSE,
  evcodes = FALSE,
  user_threshold = 0.05,
  correction_method = "fdr",
  domain_scope = "annotated",
  custom_bg = NULL,
  sources = "GO:MF",
  as_short_link = FALSE)

# Slim down GO categories using parent nodes
length(gprofiler_Male_GO_MF$result$term_id)
parent_Male_MF<-get_parent_nodes(gprofiler_Male_GO_MF$result$term_id, term_df = NULL, graph_path_df = NULL, godir = NULL)
slim_Male_MF<-parent_Male_MF %>% filter(parent_go_id %in% c("GO:0045182", "GO:0140110", "GO:0060089", "GO:0005488", "GO:0005198", "GO:0098772", "GO:0003824", "GO:0005215", "GO:0038024" ,"GO:0016209","GO:0140104") ) #slim terms
# parent_Male_MF%>% filter(child_go_id %in% slim_Male_MF$child_go_id ==FALSE)
length(unique(slim_Male_MF$parent_name))
length(unique(slim_Male_MF$child_go_id))

#CC MALE overlap
gprofiler_Male_GO_CC <-gost(
  annotLookupM$external_gene_name,
  organism = "ggallus",
  ordered_query = FALSE,
  multi_query = FALSE,
  significant = FALSE,
  exclude_iea = FALSE,
  measure_underrepresentation = FALSE,
  evcodes = FALSE,
  user_threshold = 0.05,
  correction_method = "fdr",
  domain_scope = "annotated",
  custom_bg = NULL,
  sources = "GO:CC",
  as_short_link = FALSE)

# Slim down GO categories using parent nodes
length(gprofiler_Male_GO_CC$result$term_id)
parent_Male_CC<-get_parent_nodes(gprofiler_Male_GO_CC$result$term_id, term_df = NULL, graph_path_df = NULL, godir = NULL)
slim_Male_CC<-parent_Male_CC %>% filter(parent_go_id %in% c("GO:0044456", "GO:0044425", "GO:0016020", "GO:0045202", "GO:0044422", "GO:0044421", "GO:0030054", "GO:0031974", "GO:0032991", "GO:0099080", "GO:0005576", "GO:0005623", "GO:0044464", "GO:0043226", "GO:0005622", "GO:0005737", "GO:0110165") ) #slim terms
#parent_Male_CC%>% filter(child_go_id %in% slim_Male_CC$child_go_id ==FALSE)
length(unique(slim_Male_CC$parent_name))
length(unique(slim_Male_CC$child_go_id))

# REFORMAT DATA FOR GO PLOT DE
# Add male and female column
slim_Female_BP$sex<-"Female"
slim_Male_BP$sex<-"Male"
slim_Female_MF$sex<-"Female"
slim_Male_MF$sex<-"Male"
slim_Female_CC$sex<-"Female"
slim_Male_CC$sex<-"Male"

# Combine classification tables together and summarize counts and add percentages
BP<-as.data.frame(rbind(slim_Female_BP,slim_Male_BP)%>%
  dplyr::select(child_go_id, parent_go_id, parent_name, sex) %>%
  group_by(parent_go_id, parent_name, sex) %>%
  summarise("n"=n_distinct(child_go_id)))
BP<- BP%>%
  mutate("percentage_of_total"=ifelse(sex=="Female",(n/sum(BP%>%filter(sex=="Female") %>% dplyr::select(n))*100), (n/sum(BP%>%filter(sex=="Male") %>% dplyr::select(n))*100)))
MF<-as.data.frame(rbind(slim_Female_MF,slim_Male_MF) %>%
  dplyr::select(child_go_id, parent_go_id, parent_name, sex) %>%
  group_by(parent_go_id, parent_name, sex) %>%
  summarise("n"=n_distinct(child_go_id)))
MF<-MF%>%
  mutate("percentage_of_total"=ifelse(sex=="Female",(n/sum(MF%>%filter(sex=="Female") %>% dplyr::select(n))*100), (n/sum(MF%>%filter(sex=="Male") %>% dplyr::select(n))*100)))
CC<-as.data.frame(rbind(slim_Female_CC, slim_Male_CC) %>%
  dplyr::select(child_go_id, parent_go_id, parent_name, sex) %>%
  group_by(parent_go_id, parent_name, sex) %>%
  summarise("n"=n_distinct(child_go_id)))
CC<-CC%>%
  mutate("percentage_of_total"=ifelse(sex=="Female",(n/sum(CC%>%filter(sex=="Female") %>% dplyr::select(n))*100), (n/sum(CC%>%filter(sex=="Male") %>% dplyr::select(n))*100)))

names(BP)[names(BP) == "parent_name"] <- "category"
names(MF)[names(MF) == "parent_name"] <- "category"
names(CC)[names(CC) == "parent_name"] <- "category"

## X^2 TESTS
# Statistical similarity with X^2 tests and Bonferroni correction factor of 6 and 85 

BP_Pivot<- BP %>% pivot_wider(id_cols=category, #reformat dataframe
                              names_from= sex,
                              values_from= n)
BP_Pivot$Male[is.na(BP_Pivot$Male) ] <- 0
#a<-prop.test(x = as.matrix(BP_Pivot[,2:3]),correct=FALSE)[3]  ***
#  for( i in seq_along(1:nrow(BP_Pivot))){yi<-prop.test(x = as.matrix(BP_Pivot[i,2:3]),correct=FALSE)}  ***

MF_Pivot<- MF %>% pivot_wider(id_cols=category, #reformat dataframe
                              names_from= sex,
                              values_from= n)
MF_Pivot$Male[is.na(MF_Pivot$Male) ] <- 0
MF_Pivot$Female[is.na(MF_Pivot$Female) ] <- 0

CC_Pivot<- CC %>% pivot_wider(id_cols=category, #reformat dataframe
                              names_from= sex,
                              values_from= n)
CC_Pivot$Male[is.na(CC_Pivot$Male) ] <- 0
CC_Pivot$Female[is.na(CC_Pivot$Female) ] <- 0

# Create dataframe of p-values from statistical tests
prop_test_overlap<-data.frame("type"= c("BP", "MF", "CC", rep("BP", 19),rep("MF", 9),rep("CC", 11)),
                              "category"=c(rep("all", 3),BP_Pivot$category[1:19],MF_Pivot$category[1:9],CC_Pivot$category[1:11]),
                              "pvalue"=c(prop.test(x = as.matrix(BP_Pivot[2:3]),correct=FALSE)[3]$p.value,
                                         prop.test(x = as.matrix(MF_Pivot[,2:3]),correct=FALSE)[3]$p.value,
                                         prop.test(x = as.matrix(CC_Pivot[,2:3]),correct=FALSE)[3]$p.value,
                                         prop.test(x = as.matrix(BP_Pivot[1,2:3]),correct=FALSE)[3]$p.value,
                                         prop.test(x = as.matrix(BP_Pivot[2,2:3]),correct=FALSE)[3]$p.value,
                                         prop.test(x = as.matrix(BP_Pivot[3,2:3]),correct=FALSE)[3]$p.value,
                                         prop.test(x = as.matrix(BP_Pivot[4,2:3]),correct=FALSE)[3]$p.value,
                                         prop.test(x = as.matrix(BP_Pivot[5,2:3]),correct=FALSE)[3]$p.value,
                                         prop.test(x = as.matrix(BP_Pivot[6,2:3]),correct=FALSE)[3]$p.value,
                                         prop.test(x = as.matrix(BP_Pivot[7,2:3]),correct=FALSE)[3]$p.value,
                                         prop.test(x = as.matrix(BP_Pivot[8,2:3]),correct=FALSE)[3]$p.value,
                                         prop.test(x = as.matrix(BP_Pivot[9,2:3]),correct=FALSE)[3]$p.value,
                                         prop.test(x = as.matrix(BP_Pivot[10,2:3]),correct=FALSE)[3]$p.value,
                                         prop.test(x = as.matrix(BP_Pivot[11,2:3]),correct=FALSE)[3]$p.value,
                                         prop.test(x = as.matrix(BP_Pivot[12,2:3]),correct=FALSE)[3]$p.value,
                                         prop.test(x = as.matrix(BP_Pivot[13,2:3]),correct=FALSE)[3]$p.value,
                                         prop.test(x = as.matrix(BP_Pivot[14,2:3]),correct=FALSE)[3]$p.value,
                                         prop.test(x = as.matrix(BP_Pivot[15,2:3]),correct=FALSE)[3]$p.value,
                                         prop.test(x = as.matrix(BP_Pivot[16,2:3]),correct=FALSE)[3]$p.value,
                                         prop.test(x = as.matrix(BP_Pivot[17,2:3]),correct=FALSE)[3]$p.value,
                                         prop.test(x = as.matrix(BP_Pivot[18,2:3]),correct=FALSE)[3]$p.value,
                                         prop.test(x = as.matrix(BP_Pivot[19,2:3]),correct=FALSE)[3]$p.value,
                                         prop.test(x = as.matrix(MF_Pivot[1,2:3]),correct=FALSE)[3]$p.value,
                                         prop.test(x = as.matrix(MF_Pivot[2,2:3]),correct=FALSE)[3]$p.value,
                                         prop.test(x = as.matrix(MF_Pivot[3,2:3]),correct=FALSE)[3]$p.value,
                                         prop.test(x = as.matrix(MF_Pivot[4,2:3]),correct=FALSE)[3]$p.value,
                                         prop.test(x = as.matrix(MF_Pivot[5,2:3]),correct=FALSE)[3]$p.value,
                                         prop.test(x = as.matrix(MF_Pivot[6,2:3]),correct=FALSE)[3]$p.value,
                                         prop.test(x = as.matrix(MF_Pivot[7,2:3]),correct=FALSE)[3]$p.value,
                                         prop.test(x = as.matrix(MF_Pivot[8,2:3]),correct=FALSE)[3]$p.value,
                                         prop.test(x = as.matrix(MF_Pivot[9,2:3]),correct=FALSE)[3]$p.value,
                                         prop.test(x = as.matrix(CC_Pivot[1,2:3]),correct=FALSE)[3]$p.value,
                                         prop.test(x = as.matrix(CC_Pivot[2,2:3]),correct=FALSE)[3]$p.value,
                                         prop.test(x = as.matrix(CC_Pivot[3,2:3]),correct=FALSE)[3]$p.value,
                                         prop.test(x = as.matrix(CC_Pivot[4,2:3]),correct=FALSE)[3]$p.value,
                                         prop.test(x = as.matrix(CC_Pivot[5,2:3]),correct=FALSE)[3]$p.value,
                                         prop.test(x = as.matrix(CC_Pivot[6,2:3]),correct=FALSE)[3]$p.value,
                                         prop.test(x = as.matrix(CC_Pivot[7,2:3]),correct=FALSE)[3]$p.value,
                                         prop.test(x = as.matrix(CC_Pivot[8,2:3]),correct=FALSE)[3]$p.value,
                                         prop.test(x = as.matrix(CC_Pivot[9,2:3]),correct=FALSE)[3]$p.value,
                                         prop.test(x = as.matrix(CC_Pivot[10,2:3]),correct=FALSE)[3]$p.value,
                                         prop.test(x = as.matrix(CC_Pivot[11,2:3]),correct=FALSE)[3]$p.value)) %>%
  #Sidak padj=1-(1-p)^n 
  mutate( "Sidak"=c(1-(1-pvalue[1:3])^6,1-(1-pvalue[4:42])^83)) %>%
  mutate("reject"=ifelse(Sidak<0.05, "reject", "fail to reject")   )

BP_pvalue<-prop_test_overlap%>% filter(type== "BP", category!= "all") %>% dplyr::select(category, Sidak, reject)%>% full_join(BP, by= "category")%>%
  mutate("categoryPvalue"= ifelse(Sidak<0.9995, paste0(category, " (", formatC(Sidak, format = "e", digits =2), ")", ifelse(reject=="reject",paste0("*"), paste0(""))), paste0(category)) )
MF_pvalue<-prop_test_overlap%>% filter(type== "MF", category!= "all") %>% dplyr::select(category, Sidak, reject)%>% full_join(MF, by= "category")%>%
  mutate("categoryPvalue"= ifelse(Sidak<0.9995, paste0(category, " (", formatC(Sidak, format = "e",digits = 2), ")", ifelse(reject=="reject",paste0("*"), paste0(""))), paste0(category)) )
CC_pvalue<-prop_test_overlap%>% filter(type== "CC", category!= "all") %>% dplyr::select(category, Sidak, reject)%>% full_join(CC, by= "category")%>%
  mutate("categoryPvalue"= ifelse(Sidak<0.9995, paste0(category, " (", formatC(Sidak, format = "e",digits = 2), ")", ifelse(reject=="reject",paste0("*"), paste0(""))), paste0(category)) )


# BP_pvalue<-prop_test_overlap%>% filter(type== "BP", category!= "all") %>% dplyr::select(category, Bonferroni, reject)%>% full_join(BP, by= "category")%>%
#   mutate("categoryPvalue"= paste0(category, " (", signif(Bonferroni, 3), ")", ifelse(reject=="reject",paste0("*"), paste0(""))) )
# MF_pvalue<-prop_test_overlap%>% filter(type== "MF", category!= "all") %>% dplyr::select(category, Bonferroni, reject)%>% full_join(MF, by= "category")%>%
#   mutate("categoryPvalue"= paste0(category, " (", signif(Bonferroni, 3), ")", ifelse(reject=="reject",paste0("*"), paste0(""))) )
# CC_pvalue<-prop_test_overlap%>% filter(type== "CC", category!= "all") %>% dplyr::select(category, Bonferroni, reject)%>% full_join(CC, by= "category")%>%
#   mutate("categoryPvalue"= paste0(category, " (", signif(Bonferroni, 3), ")", ifelse(reject=="reject",paste0("*"), paste0(""))) )

# GO PLOTS DE OVERLAP padj < 0.01
# Plot bar plot of percent of gene hit against total number of process hits
BP_plot<-BP_pvalue%>%
  ggplot(aes(sex, percentage_of_total, fill=categoryPvalue))+
  geom_bar(position="stack", stat="identity") +
  ylab("Percentage of Total") +
  xlab("Sex") +
  ggtitle(paste0("Biological Processes (", formatC(prop_test_overlap$Sidak, format = "e", digits =2)[1], ")*")) +
  theme_bw()+
  theme(axis.text.x = element_text(size=31),
        axis.text.y = element_text(size=31),
        axis.title.x = element_text(size=33),
        axis.title.y = element_text(size=33),
        plot.title = element_text(size=38),
        legend.text = element_text(size=33),
        legend.title = element_text(size=33))+
  labs(fill=" Category (P-Value)")+
  guides(fill=guide_legend(ncol=2))

MF_plot<-MF_pvalue%>%
  ggplot(aes(sex, percentage_of_total, fill=categoryPvalue))+
  geom_bar(position="stack", stat="identity") +
  ylab("Percentage of Total") +
  xlab("Sex") +
  ggtitle(paste0("Molecular Functions (", formatC(prop_test_overlap$Sidak, format = "e", digits =2)[2], ")")) +
  theme_bw()+
  theme(axis.text.x = element_text(size=31),
        axis.text.y = element_text(size=31),
        axis.title.x = element_text(size=33),
        axis.title.y = element_text(size=33),
        plot.title = element_text(size=38),
        legend.text = element_text(size=33),
        legend.title = element_text(size=33))  +
  labs(fill="Category (P-Value)                                                                                                                                ")+
  guides(color = guide_legend(override.aes = list(size = 15)))

CC_plot<-CC_pvalue%>%
  ggplot(aes(sex, percentage_of_total, fill=categoryPvalue))+
  geom_bar(position="stack", stat="identity") +
  ylab("Percentage of Total") +
  xlab("Sex") +
  ggtitle(paste0("Cellular Components (", formatC(prop_test_overlap$Sidak, format = "e", digits =2)[3], ")")) +
  theme_bw()+
  theme(axis.text.x = element_text(size=31),
        axis.text.y = element_text(size=31),
        axis.title.x = element_text(size=33),
        axis.title.y = element_text(size=33),
        plot.title = element_text(size=38),
        legend.text = element_text(size=33),
        legend.title = element_text(size=33))  +
  labs(fill="Category (P-Value)                                                                                                                                ")+
  guides(color = guide_legend(override.aes = list(size = 15)))


# Arrange plots together
GO_figure <-ggarrange(BP_plot, MF_plot, CC_plot,
                      widths = rep(5),
                      hjust = -0.5,
                      vjust = 1,
                      ncol = 1, nrow = 3,
                      align = "hv")

outfile.GO<-"/Users/sunshine2.0/Desktop/WORK/output/GO_bar_gProfiler_overlap_0.01.png"
ggsave(outfile.GO, 
       dpi=300, 
       device="png",
       width=72.5, 
       height=58,
       units="cm")
#_____________________________********
# Obtain gene names for each GO term
female_genes_GO_BP <-gconvert(
  gprofiler_Female_GO_BP$result$term_id,
  organism = "ggallus",
  target = "AFFY_CHICKEN",
  numeric_ns = "",
  mthreshold = Inf,
  filter_na = TRUE)
names(female_genes_GO_BP)[2]<-paste("term_id") #change column name
female_genes_go_bp <- female_genes_GO_BP %>%  #add genes to results
  dplyr::select(term_id, name) %>%
  filter(name %in% annotLookupF$external_gene_name) %>% # genes only inputted into GO not all for that term
  unique() %>%
  group_by(term_id) %>%
  mutate("genes" = paste0(name, ", ", collapse = "")) %>%
  dplyr::select(term_id, genes) %>%
  unique()%>%
  full_join(gprofiler_Female_GO_BP$result, key= "term_id")

female_genes_GO_MF <-gconvert(
  gprofiler_Female_GO_MF$result$term_id,
  organism = "ggallus",
  target = "AFFY_CHICKEN",
  numeric_ns = "",
  mthreshold = Inf,
  filter_na = TRUE)
names(female_genes_GO_MF)[2]<-paste("term_id") #change column name
female_genes_go_mf <- female_genes_GO_MF %>%  #add genes to results
  dplyr::select(term_id, name) %>%
  filter(name %in% annotLookupF$external_gene_name) %>%
  unique() %>%
  group_by(term_id) %>%
  mutate("genes" = paste0(name, ", ", collapse = "")) %>%
  dplyr::select(term_id, genes) %>%
  unique()%>%
  full_join(gprofiler_Female_GO_MF$result, key= "term_id")

female_genes_GO_CC <-gconvert(
  gprofiler_Female_GO_CC$result$term_id,
  organism = "ggallus",
  target = "AFFY_CHICKEN",
  numeric_ns = "",
  mthreshold = Inf,
  filter_na = TRUE)
names(female_genes_GO_CC)[2]<-paste("term_id") #change column name
female_genes_go_cc <- female_genes_GO_CC %>%  #add genes to results
  dplyr::select(term_id, name) %>%
  filter(name %in% annotLookupF$external_gene_name) %>%
  unique() %>%
  group_by(term_id) %>%
  mutate("genes" = paste0(name, ", ", collapse = "")) %>%
  dplyr::select(term_id, genes) %>%
  unique()%>%
  full_join(gprofiler_Female_GO_CC$result, key= "term_id")

male_genes_GO_BP <-gconvert(
  gprofiler_Male_GO_BP$result$term_id,
  organism = "ggallus",
  target = "AFFY_CHICKEN",
  numeric_ns = "",
  mthreshold = Inf,
  filter_na = TRUE)
names(male_genes_GO_BP)[2]<-paste("term_id") #change column name
male_genes_go_bp <- male_genes_GO_BP %>%  #add genes to results
  dplyr::select(term_id, name) %>%
  filter(name %in% annotLookupM$external_gene_name) %>%
  unique() %>%
  group_by(term_id) %>%
  mutate("genes" = paste0(name, ", ", collapse = "")) %>%
  dplyr::select(term_id, genes) %>%
  unique()%>%
  full_join(gprofiler_Male_GO_BP$result, key= "term_id")

male_genes_GO_MF <-gconvert(
  gprofiler_Male_GO_MF$result$term_id,
  organism = "ggallus",
  target = "AFFY_CHICKEN",
  numeric_ns = "",
  mthreshold = Inf,
  filter_na = TRUE)
names(male_genes_GO_MF)[2]<-paste("term_id") #change column name
male_genes_go_mf <- male_genes_GO_MF %>%  #add genes to results
  dplyr::select(term_id, name) %>%
  filter(name %in% annotLookupM$external_gene_name) %>%
  unique() %>%
  group_by(term_id) %>%
  mutate("genes" = paste0(name, ", ", collapse = "")) %>%
  dplyr::select(term_id, genes) %>%
  unique()%>%
  full_join(gprofiler_Male_GO_MF$result, key= "term_id")

male_genes_GO_CC <-gconvert(
  gprofiler_Male_GO_CC$result$term_id,
  organism = "ggallus",
  target = "AFFY_CHICKEN",
  numeric_ns = "",
  mthreshold = Inf,
  filter_na = TRUE)
names(male_genes_GO_CC)[2]<-paste("term_id") #change column name
male_genes_go_cc <- male_genes_GO_CC %>%  #add genes to results
  dplyr::select(term_id, name) %>%
  filter(name %in% annotLookupM$external_gene_name) %>%
  unique() %>%
  group_by(term_id) %>%
  mutate("genes" = paste0(name, ", ", collapse = "")) %>%
  dplyr::select(term_id, genes) %>%
  unique()%>%
  full_join(gprofiler_Male_GO_CC$result, key= "term_id")

# Add monthly comparison of up regulation and down regulation to GO dataframe
names(Female_DE_genes)[2]<-paste("genes")

x<-str_split(female_genes_go_bp$genes, ", ", simplify = TRUE)
x[x==""]<-NA
y <- matrix("double", nrow(female_genes_go_bp),ncol(x))
for(i in seq_along(1:ncol(x))){ z=data.frame(matrix(ncol=1,nrow=nrow(female_genes_go_bp), dimnames=list(NULL, "genes")))
z$genes<-x[,i]
y[,i]<-left_join(z, Female_DE_genes, key=genes)[,3]}
y<-as.data.frame(y, stringsasfactors=FALSE)
y<-y%>%
  mutate_all(~str_replace_na(., "")) %>%
  transform("regulated.month_comparison" = paste0(V1, V2, V3, V4, V5, V6, V7, V8, V9, V10, V11, V12, V13, V14, V15, V16, V17, V18, V19, V20, V21, V22, V23, V24, V25)) %>%
  dplyr::select(regulated.month_comparison)
female_genes_go_bp$regulated.month_comparison<-y$regulated.month_comparison

x<-str_split(female_genes_go_mf$genes, ", ", simplify = TRUE)
x[x==""]<-NA
y <- matrix("double", nrow(female_genes_go_mf),ncol(x))
for(i in seq_along(1:ncol(x))){ z=data.frame(matrix(ncol=1,nrow=nrow(female_genes_go_mf), dimnames=list(NULL, "genes")))
z$genes<-x[,i]
y[,i]<-left_join(z, Female_DE_genes, key=genes)[,3]}
y<-as.data.frame(y, stringsasfactors=FALSE)
y<-y%>%
  mutate_all(~str_replace_na(., "")) %>%
  transform("regulated.month_comparison" = paste0(V1, V2, V3, V4, V5, V6, V7, V8, V9, V10, V11, V12, V13, V14, V15, V16, V17, V18, V19, V20, V21, V22, V23, V24, V25)) %>%
  dplyr::select(regulated.month_comparison)
female_genes_go_mf$regulated.month_comparison<-y$regulated.month_comparison

x<-str_split(female_genes_go_cc$genes, ", ", simplify = TRUE)
x[x==""]<-NA
y <- matrix("double", nrow(female_genes_go_cc),ncol(x))
for(i in seq_along(1:ncol(x))){ z=data.frame(matrix(ncol=1,nrow=nrow(female_genes_go_cc), dimnames=list(NULL, "genes")))
z$genes<-x[,i]
y[,i]<-left_join(z, Female_DE_genes, key=genes)[,3]}
y<-as.data.frame(y, stringsasfactors=FALSE)
y<-y%>%
  mutate_all(~str_replace_na(., "")) %>%
  transform("regulated.month_comparison" = paste0(V1, V2, V3, V4, V5, V6, V7, V8, V9, V10, V11, V12, V13, V14, V15, V16, V17, V18, V19, V20, V21, V22, V23, V24, V25)) %>%
  dplyr::select(regulated.month_comparison)
female_genes_go_cc$regulated.month_comparison<-y$regulated.month_comparison

names(Male_DE_genes)[2]<-paste("genes")

x<-str_split(male_genes_go_bp$genes, ", ", simplify = TRUE)
x[x==""]<-NA
y <- matrix("double", nrow(male_genes_go_bp),ncol(x))
for(i in seq_along(1:ncol(x))){ z=data.frame(matrix(ncol=1,nrow=nrow(male_genes_go_bp), dimnames=list(NULL, "genes")))
z$genes<-x[,i]
y[,i]<-left_join(z, Male_DE_genes, key=genes)[,3]}
y<-as.data.frame(y, stringsasfactors=FALSE)
y<-y%>%
  mutate_all(~str_replace_na(., "")) %>%
  transform("regulated.month_comparison" = paste0(V1, V2, V3, V4, V5, V6, V7, V8, V9, V10, V11, V12, V13, V14, V15, V16, V17, V18, V19, V20, V21, V22, V23, V24, V25)) %>%
  dplyr::select(regulated.month_comparison)
male_genes_go_bp$regulated.month_comparison<-y$regulated.month_comparison

x<-str_split(male_genes_go_mf$genes, ", ", simplify = TRUE)
x[x==""]<-NA
y <- matrix("double", nrow(male_genes_go_mf),ncol(x))
for(i in seq_along(1:ncol(x))){ z=data.frame(matrix(ncol=1,nrow=nrow(male_genes_go_mf), dimnames=list(NULL, "genes")))
z$genes<-x[,i]
y[,i]<-left_join(z, Male_DE_genes, key=genes)[,3]}
y<-as.data.frame(y, stringsasfactors=FALSE)
y<-y%>%
  mutate_all(~str_replace_na(., "")) %>%
  transform("regulated.month_comparison" = paste0(V1, V2, V3, V4, V5, V6, V7, V8, V9, V10, V11, V12, V13, V14, V15, V16, V17, V18, V19, V20, V21, V22, V23, V24, V25)) %>%
  dplyr::select(regulated.month_comparison)
male_genes_go_mf$regulated.month_comparison<-y$regulated.month_comparison

x<-str_split(male_genes_go_cc$genes, ", ", simplify = TRUE)
x[x==""]<-NA
y <- matrix("double", nrow(male_genes_go_cc),ncol(x))
for(i in seq_along(1:ncol(x))){ z=data.frame(matrix(ncol=1,nrow=nrow(male_genes_go_cc), dimnames=list(NULL, "genes")))
z$genes<-x[,i]
y[,i]<-left_join(z, Male_DE_genes, key=genes)[,3]}
y<-as.data.frame(y, stringsasfactors=FALSE)
y<-y%>%
  mutate_all(~str_replace_na(., "")) %>%
  transform("regulated.month_comparison" = paste0(V1, V2, V3, V4, V5, V6, V7, V8, V9, V10, V11, V12, V13, V14, V15, V16, V17, V18, V19, V20, V21, V22, V23, V24, V25)) %>%
  dplyr::select(regulated.month_comparison)
male_genes_go_cc$regulated.month_comparison<-y$regulated.month_comparison

# Export GO dataframe
sheetss<-list("female_BP"=female_genes_go_bp, "female_MF"=female_genes_go_mf, "female_CC"= female_genes_go_cc, "male_BP"=male_genes_go_bp, "male_MF"=male_genes_go_mf, "male_CC"= male_genes_go_cc)
write.xlsx(sheetss, "/Users/sunshine2.0/Desktop/WORK/output/GO_terms_emu.xlsx")
write.xlsx(genes_males, "/Users/sunshine2.0/Desktop/WORK/output/DE_males.xlsx")
write.xlsx(genes_females, "/Users/sunshine2.0/Desktop/WORK/output/DE_females.xlsx")


# PATHWAY ANALYSIS----------------------------------------------------------------------------------------------

#OBTAIN PA DATA ON ALL EXPRESSED GENES
#female
RLE$samples$samples<-rownames(RLE$samples) #filter for female only in metadata
RLE_F_meta<-RLE$samples %>%
  filter(gender=="F") 

RLEcounts <-as.data.frame(RLE$counts)
RLEcounts$gene<-rownames(RLE$counts) #filter for female only in counts
RLE_F_counts<-RLEcounts%>%
  dplyr::select(c(RLE_F_meta$samples, gene)) %>%
  filter_all(all_vars(. > 1)) #filter for counts >1
nrow(RLE$counts)
nrow(RLE_F_counts)

#male
RLE$samples$samples<-rownames(RLE$samples) #filter for male only in metadata
RLE_M_meta<-RLE$samples %>%
  filter(gender=="M") 

RLE_M_counts<-RLEcounts%>% #filter for male only in counts
  dplyr::select(c(RLE_M_meta$samples, gene)) %>%
  filter_all(all_vars(. > 1)) #filter for counts >1
nrow(RLE$counts)
nrow(RLE_M_counts)

# remove .n from gene ID
RLE_F_counts$gene = substr(RLE_F_counts$gene,1,nchar(RLE_F_counts$gene)-2) 
RLE_M_counts$gene = substr(RLE_M_counts$gene,1,nchar(RLE_M_counts$gene)-2) 
require('biomaRt')

#get gene symbol of emsembl ID's
# mart <- useMart('ENSEMBL_MART_ENSEMBL')
# mart <- useDataset('dnovaehollandiae_gene_ensembl', mart)
mart <- useEnsembl(biomart = "ensembl", dataset = "dnovaehollandiae_gene_ensembl", version = "102")

annotLookupallF <-getBM( 
  #tried to get uniprot annotation for emu ***emu had no or few annotation for uniprot*--> used chicken uniprot annotation on panther
  mart = mart,
  attributes = c(
    'ensembl_gene_id',
    'gene_biotype',
    'external_gene_name',
    'uniprot_gn_symbol',
    'uniprot_gn_id'),
  filter = 'ensembl_gene_id',
  values = RLE_F_counts$gene,
  uniqueRows = TRUE)
write.csv(annotLookupallF$external_gene_name, "/Users/sunshine2.0/Desktop/WORK/output/annotLookup$external_gene_name_all0.01F.csv")

annotLookupallM <-getBM(
  mart = mart,
  attributes = c(
    'ensembl_gene_id',
    'gene_biotype',
    'external_gene_name',
    'uniprot_gn_symbol',
    'uniprot_gn_id'),
  filter = 'ensembl_gene_id',
  values = RLE_M_counts$gene,
  uniqueRows = TRUE)
write.csv(annotLookupallM$external_gene_name, "/Users/sunshine2.0/Desktop/WORK/output/annotLookup$external_gene_name_all0.01M.csv")

#OBTAIN PA DATA ON DE GENES
gprofiler_Female_all <-gost(
  annotLookupallF$external_gene_name,
  organism = "ggallus",
  ordered_query = FALSE,
  multi_query = FALSE,
  significant = FALSE,
  exclude_iea = FALSE,
  measure_underrepresentation = FALSE,
  evcodes = FALSE,
  user_threshold = 0.05,
  correction_method = "fdr",
  domain_scope = "annotated",
  custom_bg = NULL,
  sources = "REAC",
  as_short_link = FALSE)

gprofiler_Male_all <-gost(
  annotLookupallM$external_gene_name,
  organism = "ggallus",
  ordered_query = FALSE,
  multi_query = FALSE,
  significant = FALSE,
  exclude_iea = FALSE,
  measure_underrepresentation = FALSE,
  evcodes = FALSE,
  user_threshold = 0.05,
  correction_method = "fdr",
  domain_scope = "annotated",
  custom_bg = NULL,
  sources = "REAC",
  as_short_link = FALSE)

# Tidy result data
gprofiler_Female_all_result <-gprofiler_Female_all$result %>% 
  dplyr::select(term_name, intersection_size, term_size, term_id) %>%
  mutate("sex"="Female")
names(gprofiler_Female_all_result)[names(gprofiler_Female_all_result) == "intersection_size"] <- "intersection_size_all" #rename intersection_size and term_size with all addition
names(gprofiler_Female_all_result)[names(gprofiler_Female_all_result) == "term_size"] <- "term_size_all"
nrow(gprofiler_Female_all_result)

gprofiler_Male_all_result <-gprofiler_Male_all$result %>% 
  dplyr::select(term_name, intersection_size, term_size, term_id) %>%
  mutate("sex"="Male")
names(gprofiler_Male_all_result)[names(gprofiler_Male_all_result) == "intersection_size"] <- "intersection_size_all" #rename intersection_size and term_size with all addition
names(gprofiler_Male_all_result)[names(gprofiler_Male_all_result) == "term_size"] <- "term_size_all"
nrow(gprofiler_Male_all_result)

gprofiler_sorted_all_plot <-rbind(gprofiler_Female_all_result, gprofiler_Male_all_result) 

library("gprofiler2")
# Gene list functional enrichment analysis
gprofiler_Female <-gost(
  annotLookupF$external_gene_name,
  organism = "ggallus",
  ordered_query = FALSE,
  multi_query = FALSE,
  significant = TRUE,
  exclude_iea = FALSE,
  measure_underrepresentation = FALSE,
  evcodes = FALSE,
  user_threshold = 0.05,
  correction_method = "fdr",
  domain_scope = "annotated",
  custom_bg = NULL,
  sources = "REAC",
  as_short_link = FALSE)

gprofiler_Female_insig <-gost( #insignificant
  annotLookupF$external_gene_name,
  organism = "ggallus",
  ordered_query = FALSE,
  multi_query = FALSE,
  significant = FALSE,
  exclude_iea = FALSE,
  measure_underrepresentation = FALSE,
  evcodes = FALSE,
  user_threshold = 0.05,
  correction_method = "fdr",
  domain_scope = "annotated",
  custom_bg = NULL,
  sources = "REAC",
  as_short_link = FALSE)

gprofiler_Male <-gost(
  annotLookupM$external_gene_name,
  organism = "ggallus",
  ordered_query = FALSE,
  multi_query = FALSE,
  significant = TRUE,
  exclude_iea = FALSE,
  measure_underrepresentation = FALSE,
  evcodes = FALSE,
  user_threshold = 0.05,
  correction_method = "fdr",
  domain_scope = "annotated",
  custom_bg = NULL,
  sources = "REAC",
  as_short_link = FALSE)

gprofiler_Male_insig <-gost( #insignificant
  annotLookupM$external_gene_name,
  organism = "ggallus",
  ordered_query = FALSE,
  multi_query = FALSE,
  significant = FALSE,
  exclude_iea = FALSE,
  measure_underrepresentation = FALSE,
  evcodes = FALSE,
  user_threshold = 0.05,
  correction_method = "fdr",
  domain_scope = "annotated",
  custom_bg = NULL,
  sources = "REAC",
  as_short_link = FALSE)

# Obtain gene names for each pathway
female_genes_pathway <-gconvert(
  gprofiler_Female$result$term_id,
  organism = "ggallus",
  target = "REACTOME",
  numeric_ns = "",
  mthreshold = Inf,
  filter_na = TRUE
)

male_genes_pathway <-gconvert(
  gprofiler_Male$result$term_id,
  organism = "ggallus",
  target = "REACTOME",
  numeric_ns = "",
  mthreshold = Inf,
  filter_na = TRUE
)
names(female_genes_pathway)[2]<-paste("term_id") #change column name
names(male_genes_pathway)[2]<-paste("term_id") #change column name


###NEW-----------------
#top 5 pathway genes FEMALES
max_PA_F<-gprofiler_Female$result %>%
  filter(term_size>5)
max_PA_F<- max_PA_F[order(max_PA_F$p_value),]
top_5_PA_F<-max_PA_F[1:6,] #triglyceride metabolism is also in males

#1: Plasma lipoprotein assembly, remodeling, and clearance    
plasma_arc1<-female_genes_pathway%>%
  filter(term_id =="REAC:R-GGA-174824")
plasma_arc2<-unique(plasma_arc1$name)%>%
  as.data.frame()
colnames(plasma_arc2)<-"genes"
plasma_arc2<-plasma_arc2%>%
  mutate("pathway"=("REAC:R-GGA-174824"))%>%
  mutate("description"="Plasma lipoprotein assembly, remodeling, and clearance")%>%
  mutate("number_genes"=length(unique(plasma_arc1$name)))

#2: Assembly of active LPL and LIPC lipase complexes   
assembly1<-female_genes_pathway%>%
  filter(term_id =="REAC:R-GGA-8963889")
assembly2<-unique(assembly1$name)%>%
  as.data.frame()
colnames(assembly2)<-"genes"
assembly2<-assembly2%>%
  mutate("pathway"=("REAC:R-GGA-8963889"))%>%
  mutate("description"="Assembly of active LPL and LIPC lipase complexes ")%>%
  mutate("number_genes"=length(unique(assembly1$name)))

#3: Plasma lipoprotein remodeling   
plasma_lr1<-female_genes_pathway%>%
  filter(term_id =="REAC:R-GGA-8963899")
plasma_lr2<-unique(plasma_lr1$name)%>%
  as.data.frame()
colnames(plasma_lr2)<-"genes"
plasma_lr2<-plasma_lr2%>%
  mutate("pathway"=("REAC:R-GGA-8963899"))%>%
  mutate("description"="Plasma lipoprotein remodeling")%>%
  mutate("number_genes"=length(unique(plasma_lr1$name)))

#5: Carbohydrate metabolism
carb1<-female_genes_pathway%>%
  filter(term_id =="REAC:R-GGA-353098")
carb2<-unique(carb1$name)%>%
  as.data.frame()
colnames(carb2)<-"genes"
carb2<-carb2%>%
  mutate("pathway"=("REAC:R-GGA-353098"))%>%
  mutate("description"="Carbohydrate metabolism")%>%
  mutate("number_genes"=length(unique(carb1$name)))

#6: Glycolysis   
glycolysis1<-female_genes_pathway%>%
  filter(term_id =="REAC:R-GGA-352882")
glycolysis2<-unique(glycolysis1$name)%>%
  as.data.frame()
colnames(glycolysis2)<-"genes"
glycolysis2<-glycolysis2%>%
  mutate("pathway"=("REAC:R-GGA-352882"))%>%
  mutate("description"="Glycolysis")%>%
  mutate("number_genes"=length(unique(glycolysis1$name)))

#combine dataframes
F_PA_genes_top_5<-rbind(plasma_arc2, assembly2, plasma_lr2, carb2, glycolysis2)%>%
  mutate("sex"="female")


#top 5 pathway genes MALES
max_PA_M<-gprofiler_Male$result %>%
  filter(term_size>5)
max_PA_M<- max_PA_M[order(max_PA_M$p_value),]
top_5_PA_M<-max_PA_M[1:5,]

#1: Gluconeogenesis  
gluco1<-male_genes_pathway%>%
  filter(term_id =="REAC:R-GGA-70263")
gluco2<-unique(gluco1$name)%>%
  as.data.frame()
colnames(gluco2)<-"genes"
gluco2<-gluco2%>%
  mutate("pathway"=("REAC:R-GGA-70263"))%>%
  mutate("description"="Gluconeogenesis")%>%
  mutate("number_genes"=length(unique(gluco1$name)))

#2: Integrin cell surface interactions  
integrin1<-male_genes_pathway%>%
  filter(term_id =="REAC:R-GGA-216083")
integrin2<-unique(integrin1$name)%>%
  as.data.frame()
colnames(integrin2)<-"genes"
integrin2<-integrin2%>%
  mutate("pathway"=("REAC:R-GGA-216083"))%>%
  mutate("description"="Integrin cell surface interactions")%>%
  mutate("number_genes"=length(unique(integrin1$name)))

#3:  Keratan sulfate degradation  
keratan_d1<-male_genes_pathway%>%
  filter(term_id =="REAC:R-GGA-2022857")
keratan_d2<-unique(keratan_d1$name)%>%
  as.data.frame()
colnames(keratan_d2)<-"genes"
keratan_d2<-keratan_d2%>%
  mutate("pathway"=("REAC:R-GGA-2022857"))%>%
  mutate("description"="Keratan sulfate degradation")%>%
  mutate("number_genes"=length(unique(keratan_d1$name)))

#4:  Triglyceride metabolism 
tag_m1<-male_genes_pathway%>%
  filter(term_id =="REAC:R-GGA-8979227")
tag_m2<-unique(tag_m1$name)%>%
  as.data.frame()
colnames(tag_m2)<-"genes"
tag_m2<-tag_m2%>%
  mutate("pathway"=("REAC:R-GGA-8979227"))%>%
  mutate("description"="Triglyceride metabolism")%>%
  mutate("number_genes"=length(unique(tag_m1$name)))

#5:  Keratan sulfate biosynthesis 
keratan_b1<-male_genes_pathway%>%
  filter(term_id =="REAC:R-GGA-2022854")
keratan_b2<-unique(keratan_b1$name)%>%
  as.data.frame()
colnames(keratan_b2)<-"genes"
keratan_b2<-keratan_b2%>%
  mutate("pathway"=("REAC:R-GGA-2022854"))%>%
  mutate("description"="Keratan sulfate biosynthesis")%>%
  mutate("number_genes"=length(unique(keratan_b1$name)))

#combine dataframes
M_PA_genes_top_5<-rbind(gluco2, integrin2, keratan_d2, tag_m2, keratan_b2)%>%
  mutate("sex"="male")
PA_genes_top_5<-rbind(F_PA_genes_top_5, M_PA_genes_top_5)
PA_genes_top_5_name<-PA_genes_top_5%>%
  filter(!grepl("ENSGALG",genes))
PA_genes_top_5_id<-PA_genes_top_5%>%
  filter(grepl("ENSGALG",genes))

nrow(PA_genes_top_5)
nrow(PA_genes_top_5_name)
nrow(PA_genes_top_5_id)


#ARE GGALLUS PA GENES IN EMU
#ggallus id to entrez
#a<-getAttributes(dataset="ggallus_gene_ensembl", mart = 'ENSEMBL_MART_ENSEMBL')
mart <- useEnsembl(biomart = "ensembl", dataset = "ggallus_gene_ensembl", version = "102")
PA_genes_top_5_ggallus_id <-getBM(
  mart = mart,
  attributes = c(
    'external_gene_name',
    'ensembl_gene_id',
    'entrezgene_id'),
  filter = 'ensembl_gene_id',
  values = PA_genes_top_5_id$gene,
  uniqueRows = TRUE)

#ggallus name to entrez
mart <- useEnsembl(biomart = "ensembl", dataset = "ggallus_gene_ensembl", version = "102")
PA_genes_top_5_ggallus_name <-getBM(
  mart = mart,
  attributes = c(
    'external_gene_name',
    'ensembl_gene_id',
    'entrezgene_id'),
  filter = 'external_gene_name',
  values = PA_genes_top_5_name$gene,
  uniqueRows = TRUE)

PA_genes_top_5_ggallus<-rbind(PA_genes_top_5_ggallus_id,PA_genes_top_5_ggallus_name)

#all emu gene ids
mart <- useEnsembl(biomart = "ensembl", dataset = "ggallus_gene_ensembl", version = "102")
#a<-getAttributes(dataset="ggallus_gene_ensembl", mart = 'ENSEMBL_MART_ENSEMBL')
PA_genes_top_5_emu <-getBM(# no entrez in dataframe
  mart = mart,
  attributes = c(
    'ensembl_gene_id',
    'external_gene_name',
    'dnovaehollandiae_homolog_ensembl_gene'),
  filter = 'ensembl_gene_id',
  values = PA_genes_top_5_ggallus$ensembl_gene_id,
  uniqueRows = TRUE)

#all emu gene names
mart <- useEnsembl(biomart = "ensembl", dataset = "dnovaehollandiae_gene_ensembl", version = "102")
PA_genes_top_5_emu_names <-getBM(
  mart = mart,
  attributes = c(
    'ensembl_gene_id',
    'external_gene_name',
    'entrezgene_id'),
  filter = 'ensembl_gene_id',
  values = PA_genes_top_5_emu$dnovaehollandiae_homolog_ensembl_gene,
  uniqueRows = TRUE)
colnames(PA_genes_top_5_emu_names)<-c("dnovaehollandiae_homolog_ensembl_gene", "dnovaehollandiae_external_gene_name", 'dnovaehollandiae_entrezgene_id')

#join id dataframes
PA_genes_top_5_emu_join<-full_join(PA_genes_top_5_emu,PA_genes_top_5_emu_names, by="dnovaehollandiae_homolog_ensembl_gene")
PA_genes_top_5_join<-full_join(PA_genes_top_5_ggallus, PA_genes_top_5_emu_join,  by="ensembl_gene_id")

#join PA dataframes
colnames(PA_genes_top_5_id)[1]<-"ensembl_gene_id"
PA_genes_top_5_id_full<-full_join(PA_genes_top_5_id, PA_genes_top_5_ggallus_id, by="ensembl_gene_id")%>%
  dplyr::select(sex, pathway, description, number_genes, ensembl_gene_id )

colnames(PA_genes_top_5_name)[1]<-"external_gene_name"
PA_genes_top_5_name_full<-full_join(PA_genes_top_5_name, PA_genes_top_5_ggallus_name, by="external_gene_name")%>%
  dplyr::select(sex, pathway, description, number_genes, ensembl_gene_id )

PA_genes_top_5_full<-rbind(PA_genes_top_5_id_full, PA_genes_top_5_name_full)

#join PA and id dataframes
PA_genes_top_5_conjoined<-full_join(PA_genes_top_5_full,PA_genes_top_5_join, by="ensembl_gene_id")
nrow(PA_genes_top_5_full)
nrow(PA_genes_top_5_join)
nrow(PA_genes_top_5_conjoined)

write.xlsx(PA_genes_top_5_conjoined, "/Users/sunshine2.0/Desktop/WORK/output/Top_5_pathways_genes.xlsx")

# PA_genes_top_5_conjoined%>%
#   filter(pathway=="REAC:R-GGA-70263")%>%
#   nrow()
# 
# unique(PA_genes_top_5_conjoined$pathway)

#DUPLICATIONS: 2 emu id for one chicken id (4x , should be 138 but are 142):
# ENSGALG00000044187 (3x) for Plasma lipoprotein assembly, remodeling, and clearance &
 #Assembly of active LPL and LIPC lipase complexes & Plasma lipoprotein remodeling 
# ENSGALG00000030034 for Gluconeogenesis
#

# to explain why we cant see them in data in not observed or that the total is completely wrong because they arent even there? 
#but they could be but just arent identified yet...


###NEW---------------

# Remove pathways with less than 5 genes in it
min_PA<-gprofiler_Female$result %>%
  filter(term_size<6)

# FIND NUMBER OF APPEARANCES A GENE MAKES IN MULTIPLE PATHWAYS 
female_genes_PA <- female_genes_pathway %>% #find common #
  dplyr::select(c(term_id, name, description) ) %>%
  filter(!term_id %in% min_PA$term_id) %>% #remove pathways with less than 5 genes in it
  unique() %>%
  arrange(name) %>%
  group_by(name) %>%
  mutate("appearances_in_name" = n()) %>%
  arrange(desc(appearances_in_name)) %>%
  left_join(gprofiler_Female$result%>% dplyr::select(c(term_id, term_name)), key= term_id) %>%
  mutate("DE_overlap"= ifelse(name%in% annotLookupF$external_gene_name, paste0( "yes"), paste0("no")))

male_genes_PA <- male_genes_pathway %>% #find common #
  dplyr::select(c(term_id, name, description) ) %>%
  unique() %>%
  arrange(name) %>%
  group_by(name) %>%
  mutate("appearances_in_name" = n()) %>%
  arrange(desc(appearances_in_name)) %>%
  full_join(gprofiler_Male$result%>% dplyr::select(c(term_id, term_name)), key= term_id) %>%
  mutate("DE_overlap"= ifelse(name%in% annotLookupM$external_gene_name, paste0( "yes"), paste0("no")))

# Add up regulation and down regulation to PA sheet
names(Female_DE_genes)[2]<-"name"
female_genes_PA<-left_join(female_genes_PA, Female_DE_genes, by="name")
names(Male_DE_genes)[2]<-"name"
male_genes_PA<-left_join(male_genes_PA, Male_DE_genes, by="name")

# Export dataframe of number of appearances a genes makes in multiple pathways
sheets<-list("female_PA"=female_genes_PA, "male_PA"=male_genes_PA)
write.xlsx(sheets, "/Users/sunshine2.0/Desktop/WORK/output/Pathway_Analysis_common_genes_min.xlsx")

# Plot PA distribution of genes in multiple pathways
PA_appearances_count_F<-female_genes_PA %>%
  filter(DE_overlap == "yes") %>%
  dplyr::select(appearances_in_name, name)%>%
  unique() %>%
  group_by( appearances_in_name) %>%
  mutate("count"= n()) %>%
  dplyr::select(appearances_in_name, count) %>%
  unique() %>%
  mutate("sex"="Female")

PA_appearances_count_M<-male_genes_PA %>%
  filter(DE_overlap == "yes") %>%
  dplyr::select(appearances_in_name, name)%>%
  unique() %>%
  group_by( appearances_in_name) %>%
  mutate("count"= n()) %>%
  dplyr::select(appearances_in_name, count) %>%
  unique() %>%
  mutate("sex"="Male")

ggplot(full_join(PA_appearances_count_F, PA_appearances_count_M), aes(appearances_in_name, count, fill=sex))+
         geom_col( position =position_dodge2(width = 0.9, preserve = "single", padding = 0)) +
  theme_bw(30)+
  xlab("Number of Pathways") +
  labs(fill='Sex' ) +
  ylab("Number of Genes") 
  
ggsave("/Users/sunshine2.0/Desktop/WORK/output/PA_appearances_count.png",
       dpi=300,
       device="png",
       width=20,
       height=20,
       units="cm")


# Tidy result data for PA plot
gprofiler_Female_result <-gprofiler_Female$result %>%
  dplyr::select(term_name, intersection_size, term_size, p_value, term_id) %>%
  filter(!term_id %in% min_PA$term_id) %>%
  mutate("sex"="Female")
nrow(gprofiler_Female_result)

gprofiler_Female_result_insig <-gprofiler_Female_insig$result %>%
  dplyr::select(term_name, intersection_size, term_size, p_value, term_id) %>%
  mutate("sex"="Female")
nrow(gprofiler_Female_result_insig)

gprofiler_Male_result <-gprofiler_Male$result %>%
  dplyr::select(term_name, intersection_size, term_size, p_value, term_id) %>%
  mutate("sex"="Male")
nrow(gprofiler_Male_result)

gprofiler_Male_result_insig <-gprofiler_Male_insig$result %>%
  dplyr::select(term_name, intersection_size, term_size, p_value, term_id) %>%
  mutate("sex"="Male")
nrow(gprofiler_Male_result_insig)

#similarity in term_names between sexes
same_pathways_sexes <-gprofiler_Female_result%>%
  filter(gprofiler_Female_result$term_id %in% gprofiler_Male_result$term_id)%>%
  dplyr::select(term_name)


# MAKE PA PLOTS

#FEMALE_______________

# Add M and F datatables by row and filter only for those top pathways in females
gprofiler_Female_sorted_DE_plot <-rbind(gprofiler_Female_result, gprofiler_Male_result_insig) %>%
  filter(term_id %in% gprofiler_Female_result$term_id)
nrow(gprofiler_Female_sorted_DE_plot)

# Add rows for male data not in top female DE overlap pathways
gprofiler_Female_sorted_Male_extra<-gprofiler_Female_sorted_DE_plot %>%
  filter(gprofiler_Female_sorted_DE_plot$term_id %in% gprofiler_Male_result_insig$term_id == FALSE) %>%
  mutate(intersection_size= 0)%>%
  mutate(p_value= NA)%>%
  mutate(sex="Male")

# Combine row wise female top pathways with male extra pathways not shown significantly in DE
gprofiler_Female_sorted_DE_complete_plot<- rbind(gprofiler_Female_sorted_DE_plot, gprofiler_Female_sorted_Male_extra)

# Join DE overlap gprofiler results with all gprofiler results
gprofiler_Female_sorted_plot<-inner_join(gprofiler_Female_sorted_DE_complete_plot, gprofiler_sorted_all_plot, by=c("term_id", "term_name", "sex"))%>%
  mutate("Overlap_percent"= intersection_size/term_size*100) %>%  #add column for intersection_size %
  mutate("Not_Overlap"=intersection_size_all-intersection_size) %>% #add column for not overlap genes in pathway
  mutate("Not_Overlap_percent"=Not_Overlap/term_size*100) %>%
  mutate("Not_Observed"=term_size-Not_Overlap-intersection_size) %>%   #add column for not observed
  mutate("Not_Observed_percent"=Not_Observed/term_size*100) %>%
  mutate("totalP"=Overlap_percent+Not_Overlap_percent+Not_Observed_percent) %>%
  arrange((Overlap_percent))%>%
  arrange(sex) # arrange in order of top pathways
nrow(gprofiler_Female_sorted_plot)

gprofiler_Female_sorted_plot$term_name <- factor(gprofiler_Female_sorted_plot$term_name, levels=unique(gprofiler_Female_sorted_plot$term_name))


# Arrange percents for plotting
gprofiler_Female_sorted_plot_ready<-gprofiler_Female_sorted_plot %>%
  pivot_longer(cols      = c(Overlap_percent,Not_Overlap_percent,Not_Observed_percent),
               names_to  = "type",
               values_to = "percent")
gprofiler_Female_sorted_plot_ready$type <- factor(gprofiler_Female_sorted_plot_ready$type, levels=c("Not_Observed_percent","Not_Overlap_percent","Overlap_percent"))
gprofiler_Female_sorted_plot_ready$term_name <- factor(gprofiler_Female_sorted_plot_ready$term_name, levels=unique(gprofiler_Female_sorted_plot_ready$term_name))


# Significance X^2 testing with Sidak correction factor of 19
PA_female_pvalue<-data.frame("term_name"= gprofiler_Female_sorted_plot$term_name[1:10],
                           "pvalue"=c(
                             prop.test(x = gprofiler_Female_sorted_plot %>%
                                         dplyr::select(c(term_name, sex, Overlap_percent,Not_Overlap_percent,Not_Observed_percent)) %>%
                                         filter(term_name=="Plasma lipoprotein assembly, remodeling, and clearance") %>%
                                         dplyr::select(-sex, -term_name)%>%
                                         t(),correct=FALSE)[3]$p.value,
                             prop.test(x = gprofiler_Female_sorted_plot %>%
                                         dplyr::select(c(term_name, sex, Overlap_percent,Not_Overlap_percent,Not_Observed_percent)) %>%
                                         filter(term_name=="Glutathione synthesis and recycling") %>%
                                         dplyr::select(-sex, -term_name)%>%
                                         t(),correct=FALSE)[3]$p.value,
                             prop.test(x = gprofiler_Female_sorted_plot %>%
                                         dplyr::select(c(term_name, sex, Overlap_percent,Not_Overlap_percent,Not_Observed_percent)) %>%
                                         filter(term_name=="IL-6-type cytokine receptor ligand interactions") %>%
                                         dplyr::select(-sex, -term_name)%>%
                                         t(),correct=FALSE)[3]$p.value,
                             prop.test(x = gprofiler_Female_sorted_plot %>%
                                         dplyr::select(c(term_name, sex, Overlap_percent,Not_Overlap_percent)) %>%
                                         filter(term_name=="Triglyceride catabolism") %>%
                                         dplyr::select(-sex, -term_name)%>%
                                         t(),correct=FALSE)[3]$p.value,
                             prop.test(x = gprofiler_Female_sorted_plot %>%
                                         dplyr::select(c(term_name, sex, Overlap_percent,Not_Overlap_percent,Not_Observed_percent)) %>%
                                         filter(term_name=="Triglyceride metabolism") %>%
                                         dplyr::select(-sex, -term_name)%>%
                                         t(),correct=FALSE)[3]$p.value,
                             prop.test(x = gprofiler_Female_sorted_plot %>%
                                         dplyr::select(c(term_name, sex, Overlap_percent,Not_Overlap_percent,Not_Observed_percent)) %>%
                                         filter(term_name=="Triglyceride biosynthesis") %>%
                                         dplyr::select(-sex, -term_name)%>%
                                         t(),correct=FALSE)[3]$p.value,
                             prop.test(x = gprofiler_Female_sorted_plot %>%
                                         dplyr::select(c(term_name, sex, Overlap_percent,Not_Overlap_percent,Not_Observed_percent)) %>%
                                         filter(term_name=="Plasma lipoprotein remodeling") %>%
                                         dplyr::select(-sex, -term_name)%>%
                                         t(),correct=FALSE)[3]$p.value,
                             prop.test(x = gprofiler_Female_sorted_plot %>%
                                         dplyr::select(c(term_name, sex, Overlap_percent,Not_Overlap_percent,Not_Observed_percent)) %>%
                                         filter(term_name=="Carbohydrate metabolism") %>%
                                         dplyr::select(-sex, -term_name)%>%
                                         t(),correct=FALSE)[3]$p.value,
                             prop.test(x = gprofiler_Female_sorted_plot %>%
                                         dplyr::select(c(term_name, sex, Overlap_percent,Not_Overlap_percent,Not_Observed_percent)) %>%
                                         filter(term_name=="Glycolysis") %>%
                                         dplyr::select(-sex, -term_name)%>%
                                         t(),correct=FALSE)[3]$p.value,
                             prop.test(x = gprofiler_Female_sorted_plot %>%
                                         dplyr::select(c(term_name, sex, Overlap_percent,Not_Overlap_percent,Not_Observed_percent)) %>%
                                         filter(term_name=="Assembly of active LPL and LIPC lipase complexes") %>%
                                         dplyr::select(-sex, -term_name)%>%
                                         t(),correct=FALSE)[3]$p.value))
                            
#Sidak padj=1-(1-p)^n 
PA_female_pvalue<-PA_female_pvalue%>% mutate( "Sidak"=1-(1-pvalue)^19) %>%
  mutate("reject"=ifelse(Sidak<0.05, "reject", "fail to reject"))

gprofiler_Female_sorted_plot_ready_pvalue<-PA_female_pvalue %>% full_join(gprofiler_Female_sorted_plot_ready, by= "term_name")%>%
  mutate("term_namePvalue"= ifelse(Sidak!=1, paste0(term_name, " (", formatC(Sidak, format = "e", digits =2), ")", ifelse(reject=="reject",paste0("*"), paste0("")) ), paste0(term_name)))



# Add '#' for same pathways for both sexes
gprofiler_Female_sorted_plot_ready_pvalue$term_namePvalue<-
  ifelse(gprofiler_Female_sorted_plot_ready_pvalue$term_name %in% same_pathways_sexes$term_name, paste0("# ", gprofiler_Female_sorted_plot_ready_pvalue$term_namePvalue), paste0("", gprofiler_Female_sorted_plot_ready_pvalue$term_namePvalue))
gprofiler_Female_sorted_plot_ready_pvalue$term_namePvalue <- factor(gprofiler_Female_sorted_plot_ready_pvalue$term_namePvalue, levels=unique(gprofiler_Female_sorted_plot_ready_pvalue$term_namePvalue))


# # Overlap f vs M all type
# prop.test(x = as.matrix(gprofiler_Female_sorted_plot %>%
#                           dplyr::select(c( term_name, sex, Overlap_percent)) %>%
#                           pivot_wider(id_cols=term_name,
#                                       names_from=sex,
#                                       values_from=Overlap_percent)%>%
#                           dplyr::select( Male, Female))
#           ,correct=FALSE)[3]$p.value


# PLOT DE FEMALE ENRICHED PATHWAYS 
#https://stackoverflow.com/questions/18265941/two-horizontal-bar-charts-with-shared-axis-in-ggplot2-similar-to-population-pyr
library(grid)
library('ggtext')
g.mid<-ggplot(gprofiler_Female_sorted_plot_ready_pvalue,aes(x=1,y=term_namePvalue))+
  geom_text(aes(label=term_namePvalue), size=7.5)+
  geom_segment(aes(x=0.94,xend=0.943,yend=term_namePvalue))+
  geom_segment(aes(x=1.057,xend=1.06,yend=term_namePvalue))+
  geom_text(aes(x=1.054, label=term_size), size=7.5)+
  scale_x_continuous(expand=c(0,0),limits=c(0.94,1.06))+
  ylab(NULL)+
  xlab(" ")+
  ggtitle(" ")+
  theme(axis.title.y=element_blank(),
        panel.grid=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        panel.background=element_blank(),
        axis.text.x=element_text(color=NA),
        axis.ticks.x=element_line(color=NA),
        axis.text = element_text(size=24),
        axis.title.x= element_text(size=25),
        plot.title = element_text(size = 27, face = "bold"),
        plot.margin = unit(c(1,-1,1,-1), "mm"))+
  coord_cartesian(clip = "off") + #keep text from being cut off
  annotation_custom(grid::textGrob("Pathway (P-Value)", gp = gpar(fontsize = 20, fontface="bold")), 
                    xmin = .94,  ymin = 10.45) +
  annotation_custom(grid::textGrob("of Genes", gp = gpar(fontsize = 20, fontface="bold")), 
                    xmin = 1.044,  ymin = 10.45) +
  annotation_custom(grid::textGrob("Number", gp = gpar(fontsize = 20, fontface="bold")), 
                    xmin = 1.044,  ymin = 11.1) 
  
  
   # annotate("text", x = 1.026, y = 10.45, size = 7,
          #   label = "Pathway                            Number of Genes", fontface = "bold")
   # annotate("text", x = 1.05, y = 11.5, size = 10,
   #          label = "of Genes")
   
cols <- c("Overlap_percent" = "#ff6666", "Not_Overlap_percent" = "#6666ff", "Not_Observed_percent" = "#818181")

g1 <- ggplot(data = gprofiler_Female_sorted_plot_ready %>% filter(sex=="Female"), aes(x = percent, y = term_name)) +
  geom_bar(stat = "identity",aes(fill=type)) + ggtitle("Female") +
  xlab("Percentage of Total")+
  theme( axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.x= element_text(size=25),
        legend.position = "none",
        plot.margin = unit(c(1,-1,1,0), "mm"),
        axis.text = element_text(size=24),
        plot.title = element_text(size = 29, face = "bold"),
        legend.text=element_text(size=24)) +
  scale_x_reverse(labels = function(x) paste0(x, "%"), limits=c(103,0))+ # Add percent sign
  scale_fill_manual(values = cols , name = "", labels = c("Not Observed","Not Differentially Expressed","Differentially Expressed" ))+
  guides(color = guide_legend(override.aes = list(size = 3)))
  

g2 <- ggplot(data = gprofiler_Female_sorted_plot_ready%>% filter(sex=="Male"), aes(x = percent, y = term_name)) +
  geom_bar(stat = "identity",aes(fill=type)) + ggtitle("Male") +
  xlab("Percentage of Total")+
  theme( axis.title.y = element_blank(),
        axis.text.y = element_blank(), axis.ticks.y = element_blank(), plot.margin = unit(c(1,0,1,-1), "mm"))+
 # scale_fill_discrete(labels = c("Not Observed","Not Differentially Expressed", "Differentially Expressed"))+
  theme(legend.title = element_blank(),
        legend.direction = "vertical",
        axis.text = element_text(size=24),
        legend.text=element_text(size=24),
        axis.title.x= element_text(size=25),
        plot.title = element_text(size = 29, face = "bold"))+
  scale_x_continuous(labels = function(x) paste0(x, "%"), limits=c(0,103))+
  scale_fill_manual(values = cols, drop=FALSE, name = "", labels = c("Not Observed","Not Differentially Expressed","Differentially Expressed" ))

# # Multiple plot function
# #
# # ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# # - cols:   Number of columns in layout
# # - layout: A matrix specifying the layout. If present, 'cols' is ignored.
# #
# # If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# # then plot 1 will go in the upper left, 2 will go in the upper right, and
# # 3 will go all the way across the bottom.
# #http://www.cookbook-r.com/Graphs/Multiple_graphs_on_one_page_(ggplot2)/
# multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
#   library(grid)
# 
#   # Make a list from the ... arguments and plotlist
#   plots <- c(list(...), plotlist)
# 
#   numPlots = length(plots)
# 
#   # If layout is NULL, then use 'cols' to determine layout
#   if (is.null(layout)) {
#     # Make the panel
#     # ncol: Number of columns of plots
#     # nrow: Number of rows needed, calculated from # of cols
#     layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
#                      ncol = cols, nrow = ceiling(numPlots/cols))
#   }
# 
#   if (numPlots==1) {
#     print(plots[[1]])
# 
#   } else {
#     # Set up the page
#     grid.newpage()
#     pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
# 
#     # Make each plot, in the correct location
#     for (i in 1:numPlots) {
#       # Get the i,j matrix positions of the regions that contain this subplot
#       matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
# 
#       print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
#                                       layout.pos.col = matchidx$col))
#     }
#   }
# }

ggarrange(g1,g.mid, g2,  ncol = 3, nrow = 1,common.legend = TRUE,widths=c(1.5,2.5,1.5),
          legend = "bottom")


outfile.999<-"/Users/sunshine2.0/Desktop/WORK/output/pathway_analysis_gProfiler_female_top_overlap_0.01.png"

ggsave(outfile.999, 
       dpi=300, 
       device="png",
       width=23, 
       height=10,
       units="in")

#MALE_____________

# Add M and F datatables by row and filter only for those top pathways in males
gprofiler_Male_sorted_DE_plot <-rbind(gprofiler_Male_result, gprofiler_Female_result_insig) %>%
  filter(term_id %in% gprofiler_Male_result$term_id)
nrow(gprofiler_Male_sorted_DE_plot)

# Add rows for female data not in top male DE overlap pathways
gprofiler_Male_sorted_Female_extra<-gprofiler_Male_sorted_DE_plot %>%
  filter(gprofiler_Male_sorted_DE_plot$term_id %in% gprofiler_Female_result_insig$term_id == FALSE) %>%
  mutate(intersection_size= 0)%>%
  mutate(p_value= NA)%>%
  mutate(sex="Female")

# Combine row wise male top pathways with female extra pathways not shown significantly in DE
gprofiler_Male_sorted_DE_complete_plot<- rbind(gprofiler_Male_sorted_DE_plot, gprofiler_Male_sorted_Female_extra)

# Join DE overlap gprofiler results with all gprofiler results
gprofiler_Male_sorted_plot<-inner_join(gprofiler_Male_sorted_DE_complete_plot, gprofiler_sorted_all_plot, by=c("term_id", "term_name", "sex"))%>%
  mutate("Overlap_percent"= intersection_size/term_size*100) %>%  #add column for intersection_size %
  mutate("Not_Overlap"=intersection_size_all-intersection_size) %>% #add column for not overlap genes in pathway
  mutate("Not_Overlap_percent"=Not_Overlap/term_size*100) %>%
  mutate("Not_Observed"=term_size-Not_Overlap-intersection_size) %>%   #add column for not observed
  mutate("Not_Observed_percent"=Not_Observed/term_size*100) %>%
  mutate("totalP"=Overlap_percent+Not_Overlap_percent+Not_Observed_percent) %>%
  arrange((Overlap_percent))%>%
  arrange(desc(sex)) # arrange in order of top pathways
nrow(gprofiler_Male_sorted_plot)

gprofiler_Male_sorted_plot$term_name <- factor(gprofiler_Male_sorted_plot$term_name, levels=unique(gprofiler_Male_sorted_plot$term_name))

# Arrange percents for plotting
gprofiler_Male_sorted_plot_ready<-gprofiler_Male_sorted_plot %>%
  pivot_longer(cols      = c(Overlap_percent,Not_Overlap_percent,Not_Observed_percent),
               names_to  = "type",
               values_to = "percent")
gprofiler_Male_sorted_plot_ready$type <- factor(gprofiler_Male_sorted_plot_ready$type, levels=c("Not_Observed_percent","Not_Overlap_percent","Overlap_percent"))
gprofiler_Male_sorted_plot$term_name <- factor(gprofiler_Male_sorted_plot$term_name, levels=unique(gprofiler_Male_sorted_plot$term_name))


# Significance X^2 testing with Sidak correction factor of 19
PA_male_pvalue<-data.frame("term_name"= gprofiler_Male_sorted_plot$term_name[1:9],
                           "pvalue"=c(
prop.test(x = gprofiler_Male_sorted_plot %>%
            dplyr::select(c(term_name, sex, Overlap_percent,Not_Overlap_percent,Not_Observed_percent)) %>%
            filter(term_name=="Metabolism of carbohydrates") %>%
            dplyr::select(-sex, -term_name)%>%
            t(),correct=FALSE)[3]$p.value,
prop.test(x = gprofiler_Male_sorted_plot %>%
            dplyr::select(c(term_name, sex, Overlap_percent,Not_Overlap_percent,Not_Observed_percent)) %>%
            filter(term_name=="Glucose metabolism") %>%
            dplyr::select(-sex, -term_name)%>%
            t(),correct=FALSE)[3]$p.value,
prop.test(x = gprofiler_Male_sorted_plot %>%
            dplyr::select(c(term_name, sex, Overlap_percent,Not_Overlap_percent,Not_Observed_percent)) %>%
            filter(term_name=="ROS and RNS production in phagocytes") %>%
            dplyr::select(-sex, -term_name)%>%
            t(),correct=FALSE)[3]$p.value,
prop.test(x = gprofiler_Male_sorted_plot %>%
            dplyr::select(c(term_name, sex, Overlap_percent,Not_Overlap_percent,Not_Observed_percent)) %>%
            filter(term_name=="Keratan sulfate/keratin metabolism") %>%
            dplyr::select(-sex, -term_name)%>%
            t(),correct=FALSE)[3]$p.value,
prop.test(x = gprofiler_Male_sorted_plot %>%
            dplyr::select(c(term_name, sex, Overlap_percent,Not_Overlap_percent,Not_Observed_percent)) %>%
            filter(term_name=="Keratan sulfate biosynthesis") %>%
            dplyr::select(-sex, -term_name)%>%
            t(),correct=FALSE)[3]$p.value,
prop.test(x = gprofiler_Male_sorted_plot %>%
            dplyr::select(c(term_name, sex, Overlap_percent,Not_Overlap_percent,Not_Observed_percent)) %>%
            filter(term_name=="Gluconeogenesis") %>%
            dplyr::select(-sex, -term_name)%>%
            t(),correct=FALSE)[3]$p.value,
prop.test(x = gprofiler_Male_sorted_plot %>%
            dplyr::select(c(term_name, sex, Overlap_percent,Not_Overlap_percent,Not_Observed_percent)) %>%
            filter(term_name=="Integrin cell surface interactions") %>%
            dplyr::select(-sex, -term_name)%>%
            t(),correct=FALSE)[3]$p.value,
prop.test(x = gprofiler_Male_sorted_plot %>%
            dplyr::select(c(term_name, sex, Overlap_percent,Not_Overlap_percent,Not_Observed_percent)) %>%
            filter(term_name=="Triglyceride metabolism") %>%
            dplyr::select(-sex, -term_name)%>%
            t(),correct=FALSE)[3]$p.value,
prop.test(x = gprofiler_Male_sorted_plot %>%
            dplyr::select(c(term_name, sex, Overlap_percent,Not_Overlap_percent,Not_Observed_percent)) %>%
            filter(term_name=="Keratan sulfate degradation") %>%
            dplyr::select(-sex, -term_name)%>%
            t(),correct=FALSE)[3]$p.value))%>%
  #Sidak padj=1-(1-p)^n 
mutate( "Sidak"=1-(1-pvalue)^19)
PA_male_pvalue<-PA_male_pvalue%>%
  mutate("reject"=ifelse(Sidak<0.05, "reject", "fail to reject"))

gprofiler_Male_sorted_plot_ready_pvalue<-PA_male_pvalue %>% full_join(gprofiler_Male_sorted_plot_ready, by= "term_name")%>%
  mutate("term_namePvalue"= ifelse(Sidak<0.999, paste0(term_name, " (", formatC(Sidak, format = "e", digits =2), ")", ifelse(reject=="reject",paste0("*"), paste0("")) ), paste0(term_name)))

#add extra *
gprofiler_Male_sorted_plot_ready_pvalue$term_namePvalue<-
  ifelse(gprofiler_Male_sorted_plot_ready_pvalue$term_namePvalue =="Integrin cell surface interactions (1.95e-02)*", paste0( gprofiler_Male_sorted_plot_ready_pvalue$term_namePvalue,"*"), paste0( gprofiler_Male_sorted_plot_ready_pvalue$term_namePvalue))

# Add '#' for same pathways for both sexes
gprofiler_Male_sorted_plot_ready_pvalue$term_namePvalue<-
  ifelse(gprofiler_Male_sorted_plot_ready_pvalue$term_name %in% same_pathways_sexes$term_name, paste0("# ", gprofiler_Male_sorted_plot_ready_pvalue$term_namePvalue), paste0("", gprofiler_Male_sorted_plot_ready_pvalue$term_namePvalue))
gprofiler_Male_sorted_plot_ready_pvalue$term_namePvalue <- factor(gprofiler_Male_sorted_plot_ready_pvalue$term_namePvalue, levels=unique(gprofiler_Male_sorted_plot_ready_pvalue$term_namePvalue))
gprofiler_Male_sorted_plot_ready_pvalue$type <- factor(gprofiler_Male_sorted_plot_ready_pvalue$type, levels=c("Overlap_percent", "Not_Overlap_percent", "Not_Observed_percent"))
gprofiler_Male_sorted_plot_ready$type <- factor(gprofiler_Male_sorted_plot_ready$type, levels=c("Not_Observed_percent", "Not_Overlap_percent", "Overlap_percent"))

# #overlap f vs M all type
# prop.test(x = as.matrix(gprofiler_Male_sorted_plot %>%
#             dplyr::select(c( term_name, sex, Overlap_percent)) %>%
#             pivot_wider(id_cols=term_name,
#                            names_from=sex,
#                            values_from=Overlap_percent)%>%
#             dplyr::select( Male, Female))
#             ,correct=FALSE)[3]$p.value


# PLOT DE MALE ENRICHED PATHWAYS 
#https://stackoverflow.com/questions/18265941/two-horizontal-bar-charts-with-shared-axis-in-ggplot2-similar-to-population-pyr
library(grid)
g.mid<-ggplot(gprofiler_Male_sorted_plot_ready_pvalue,aes(x=1,y=term_namePvalue))+geom_text(aes(label=term_namePvalue),size=8)+
  geom_segment(aes(x=0.94,xend=0.95,yend=term_namePvalue))+
  geom_segment(aes(x=1.05,xend=1.06,yend=term_namePvalue))+
  geom_text(aes(x=1.045, label=term_size), size=8)+
  ggtitle("")+
  ylab(NULL)+
  xlab(" ")+
  ggtitle(" ")+
  scale_x_continuous(expand=c(0,0),limits=c(0.94,1.06))+
  theme(axis.title=element_blank(),
        panel.grid=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        panel.background=element_blank(),
        axis.text.x=element_text(color=NA),
        axis.ticks.x=element_line(color=NA),
        plot.margin = unit(c(1,-1,1,-1), "mm"),
        axis.text = element_text(size=24),
        axis.title.x= element_text(size=25),
        plot.title = element_text(size = 27, face = "bold"))+
  coord_cartesian(clip = "off") + #keep text from being cut off
  annotation_custom(grid::textGrob("Pathway (P-Value)", gp = gpar(fontsize = 21, fontface="bold")), 
                    xmin = .94,  ymin = 9.45) +
  annotation_custom(grid::textGrob("of Genes", gp = gpar(fontsize = 21, fontface="bold")), 
                    xmin = 1.03,  ymin = 9.45) +
  annotation_custom(grid::textGrob("Number", gp = gpar(fontsize = 21, fontface="bold")), 
                    xmin = 1.03,  ymin = 10.1) 
  # annotate("text", x = 0.9809, y = 8, size = 10,
  #          label = "*",
  #          fontface = "bold")+
  # annotate("text", x = 1.0271, y = 9, size = 10,
  #          label = "*",
  #          fontface = "bold")+
  # annotate("text", x = 1.0281, y = 7, size = 10,
  #          label = "*",
  #          fontface = "bold")
cols <- c("Overlap_percent" = "#ff6666", "Not_Overlap_percent" = "#6666ff", "Not_Observed_percent" = "#818181")

g1 <- ggplot(data = gprofiler_Male_sorted_plot_ready%>% filter(sex=="Female"), aes(x = percent, y = term_name)) +
  geom_bar(stat = "identity",aes(fill=type)) + ggtitle("Female") +xlab("Percentage of Total")+
  theme( axis.title.y = element_blank(),
        axis.text.y = element_blank(), axis.ticks.y = element_blank(), plot.margin = unit(c(1,-1,1,0), "mm"))+
  theme(legend.position = "none",
        axis.text = element_text(size=24),
        axis.title.x= element_text(size=25),
        plot.title = element_text(size = 29, face = "bold"),
        legend.text=element_text(size=24))+
  scale_x_reverse(labels = function(x) paste0(x, "%"), limits=c(103,0))+ # Add percent sign
  scale_fill_manual(values = cols , name = "", labels = c("Not Observed","Not Differentially Expressed","Differentially Expressed" ))+
  guides(color = guide_legend(override.aes = list(size = 3)))


g2 <-   ggplot(data = gprofiler_Male_sorted_plot_ready %>% filter(sex=="Male"), aes(x = percent, y = term_name)) +
  geom_bar(stat = "identity",aes(fill=type)) + ggtitle("Male") +xlab("Percentage of Total")+
  #scale_fill_discrete(labels = c("Differentially Expressed","Not Differentially Expressed","Not Observed                                                                                                                                                                                                                                              " ))+
  theme(
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.x= element_text(size=25),

        legend.title = element_blank(),
        legend.direction = "vertical",
        legend.text=element_text(size=24),
        plot.margin = unit(c(1,0,1,-1), "mm"),
        axis.text = element_text(size=24),
        plot.title = element_text(size = 29, face = "bold")) +
  scale_x_continuous(labels = function(x) paste0(x, "%"), limits=c(0,103))+
  scale_fill_manual(values = cols, drop=FALSE, name = "", labels = c("Not Observed","Not Differentially Expressed","Differentially Expressed" ))


ggarrange(g1,g.mid, g2,  ncol = 3, nrow = 1,common.legend = TRUE,widths=c(1.5,2.5,1.5),
                    legend = "bottom")
outfile.999<-"/Users/sunshine2.0/Desktop/WORK/output/pathway_analysis_gProfiler_male_top_overlap_0.01.png"

ggsave(outfile.999, 
       dpi=300, 
       device="png",
       width=23, 
       height=10,
       units="in")

























#HEIRARICHAL CLUSTERING
# extract 8 sample clusters and 5 gene clusters -----------------------------------------------------------------------

cluster_samples <- cutree(pheatmap0.overlap01$tree_col, k = 8)
cluster_genes <- cutree(pheatmap0.overlap01$tree_row, k = 5)
cluster_samples<-as.data.frame(cluster_samples)
cluster_genes<-as.data.frame(cluster_genes)
head(cluster_samples)
head(cluster_genes)

# save 5 gene clusters for GO of each cluster
cluster_genes_save<-cluster_genes
cluster_genes_save$genes<- rownames(cluster_genes_save) # rownames into new column
cluster_genes_save <-cluster_genes_save %>% arrange(cluster_genes) #order by gene cluster
write.csv(cluster_genes_save, "/Users/sunshine2.0/Desktop/WORK/output/cluster_genes_overlap_0.01.csv")

# Add cluster assignments in log2tpm dataframe
colnames(topDat)==rownames(cluster_samples)
topDat_clusters <- rbind(topDat, "cluster_samples"=(cluster_samples$cluster_samples)) # add row with sample cluster numbers
tail(topDat_clusters)
rownames(topDat_clusters)==rownames(cluster_genes)
topDat_clusters <- cbind(topDat_clusters, "cluster_genes"=c(cluster_genes$cluster_genes, NA)) # add column with sample cluster number
head(topDat_clusters)
topDat_clusters$genes <- rownames(topDat_clusters) # genes ID new column
topDat_clusters <- rbind(topDat_clusters, "samples"=colnames(topDat_clusters)) # samples ID new column
topDat_clusters_arrange<-topDat_clusters %>% arrange(cluster_genes) %>% t() %>% as.data.frame() %>% arrange(V105) %>% t()# arrange by cluster order

# save values for excel to calculate median expression values in log2tpm
write.csv(topDat_clusters_arrange, "/Users/sunshine2.0/Desktop/WORK/output/clusters_genes_samples_overlap_0.01.csv")

#----------------------------------------------------------------------------------------------------------------------

# Get list of DE genes less than adjusted p-value cut off of 0.05 and 1 < LFC for heatmap clustering (upregulated)

length(unique(c(DEgenelist1.up05, DEgenelist2.up05, DEgenelist3.up05, DEgenelist4.up05, DEgenelist5.up05, DEgenelist6.up05, DEgenelist7.up05, DEgenelist8.up05, DEgenelist9.up05)))
DElist0.up05<-unique(c(DEgenelist1.up05, DEgenelist2.up05, DEgenelist3.up05, DEgenelist4.up05, DEgenelist5.up05, DEgenelist6.up05, DEgenelist7.up05, DEgenelist8.up05, DEgenelist9.up05))

# Get list of DE genes less than adjusted p-value cut off of 0.05 and LFC < -1 for heatmap clustering (downregulated)

length(unique(c(DEgenelist1.down05, DEgenelist2.down05, DEgenelist3.down05, DEgenelist4.down05, DEgenelist5.down05, DEgenelist6.down05, DEgenelist7.down05, DEgenelist8.down05, DEgenelist9.down05)))
DElist0.down05<-unique(c(DEgenelist1.down05, DEgenelist2.down05, DEgenelist3.down05, DEgenelist4.down05, DEgenelist5.down05, DEgenelist6.down05, DEgenelist7.down05, DEgenelist8.down05, DEgenelist9.down05))

# Get list of overlap of upregulated and downregulated genes 0.05

length(DElist0.up05) +length(DElist0.down05)-length(DElist0.05)
up0.05<-as.data.frame(DElist0.up05, stringsAsFactors=FALSE) %>%
  rename( DElist0.up05= "x")
down0.05<-as.data.frame("x"=DElist0.down05, stringsAsFactors=FALSE) %>%
  rename( DElist0.down05= "x")
overlap0.05<-up0.05 %>%
  inner_join(down0.05, by= "x")
# Scaled data of topGenes
topDat <- log2tpm[rownames(log2tpm)%in% overlap0.05$x, ]
nrow(topDat )

# Heatmap overlap of upregulated and downregulated genes padj < 0.05

pheatmap0.overlap05<-pheatmap(topDat, cluster_rows = TRUE, scale = "none", clustering_method = "average",
                              clustering_distance_cols = "euclidean", clustering_distance_rows = "euclidean",
                              annotation = samples[, c("month", "gender","grp")], show_rownames = FALSE)

save_pheatmap_pdf(pheatmap0.overlap05, "/Users/sunshine2.0/Desktop/WORK/output/heatmap0.05overlap.pdf")



# Get list of 5 most up-regulated and down-regulated DE genes

DEsubsetAll<-rbind(DEsubset1,DEsubset2,DEsubset3,DEsubset4,DEsubset5,DEsubset6,DEsubset7,DEsubset8,DEsubset9)
DEsubsetAll<-DEsubsetAll %>%
  mutate("absLFC"=abs(log2FoldChange)) %>%
  arrange(desc(log2FoldChange))
top10DElist<-unique(DEsubsetAll$gene)
top10DElist<-top10DElist[c(1:5,(length(top10DElist)-4):length(top10DElist))]

# Scaled data of topGenes
topDat <- log2tpm[rownames(log2tpm)%in% top10DElist, ]
nrow(topDat )

# Heatmap 5 most up-regulated and down-regulated padj < 0.05

pheatmap0.top01<-pheatmap(topDat, cluster_rows = TRUE, scale = "none", clustering_method = "average",
                           clustering_distance_cols = "euclidean", clustering_distance_rows = "euclidean",
                           annotation = samples[, c("month", "gender","grp")], show_rownames = FALSE)

save_pheatmap_pdf(pheatmap0.top01, "/Users/sunshine2.0/Desktop/WORK/output/heatmap0.05top10.pdf")

# Get 10 most DE gene names by 5 highest and 5 lowest log2FC with adjusted p-value < 0.05 (|highest log2FC| are all downregulated genes) and 1 < LFC < -1

# Subset log2tpm to those top genes
top10DE <- subset(log2tpm, rownames(log2tpm) %in% top10DElist)
top10DE <-log2tpm[rownames(log2tpm)%in% top10DElist, ]
top10DE<-transform(top10DE, "regulated"=case_when(rownames(top10DE) %in% top10DElist[1:5] ~ paste("up"),
                                       rownames(top10DE) %in% top10DElist[(length(top10DElist)-5):length(top10DElist)] ~ paste("down")))
top10DE<-mutate(top10DE, "gene" = rownames(top10DE))

# Edit meta data month to full name
samples<-transform(samples,month= case_when(month=="A" ~ "April",
                 month=="J" ~ "June",
                 month=="N" ~ "Nov"))

# Prep dataframe for plotting
samples<-transform(samples, "samples"=rownames(samples))
top10DEcomb<-ColumnComb(DataFrame=top10DE,MetaData=samples,
                       ColNameCol=6,NewColName="Log2.TPM")

# Plot top 10 differentially expressed
class(top10DEcomb$Emu)
top10DEcomb$month<-factor(top10DEcomb$month,
                          levels=unique(top10DEcomb$month))
top10DEcomb<-top10DEcomb%>%
  arrange(regulated)
top10DEcomb$regulated<-factor(top10DEcomb$regulated,
                          levels=unique(top10DEcomb$regulated))
top10DEcomb$gene<-factor(top10DEcomb$gene,
                              levels=unique(top10DEcomb$gene))

  ggplot(top10DEcomb, aes(month,Log2.TPM)) +
    geom_point(aes(group=Emu, color=Emu)) +
    geom_line(aes(group=Emu, color=Emu)) +
    facet_grid(~gene) +
    theme_bw() +
    ggtitle("Top 5 Most Downregulated and Top 5 Most Upregulated Genes")+
    xlab('Month') +
    scale_color_manual(values=c("#46A50C", "#009345", "#009391", "#003193","#920CA5", "#D0159D","#D01534", "#A5620C"))

  outfiletop10<-"/Users/sunshine2.0/Desktop/WORK/output/top5updown_regulated.png"

  ggsave(outfiletop10,
         dpi=300,
         device="png",
         width=45,
         height=20,
         units="cm")

