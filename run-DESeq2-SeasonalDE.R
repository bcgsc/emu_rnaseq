library(tximport)
library(tidyverse)
library(rhdf5)
library(DESeq2)

## TO RUN DIFFERENTIAL EXPRESSION ANALYSIS (SEASONAL COMPARISONS)

# Read in metadata file
# samples <- read_csv("/path/to/metadata.csv")
samples <- read.table("/Users/sunshine2.0/Desktop/WORK/data/emu.SampleInfo", head = FALSE, sep = "", col.names = c("sample","monthEmu","gender"))
samples <- transform(samples, month = substr(monthEmu, 1, 1), Emu = substr(monthEmu, 2, 2))

# Read in the files where your abundance is (quants from Salmon)
# files should be a vector containing all file paths
files <- scan("/Users/sunshine2.0/Desktop/WORK/data/isoform_abundance_files.txt", character(), quote="")

# Name the paths in the vector after which treatment they correspond to
# You might have to change this line depending on your column names in metadata
names(files) <- samples$sample

# Read in a csv that maps isoforms to their gene
tx2gene <- read.table("/Users/sunshine2.0/Desktop/WORK/data/gene_map.tsv", sep="\t")

# Import the abundance from Salmon
salmon_quant <- tximport(files,type="salmon",tx2gene=tx2gene) #, ignoreTxVersion = TRUE, ignoreAfterBar = TRUE)

# Change rownames to correspond with the column names
rownames(samples) <- colnames(salmon_quant$counts)

# Use DESeq2
ddsAll <- DESeqDataSetFromTximport(salmon_quant, samples, ~month)

# Adjusted p-value threshold default, can be modified (default alpha=0.1)
padj_cutoff<-0.05
padj_cutoff1<-0.01

# Log fold change threshold default, can be modified
lfc <-0
lfc_cutoff <-1

## PLOT VOLCANO PLOTS For each condition vs control
# JUNE VS APRIL

# Subset the relevant columns from the full dataset and choose the control condition so log2 fold changes are calculated as treatment over control
dds <- ddsAll[ , (ddsAll$month == "A")|(ddsAll$month == "J") ]
dds$month <- droplevels( dds$month )
dds$month <- relevel(dds$month, ref = "A")
# check if subsetted correctly
#as.data.frame( colData(dds) )

dds <- DESeq(dds)
res <- results(dds, alpha=padj_cutoff)

# Number of adjusted p-values less than padj_cutoff
sum(res$padj < padj_cutoff, na.rm=TRUE)
# summary (default lfcThreshold = 0)
summary(res)
sum(res$baseMean ==0)

# Convert to dataframe
df <- as.data.frame(res)

# Make summary table similar to summary(res)
total<- nrow(df)-sum(df$baseMean ==0, na.rm=TRUE)
lowCountsSubset<-df[ ((!is.na(df$pvalue))&(is.na(df$padj))),]
dfsummary1 <- data.frame("LFCup" = c(round(sum(df$log2FoldChange>0& df$padj<padj_cutoff, na.rm=TRUE), digits=0),round(sum(df$log2FoldChange>0& df$padj<padj_cutoff, na.rm=TRUE)/total*100, digits = 1)), 
                      "LFCdown" = c(round(sum(df$log2FoldChange<0& df$padj<padj_cutoff, na.rm=TRUE), digits = 0),round(sum(df$log2FoldChange<0& df$padj<padj_cutoff, na.rm=TRUE)/total*100, digits = 1)), 
                      "outliers"= c(round(sum(!is.na(df$log2FoldChange)&is.na(df$pvalue)), digits = 0),round(sum(!is.na(df$log2FoldChange)&is.na(df$pvalue))/round(sum(is.na(df$padj)), digits = 0)*100, digits = 1)),
                      "lowCounts"= c(round(sum(!is.na(df$pvalue)&is.na(df$padj)), digits = 0),round(sum(!is.na(df$pvalue)&is.na(df$padj))/round(sum(is.na(df$padj)), digits = 0)*100, digits = 1)),
                      "baseMean"= c(round(sum(df$baseMean ==0, na.rm=TRUE), digits = 0),round(sum(df$baseMean ==0, na.rm=TRUE)/round(sum(is.na(df$padj)), digits = 0)*100, digits = 1)),
                      "total"= c(round(total, digits = 0), round(sum((sum(df$log2FoldChange>0& df$padj<padj_cutoff, na.rm=TRUE)/total*100),(sum(df$log2FoldChange<0& df$padj<padj_cutoff, na.rm=TRUE)/total*100),(sum(!is.na(df$log2FoldChange)&is.na(df$pvalue))/round(sum(is.na(df$padj)), digits = 0)*100),(sum(!is.na(df$pvalue)&is.na(df$padj))/round(sum(is.na(df$padj)), digits = 0)*100),(sum(df$baseMean ==0, na.rm=TRUE)/round(sum(is.na(df$padj)), digits = 0)*100), digits = 1))),
                      "treatment" = "MF_JvsA",
                      "maxlowCounts"=c(round(max(lowCountsSubset$baseMean), digits=0),NA))
                    
# Make a data table of df statistics
dfstats1<-data.frame("#NA"=c(round(sum(is.na(df$padj)), digits = 0),NA),
                     "#sig"=c(sum(df$padj < padj_cutoff, na.rm=TRUE),round(sum(df$padj < padj_cutoff, na.rm=TRUE)/(nrow(df)-sum(is.na(df$padj)))*100, digits = 1)),
                     "#nonsig"= c(sum(df$padj > padj_cutoff | df$padj == padj_cutoff, na.rm=TRUE),round(sum(df$padj > padj_cutoff | df$padj == padj_cutoff, na.rm=TRUE)/(nrow(df)-sum(is.na(df$padj)))*100, digits = 1)),
                     "total"= c(nrow(df), NA),
                     "treatment" = "MF_JvsA")

# Add -log(pvalue) column, and padj category 
df_log <- mutate(df, logp=-log10(pvalue), padjcat=padj,
                 lfccat=log2FoldChange)

df_log <- transform(df_log, padjcat=case_when(padj < padj_cutoff & padj >= padj_cutoff1 ~ paste(padj_cutoff1,"< P-Value <", padj_cutoff),
                                              padj >= padj_cutoff ~ paste("P-Value >=", padj_cutoff),
                                              padj < padj_cutoff1 ~ paste("P-Value <", padj_cutoff1)))

df_log<- transform(df_log, lfccat= case_when(lfccat > lfc_cutoff ~ paste("LFC > +", lfc_cutoff), 
                                             lfccat < -lfc_cutoff ~ paste("LFC < -", lfc_cutoff),
                                             lfccat <= lfc_cutoff & lfccat >= (-lfc_cutoff) ~ paste("-",lfc_cutoff," =< LFC =< +", lfc_cutoff)))

df_log <- mutate(df_log, cat=paste(padjcat,";",lfccat))
df_log$cat[is.na(df_log$padjcat)] <- NA
df_log$cat[df_log$lfccat == paste("-",lfc_cutoff," =< LFC =< +", lfc_cutoff)] <- NA
df_log$cat[df_log$padjcat == paste("P-Value >=", padj_cutoff)] <- NA

upcat<-df_log %>%
  filter(cat == paste(padj_cutoff1,"< P-Value <", padj_cutoff,";","LFC > +", lfc_cutoff) ) %>%
  nrow()

upperupcat<-df_log %>%
  filter(cat == paste("P-Value <", padj_cutoff1,";","LFC > +", lfc_cutoff) ) %>%
  nrow()

downcat<- df_log %>%
  filter(cat == paste(padj_cutoff1,"< P-Value <", padj_cutoff,";","LFC < -", lfc_cutoff) ) %>%
  nrow()

upperdowncat<- df_log %>%
  filter(cat == paste("P-Value <", padj_cutoff1,";","LFC < -", lfc_cutoff) ) %>%
  nrow()

# Subset significant genes into list of gene names only
df_log <- mutate(df_log, 'gene'=rownames(df))
DEsubset1 <- df_log[  !is.na(df_log$cat) , ]
DEgenelist1 <- DEsubset1$gene #all significant genes meeting P-Value cutoff <0.05 

DEsubset.01 <- df_log %>%
  filter(cat== paste("P-Value <", padj_cutoff1,";","LFC < -", lfc_cutoff) | cat== paste("P-Value <", padj_cutoff1,";","LFC > +", lfc_cutoff) )
DEgenelist1.01 <- DEsubset.01$gene #all significant genes meeting padj cutoff <0.01

DEsubset.up01 <- df_log %>%
  filter(cat== paste("P-Value <", padj_cutoff1,";","LFC > +", lfc_cutoff) )
DEgenelist1.up01 <- DEsubset.up01$gene
length(DEgenelist1.up01) #all up-regulated genes meeting padj cutoff <0.01

DEsubset.down05 <- df_log %>%
  filter(cat== paste("P-Value <", padj_cutoff,";","LFC < -", lfc_cutoff) )
DEgenelist1.down01 <- DEsubset.down01$gene
length(DEgenelist1.down01) #all down-regulated genes meeting padj cutoff <0.01

DEsubset.up05 <-  df_log[  !is.na(df_log$cat) , ]
DEsubset.up05 <- DEsubset.up05 %>%
  filter(lfccat== paste("LFC > +", lfc_cutoff) )
DEgenelist1.up05 <- DEsubset.up05$gene
length(DEgenelist1.up05) #all up-regulated genes meeting padj cutoff <0.05

DEsubset.down05 <-  df_log[  !is.na(df_log$cat) , ]
DEsubset.down05 <- DEsubset.down05 %>%
  filter(lfccat== paste("LFC < -", lfc_cutoff) )
DEgenelist1.down05 <- DEsubset.down05$gene
length(DEgenelist1.down05) #all down-regulated genes meeting padj cutoff <0.05

# Volacno plot font sizes
axisfont <- 14
labelfont <- 5
allfont <- 18

# Reassign rownames
rownames(df_log) <- rownames(df)

# Convert padjcat to factor
df_log$padjcat <- as.factor(df_log$padjcat)

# Replace NA with Not Significant or Substantial
df_log<-df_log %>%
  mutate(cat = if_else(is.na(cat), "Not Significant or Substantial", cat))
df_log$cat<-factor(df_log$cat,
                   levels=c("P-Value < 0.01 ; LFC < - 1" ,"P-Value < 0.01 ; LFC > + 1", "0.01 < P-Value < 0.05 ; LFC < - 1" , "0.01 < P-Value < 0.05 ; LFC > + 1", "Not Significant or Substantial"))
cols <- c("0.01 < P-Value < 0.05 ; LFC < - 1" = "#00ace6", "0.01 < P-Value < 0.05 ; LFC > + 1" = "#ff66cc", "P-Value < 0.01 ; LFC < - 1" = "#0000ff", "P-Value < 0.01 ; LFC > + 1" = "#ff0000", "Not Significant or Substantial"="#818181")

# Plot discrete volcano plot
p1<-ggplot(df_log, aes(x=log2FoldChange, y=logp, color=cat)) + 
  geom_point(na.rm = TRUE,alpha=0.5) +
  geom_vline(xintercept=c(-1,1),linetype="dotted") + 
  geom_hline(yintercept=c(-log10(padj_cutoff/10),-log10(padj_cutoff1/36)) ,linetype="dotted") + 
  theme(text = element_text(size=allfont),
        axis.text = element_text(size=axisfont)) + 
  ggtitle("June vs. April") +
  ylab(expression(paste("\u2013", log[10], "(P-Value)"))) + 
  theme_bw()+
  xlab(expression(paste(log[2],"(FoldChange)"))) +
  labs(color='Genes' ) +
  xlim(-8,8)+
  scale_colour_manual(values = cols, drop=FALSE, name = "Genes", labels = c("P-Value < 0.01 ; LFC < - 1", "P-Value < 0.01 ; LFC > + 1", "0.01 \u2264 P-Value < 0.05 ; LFC < - 1", "0.01 \u2264 P-Value < 0.05 ; LFC > + 1",  "Not Significant or Substantial"))+
  guides(color = guide_legend(override.aes = list(size = 3)))

ylim<-ggplot_build(p1)$layout$panel_scales_y[[1]]$range$range[2]

p1<-p1  + geom_text(x=6.25, y=ylim-2.45*(ylim/8.25), label= upcat, color='#ff66cc') 
p1<-p1 + geom_text(x=-6.25, y=ylim-2.45*(ylim/8.25), label= downcat, color='#00ace6') 
p1<-p1 + geom_text(x=6.25, y=ylim-.7*(ylim/8.25), label= upperupcat, color='#ff0000') 
p1<-p1 + geom_text(x=-6.25, y=ylim-.7*(ylim/8.25), label= upperdowncat, color='#0000ff') 

outfile1<-"/Users/sunshine2.0/Desktop/WORK/output/volcanoJvsA_"
outfile2<-".png"
outfile<-paste(outfile1,padj_cutoff,outfile2, sep="")

ggsave(outfile, 
       dpi=300, 
       device="png",
       width=23, 
       height=14,
       units="cm")

# Plot continuous volcano plot
# ggplot(df_log, aes(x=log2FoldChange, y=logp, color=padj)) +
#   # geom_point(na.rm = TRUE, alpha=0.3) +
#   geom_vline(xintercept=lfc, linetype="dotted") +
#   theme(text = element_text(size=allfont),
#         axis.text = element_text(size=axisfont)) +
#   ggtitle("June vs. April") +
#   ylab("\u2013log10(P-Value)") +
#   xlab("log2(FoldChange)")  +
#   labs(color='Adj P-Value' )
# 
# outfile1 <- "/home/krwright/Downloads/volcano_continuousJvsA.png"
# 
# ggsave(outfile1,
#        dpi=300,
#        device="png",
#        width=23,
#        height=14,
#        units="cm")

# NOVEMBER VS JUNE 

# Subset the relevant columns from the full dataset and choose the control condition so log2 fold changes are calculated as treatment over control
dds <- ddsAll[ , (ddsAll$month == "N")|(ddsAll$month == "J") ]
dds$month <- droplevels( dds$month )
dds$month <- relevel(dds$month, ref = "J")
# check if subsetted correctly
#as.data.frame( colData(dds) )

dds <- DESeq(dds)
res <- results(dds, alpha=padj_cutoff)

# Number of adjusted p-values less than 0.1
sum(res$padj < padj_cutoff, na.rm=TRUE)
# summary (default lfcThreshold = 0)
summary(res)
sum(res$baseMean ==0)

# Convert to dataframe
df <- as.data.frame(res)

# Make summary table similar to summary(res)
total<- nrow(df)-sum(df$baseMean ==0, na.rm=TRUE)
lowCountsSubset<-df[ ((!is.na(df$pvalue))&(is.na(df$padj))),]
dfsummary2 <-  data.frame("LFCup" = c(round(sum(df$log2FoldChange>0& df$padj<padj_cutoff, na.rm=TRUE), digits=0),round(sum(df$log2FoldChange>0& df$padj<padj_cutoff, na.rm=TRUE)/total*100, digits = 1)), 
                          "LFCdown" = c(round(sum(df$log2FoldChange<0& df$padj<padj_cutoff, na.rm=TRUE), digits = 0),round(sum(df$log2FoldChange<0& df$padj<padj_cutoff, na.rm=TRUE)/total*100, digits = 1)), 
                          "outliers"= c(round(sum(!is.na(df$log2FoldChange)&is.na(df$pvalue)), digits = 0),round(sum(!is.na(df$log2FoldChange)&is.na(df$pvalue))/round(sum(is.na(df$padj)), digits = 0)*100, digits = 1)),
                          "lowCounts"= c(round(sum(!is.na(df$pvalue)&is.na(df$padj)), digits = 0),round(sum(!is.na(df$pvalue)&is.na(df$padj))/round(sum(is.na(df$padj)), digits = 0)*100, digits = 1)),
                          "baseMean"= c(round(sum(df$baseMean ==0, na.rm=TRUE), digits = 0),round(sum(df$baseMean ==0, na.rm=TRUE)/round(sum(is.na(df$padj)), digits = 0)*100, digits = 1)),
                          "total"= c(round(total, digits = 0), round(sum((sum(df$log2FoldChange>0& df$padj<padj_cutoff, na.rm=TRUE)/total*100),(sum(df$log2FoldChange<0& df$padj<padj_cutoff, na.rm=TRUE)/total*100),(sum(!is.na(df$log2FoldChange)&is.na(df$pvalue))/round(sum(is.na(df$padj)), digits = 0)*100),(sum(!is.na(df$pvalue)&is.na(df$padj))/round(sum(is.na(df$padj)), digits = 0)*100),(sum(df$baseMean ==0, na.rm=TRUE)/round(sum(is.na(df$padj)), digits = 0)*100), digits = 1))),
                          "treatment" = "MF_NvsJ",
                          "maxlowCounts"=c(round(max(lowCountsSubset$baseMean), digits=0),NA))
   
# Make a data table of df statistics
dfstats2<-data.frame("#NA"=c(round(sum(is.na(df$padj)), digits = 0),NA),
                     "#sig"=c(sum(df$padj < padj_cutoff, na.rm=TRUE),round(sum(df$padj < padj_cutoff, na.rm=TRUE)/(nrow(df)-sum(is.na(df$padj)))*100, digits = 1)),
                     "#nonsig"= c(sum(df$padj > padj_cutoff | df$padj == padj_cutoff, na.rm=TRUE),round(sum(df$padj > padj_cutoff | df$padj == padj_cutoff, na.rm=TRUE)/(nrow(df)-sum(is.na(df$padj)))*100, digits = 1)),
                     "total"= c(nrow(df), NA),
                     "treatment" = "MF_NvsJ")

# Add -log(pvalue) column, and padj category 
df_log <- mutate(df, logp=-log10(pvalue), padjcat=padj,
                 lfccat=log2FoldChange)

df_log <- transform(df_log, padjcat=case_when(padj < padj_cutoff & padj >= padj_cutoff1 ~ paste(padj_cutoff1,"< padj <", padj_cutoff),
                                              padj >= padj_cutoff ~ paste("padj >=", padj_cutoff),
                                              padj < padj_cutoff1 ~ paste("padj <", padj_cutoff1)))

df_log<- transform(df_log, lfccat= case_when(lfccat > lfc_cutoff ~ paste("lfc > +", lfc_cutoff), 
                                             lfccat < -lfc_cutoff ~ paste("lfc < -", lfc_cutoff),
                                             lfccat <= lfc_cutoff & lfccat >= (-lfc_cutoff) ~ paste("-",lfc_cutoff," =< lfc =< +", lfc_cutoff)))

df_log <- mutate(df_log, cat=paste(padjcat,";",lfccat))
df_log$cat[is.na(df_log$padjcat)] <- NA
df_log$cat[df_log$lfccat == paste("-",lfc_cutoff," =< lfc =< +", lfc_cutoff)] <- NA
df_log$cat[df_log$padjcat == paste("padj >=", padj_cutoff)] <- NA

upcat<-df_log %>%
  filter(cat == paste(padj_cutoff1,"< padj <", padj_cutoff,";","lfc > +", lfc_cutoff) ) %>%
  nrow()

upperupcat<-df_log %>%
  filter(cat == paste("padj <", padj_cutoff1,";","lfc > +", lfc_cutoff) ) %>%
  nrow()

downcat<- df_log %>%
  filter(cat == paste(padj_cutoff1,"< padj <", padj_cutoff,";","lfc < -", lfc_cutoff) ) %>%
  nrow()

upperdowncat<- df_log %>%
  filter(cat == paste("padj <", padj_cutoff1,";","lfc < -", lfc_cutoff) ) %>%
  nrow()

# Subset significant genes into list of gene names only
df_log <- mutate(df_log, 'gene'=rownames(df))
DEsubset2 <- df_log[  !is.na(df_log$cat) , ]
DEgenelist2 <- DEsubset2$gene #all significant genes meeting padj cutoff <0.05 

DEsubset.01 <- df_log %>%
  filter(cat== paste("padj <", padj_cutoff1,";","lfc < -", lfc_cutoff) | cat== paste("padj <", padj_cutoff1,";","lfc > +", lfc_cutoff) )
DEgenelist2.01 <- DEsubset.01$gene #all significant genes meeting padj cutoff <0.01

DEsubset.up01 <- df_log %>%
  filter(cat== paste("padj <", padj_cutoff1,";","lfc > +", lfc_cutoff) )
DEgenelist2.up01 <- DEsubset.up01$gene #all up-regulated genes meeting padj cutoff <0.01

DEsubset.down01 <- df_log %>%
  filter(cat== paste("padj <", padj_cutoff1,";","lfc < -", lfc_cutoff) )
DEgenelist2.down01 <- DEsubset.down01$gene #all down-regulated genes meeting padj cutoff <0.01

DEsubset.up05 <-  df_log[  !is.na(df_log$cat) , ]
DEsubset.up05 <- DEsubset.up05 %>%
  filter(lfccat== paste("lfc > +", lfc_cutoff) )
DEgenelist2.up05 <- DEsubset.up05$gene #all up-regulated genes meeting padj cutoff <0.05

DEsubset.down05 <-  df_log[  !is.na(df_log$cat) , ]
DEsubset.down05 <- DEsubset.down05 %>%
  filter(lfccat== paste("lfc < -", lfc_cutoff) )
DEgenelist2.down05 <- DEsubset.down05$gene #all down-regulated genes meeting padj cutoff <0.05

# Volcano plot font sizes
axisfont <- 14
labelfont <- 5
allfont <- 18

# Reassign rownames
rownames(df_log) <- rownames(df)

# Convert padjcat to factor
df_log$padjcat <- as.factor(df_log$padjcat)

# Replace NA with Not Significant or Substantial
df_log<-df_log %>%
  mutate(cat = if_else(is.na(cat), "Not Significant or Substantial", cat))
df_log$cat<-factor(df_log$cat,
                   levels=c("padj < 0.01 ; lfc < - 1" ,"padj < 0.01 ; lfc > + 1", "0.01 < padj < 0.05 ; lfc < - 1" , "0.01 < padj < 0.05 ; lfc > + 1", "Not Significant or Substantial"))
cols <- c("0.01 < padj < 0.05 ; lfc < - 1" = "#00ace6", "0.01 < padj < 0.05 ; lfc > + 1" = "#ff66cc", "padj < 0.01 ; lfc < - 1" = "#0000ff", "padj < 0.01 ; lfc > + 1" = "#ff0000", "Not Significant or Substantial"="#818181")

# Plot discrete volcano plot
p2<-ggplot(df_log, aes(x=log2FoldChange, y=logp, color=cat)) + 
  geom_point(na.rm = TRUE,alpha=0.5) +
  geom_vline(xintercept=c(-1,1),linetype="dotted") + 
  geom_hline(yintercept=c(2.25,3.34) ,linetype="dotted") + 
  theme(text = element_text(size=allfont),
        axis.text = element_text(size=axisfont)) + 
  ggtitle("Nov vs. June") +
  ylab(expression(paste("\u2013", log[10], "(P-Value)"))) + 
  theme_bw()+
  xlab(expression(paste(log[2],"(FoldChange)"))) +
  labs(color='Genes' ) +
  xlim(-8,8)+
  scale_colour_manual(values = cols, drop=FALSE, name = "Genes", labels = c("P-Value < 0.01 ; LFC < - 1", "P-Value < 0.01 ; LFC > + 1", "0.01 \u2264 P-Value < 0.05 ; LFC < - 1", "0.01 \u2264 P-Value < 0.05 ; LFC > + 1",  "Not Significant or Substantial"))+
  guides(color = guide_legend(override.aes = list(size = 3)))

ylim<-ggplot_build(p2)$layout$panel_scales_y[[1]]$range$range[2]

p2<-p2  + geom_text(x=6.25, y=ylim-2.45*(ylim/8.25), label= upcat, color='#ff66cc') 
p2<-p2 + geom_text(x=-6.25, y=ylim-2.45*(ylim/8.25), label= downcat, color='#00ace6') 
p2<-p2 + geom_text(x=6.25, y=ylim-.7*(ylim/8.25), label= upperupcat, color='#ff0000') 
p2<-p2 + geom_text(x=-6.25, y=ylim-.7*(ylim/8.25), label= upperdowncat, color='#0000ff') 


outfile1<-"/Users/sunshine2.0/Desktop/WORK/output/volcanoNvsJ_"
outfile2<-".png"
outfile<-paste(outfile1,padj_cutoff,outfile2, sep="")

ggsave(outfile, 
       dpi=300, 
       device="png",
       width=23, 
       height=14,
       units="cm")
 
# # Plot continuous volcano plot
# ggplot(df_log, aes(x=log2FoldChange, y=logp, color=padj)) + 
#   geom_point(na.rm = TRUE, alpha=0.3) + 
#   geom_vline(xintercept=lfc, linetype="dotted") + 
#   theme(text = element_text(size=allfont),
#         axis.text = element_text(size=axisfont)) + 
#   ggtitle("November vs. June") +
#   ylab("\u2013log10(P-Value)") + 
#   xlab("log2(FoldChange)")  +
#   labs(color='Adj P-Value' ) +
#   xlim(-8,8)
# 
# outfile1 <- "/Users/sunshine2.0/Desktop/WORK/output/volcano_continuousNvsJ.png"  
# 
# ggsave(outfile1, 
#        dpi=300, 
#        device="png", 
#        width=23, 
#        height=14, 
#        units="cm")

# APRIL VS NOVEMBER

# Subset the relevant columns from the full dataset and choose the control condition so log2 fold changes are calculated as treatment over control
dds <- ddsAll[ , (ddsAll$month == "N")|(ddsAll$month == "A") ]
dds$month <- droplevels( dds$month )
dds$month <- relevel(dds$month, ref = "N")
# check if subsetted correctly
#as.data.frame( colData(dds) )

dds <- DESeq(dds)
res <- results(dds, alpha=padj_cutoff)

# Number of adjusted p-values less than 0.1
sum(res$padj < padj_cutoff, na.rm=TRUE)
# summary (default lfcThreshold = 0)
summary(res)
sum(res$baseMean ==0)

# Convert to dataframe
df <- as.data.frame(res)

# Make summary table similar to summary(res)
total<- nrow(df)-sum(df$baseMean ==0, na.rm=TRUE)
lowCountsSubset<-df[ ((!is.na(df$pvalue))&(is.na(df$padj))),]
dfsummary3 <-  data.frame("LFCup" = c(round(sum(df$log2FoldChange>0& df$padj<padj_cutoff, na.rm=TRUE), digits=0),round(sum(df$log2FoldChange>0& df$padj<padj_cutoff, na.rm=TRUE)/total*100, digits = 1)), 
                          "LFCdown" = c(round(sum(df$log2FoldChange<0& df$padj<padj_cutoff, na.rm=TRUE), digits = 0),round(sum(df$log2FoldChange<0& df$padj<padj_cutoff, na.rm=TRUE)/total*100, digits = 1)), 
                          "outliers"= c(round(sum(!is.na(df$log2FoldChange)&is.na(df$pvalue)), digits = 0),round(sum(!is.na(df$log2FoldChange)&is.na(df$pvalue))/round(sum(is.na(df$padj)), digits = 0)*100, digits = 1)),
                          "lowCounts"= c(round(sum(!is.na(df$pvalue)&is.na(df$padj)), digits = 0),round(sum(!is.na(df$pvalue)&is.na(df$padj))/round(sum(is.na(df$padj)), digits = 0)*100, digits = 1)),
                          "baseMean"= c(round(sum(df$baseMean ==0, na.rm=TRUE), digits = 0),round(sum(df$baseMean ==0, na.rm=TRUE)/round(sum(is.na(df$padj)), digits = 0)*100, digits = 1)),
                          "total"= c(round(total, digits = 0), round(sum((sum(df$log2FoldChange>0& df$padj<padj_cutoff, na.rm=TRUE)/total*100),(sum(df$log2FoldChange<0& df$padj<padj_cutoff, na.rm=TRUE)/total*100),(sum(!is.na(df$log2FoldChange)&is.na(df$pvalue))/round(sum(is.na(df$padj)), digits = 0)*100),(sum(!is.na(df$pvalue)&is.na(df$padj))/round(sum(is.na(df$padj)), digits = 0)*100),(sum(df$baseMean ==0, na.rm=TRUE)/round(sum(is.na(df$padj)), digits = 0)*100), digits = 1))),
                          "treatment" = "MF_NvsA",
                          "maxlowCounts"=c(round(max(lowCountsSubset$baseMean), digits=0),NA))

# Make a data table of df statistics
dfstats3<-data.frame("#NA"=c(round(sum(is.na(df$padj)), digits = 0),NA),
                     "#sig"=c(sum(df$padj < padj_cutoff, na.rm=TRUE),round(sum(df$padj < padj_cutoff, na.rm=TRUE)/(nrow(df)-sum(is.na(df$padj)))*100, digits = 1)),
                     "#nonsig"= c(sum(df$padj > padj_cutoff | df$padj == padj_cutoff, na.rm=TRUE),round(sum(df$padj > padj_cutoff | df$padj == padj_cutoff, na.rm=TRUE)/(nrow(df)-sum(is.na(df$padj)))*100, digits = 1)),
                     "total"= c(nrow(df), NA),
                     "treatment" = "MF_NvsA")

# Add -log(pvalue) column, and padj category 
df_log <- mutate(df, logp=-log10(pvalue), padjcat=padj,
                 lfccat=log2FoldChange)

df_log <- transform(df_log, padjcat=case_when(padj < padj_cutoff & padj >= padj_cutoff1 ~ paste(padj_cutoff1,"< padj <", padj_cutoff),
                                              padj >= padj_cutoff ~ paste("padj >=", padj_cutoff),
                                              padj < padj_cutoff1 ~ paste("padj <", padj_cutoff1)))

df_log<- transform(df_log, lfccat= case_when(lfccat > lfc_cutoff ~ paste("lfc > +", lfc_cutoff), 
                                             lfccat < -lfc_cutoff ~ paste("lfc < -", lfc_cutoff),
                                             lfccat <= lfc_cutoff & lfccat >= (-lfc_cutoff) ~ paste("-",lfc_cutoff," =< lfc =< +", lfc_cutoff)))

df_log <- mutate(df_log, cat=paste(padjcat,";",lfccat))
df_log$cat[is.na(df_log$padjcat)] <- NA
df_log$cat[df_log$lfccat == paste("-",lfc_cutoff," =< lfc =< +", lfc_cutoff)] <- NA
df_log$cat[df_log$padjcat == paste("padj >=", padj_cutoff)] <- NA

upcat<-df_log %>%
  filter(cat == paste(padj_cutoff1,"< padj <", padj_cutoff,";","lfc > +", lfc_cutoff) ) %>%
  nrow()

upperupcat<-df_log %>%
  filter(cat == paste("padj <", padj_cutoff1,";","lfc > +", lfc_cutoff) ) %>%
  nrow()

downcat<- df_log %>%
  filter(cat == paste(padj_cutoff1,"< padj <", padj_cutoff,";","lfc < -", lfc_cutoff) ) %>%
  nrow()

upperdowncat<- df_log %>%
  filter(cat == paste("padj <", padj_cutoff1,";","lfc < -", lfc_cutoff) ) %>%
  nrow()

# Subset significant genes into list of gene names only
df_log <- mutate(df_log, 'gene'=rownames(df))
DEsubset3 <- df_log[  !is.na(df_log$cat) , ]
DEgenelist3 <- DEsubset3$gene #all significant genes meeting padj cutoff <0.05

DEsubset.01 <- df_log %>%
  filter(cat== paste("padj <", padj_cutoff1,";","lfc < -", lfc_cutoff) | cat== paste("padj <", padj_cutoff1,";","lfc > +", lfc_cutoff) )
DEgenelist3.01 <- DEsubset.01$gene #all significant genes meeting padj cutoff <0.01

DEsubset.up01 <- df_log %>%
  filter(cat== paste("padj <", padj_cutoff1,";","lfc > +", lfc_cutoff) )
DEgenelist3.up01 <- DEsubset.up01$gene #all up-regulated genes meeting padj cutoff <0.01

DEsubset.down01 <- df_log %>%
  filter(cat== paste("padj <", padj_cutoff1,";","lfc < -", lfc_cutoff) )
DEgenelist3.down01 <- DEsubset.down01$gene #all down-regulated genes meeting padj cutoff <0.01

DEsubset.up05 <-  df_log[  !is.na(df_log$cat) , ]
DEsubset.up05 <- DEsubset.up05 %>%
  filter(lfccat== paste("lfc > +", lfc_cutoff) )
DEgenelist3.up05 <- DEsubset.up05$gene #all up-regulated genes meeting padj cutoff <0.05

DEsubset.down05 <-  df_log[  !is.na(df_log$cat) , ]
DEsubset.down05 <- DEsubset.down05 %>%
  filter(lfccat== paste("lfc < -", lfc_cutoff) )
DEgenelist3.down05 <- DEsubset.down05$gene #all down-regulated genes meeting padj cutoff <0.05

# Volcano plot font sizes
axisfont <- 14
labelfont <- 5
allfont <- 18

# Reassign rownames
rownames(df_log) <- rownames(df)

# Convert padjcat to factor
df_log$padjcat <- as.factor(df_log$padjcat)

# Replace NA with Not Significant or Substantial
df_log<-df_log %>%
  mutate(cat = if_else(is.na(cat), "Not Significant or Substantial", cat))
df_log$cat<-factor(df_log$cat,
                   levels=c("padj < 0.01 ; lfc < - 1" ,"padj < 0.01 ; lfc > + 1", "0.01 < padj < 0.05 ; lfc < - 1" , "0.01 < padj < 0.05 ; lfc > + 1", "Not Significant or Substantial"))
cols <- c("0.01 < padj < 0.05 ; lfc < - 1" = "#00ace6", "0.01 < padj < 0.05 ; lfc > + 1" = "#ff66cc", "padj < 0.01 ; lfc < - 1" = "#0000ff", "padj < 0.01 ; lfc > + 1" = "#ff0000", "Not Significant or Substantial"="#818181")

# Plot discrete volcano plot
p3<-ggplot(df_log, aes(x=log2FoldChange, y=logp, color=cat)) + 
  geom_point(na.rm = TRUE,alpha=0.5) +
  geom_vline(xintercept=c(-1,1),linetype="dotted") + 
  geom_hline(yintercept=c(2.1,3.37) ,linetype="dotted") + 
  theme(text = element_text(size=allfont),
        axis.text = element_text(size=axisfont)) + 
  ggtitle("April vs. Nov") +
  ylab(expression(paste("\u2013", log[10], "(P-Value)"))) + 
  theme_bw()+
  xlab(expression(paste(log[2],"(FoldChange)"))) +
  labs(color='Genes' ) +
  xlim(-8,8)+
  scale_colour_manual(values = cols, drop=FALSE, name = "Genes", labels = c("P-Value < 0.01 ; LFC < - 1", "P-Value < 0.01 ; LFC > + 1", "0.01 \u2264 P-Value < 0.05 ; LFC < - 1", "0.01 \u2264 P-Value < 0.05 ; LFC > + 1",  "Not Significant or Substantial"))+
  guides(color = guide_legend(override.aes = list(size = 3)))

ylim<-ggplot_build(p3)$layout$panel_scales_y[[1]]$range$range[2]

p3<-p3  + geom_text(x=6.25, y=ylim-2.45*(ylim/8.25), label= upcat, color='#ff66cc') 
p3<-p3 + geom_text(x=-6.25, y=ylim-2.45*(ylim/8.25), label= downcat, color='#00ace6') 
p3<-p3 + geom_text(x=6.25, y=ylim-.7*(ylim/8.25), label= upperupcat, color='#ff0000') 
p3<-p3 + geom_text(x=-6.25, y=ylim-.7*(ylim/8.25), label= upperdowncat, color='#0000ff') 

outfile1<-"/Users/sunshine2.0/Desktop/WORK/output/volcanoAvaN_"
outfile2<-".png"
outfile<-paste(outfile1,padj_cutoff,outfile2, sep="")

ggsave(outfile, 
       dpi=300, 
       device="png",
       width=23, 
       height=14,
       units="cm")

# # Plot continuous volano plot
# ggplot(df_log, aes(x=log2FoldChange, y=logp, color=padj)) + 
#   geom_point(na.rm = TRUE, alpha=0.3) + 
#   geom_vline(xintercept=lfc, linetype="dotted") + 
#   theme(text = element_text(size=allfont),
#         axis.text = element_text(size=axisfont)) + 
#   ggtitle("November vs. April") +
#   ylab("\u2013log10(P-Value)") + 
#   xlab("log2(FoldChange)")  +
#   labs(color='Adj P-Value' )
# 
# outfile1 <- "/home/krwright/Downloads/volcano_continuousNvsA.png"  
# 
# ggsave(outfile1, 
#        dpi=300, 
#        device="png", 
#        width=23, 
#        height=14, 
#        units="cm")


#_________________________________________ FEMALE SUBSETS_________________________________________________________

# JUNE VS APRIL

# Subset the relevant columns from the full dataset and choose the control condition so log2 fold changes are calculated as treatment over control
ddsMF <- ddsAll[ , (ddsAll$month == "A")|(ddsAll$month == "J") ]
dds <- ddsMF[ , (ddsMF$gender == "F") ]
dds$month <- droplevels( dds$month )
dds$month <- relevel(dds$month, ref = "A")
# check if subsetted correctly
#as.data.frame( colData(dds) )

dds <- DESeq(dds)
res <- results(dds, alpha=padj_cutoff)

# Number of adjusted p-values less than 0.1
sum(res$padj < padj_cutoff, na.rm=TRUE)
# summary (default lfcThreshold = 0)
summary(res)
sum(res$baseMean ==0)
#sum(is.na(res$padj))

# Convert to dataframe
df <- as.data.frame(res)

# Make summary table similar to summary(res)
total<- nrow(df)-sum(df$baseMean ==0, na.rm=TRUE)
lowCountsSubset<-df[ ((!is.na(df$pvalue))&(is.na(df$padj))),]
dfsummary4 <-  data.frame("LFCup" = c(round(sum(df$log2FoldChange>0& df$padj<padj_cutoff, na.rm=TRUE), digits=0),round(sum(df$log2FoldChange>0& df$padj<padj_cutoff, na.rm=TRUE)/total*100, digits = 1)), 
                          "LFCdown" = c(round(sum(df$log2FoldChange<0& df$padj<padj_cutoff, na.rm=TRUE), digits = 0),round(sum(df$log2FoldChange<0& df$padj<padj_cutoff, na.rm=TRUE)/total*100, digits = 1)), 
                          "outliers"= c(round(sum(!is.na(df$log2FoldChange)&is.na(df$pvalue)), digits = 0),round(sum(!is.na(df$log2FoldChange)&is.na(df$pvalue))/round(sum(is.na(df$padj)), digits = 0)*100, digits = 1)),
                          "lowCounts"= c(round(sum(!is.na(df$pvalue)&is.na(df$padj)), digits = 0),round(sum(!is.na(df$pvalue)&is.na(df$padj))/round(sum(is.na(df$padj)), digits = 0)*100, digits = 1)),
                          "baseMean"= c(round(sum(df$baseMean ==0, na.rm=TRUE), digits = 0),round(sum(df$baseMean ==0, na.rm=TRUE)/round(sum(is.na(df$padj)), digits = 0)*100, digits = 1)),
                          "total"= c(round(total, digits = 0), round(sum((sum(df$log2FoldChange>0& df$padj<padj_cutoff, na.rm=TRUE)/total*100),(sum(df$log2FoldChange<0& df$padj<padj_cutoff, na.rm=TRUE)/total*100),(sum(!is.na(df$log2FoldChange)&is.na(df$pvalue))/round(sum(is.na(df$padj)), digits = 0)*100),(sum(!is.na(df$pvalue)&is.na(df$padj))/round(sum(is.na(df$padj)), digits = 0)*100),(sum(df$baseMean ==0, na.rm=TRUE)/round(sum(is.na(df$padj)), digits = 0)*100), digits = 1))),
                          "treatment" = "F_JvsA",
                          "maxlowCounts"=c(round(max(lowCountsSubset$baseMean), digits=0),NA))

# Make a data table of df statistics
dfstats4<-data.frame("#NA"=c(round(sum(is.na(df$padj)), digits = 0),NA),
                     "#sig"=c(sum(df$padj < padj_cutoff, na.rm=TRUE),round(sum(df$padj < padj_cutoff, na.rm=TRUE)/(nrow(df)-sum(is.na(df$padj)))*100, digits = 1)),
                     "#nonsig"= c(sum(df$padj > padj_cutoff | df$padj == padj_cutoff, na.rm=TRUE),round(sum(df$padj > padj_cutoff | df$padj == padj_cutoff, na.rm=TRUE)/(nrow(df)-sum(is.na(df$padj)))*100, digits = 1)),
                     "total"= c(nrow(df), NA),
                     "treatment" = "F_JvsA")

# Add -log(pvalue) column, and padj category 
df_log <- mutate(df, logp=-log10(pvalue), padjcat=padj,
                 lfccat=log2FoldChange)

df_log <- transform(df_log, padjcat=case_when(padj < padj_cutoff & padj >= padj_cutoff1 ~ paste(padj_cutoff1,"< padj <", padj_cutoff),
                                              padj >= padj_cutoff ~ paste("padj >=", padj_cutoff),
                                              padj < padj_cutoff1 ~ paste("padj <", padj_cutoff1)))

df_log<- transform(df_log, lfccat= case_when(lfccat > lfc_cutoff ~ paste("lfc > +", lfc_cutoff), 
                                             lfccat < -lfc_cutoff ~ paste("lfc < -", lfc_cutoff),
                                             lfccat <= lfc_cutoff & lfccat >= (-lfc_cutoff) ~ paste("-",lfc_cutoff," =< lfc =< +", lfc_cutoff)))

df_log <- mutate(df_log, cat=paste(padjcat,";",lfccat))
df_log$cat[is.na(df_log$padjcat)] <- NA
df_log$cat[df_log$lfccat == paste("-",lfc_cutoff," =< lfc =< +", lfc_cutoff)] <- NA
df_log$cat[df_log$padjcat == paste("padj >=", padj_cutoff)] <- NA

upcat<-df_log %>%
  filter(cat == paste(padj_cutoff1,"< padj <", padj_cutoff,";","lfc > +", lfc_cutoff) ) %>%
  nrow()

upperupcat<-df_log %>%
  filter(cat == paste("padj <", padj_cutoff1,";","lfc > +", lfc_cutoff) ) %>%
  nrow()

downcat<- df_log %>%
  filter(cat == paste(padj_cutoff1,"< padj <", padj_cutoff,";","lfc < -", lfc_cutoff) ) %>%
  nrow()

upperdowncat<- df_log %>%
  filter(cat == paste("padj <", padj_cutoff1,";","lfc < -", lfc_cutoff) ) %>%
  nrow()

# Subset significant genes into list of gene names only
df_log <- mutate(df_log, 'gene'=rownames(df))
DEsubset4 <- df_log[  !is.na(df_log$cat) , ]
DEgenelist4 <- DEsubset4$gene #all significant genes meeting padj cutoff <0.05

DEsubset.01 <- df_log %>%
  filter(cat== paste("padj <", padj_cutoff1,";","lfc < -", lfc_cutoff) | cat== paste("padj <", padj_cutoff1,";","lfc > +", lfc_cutoff) )
DEgenelist4.01 <- DEsubset.01$gene #all significant genes meeting padj cutoff <0.01

DEgenelist4.01p<-df_log%>% # Obtain list of DE0.01 genes with padj
  dplyr::select(c(gene, padj)) %>%
  filter(gene %in% DEgenelist4.01) %>%
  mutate("month.comparison"="June vs. April")

DEsubset.up01 <- df_log %>%
  filter(cat== paste("padj <", padj_cutoff1,";","lfc > +", lfc_cutoff) )
DEgenelist4.up01 <- DEsubset.up01$gene #all up-regulated genes meeting padj cutoff <0.01

DEsubset.down01 <- df_log %>%
  filter(cat== paste("padj <", padj_cutoff1,";","lfc < -", lfc_cutoff) )
DEgenelist4.down01 <- DEsubset.down01$gene #all down-regulated genes meeting padj cutoff <0.01

DEsubset.up05 <-  df_log[  !is.na(df_log$cat) , ]
DEsubset.up05 <- DEsubset.up05 %>%
  filter(lfccat== paste("lfc > +", lfc_cutoff) )
DEgenelist4.up05 <- DEsubset.up05$gene #all up-regulated genes meeting padj cutoff <0.05

DEsubset.down05 <-  df_log[  !is.na(df_log$cat) , ]
DEsubset.down05 <- DEsubset.down05 %>%
  filter(lfccat== paste("lfc < -", lfc_cutoff) )
DEgenelist4.down05 <- DEsubset.down05$gene #all down-regulated genes meeting padj cutoff <0.05

DE_padj4<-df_log%>% # Obtain list of DE0.01 genes with padj
  dplyr::select(c(gene, padj)) %>%
  mutate("month.comparison"="June vs. April")

df_log_4<-df_log

# Volcano plot font sizes
axisfont <- 14
labelfont <- 5
allfont <- 18

# Reassign rownames
rownames(df_log) <- rownames(df)

# Convert padjcat to factor
df_log$padjcat <- as.factor(df_log$padjcat)

# Replace NA with Not Significant or Substantial
df_log<-df_log %>%
  mutate(cat = if_else(is.na(cat), "Not Significant or Substantial", cat))
df_log$cat<-factor(df_log$cat,
                   levels=c("padj < 0.01 ; lfc < - 1" ,"padj < 0.01 ; lfc > + 1", "0.01 < padj < 0.05 ; lfc < - 1" , "0.01 < padj < 0.05 ; lfc > + 1", "Not Significant or Substantial"))
cols <- c("0.01 < padj < 0.05 ; lfc < - 1" = "#00ace6", "0.01 < padj < 0.05 ; lfc > + 1" = "#ff66cc", "padj < 0.01 ; lfc < - 1" = "#0000ff", "padj < 0.01 ; lfc > + 1" = "#ff0000", "Not Significant or Substantial"="#818181")

# Plot discrete volcano plot
p4<-ggplot(df_log, aes(x=log2FoldChange, y=logp, color=cat)) + 
  geom_point(na.rm = TRUE,alpha=0.5) +
  geom_vline(xintercept=c(-1,1),linetype="dotted") + 
  geom_hline(yintercept=c(4.7,6.2) ,linetype="dotted") + 
  theme(text = element_text(size=allfont),
  axis.text = element_text(size=axisfont)) + 
  ggtitle("Female: June vs. April") +
  ylab(expression(paste("\u2013", log[10], "(P-Value)"))) + 
  theme_bw()+
  xlab(expression(paste(log[2],"(FoldChange)"))) +
  labs(color='Genes' ) +
  xlim(-8,8) +
  scale_colour_manual(values = cols, drop=FALSE, name = "Genes", labels = c("P-Value < 0.01 ; LFC < - 1", "P-Value < 0.01 ; LFC > + 1", "0.01 \u2264 P-Value < 0.05 ; LFC < - 1", "0.01 \u2264 P-Value < 0.05 ; LFC > + 1",  "Not Significant or Substantial"))+
  guides(color = guide_legend(override.aes = list(size = 3)))

ylim<-ggplot_build(p4)$layout$panel_scales_y[[1]]$range$range[2]

p4<-p4  + geom_text(x=6.25, y=ylim-2.45*(ylim/8.25), label= upcat, color='#ff66cc') 
p4<-p4 + geom_text(x=-6.25, y=ylim-2.45*(ylim/8.25), label= downcat, color='#00ace6') 
p4<-p4 + geom_text(x=6.25, y=ylim-.7*(ylim/8.25), label= upperupcat, color='#ff0000') 
p4<-p4 + geom_text(x=-6.25, y=ylim-.7*(ylim/8.25), label= upperdowncat, color='#0000ff') 

outfile1<-"/Users/sunshine2.0/Desktop/WORK/output/volcanoFemaleJvsA_"
outfile2<-".png"
outfile<-paste(outfile1,padj_cutoff,outfile2, sep="")

ggsave(outfile, 
       dpi=300, 
       device="png",
       width=23, 
       height=14,
       units="cm")

# # Plot continuous volcano plot
# ggplot(df_log, aes(x=log2FoldChange, y=logp, color=padj)) + 
#   geom_point(na.rm = TRUE, alpha=0.3) + 
#   geom_vline(xintercept=lfc, linetype="dotted") + 
#   theme(text = element_text(size=allfont),
#         axis.text = element_text(size=axisfont)) + 
#   ggtitle("Female: June vs. April") +
#   ylab("\u2013log10(P-Value)") + 
#   xlab("log2(FoldChange)")  +
#   labs(color='Adj P-Value' ) +
#   xlim(-8,8)
# 
# outfile1 <- "/home/krwright/Downloads/volcano_continuousFemaleJvsA.png"  
# 
# ggsave(outfile1, 
#        dpi=300, 
#        device="png", 
#        width=23, 
#        height=14, 
#        units="cm")

# NOVEMBER VS JUNE 

# Subset the relevant columns from the full dataset and choose the control condition so log2 fold changes are calculated as treatment over control
ddsMF <- ddsAll[ , (ddsAll$month == "N")|(ddsAll$month == "J") ]
dds <- ddsMF[ , (ddsMF$gender == "F") ]
dds$month <- droplevels( dds$month )
dds$month <- relevel(dds$month, ref = "J")
# check if subsetted correctly
#as.data.frame( colData(dds) )

dds <- DESeq(dds)
res <- results(dds, alpha=padj_cutoff)

# Number of adjusted p-values less than 0.1
sum(res$padj < padj_cutoff, na.rm=TRUE)
# summary (default lfcThreshold = 0)
summary(res)
sum(res$baseMean ==0)

# Convert to dataframe
df <- as.data.frame(res)

# Make summary table similar to summary(res)
total<- nrow(df)-sum(df$baseMean ==0, na.rm=TRUE)
lowCountsSubset<-df[ ((!is.na(df$pvalue))&(is.na(df$padj))),]
dfsummary5 <-  data.frame("LFCup" = c(round(sum(df$log2FoldChange>0& df$padj<padj_cutoff, na.rm=TRUE), digits=0),round(sum(df$log2FoldChange>0& df$padj<padj_cutoff, na.rm=TRUE)/total*100, digits = 1)), 
                          "LFCdown" = c(round(sum(df$log2FoldChange<0& df$padj<padj_cutoff, na.rm=TRUE), digits = 0),round(sum(df$log2FoldChange<0& df$padj<padj_cutoff, na.rm=TRUE)/total*100, digits = 1)), 
                          "outliers"= c(round(sum(!is.na(df$log2FoldChange)&is.na(df$pvalue)), digits = 0),round(sum(!is.na(df$log2FoldChange)&is.na(df$pvalue))/round(sum(is.na(df$padj)), digits = 0)*100, digits = 1)),
                          "lowCounts"= c(round(sum(!is.na(df$pvalue)&is.na(df$padj)), digits = 0),round(sum(!is.na(df$pvalue)&is.na(df$padj))/round(sum(is.na(df$padj)), digits = 0)*100, digits = 1)),
                          "baseMean"= c(round(sum(df$baseMean ==0, na.rm=TRUE), digits = 0),round(sum(df$baseMean ==0, na.rm=TRUE)/round(sum(is.na(df$padj)), digits = 0)*100, digits = 1)),
                          "total"= c(round(total, digits = 0), round(sum((sum(df$log2FoldChange>0& df$padj<padj_cutoff, na.rm=TRUE)/total*100),(sum(df$log2FoldChange<0& df$padj<padj_cutoff, na.rm=TRUE)/total*100),(sum(!is.na(df$log2FoldChange)&is.na(df$pvalue))/round(sum(is.na(df$padj)), digits = 0)*100),(sum(!is.na(df$pvalue)&is.na(df$padj))/round(sum(is.na(df$padj)), digits = 0)*100),(sum(df$baseMean ==0, na.rm=TRUE)/round(sum(is.na(df$padj)), digits = 0)*100), digits = 1))),
                          "treatment" = "F_NvsJ",
                          "maxlowCounts"=c(round(max(lowCountsSubset$baseMean), digits=0),NA))

# Make a data table of df statistics
dfstats5<-data.frame("#NA"=c(round(sum(is.na(df$padj)), digits = 0),NA),
                     "#sig"=c(sum(df$padj < padj_cutoff, na.rm=TRUE),round(sum(df$padj < padj_cutoff, na.rm=TRUE)/(nrow(df)-sum(is.na(df$padj)))*100, digits = 1)),
                     "#nonsig"= c(sum(df$padj > padj_cutoff | df$padj == padj_cutoff, na.rm=TRUE),round(sum(df$padj > padj_cutoff | df$padj == padj_cutoff, na.rm=TRUE)/(nrow(df)-sum(is.na(df$padj)))*100, digits = 1)),
                     "total"= c(nrow(df), NA),
                     "treatment" = "F_NvsJ")

# Add -log(pvalue) column, and padj category 
df_log <- mutate(df, logp=-log10(pvalue), padjcat=padj,
                 lfccat=log2FoldChange)

df_log <- transform(df_log, padjcat=case_when(padj < padj_cutoff & padj >= padj_cutoff1 ~ paste(padj_cutoff1,"< padj <", padj_cutoff),
                                              padj >= padj_cutoff ~ paste("padj >=", padj_cutoff),
                                              padj < padj_cutoff1 ~ paste("padj <", padj_cutoff1)))

df_log<- transform(df_log, lfccat= case_when(lfccat > lfc_cutoff ~ paste("lfc > +", lfc_cutoff), 
                                             lfccat < -lfc_cutoff ~ paste("lfc < -", lfc_cutoff),
                                             lfccat <= lfc_cutoff & lfccat >= (-lfc_cutoff) ~ paste("-",lfc_cutoff," =< lfc =< +", lfc_cutoff)))

df_log <- mutate(df_log, cat=paste(padjcat,";",lfccat))
df_log$cat[is.na(df_log$padjcat)] <- NA
df_log$cat[df_log$lfccat == paste("-",lfc_cutoff," =< lfc =< +", lfc_cutoff)] <- NA
df_log$cat[df_log$padjcat == paste("padj >=", padj_cutoff)] <- NA

upcat<-df_log %>%
  filter(cat == paste(padj_cutoff1,"< padj <", padj_cutoff,";","lfc > +", lfc_cutoff) ) %>%
  nrow()

upperupcat<-df_log %>%
  filter(cat == paste("padj <", padj_cutoff1,";","lfc > +", lfc_cutoff) ) %>%
  nrow()

downcat<- df_log %>%
  filter(cat == paste(padj_cutoff1,"< padj <", padj_cutoff,";","lfc < -", lfc_cutoff) ) %>%
  nrow()

upperdowncat<- df_log %>%
  filter(cat == paste("padj <", padj_cutoff1,";","lfc < -", lfc_cutoff) ) %>%
  nrow()

# Subset significant genes into list of gene names only
df_log <- mutate(df_log, 'gene'=rownames(df))
DEsubset5 <- df_log[  !is.na(df_log$cat) , ]
DEgenelist5 <- DEsubset5$gene #all significant genes meeting padj cutoff <0.05

DEsubset.01 <- df_log %>%
  filter(cat== paste("padj <", padj_cutoff1,";","lfc < -", lfc_cutoff) | cat== paste("padj <", padj_cutoff1,";","lfc > +", lfc_cutoff) )
DEgenelist5.01 <- DEsubset.01$gene #all significant genes meeting padj cutoff <0.01

DEgenelist5.01p<-df_log%>% # Obtain list of DE0.01 genes with padj
  dplyr::select(c(gene, padj)) %>%
  filter(gene %in% DEgenelist5.01) %>%
  mutate("month.comparison"="Nov vs. June")

DEsubset.up01 <- df_log %>%
  filter(cat== paste("padj <", padj_cutoff1,";","lfc > +", lfc_cutoff) )
DEgenelist5.up01 <- DEsubset.up01$gene #all up-regulated genes meeting padj cutoff <0.01

DEsubset.down01 <- df_log %>%
  filter(cat== paste("padj <", padj_cutoff1,";","lfc < -", lfc_cutoff) )
DEgenelist5.down01 <- DEsubset.down01$gene #all down-regulated genes meeting padj cutoff <0.01

DEsubset.up05 <-  df_log[  !is.na(df_log$cat) , ]
DEsubset.up05 <- DEsubset.up05 %>%
  filter(lfccat== paste("lfc > +", lfc_cutoff) )
DEgenelist5.up05 <- DEsubset.up05$gene #all up-regulated genes meeting padj cutoff <0.05

DEsubset.down05 <-  df_log[  !is.na(df_log$cat) , ]
DEsubset.down05 <- DEsubset.down05 %>%
  filter(lfccat== paste("lfc < -", lfc_cutoff) )
DEgenelist5.down05 <- DEsubset.down05$gene #all down-regulated genes meeting padj cutoff <0.05

DE_padj5<-df_log%>% # Obtain list of DE0.01 genes with padj
  dplyr::select(c(gene, padj)) %>%
  mutate("month.comparison"="Nov vs. June")

df_log_5<-df_log

# Volcano plot font sizes
axisfont <- 14
labelfont <- 5
allfont <- 18

# Reassign rownames
rownames(df_log) <- rownames(df)

# Convert padjcat to factor
df_log$padjcat <- as.factor(df_log$padjcat)

# Replace NA with Not Significant or Substantial
df_log<-df_log %>%
  mutate(cat = if_else(is.na(cat), "Not Significant or Substantial", cat))
df_log$cat<-factor(df_log$cat,
                   levels=c("padj < 0.01 ; lfc < - 1" ,"padj < 0.01 ; lfc > + 1", "0.01 < padj < 0.05 ; lfc < - 1" , "0.01 < padj < 0.05 ; lfc > + 1", "Not Significant or Substantial"))
cols <- c("0.01 < padj < 0.05 ; lfc < - 1" = "#00ace6", "0.01 < padj < 0.05 ; lfc > + 1" = "#ff66cc", "padj < 0.01 ; lfc < - 1" = "#0000ff", "padj < 0.01 ; lfc > + 1" = "#ff0000", "Not Significant or Substantial"="#818181")

# Plot discrete volcano plot
p5<-ggplot(df_log, aes(x=log2FoldChange, y=logp, color=cat)) + 
  geom_point(na.rm = TRUE,alpha=0.5) +
  geom_vline(xintercept=c(-1,1),linetype="dotted") + 
  geom_hline(yintercept=c(2.83,4) ,linetype="dotted") + 
  theme(text = element_text(size=allfont),
        axis.text = element_text(size=axisfont)) + 
  ggtitle("Female: Nov vs. June") +
  ylab(expression(paste("\u2013", log[10], "(P-Value)"))) + 
  theme_bw()+
  xlab(expression(paste(log[2],"(FoldChange)"))) +
  labs(color='Genes' ) +
  xlim(-8,8)+
  scale_colour_manual(values = cols, drop=FALSE, name = "Genes", labels = c("P-Value < 0.01 ; LFC < - 1", "P-Value < 0.01 ; LFC > + 1", "0.01 \u2264 P-Value < 0.05 ; LFC < - 1", "0.01 \u2264 P-Value < 0.05 ; LFC > + 1",  "Not Significant or Substantial"))+
  guides(color = guide_legend(override.aes = list(size = 3)))

ylim<-ggplot_build(p5)$layout$panel_scales_y[[1]]$range$range[2]

p5<-p5  + geom_text(x=6.25, y=ylim-2.45*(ylim/8.25), label= upcat, color='#ff66cc') 
p5<-p5 + geom_text(x=-6.25, y=ylim-2.45*(ylim/8.25), label= downcat, color='#00ace6') 
p5<-p5 + geom_text(x=6.25, y=ylim-.7*(ylim/8.25), label= upperupcat, color='#ff0000') 
p5<-p5 + geom_text(x=-6.25, y=ylim-.7*(ylim/8.25), label= upperdowncat, color='#0000ff') 

outfile1<-"/Users/sunshine2.0/Desktop/WORK/output/volcanoFemaleNvsJ_"
outfile2<-".png"
outfile<-paste(outfile1,padj_cutoff,outfile2, sep="")

ggsave(outfile, 
       dpi=300, 
       device="png",
       width=23, 
       height=14,
       units="cm")

# # Plot continuous volcano plot
# ggplot(df_log, aes(x=log2FoldChange, y=logp, color=padj)) + 
#   geom_point(na.rm = TRUE, alpha=0.3) + 
#   geom_vline(xintercept=lfc, linetype="dotted") + 
#   theme(text = element_text(size=allfont),
#         axis.text = element_text(size=axisfont)) + 
#   ggtitle("Female: November vs. June") +
#   ylab("\u2013log10(P-Value)") + 
#   xlab("log2(FoldChange)")  +
#   labs(color='Adj P-Value' ) +
#   xlim(-8,8)
# 
# outfile1 <- "/home/krwright/Downloads/volcano_continuousFemaleNvsJ.png"  
# 
# ggsave(outfile1, 
#        dpi=300, 
#        device="png", 
#        width=23, 
#        height=14, 
#        units="cm")

# APRIL VS NOVEMBER

# Subset the relevant columns from the full dataset and choose the control condition so log2 fold changes are calculated as treatment over control
ddsMF <- ddsAll[ , (ddsAll$month == "N")|(ddsAll$month == "A") ]
dds <- ddsMF[ , (ddsMF$gender == "F") ]
dds$month <- droplevels( dds$month )
dds$month <- relevel(dds$month, ref = "N")
# check if subsetted correctly
#as.data.frame( colData(dds) )

dds <- DESeq(dds)
res <- results(dds, alpha=padj_cutoff)

# Number of adjusted p-values less than 0.1
sum(res$padj < padj_cutoff, na.rm=TRUE)
# summary (default lfcThreshold = 0)
summary(res)
sum(res$baseMean ==0)

# Convert to dataframe
df <- as.data.frame(res)

# Make summary table similar to summary(res)
total<- nrow(df)-sum(df$baseMean ==0, na.rm=TRUE)
lowCountsSubset<-df[ ((!is.na(df$pvalue))&(is.na(df$padj))),]
dfsummary6 <-  data.frame("LFCup" = c(round(sum(df$log2FoldChange>0& df$padj<padj_cutoff, na.rm=TRUE), digits=0),round(sum(df$log2FoldChange>0& df$padj<padj_cutoff, na.rm=TRUE)/total*100, digits = 1)), 
                          "LFCdown" = c(round(sum(df$log2FoldChange<0& df$padj<padj_cutoff, na.rm=TRUE), digits = 0),round(sum(df$log2FoldChange<0& df$padj<padj_cutoff, na.rm=TRUE)/total*100, digits = 1)), 
                          "outliers"= c(round(sum(!is.na(df$log2FoldChange)&is.na(df$pvalue)), digits = 0),round(sum(!is.na(df$log2FoldChange)&is.na(df$pvalue))/round(sum(is.na(df$padj)), digits = 0)*100, digits = 1)),
                          "lowCounts"= c(round(sum(!is.na(df$pvalue)&is.na(df$padj)), digits = 0),round(sum(!is.na(df$pvalue)&is.na(df$padj))/round(sum(is.na(df$padj)), digits = 0)*100, digits = 1)),
                          "baseMean"= c(round(sum(df$baseMean ==0, na.rm=TRUE), digits = 0),round(sum(df$baseMean ==0, na.rm=TRUE)/round(sum(is.na(df$padj)), digits = 0)*100, digits = 1)),
                          "total"= c(round(total, digits = 0), round(sum((sum(df$log2FoldChange>0& df$padj<padj_cutoff, na.rm=TRUE)/total*100),(sum(df$log2FoldChange<0& df$padj<padj_cutoff, na.rm=TRUE)/total*100),(sum(!is.na(df$log2FoldChange)&is.na(df$pvalue))/round(sum(is.na(df$padj)), digits = 0)*100),(sum(!is.na(df$pvalue)&is.na(df$padj))/round(sum(is.na(df$padj)), digits = 0)*100),(sum(df$baseMean ==0, na.rm=TRUE)/round(sum(is.na(df$padj)), digits = 0)*100), digits = 1))),
                          "treatment" = "F_NvsA",
                          "maxlowCounts"=c(round(max(lowCountsSubset$baseMean), digits=0),NA))

# Make a data table of df statistics
dfstats6<-data.frame("#NA"=c(round(sum(is.na(df$padj)), digits = 0),NA),
                     "#sig"=c(sum(df$padj < padj_cutoff, na.rm=TRUE),round(sum(df$padj < padj_cutoff, na.rm=TRUE)/(nrow(df)-sum(is.na(df$padj)))*100, digits = 1)),
                     "#nonsig"= c(sum(df$padj > padj_cutoff | df$padj == padj_cutoff, na.rm=TRUE),round(sum(df$padj > padj_cutoff | df$padj == padj_cutoff, na.rm=TRUE)/(nrow(df)-sum(is.na(df$padj)))*100, digits = 1)),
                     "total"= c(nrow(df), NA),
                     "treatment" = "F_NvsA")

# Add -log(pvalue) column, and padj category 
df_log <- mutate(df, logp=-log10(pvalue), padjcat=padj,
                 lfccat=log2FoldChange)

df_log <- transform(df_log, padjcat=case_when(padj < padj_cutoff & padj >= padj_cutoff1 ~ paste(padj_cutoff1,"< padj <", padj_cutoff),
                                              padj >= padj_cutoff ~ paste("padj >=", padj_cutoff),
                                              padj < padj_cutoff1 ~ paste("padj <", padj_cutoff1)))

df_log<- transform(df_log, lfccat= case_when(lfccat > lfc_cutoff ~ paste("lfc > +", lfc_cutoff), 
                                             lfccat < -lfc_cutoff ~ paste("lfc < -", lfc_cutoff),
                                             lfccat <= lfc_cutoff & lfccat >= (-lfc_cutoff) ~ paste("-",lfc_cutoff," =< lfc =< +", lfc_cutoff)))

df_log <- mutate(df_log, cat=paste(padjcat,";",lfccat))
df_log$cat[is.na(df_log$padjcat)] <- NA
df_log$cat[df_log$lfccat == paste("-",lfc_cutoff," =< lfc =< +", lfc_cutoff)] <- NA
df_log$cat[df_log$padjcat == paste("padj >=", padj_cutoff)] <- NA

upcat<-df_log %>%
  filter(cat == paste(padj_cutoff1,"< padj <", padj_cutoff,";","lfc > +", lfc_cutoff) ) %>%
  nrow()

upperupcat<-df_log %>%
  filter(cat == paste("padj <", padj_cutoff1,";","lfc > +", lfc_cutoff) ) %>%
  nrow()

downcat<- df_log %>%
  filter(cat == paste(padj_cutoff1,"< padj <", padj_cutoff,";","lfc < -", lfc_cutoff) ) %>%
  nrow()

upperdowncat<- df_log %>%
  filter(cat == paste("padj <", padj_cutoff1,";","lfc < -", lfc_cutoff) ) %>%
  nrow()

# Subset significant genes into list of gene names only
df_log <- mutate(df_log, 'gene'=rownames(df))
DEsubset6 <- df_log[  !is.na(df_log$cat) , ]
DEgenelist6 <- DEsubset6$gene #all significant genes meeting padj cutoff <0.05

DEsubset.01 <- df_log %>%
  filter(cat== paste("padj <", padj_cutoff1,";","lfc < -", lfc_cutoff) | cat== paste("padj <", padj_cutoff1,";","lfc > +", lfc_cutoff) )
DEgenelist6.01 <- DEsubset.01$gene #all significant genes meeting padj cutoff <0.01

DEgenelist6.01p<-df_log%>% # Obtain list of DE0.01 genes with padj
  dplyr::select(c(gene, padj)) %>%
  filter(gene %in% DEgenelist6.01) %>%
  mutate("month.comparison"="April vs. Nov")

DEsubset.up01 <- df_log %>%
  filter(cat== paste("padj <", padj_cutoff1,";","lfc > +", lfc_cutoff) )
DEgenelist6.up01 <- DEsubset.up01$gene #all up-regulated genes meeting padj cutoff <0.01

DEsubset.down01 <- df_log %>%
  filter(cat== paste("padj <", padj_cutoff1,";","lfc < -", lfc_cutoff) )
DEgenelist6.down01 <- DEsubset.down01$gene #all down-regulated genes meeting padj cutoff <0.01

DEsubset.up05 <-  df_log[  !is.na(df_log$cat) , ]
DEsubset.up05 <- DEsubset.up05 %>%
  filter(lfccat== paste("lfc > +", lfc_cutoff) )
DEgenelist6.up05 <- DEsubset.up05$gene #all up-regulated genes meeting padj cutoff <0.05

DEsubset.down05 <-  df_log[  !is.na(df_log$cat) , ]
DEsubset.down05 <- DEsubset.down05 %>%
  filter(lfccat== paste("lfc < -", lfc_cutoff) )
DEgenelist6.down05 <- DEsubset.down05$gene #all down-regulated genes meeting padj cutoff <0.05

DE_padj6<-df_log%>% # Obtain list of DE0.01 genes with padj
  dplyr::select(c(gene, padj)) %>%
  mutate("month.comparison"="April vs. Nov")

df_log_6<-df_log

# Volcano plot font sizes
axisfont <- 14
labelfont <- 5
allfont <- 18

# Reassign rownames
rownames(df_log) <- rownames(df)

# Convert padjcat to factor
df_log$padjcat <- as.factor(df_log$padjcat)

# Replace NA with Not Significant or Substantial
df_log<-df_log %>%
  mutate(cat = if_else(is.na(cat), "Not Significant or Substantial", cat))
df_log$cat<-factor(df_log$cat,
                   levels=c("padj < 0.01 ; lfc < - 1" ,"padj < 0.01 ; lfc > + 1", "0.01 < padj < 0.05 ; lfc < - 1" , "0.01 < padj < 0.05 ; lfc > + 1", "Not Significant or Substantial"))
cols <- c("0.01 < padj < 0.05 ; lfc < - 1" = "#00ace6", "0.01 < padj < 0.05 ; lfc > + 1" = "#ff66cc", "padj < 0.01 ; lfc < - 1" = "#0000ff", "padj < 0.01 ; lfc > + 1" = "#ff0000", "Not Significant or Substantial"="#818181")

# Plot discrete volcano plot
p6<-ggplot(df_log, aes(x=log2FoldChange, y=logp, color=cat)) + 
  geom_point(na.rm = TRUE,alpha=0.5) +
  geom_vline(xintercept=c(-1,1),linetype="dotted") + 
  geom_hline(yintercept=c(2.88,3.88) ,linetype="dotted") + 
  theme(text = element_text(size=allfont),
        axis.text = element_text(size=axisfont)) + 
  ggtitle("Female: April vs. Nov") +
  ylab(expression(paste("\u2013", log[10], "(P-Value)"))) + 
  theme_bw()+
  xlab(expression(paste(log[2],"(FoldChange)"))) +
  labs(color='Genes' ) +
  xlim(-8,8)+
  scale_colour_manual(values = cols, drop=FALSE, name = "Genes", labels = c("P-Value < 0.01 ; LFC < - 1", "P-Value < 0.01 ; LFC > + 1", "0.01 \u2264 P-Value < 0.05 ; LFC < - 1", "0.01 \u2264 P-Value < 0.05 ; LFC > + 1",  "Not Significant or Substantial"))+
  guides(color = guide_legend(override.aes = list(size = 3)))

ylim<-ggplot_build(p6)$layout$panel_scales_y[[1]]$range$range[2]

p6<-p6  + geom_text(x=6.25, y=ylim-2.45*(ylim/8.25), label= upcat, color='#ff66cc') 
p6<-p6 + geom_text(x=-6.25, y=ylim-2.45*(ylim/8.25), label= downcat, color='#00ace6') 
p6<-p6 + geom_text(x=6.25, y=ylim-.7*(ylim/8.25), label= upperupcat, color='#ff0000') 
p6<-p6 + geom_text(x=-6.25, y=ylim-.7*(ylim/8.25), label= upperdowncat, color='#0000ff') 

outfile1<-"/Users/sunshine2.0/Desktop/WORK/output/volcanoFemaleAvaN_"
outfile2<-".png"
outfile<-paste(outfile1,padj_cutoff,outfile2, sep="")

ggsave(outfile, 
       dpi=300, 
       device="png",
       width=23, 
       height=14,
       units="cm")

# 
# # Plot continuous volacno plot
# ggplot(df_log, aes(x=log2FoldChange, y=logp, color=padj)) + 
#   geom_point(na.rm = TRUE, alpha=0.3) + 
#   geom_vline(xintercept=lfc, linetype="dotted") + 
#   theme(text = element_text(size=allfont),
#         axis.text = element_text(size=axisfont)) + 
#   ggtitle("Female: November vs. April") +
#   ylab("\u2013log10(P-Value)") + 
#   xlab("log2(FoldChange)")  +
#   labs(color='Adj P-Value' ) +
#   xlim(-10,10)
# 
# outfile1 <- "/home/krwright/Downloads/volcano_continuousFemaleNvsA.png"  
# 
# ggsave(outfile1, 
#        dpi=300, 
#        device="png", 
#        width=23, 
#        height=14, 
#        units="cm")



#_________________________________________ MALE SUBSETS_________________________________________________________

# JUNE VS APRIL

# Subset the relevant columns from the full dataset and choose the control condition so log2 fold changes are calculated as treatment over control
ddsMF <- ddsAll[ , (ddsAll$month == "A")|(ddsAll$month == "J") ]
dds <- ddsMF[ , (ddsMF$gender == "M") ]
dds$month <- droplevels( dds$month )
dds$month <- relevel(dds$month, ref = "A")
# check if subsetted correctly
#as.data.frame( colData(dds) )

dds <- DESeq(dds)
res <- results(dds, alpha=padj_cutoff)

# Number of adjusted p-values less than 0.1
sum(res$padj < padj_cutoff, na.rm=TRUE)
# summary (default lfcThreshold = 0)
summary(res)
sum(res$baseMean ==0)

# Convert to dataframe
df <- as.data.frame(res)

# Make summary table similar to summary(res)
total<- nrow(df)-sum(df$baseMean ==0, na.rm=TRUE)
lowCountsSubset<-df[ ((!is.na(df$pvalue))&(is.na(df$padj))),]
dfsummary7 <-  data.frame("LFCup" = c(round(sum(df$log2FoldChange>0& df$padj<padj_cutoff, na.rm=TRUE), digits=0),round(sum(df$log2FoldChange>0& df$padj<padj_cutoff, na.rm=TRUE)/total*100, digits = 1)), 
                          "LFCdown" = c(round(sum(df$log2FoldChange<0& df$padj<padj_cutoff, na.rm=TRUE), digits = 0),round(sum(df$log2FoldChange<0& df$padj<padj_cutoff, na.rm=TRUE)/total*100, digits = 1)), 
                          "outliers"= c(round(sum(!is.na(df$log2FoldChange)&is.na(df$pvalue)), digits = 0),round(sum(!is.na(df$log2FoldChange)&is.na(df$pvalue))/round(sum(is.na(df$padj)), digits = 0)*100, digits = 1)),
                          "lowCounts"= c(round(sum(!is.na(df$pvalue)&is.na(df$padj)), digits = 0),round(sum(!is.na(df$pvalue)&is.na(df$padj))/round(sum(is.na(df$padj)), digits = 0)*100, digits = 1)),
                          "baseMean"= c(round(sum(df$baseMean ==0, na.rm=TRUE), digits = 0),round(sum(df$baseMean ==0, na.rm=TRUE)/round(sum(is.na(df$padj)), digits = 0)*100, digits = 1)),
                          "total"= c(round(total, digits = 0), round(sum((sum(df$log2FoldChange>0& df$padj<padj_cutoff, na.rm=TRUE)/total*100),(sum(df$log2FoldChange<0& df$padj<padj_cutoff, na.rm=TRUE)/total*100),(sum(!is.na(df$log2FoldChange)&is.na(df$pvalue))/round(sum(is.na(df$padj)), digits = 0)*100),(sum(!is.na(df$pvalue)&is.na(df$padj))/round(sum(is.na(df$padj)), digits = 0)*100),(sum(df$baseMean ==0, na.rm=TRUE)/round(sum(is.na(df$padj)), digits = 0)*100), digits = 1))),
                          "treatment" = "M_JvsA",
                          "maxlowCounts"=c(round(max(lowCountsSubset$baseMean), digits=0),NA))

# Make a data table of df statistics
dfstats7<-data.frame("#NA"=c(round(sum(is.na(df$padj)), digits = 0),NA),
                     "#sig"=c(sum(df$padj < padj_cutoff, na.rm=TRUE),round(sum(df$padj < padj_cutoff, na.rm=TRUE)/(nrow(df)-sum(is.na(df$padj)))*100, digits = 1)),
                     "#nonsig"= c(sum(df$padj > padj_cutoff | df$padj == padj_cutoff, na.rm=TRUE),round(sum(df$padj > padj_cutoff | df$padj == padj_cutoff, na.rm=TRUE)/(nrow(df)-sum(is.na(df$padj)))*100, digits = 1)),
                     "total"= c(nrow(df), NA),
                     "treatment" = "M_JvsA")

# Add -log(pvalue) column, and padj category 
df_log <- mutate(df, logp=-log10(pvalue), padjcat=padj,
                 lfccat=log2FoldChange)

df_log <- transform(df_log, padjcat=case_when(padj < padj_cutoff & padj >= padj_cutoff1 ~ paste(padj_cutoff1,"< padj <", padj_cutoff),
                                              padj >= padj_cutoff ~ paste("padj >=", padj_cutoff),
                                              padj < padj_cutoff1 ~ paste("padj <", padj_cutoff1)))

df_log<- transform(df_log, lfccat= case_when(lfccat > lfc_cutoff ~ paste("lfc > +", lfc_cutoff), 
                                             lfccat < -lfc_cutoff ~ paste("lfc < -", lfc_cutoff),
                                             lfccat <= lfc_cutoff & lfccat >= (-lfc_cutoff) ~ paste("-",lfc_cutoff," =< lfc =< +", lfc_cutoff)))

df_log <- mutate(df_log, cat=paste(padjcat,";",lfccat))
df_log$cat[is.na(df_log$padjcat)] <- NA
df_log$cat[df_log$lfccat == paste("-",lfc_cutoff," =< lfc =< +", lfc_cutoff)] <- NA
df_log$cat[df_log$padjcat == paste("padj >=", padj_cutoff)] <- NA

upcat<-df_log %>%
  filter(cat == paste(padj_cutoff1,"< padj <", padj_cutoff,";","lfc > +", lfc_cutoff) ) %>%
  nrow()

upperupcat<-df_log %>%
  filter(cat == paste("padj <", padj_cutoff1,";","lfc > +", lfc_cutoff) ) %>%
  nrow()

downcat<- df_log %>%
  filter(cat == paste(padj_cutoff1,"< padj <", padj_cutoff,";","lfc < -", lfc_cutoff) ) %>%
  nrow()

upperdowncat<- df_log %>%
  filter(cat == paste("padj <", padj_cutoff1,";","lfc < -", lfc_cutoff) ) %>%
  nrow()

# Subset significant genes into list of gene names only
df_log <- mutate(df_log, 'gene'=rownames(df))
DEsubset7 <- df_log[  !is.na(df_log$cat) , ]
DEgenelist7 <- DEsubset7$gene #all significant genes meeting padj cutoff <0.05

DEsubset.01 <- df_log %>%
  filter(cat== paste("padj <", padj_cutoff1,";","lfc < -", lfc_cutoff) | cat== paste("padj <", padj_cutoff1,";","lfc > +", lfc_cutoff) )
DEgenelist7.01 <- DEsubset.01$gene #all significant genes meeting padj cutoff <0.01

DEgenelist7.01p<-df_log%>% # Obtain list of DE0.01 genes with padj
  dplyr::select(c(gene, padj)) %>%
  filter(gene %in% DEgenelist7.01) %>%
  mutate("month.comparison"="June vs. April")

DEsubset.up01 <- df_log %>%
  filter(cat== paste("padj <", padj_cutoff1,";","lfc > +", lfc_cutoff) )
DEgenelist7.up01 <- DEsubset.up01$gene #all up-regulated genes meeting padj cutoff <0.01

DEsubset.down01 <- df_log %>%
  filter(cat== paste("padj <", padj_cutoff1,";","lfc < -", lfc_cutoff) )
DEgenelist7.down01 <- DEsubset.down01$gene #all down-regulated genes meeting padj cutoff <0.01

DEsubset.up05 <-  df_log[  !is.na(df_log$cat) , ]
DEsubset.up05 <- DEsubset.up05 %>%
  filter(lfccat== paste("lfc > +", lfc_cutoff) )
DEgenelist7.up05 <- DEsubset.up05$gene #all up-regulated genes meeting padj cutoff <0.05

DEsubset.down05 <-  df_log[  !is.na(df_log$cat) , ]
DEsubset.down05 <- DEsubset.down05 %>%
  filter(lfccat== paste("lfc < -", lfc_cutoff) )
DEgenelist7.down05 <- DEsubset.down05$gene #all down-regulated genes meeting padj cutoff <0.05

DE_padj7<-df_log%>% # Obtain list of DE0.01 genes with padj
  dplyr::select(c(gene, padj)) %>%
  mutate("month.comparison"="June vs. April")

df_log_7<-df_log

# Volcano plot font sizes
axisfont <- 14
labelfont <- 5
allfont <- 18

# Reassign rownames
rownames(df_log) <- rownames(df)

# Convert padjcat to factor
df_log$padjcat <- as.factor(df_log$padjcat)

# Replace NA with Not Significant or Substantial
df_log<-df_log %>%
  mutate(cat = if_else(is.na(cat), "Not Significant or Substantial", cat))
df_log$cat<-factor(df_log$cat,
                   levels=c("padj < 0.01 ; lfc < - 1" ,"padj < 0.01 ; lfc > + 1", "0.01 < padj < 0.05 ; lfc < - 1" , "0.01 < padj < 0.05 ; lfc > + 1", "Not Significant or Substantial"))
cols <- c("0.01 < padj < 0.05 ; lfc < - 1" = "#00ace6", "0.01 < padj < 0.05 ; lfc > + 1" = "#ff66cc", "padj < 0.01 ; lfc < - 1" = "#0000ff", "padj < 0.01 ; lfc > + 1" = "#ff0000", "Not Significant or Substantial"="#818181")

# Plot discrete volcano plot
p7<-ggplot(df_log, aes(x=log2FoldChange, y=logp, color=cat)) + 
  geom_point(na.rm = TRUE,alpha=0.5) +
  geom_vline(xintercept=c(-1,1),linetype="dotted") + 
  geom_hline(yintercept=c(2.43,3.68) ,linetype="dotted") + 
  theme(text = element_text(size=allfont),
        axis.text = element_text(size=axisfont)) + 
  ggtitle("Male: June vs. April") +
  ylab(expression(paste("\u2013", log[10], "(P-Value)"))) + 
  theme_bw()+
  xlab(expression(paste(log[2],"(FoldChange)"))) +
  labs(color='Genes' ) +
  xlim(-8,8)+
  scale_colour_manual(values = cols, drop=FALSE, name = "Genes", labels = c("P-Value < 0.01 ; LFC < - 1", "P-Value < 0.01 ; LFC > + 1", "0.01 \u2264 P-Value < 0.05 ; LFC < - 1", "0.01 \u2264 P-Value < 0.05 ; LFC > + 1",  "Not Significant or Substantial"))+
  guides(color = guide_legend(override.aes = list(size = 3)))

ylim<-ggplot_build(p7)$layout$panel_scales_y[[1]]$range$range[2]

p7<-p7  + geom_text(x=6.25, y=ylim-2.45*(ylim/8.25), label= upcat, color='#ff66cc') 
p7<-p7 + geom_text(x=-6.25, y=ylim-2.45*(ylim/8.25), label= downcat, color='#00ace6') 
p7<-p7 + geom_text(x=6.25, y=ylim-.7*(ylim/8.25), label= upperupcat, color='#ff0000') 
p7<-p7 + geom_text(x=-6.25, y=ylim-.7*(ylim/8.25), label= upperdowncat, color='#0000ff') 

  outfile1<-"/Users/sunshine2.0/Desktop/WORK/output/volcanoMaleJvsA_"
  outfile2<-".png"
  outfile<-paste(outfile1,padj_cutoff,outfile2, sep="")

ggsave(outfile, 
       dpi=300, 
       device="png",
       width=23, 
       height=14,
       units="cm")

# Plot continuous volcano plot
# ggplot(df_log, aes(x=log2FoldChange, y=logp, color=padj)) + 
#   geom_point(na.rm = TRUE, alpha=0.3) + 
#   geom_vline(xintercept=lfc, linetype="dotted") + 
#   theme(text = element_text(size=allfont),
#         axis.text = element_text(size=axisfont)) + 
#   ggtitle("Male: June vs. April") +
#   ylab("\u2013log10(P-Value)") + 
#   xlab("log2(FoldChange)")  +
#   labs(color='Adj P-Value' )
# 
# outfile1 <- "/home/krwright/Downloads/volcano_continuousMaleJvsA.png"  
# 
# ggsave(outfile1, 
#        dpi=300, 
#        device="png", 
#        width=23, 
#        height=14, 
#        units="cm")

# NOVEMBER VS JUNE 

# Subset the relevant columns from the full dataset and choose the control condition so log2 fold changes are calculated as treatment over control
ddsMF <- ddsAll[ , (ddsAll$month == "N")|(ddsAll$month == "J") ]
dds <- ddsMF[ , (ddsMF$gender == "M") ]
dds$month <- droplevels( dds$month )
dds$month <- relevel(dds$month, ref = "J")
# check if subsetted correctly
#as.data.frame( colData(dds) )

dds <- DESeq(dds)
res <- results(dds, alpha=padj_cutoff)

# Number of adjusted p-values less than 0.1
sum(res$padj < padj_cutoff, na.rm=TRUE)
# summary (default lfcThreshold = 0)
summary(res)
sum(res$baseMean ==0)

# Convert to dataframe
df <- as.data.frame(res)

# Make summary table similar to summary(res)
total<- nrow(df)-sum(df$baseMean ==0, na.rm=TRUE)
lowCountsSubset<-df[ ((!is.na(df$pvalue))&(is.na(df$padj))),]
dfsummary8 <-  data.frame("LFCup" = c(round(sum(df$log2FoldChange>0& df$padj<padj_cutoff, na.rm=TRUE), digits=0),round(sum(df$log2FoldChange>0& df$padj<padj_cutoff, na.rm=TRUE)/total*100, digits = 1)), 
                          "LFCdown" = c(round(sum(df$log2FoldChange<0& df$padj<padj_cutoff, na.rm=TRUE), digits = 0),round(sum(df$log2FoldChange<0& df$padj<padj_cutoff, na.rm=TRUE)/total*100, digits = 1)), 
                          "outliers"= c(round(sum(!is.na(df$log2FoldChange)&is.na(df$pvalue)), digits = 0),round(sum(!is.na(df$log2FoldChange)&is.na(df$pvalue))/round(sum(is.na(df$padj)), digits = 0)*100, digits = 1)),
                          "lowCounts"= c(round(sum(!is.na(df$pvalue)&is.na(df$padj)), digits = 0),round(sum(!is.na(df$pvalue)&is.na(df$padj))/round(sum(is.na(df$padj)), digits = 0)*100, digits = 1)),
                          "baseMean"= c(round(sum(df$baseMean ==0, na.rm=TRUE), digits = 0),round(sum(df$baseMean ==0, na.rm=TRUE)/round(sum(is.na(df$padj)), digits = 0)*100, digits = 1)),
                          "total"= c(round(total, digits = 0), round(sum((sum(df$log2FoldChange>0& df$padj<padj_cutoff, na.rm=TRUE)/total*100),(sum(df$log2FoldChange<0& df$padj<padj_cutoff, na.rm=TRUE)/total*100),(sum(!is.na(df$log2FoldChange)&is.na(df$pvalue))/round(sum(is.na(df$padj)), digits = 0)*100),(sum(!is.na(df$pvalue)&is.na(df$padj))/round(sum(is.na(df$padj)), digits = 0)*100),(sum(df$baseMean ==0, na.rm=TRUE)/round(sum(is.na(df$padj)), digits = 0)*100), digits = 1))),
                          "treatment" = "M_NvsJ",
                          "maxlowCounts"=c(round(max(lowCountsSubset$baseMean), digits=0),NA))

# Make a data table of df statistics
dfstats8<-data.frame("#NA"=c(round(sum(is.na(df$padj)), digits = 0),NA),
                     "#sig"=c(sum(df$padj < padj_cutoff, na.rm=TRUE),round(sum(df$padj < padj_cutoff, na.rm=TRUE)/(nrow(df)-sum(is.na(df$padj)))*100, digits = 1)),
                     "#nonsig"= c(sum(df$padj > padj_cutoff | df$padj == padj_cutoff, na.rm=TRUE),round(sum(df$padj > padj_cutoff | df$padj == padj_cutoff, na.rm=TRUE)/(nrow(df)-sum(is.na(df$padj)))*100, digits = 1)),
                     "total"= c(nrow(df), NA),
                     "treatment" = "M_NvsJ")

# Add -log(pvalue) column, and padj category 
df_log <- mutate(df, logp=-log10(pvalue), padjcat=padj,
                 lfccat=log2FoldChange)

df_log <- transform(df_log, padjcat=case_when(padj < padj_cutoff & padj >= padj_cutoff1 ~ paste(padj_cutoff1,"< padj <", padj_cutoff),
                                              padj >= padj_cutoff ~ paste("padj >=", padj_cutoff),
                                              padj < padj_cutoff1 ~ paste("padj <", padj_cutoff1)))

df_log<- transform(df_log, lfccat= case_when(lfccat > lfc_cutoff ~ paste("lfc > +", lfc_cutoff), 
                                             lfccat < -lfc_cutoff ~ paste("lfc < -", lfc_cutoff),
                                             lfccat <= lfc_cutoff & lfccat >= (-lfc_cutoff) ~ paste("-",lfc_cutoff," =< lfc =< +", lfc_cutoff)))

df_log <- mutate(df_log, cat=paste(padjcat,";",lfccat))
df_log$cat[is.na(df_log$padjcat)] <- NA
df_log$cat[df_log$lfccat == paste("-",lfc_cutoff," =< lfc =< +", lfc_cutoff)] <- NA
df_log$cat[df_log$padjcat == paste("padj >=", padj_cutoff)] <- NA

upcat<-df_log %>%
  filter(cat == paste(padj_cutoff1,"< padj <", padj_cutoff,";","lfc > +", lfc_cutoff) ) %>%
  nrow()

upperupcat<-df_log %>%
  filter(cat == paste("padj <", padj_cutoff1,";","lfc > +", lfc_cutoff) ) %>%
  nrow()

downcat<- df_log %>%
  filter(cat == paste(padj_cutoff1,"< padj <", padj_cutoff,";","lfc < -", lfc_cutoff) ) %>%
  nrow()

upperdowncat<- df_log %>%
  filter(cat == paste("padj <", padj_cutoff1,";","lfc < -", lfc_cutoff) ) %>%
  nrow()

# Subset significant genes into list of gene names only
df_log <- mutate(df_log, 'gene'=rownames(df))
DEsubset8 <- df_log[  !is.na(df_log$cat) , ]
DEgenelist8 <- DEsubset8$gene #all significant genes meeting padj cutoff <0.05

DEsubset.01 <- df_log %>%
  filter(cat== paste("padj <", padj_cutoff1,";","lfc < -", lfc_cutoff) | cat== paste("padj <", padj_cutoff1,";","lfc > +", lfc_cutoff) )
DEgenelist8.01 <- DEsubset.01$gene #all significant genes meeting padj cutoff <0.01

DEgenelist8.01p<-df_log%>% # Obtain list of DE0.01 genes with padj
  dplyr::select(c(gene, padj)) %>%
  filter(gene %in% DEgenelist8.01) %>%
  mutate("month.comparison"="Nov vs. June")

DEsubset.up01 <- df_log %>%
  filter(cat== paste("padj <", padj_cutoff1,";","lfc > +", lfc_cutoff) )
DEgenelist8.up01 <- DEsubset.up01$gene #all up-regulated genes meeting padj cutoff <0.01

DEsubset.down01 <- df_log %>%
  filter(cat== paste("padj <", padj_cutoff1,";","lfc < -", lfc_cutoff) )
DEgenelist8.down01 <- DEsubset.down01$gene #all down-regulated genes meeting padj cutoff <0.01

DEsubset.up05 <-  df_log[  !is.na(df_log$cat) , ]
DEsubset.up05 <- DEsubset.up05 %>%
  filter(lfccat== paste("lfc > +", lfc_cutoff) )
DEgenelist8.up05 <- DEsubset.up05$gene #all up-regulated genes meeting padj cutoff <0.05

DEsubset.down05 <-  df_log[  !is.na(df_log$cat) , ]
DEsubset.down05 <- DEsubset.down05 %>%
  filter(lfccat== paste("lfc < -", lfc_cutoff) )
DEgenelist8.down05 <- DEsubset.down05$gene #all down-regulated genes meeting padj cutoff <0.05

DE_padj8<-df_log%>% # Obtain list of DE0.01 genes with padj
  dplyr::select(c(gene, padj)) %>%
  mutate("month.comparison"="Nov vs. June")

df_log_8<-df_log

# Volcano plot font sizes
axisfont <- 14
labelfont <- 5
allfont <- 18

# Reassign rownames
rownames(df_log) <- rownames(df)

# Convert padjcat to factor
df_log$padjcat <- as.factor(df_log$padjcat)

# Replace NA with Not Significant or Substantial
df_log<-df_log %>%
  mutate(cat = if_else(is.na(cat), "Not Significant or Substantial", cat))
df_log$cat<-factor(df_log$cat,
                   levels=c("padj < 0.01 ; lfc < - 1" ,"padj < 0.01 ; lfc > + 1", "0.01 < padj < 0.05 ; lfc < - 1" , "0.01 < padj < 0.05 ; lfc > + 1", "Not Significant or Substantial"))
cols <- c("0.01 < padj < 0.05 ; lfc < - 1" = "#00ace6", "0.01 < padj < 0.05 ; lfc > + 1" = "#ff66cc", "padj < 0.01 ; lfc < - 1" = "#0000ff", "padj < 0.01 ; lfc > + 1" = "#ff0000", "Not Significant or Substantial"="#818181")

# Plot discrete volcano plot
p8<-ggplot(df_log, aes(x=log2FoldChange, y=logp, color=cat)) + 
  geom_point(na.rm = TRUE,alpha=0.5) +
  geom_vline(xintercept=c(-1,1),linetype="dotted") + 
  geom_hline(yintercept=c(2.58,3.78) ,linetype="dotted") + 
  theme(text = element_text(size=allfont),
        axis.text = element_text(size=axisfont)) + 
  ggtitle("Male: Nov vs. June") +
  ylab(expression(paste("\u2013", log[10], "(P-Value)"))) + 
  theme_bw()+
  xlab(expression(paste(log[2],"(FoldChange)"))) +
  labs(color='Genes' ) +
  xlim(-8,8)+
  scale_colour_manual(values = cols, drop=FALSE, name = "Genes", labels = c("P-Value < 0.01 ; LFC < - 1", "P-Value < 0.01 ; LFC > + 1", "0.01 \u2264 P-Value < 0.05 ; LFC < - 1", "0.01 \u2264 P-Value < 0.05 ; LFC > + 1",  "Not Significant or Substantial"))+
  guides(color = guide_legend(override.aes = list(size = 3)))

ylim<-ggplot_build(p8)$layout$panel_scales_y[[1]]$range$range[2]

p8<-p8  + geom_text(x=6.25, y=ylim-2.45*(ylim/8.25), label= upcat, color='#ff66cc') 
p8<-p8 + geom_text(x=-6.25, y=ylim-2.45*(ylim/8.25), label= downcat, color='#00ace6') 
p8<-p8 + geom_text(x=6.25, y=ylim-.7*(ylim/8.25), label= upperupcat, color='#ff0000') 
p8<-p8 + geom_text(x=-6.25, y=ylim-.7*(ylim/8.25), label= upperdowncat, color='#0000ff') 

outfile1<-"/Users/sunshine2.0/Desktop/WORK/output/volcanoMaleNvsJ_"
outfile2<-".png"
outfile<-paste(outfile1,padj_cutoff,outfile2, sep="")

ggsave(outfile, 
       dpi=300, 
       device="png",
       width=23, 
       height=14,
       units="cm")

# # Plot continuous volcano plot
# ggplot(df_log, aes(x=log2FoldChange, y=logp, color=padj)) + 
#   geom_point(na.rm = TRUE, alpha=0.3) + 
#   geom_vline(xintercept=lfc, linetype="dotted") + 
#   theme(text = element_text(size=allfont),
#         axis.text = element_text(size=axisfont)) + 
#   ggtitle("Male: November vs. June") +
#   ylab("\u2013log10(P-Value)") + 
#   xlab("log2(FoldChange)")  +
#   labs(color='Adj P-Value' ) +
#   xlim(-8,8)
# 
# outfile1 <- "/home/krwright/Downloads/volcano_continuousMaleNvsJ.png"  
# 
# ggsave(outfile1, 
#        dpi=300, 
#        device="png", 
#        width=23, 
#        height=14, 
#        units="cm")

# APRIL vs. NOVEMBER 

# Subset the relevant columns from the full dataset and choose the control condition so log2 fold changes are calculated as treatment over control
ddsMF <- ddsAll[ , (ddsAll$month == "N")|(ddsAll$month == "A") ]
dds <- ddsMF[ , (ddsMF$gender == "M") ]
dds$month <- droplevels( dds$month )
dds$month <- relevel(dds$month, ref = "N")
# check if subsetted correctly
#as.data.frame( colData(dds) )

dds <- DESeq(dds)
res <- results(dds, alpha=padj_cutoff)

# Number of adjusted p-values less than 0.1
sum(res$padj < padj_cutoff, na.rm=TRUE)
# summary (default lfcThreshold = 0)
summary(res)
sum(res$baseMean ==0)

# Convert to dataframe
df <- as.data.frame(res)

# Make summary table similar to summary(res)
total<- nrow(df)-sum(df$baseMean ==0, na.rm=TRUE)
lowCountsSubset<-df[ ((!is.na(df$pvalue))&(is.na(df$padj))),]
dfsummary9 <-  data.frame("LFCup" = c(round(sum(df$log2FoldChange>0& df$padj<padj_cutoff, na.rm=TRUE), digits=0),round(sum(df$log2FoldChange>0& df$padj<padj_cutoff, na.rm=TRUE)/total*100, digits = 1)), 
                          "LFCdown" = c(round(sum(df$log2FoldChange<0& df$padj<padj_cutoff, na.rm=TRUE), digits = 0),round(sum(df$log2FoldChange<0& df$padj<padj_cutoff, na.rm=TRUE)/total*100, digits = 1)), 
                          "outliers"= c(round(sum(!is.na(df$log2FoldChange)&is.na(df$pvalue)), digits = 0),round(sum(!is.na(df$log2FoldChange)&is.na(df$pvalue))/round(sum(is.na(df$padj)), digits = 0)*100, digits = 1)),
                          "lowCounts"= c(round(sum(!is.na(df$pvalue)&is.na(df$padj)), digits = 0),round(sum(!is.na(df$pvalue)&is.na(df$padj))/round(sum(is.na(df$padj)), digits = 0)*100, digits = 1)),
                          "baseMean"= c(round(sum(df$baseMean ==0, na.rm=TRUE), digits = 0),round(sum(df$baseMean ==0, na.rm=TRUE)/round(sum(is.na(df$padj)), digits = 0)*100, digits = 1)),
                          "total"= c(round(total, digits = 0), round(sum((sum(df$log2FoldChange>0& df$padj<padj_cutoff, na.rm=TRUE)/total*100),(sum(df$log2FoldChange<0& df$padj<padj_cutoff, na.rm=TRUE)/total*100),(sum(!is.na(df$log2FoldChange)&is.na(df$pvalue))/round(sum(is.na(df$padj)), digits = 0)*100),(sum(!is.na(df$pvalue)&is.na(df$padj))/round(sum(is.na(df$padj)), digits = 0)*100),(sum(df$baseMean ==0, na.rm=TRUE)/round(sum(is.na(df$padj)), digits = 0)*100), digits = 1))),
                          "treatment" = "M_NvsA",
                          "maxlowCounts"=c(round(max(lowCountsSubset$baseMean), digits=0),NA))

# Make a data table of df statistics
dfstats9<-data.frame("#NA"=c(round(sum(is.na(df$padj)), digits = 0),NA),
                     "#sig"=c(sum(df$padj < padj_cutoff, na.rm=TRUE),round(sum(df$padj < padj_cutoff, na.rm=TRUE)/(nrow(df)-sum(is.na(df$padj)))*100, digits = 1)),
                     "#nonsig"= c(sum(df$padj > padj_cutoff | df$padj == padj_cutoff, na.rm=TRUE),round(sum(df$padj > padj_cutoff | df$padj == padj_cutoff, na.rm=TRUE)/(nrow(df)-sum(is.na(df$padj)))*100, digits = 1)),
                     "total"= c(nrow(df), NA),
                     "treatment" = "M_NvsA")

# Add -log(pvalue) column, and padj category 
df_log <- mutate(df, logp=-log10(pvalue), padjcat=padj,
                 lfccat=log2FoldChange)

df_log <- transform(df_log, padjcat=case_when(padj < padj_cutoff & padj >= padj_cutoff1 ~ paste(padj_cutoff1,"< padj <", padj_cutoff),
                                              padj >= padj_cutoff ~ paste("padj >=", padj_cutoff),
                                              padj < padj_cutoff1 ~ paste("padj <", padj_cutoff1)))

df_log<- transform(df_log, lfccat= case_when(lfccat > lfc_cutoff ~ paste("lfc > +", lfc_cutoff), 
                                             lfccat < -lfc_cutoff ~ paste("lfc < -", lfc_cutoff),
                                             lfccat <= lfc_cutoff & lfccat >= (-lfc_cutoff) ~ paste("-",lfc_cutoff," =< lfc =< +", lfc_cutoff)))

df_log <- mutate(df_log, cat=paste(padjcat,";",lfccat))
df_log$cat[is.na(df_log$padjcat)] <- NA
df_log$cat[df_log$lfccat == paste("-",lfc_cutoff," =< lfc =< +", lfc_cutoff)] <- NA
df_log$cat[df_log$padjcat == paste("padj >=", padj_cutoff)] <- NA

upcat<-df_log %>%
  filter(cat == paste(padj_cutoff1,"< padj <", padj_cutoff,";","lfc > +", lfc_cutoff) ) %>%
  nrow()

upperupcat<-df_log %>%
  filter(cat == paste("padj <", padj_cutoff1,";","lfc > +", lfc_cutoff) ) %>%
  nrow()

downcat<- df_log %>%
  filter(cat == paste(padj_cutoff1,"< padj <", padj_cutoff,";","lfc < -", lfc_cutoff) ) %>%
  nrow()

upperdowncat<- df_log %>%
  filter(cat == paste("padj <", padj_cutoff1,";","lfc < -", lfc_cutoff) ) %>%
  nrow()

# Subset significant genes into list of gene names only
df_log <- mutate(df_log, 'gene'=rownames(df))
DEsubset9 <- df_log[  !is.na(df_log$cat) , ]
DEgenelist9 <- DEsubset9$gene #all significant genes meeting padj cutoff <0.05

DEsubset.01 <- df_log %>%
  filter(cat== paste("padj <", padj_cutoff1,";","lfc < -", lfc_cutoff) | cat== paste("padj <", padj_cutoff1,";","lfc > +", lfc_cutoff) )
DEgenelist9.01 <- DEsubset.01$gene #all significant genes meeting padj cutoff <0.01

DEsubset.up01 <- df_log %>%
  filter(cat== paste("padj <", padj_cutoff1,";","lfc > +", lfc_cutoff) )
DEgenelist9.up01 <- DEsubset.up01$gene #all up-regulated genes meeting padj cutoff <0.01

DEgenelist9.01p<-df_log%>% # Obtain list of DE0.01 genes with padj
  dplyr::select(c(gene, padj)) %>%
  filter(gene %in% DEgenelist9.01) %>%
  mutate("month.comparison"="April vs. Nov")

DEsubset.down01 <- df_log %>%
  filter(cat== paste("padj <", padj_cutoff1,";","lfc < -", lfc_cutoff) )
DEgenelist9.down01 <- DEsubset.down01$gene #all down-regulated genes meeting padj cutoff <0.01

DEsubset.up05 <-  df_log[  !is.na(df_log$cat) , ]
DEsubset.up05 <- DEsubset.up05 %>%
  filter(lfccat== paste("lfc > +", lfc_cutoff) )
DEgenelist9.up05 <- DEsubset.up05$gene #all up-regulated genes meeting padj cutoff <0.05

DEsubset.down05 <-  df_log[  !is.na(df_log$cat) , ]
DEsubset.down05 <- DEsubset.down05 %>%
  filter(lfccat== paste("lfc < -", lfc_cutoff) )
DEgenelist9.down05 <- DEsubset.down05$gene #all down-regulated genes meeting padj cutoff <0.05

DE_padj9<-df_log%>% # Obtain list of DE0.01 genes with padj
  dplyr::select(c(gene, padj)) %>%
  mutate("month.comparison"="April vs. Nov")

df_log_9<-df_log

# Volcano plot font sizes
axisfont <- 14
labelfont <- 5
allfont <- 18

# Reassign rownames
rownames(df_log) <- rownames(df)

# Convert padjcat to factor
df_log$padjcat <- as.factor(df_log$padjcat)

# Replace NA with Not Significant or Substantial
df_log<-df_log %>%
  mutate(cat = if_else(is.na(cat), "Not Significant or Substantial", cat))
df_log$cat<-factor(df_log$cat,
                   levels=c("padj < 0.01 ; lfc < - 1" ,"padj < 0.01 ; lfc > + 1", "0.01 < padj < 0.05 ; lfc < - 1" , "0.01 < padj < 0.05 ; lfc > + 1", "Not Significant or Substantial"))
cols <- c("0.01 < padj < 0.05 ; lfc < - 1" = "#00ace6", "0.01 < padj < 0.05 ; lfc > + 1" = "#ff66cc", "padj < 0.01 ; lfc < - 1" = "#0000ff", "padj < 0.01 ; lfc > + 1" = "#ff0000", "Not Significant or Substantial"="#818181")

# Plot discrete volcano plot
p9<-ggplot(df_log, aes(x=log2FoldChange, y=logp, color=cat)) + 
  geom_point(na.rm = TRUE,alpha=0.5) +
  geom_vline(xintercept=c(-1,1),linetype="dotted") + 
  geom_hline(yintercept=c(2.72, 3.9) ,linetype="dotted") + 
  theme(text = element_text(size=allfont),
        axis.text = element_text(size=axisfont)) + 
  ggtitle("Male: April vs. Nov") +
  ylab(expression(paste("\u2013", log[10], "(P-Value)"))) + 
  theme_bw()+
  xlab(expression(paste(log[2],"(FoldChange)"))) +
  labs(color='Genes' ) +
  xlim(-8,8)+
  scale_colour_manual(values = cols, drop=FALSE, name = "Genes", labels = c("P-Value < 0.01 ; LFC < - 1", "P-Value < 0.01 ; LFC > + 1", "0.01 \u2264 P-Value < 0.05 ; LFC < - 1", "0.01 \u2264 P-Value < 0.05 ; LFC > + 1",  "Not Significant or Substantial"))+
  guides(color = guide_legend(override.aes = list(size = 3)))

ylim<-ggplot_build(p9)$layout$panel_scales_y[[1]]$range$range[2]

p9<-p9  + geom_text(x=6.25, y=ylim-2.45*(ylim/8.25), label= upcat, color='#ff66cc') 
p9<-p9 + geom_text(x=-6.25, y=ylim-2.45*(ylim/8.25), label= downcat, color='#00ace6') 
p9<-p9 + geom_text(x=6.25, y=ylim-.7*(ylim/8.25), label= upperupcat, color='#ff0000') 
p9<-p9 + geom_text(x=-6.25, y=ylim-.7*(ylim/8.25), label= upperdowncat, color='#0000ff') 

outfile1<-"/Users/sunshine2.0/Desktop/WORK/output/volcanoMaleAvaN_"
outfile2<-".png"
outfile<-paste(outfile1,padj_cutoff,outfile2, sep="")

ggsave(outfile, 
       dpi=300, 
       device="png",
       width=23, 
       height=14,
       units="cm")

# # Plot continuous volcano plot
# ggplot(df_log, aes(x=log2FoldChange, y=logp, color=padj)) + 
#   geom_point(na.rm = TRUE, alpha=0.3) + 
#   geom_vline(xintercept=lfc, linetype="dotted") + 
#   theme(text = element_text(size=allfont),
#         axis.text = element_text(size=axisfont)) + 
#   ggtitle("Male: November vs. April") +
#   ylab("\u2013log10(P-Value)") + 
#   xlab("log2(FoldChange)")  +
#   labs(color='Adj P-Value' ) +
#   xlim(-10,10)
# 
# outfile1 <- "/home/krwright/Downloads/volcano_continuousMaleNvsA.png"  
# 
# ggsave(outfile1, 
#        dpi=300, 
#        device="png", 
#        width=23, 
#        height=14, 
#        units="cm")

## EXPORT SUMMARY, STATS, AND GENE LISTS SHEETS
# Name outfiles 
outfile1.1<-"/Users/sunshine2.0/Desktop/WORK/output/dfsummary_emu_"
outfile2.1<-".xlsx"
outfile.1<-paste(outfile1.1,padj_cutoff,outfile2.1, sep="")

outfile1.2<-"/Users/sunshine2.0/Desktop/WORK/output/dfstats_emu_"
outfile2.2<-".xlsx"
outfile.2<-paste(outfile1.2,padj_cutoff,outfile2.2, sep="")

# gene list outfiles padj < 0.05
outfile1.3<-"/Users/sunshine2.0/Desktop/WORK/output/genelist_JvaA_"
outfile2.3<-".csv"
outfile.3<-paste(outfile1.3,padj_cutoff,outfile2.3, sep="")

outfile1.4<-"/Users/sunshine2.0/Desktop/WORK/output/genelist_NvaJ_"
outfile2.4<-".csv"
outfile.4<-paste(outfile1.4,padj_cutoff,outfile2.4, sep="")

outfile1.5<-"/Users/sunshine2.0/Desktop/WORK/output/genelist_NvaA_"
outfile2.5<-".csv"
outfile.5<-paste(outfile1.5,padj_cutoff,outfile2.5, sep="")

outfile1.6<-"/Users/sunshine2.0/Desktop/WORK/output/genelist_FemaleJvaA_"
outfile2.6<-".csv"
outfile.6<-paste(outfile1.6,padj_cutoff,outfile2.6, sep="")

outfile1.7<-"/Users/sunshine2.0/Desktop/WORK/output/genelist_FemaleNvaJ_"
outfile2.7<-".csv"
outfile.7<-paste(outfile1.7,padj_cutoff,outfile2.7, sep="")

outfile1.8<-"/Users/sunshine2.0/Desktop/WORK/output/genelist_FemaleNvaA_"
outfile2.8<-".csv"
outfile.8<-paste(outfile1.8,padj_cutoff,outfile2.8, sep="")

outfile1.9<-"/Users/sunshine2.0/Desktop/WORK/output/genelist_MaleJvaA_"
outfile2.9<-".csv"
outfile.9<-paste(outfile1.9,padj_cutoff,outfile2.9, sep="")

outfile1.10<-"/Users/sunshine2.0/Desktop/WORK/output/genelist_MaleNvaJ_"
outfile2.10<-".csv"
outfile.10<-paste(outfile1.10,padj_cutoff,outfile2.10, sep="")

outfile1.11<-"/Users/sunshine2.0/Desktop/WORK/output/genelist_MaleNvaA_"
outfile2.11<-".csv"
outfile.11<-paste(outfile1.11,padj_cutoff,outfile2.11, sep="")

# gene list outfiles padj < 0.01
outfile.33<-paste(outfile1.3,padj_cutoff1,outfile2.3, sep="")

outfile.44<-paste(outfile1.4,padj_cutoff1,outfile2.4, sep="")

outfile.55<-paste(outfile1.5,padj_cutoff1,outfile2.5, sep="")

outfile.66<-paste(outfile1.6,padj_cutoff1,outfile2.6, sep="")

outfile.77<-paste(outfile1.7,padj_cutoff1,outfile2.7, sep="")

outfile.88<-paste(outfile1.8,padj_cutoff1,outfile2.8, sep="")

outfile.99<-paste(outfile1.9,padj_cutoff1,outfile2.9, sep="")

outfile.1010<-paste(outfile1.10,padj_cutoff1,outfile2.10, sep="")

outfile.1111<-paste(outfile1.11,padj_cutoff1,outfile2.11, sep="")

# Combine and export excel sheet for dfsummaries padj < 0.05
library(xlsx)
dfsummary <- rbind(dfsummary1,dfsummary2,dfsummary3,dfsummary4,dfsummary5,dfsummary6,dfsummary7,dfsummary8,dfsummary9)
dfsummary <- cbind(dfsummary,"numberType"=c("#","%","#","%","#","%","#","%","#","%","#","%","#","%","#","%","#","%" ))
write.xlsx(dfsummary, outfile.1)

# Combine and export excel sheet for dfstats padj < 0.05
dfstats <- rbind(dfstats1,dfstats2,dfstats3,dfstats4,dfstats5,dfstats6,dfstats7,dfstats8,dfstats9)
dfstats <- cbind(dfstats,"numberType"=c("#","%","#","%","#","%","#","%","#","%","#","%","#","%","#","%","#","%" ))
write.xlsx(dfstats, outfile.2)

# Export excel sheet for DE gene lists of padj < 0.05
write.table(DEgenelist1, outfile.3, sep=",", row.names = FALSE, col.names = FALSE)
write.table(DEgenelist2, outfile.4, sep=",", row.names = FALSE, col.names = FALSE)
write.table(DEgenelist3, outfile.5, sep=",", row.names = FALSE, col.names = FALSE)
write.table(DEgenelist4, outfile.6, sep=",", row.names = FALSE, col.names = FALSE)
write.table(DEgenelist5, outfile.7, sep=",", row.names = FALSE, col.names = FALSE)
write.table(DEgenelist6, outfile.8, sep=",", row.names = FALSE, col.names = FALSE)
write.table(DEgenelist7, outfile.9, sep=",", row.names = FALSE, col.names = FALSE)
write.table(DEgenelist8, outfile.10, sep=",", row.names = FALSE, col.names = FALSE)
write.table(DEgenelist9, outfile.11, sep=",", row.names = FALSE, col.names = FALSE)

# Export excel sheet for DE gene lists of padj < 0.01
write.table(DEgenelist1.01, outfile.33, sep=",", row.names = FALSE, col.names = FALSE)
write.table(DEgenelist2.01, outfile.44, sep=",", row.names = FALSE, col.names = FALSE)
write.table(DEgenelist3.01, outfile.55, sep=",", row.names = FALSE, col.names = FALSE)
write.table(DEgenelist4.01, outfile.66, sep=",", row.names = FALSE, col.names = FALSE)
write.table(DEgenelist5.01, outfile.77, sep=",", row.names = FALSE, col.names = FALSE)
write.table(DEgenelist6.01, outfile.88, sep=",", row.names = FALSE, col.names = FALSE)
write.table(DEgenelist7.01, outfile.99, sep=",", row.names = FALSE, col.names = FALSE)
write.table(DEgenelist8.01, outfile.1010, sep=",", row.names = FALSE, col.names = FALSE)
write.table(DEgenelist9.01, outfile.1111, sep=",", row.names = FALSE, col.names = FALSE)

# Plot number of genes in each treatment
dfstatsSub <- dfstats[dfstats$numberType == "#",c(1:3,5)]
dfstatsSub1 <- separate(dfstatsSub,treatment, c("gender", "months"), sep="_")
dfstatsNA <- dfstatsSub1[,c(1,4,5)]
dfstatsNA$isSig <- "NA"
colnames(dfstatsNA)[1]<- "geneCount"
dfstatssig <- dfstatsSub1[,c(2,4,5)]
dfstatssig$isSig <- "significant"
colnames(dfstatssig)[1]<- "geneCount"
dfstatsnonsig <-dfstatsSub1[,c(3,4,5)]
dfstatsnonsig$isSig <- "non-significant"
colnames(dfstatsnonsig)[1]<- "geneCount"
dfstatsPlot <- rbind(dfstatssig,dfstatsnonsig,dfstatsNA)

dfstatsPlot$gender_alpha <-paste(dfstatsPlot$gender,"_", padj_cutoff, sep="")
dfstatsPlot$gender_alpha <-factor(dfstatsPlot$gender_alpha, levels= unique(dfstatsPlot$gender_alpha))
dfstatsPlot$isSig <-factor(dfstatsPlot$isSig, levels= unique(dfstatsPlot$isSig))
dfstatsPlot$months <-factor(dfstatsPlot$months, levels= unique(dfstatsPlot$months))
dfstatsPlot$gender <-factor(dfstatsPlot$gender, levels= unique(dfstatsPlot$gender))

ggplot(dfstatsPlot, aes(x=months,y=geneCount, fill= isSig))+
  geom_col(aes(group = isSig), position = "dodge")+
  facet_grid(.~gender_alpha)+
  scale_y_continuous(trans='log10')+
  scale_fill_manual(values=c("lightcoral","cyan3","azure4"))+
  theme(text=element_text(size=20))

outfile1<-"/Users/sunshine2.0/Desktop/WORK/output/geneCountbar_"
outfile2<-".png"
outfile<-paste(outfile1,padj_cutoff,outfile2, sep="")

ggsave(outfile, 
       dpi=300, 
       device="png",
       width=29, 
       height=14,
       units="cm")

# Plot % of DE genes in each treatment
dfstatsSub.1 <- dfstats[dfstats$numberType == "%",c(2,5)]
dfstatsPlot.1 <- separate(dfstatsSub.1,treatment, c("gender", "months"), sep="_")

dfstatsPlot.1$gender_alpha <-paste(dfstatsPlot.1$gender,"_", padj_cutoff, sep="")
dfstatsPlot.1$gender_alpha <-factor(dfstatsPlot.1$gender_alpha, levels= unique(dfstatsPlot.1$gender_alpha))
dfstatsPlot.1$months <-factor(dfstatsPlot.1$months, levels= unique(dfstatsPlot.1$months))
dfstatsPlot.1$gender <-factor(dfstatsPlot.1$gender, levels= unique(dfstatsPlot.1$gender))

ggplot(dfstatsPlot.1, aes(x=months,y=X.sig))+
  geom_col( fill= "lightcoral")+
  facet_grid(.~gender_alpha)+
  ylab("significant DE genes%")+
  theme(text=element_text(size=20))+
  geom_hline(yintercept=padj_cutoff*100, linetype="dashed")

outfile1.99<-"/Users/sunshine2.0/Desktop/WORK/output/DEgenePercent_"
outfile2.99<-".png"
outfile.99<-paste(outfile1.99,padj_cutoff,outfile2.99, sep="")

ggsave(outfile.99, 
       dpi=300, 
       device="png",
       width=29, 
       height=14,
       units="cm")


## ARRANGE ALL VOLCANO PLOTS INTO ONE FIGURE
library("ggpubr")

figure <-ggarrange(p1 + rremove("xylab")+ labs(title = NULL), p2 + rremove("xylab")+ labs(title = NULL), p3 +  rremove("xylab")+ labs(title = NULL),p4 + rremove("xylab")+ labs(title = NULL), p5 + rremove("xylab")+ labs(title = NULL), p6 + rremove("xylab")+ labs(title = NULL), p7 + rremove("xylab")+ labs(title = NULL), p8 + rremove("xylab")+ labs(title = NULL), p9 + rremove("xylab")+ labs(title = NULL),
          common.legend = TRUE,
          legend = "right",
          widths = rep(5), 
          hjust = -0.5,
          vjust = 1,
          ncol = 3, nrow = 3,
          align = "hv") 
  
  figure1 <-annotate_figure(figure,
                            bottom = text_grob(expression(paste(log[2],"(Fold Change)   ")), hjust=1.15),
                            left = text_grob(expression(paste("\u2013", log[10], "(P-Value)")),hjust=0.35, rot = 90))
# figure1 <-annotate_figure(figure,
#                 bottom = text_grob("log2(Fold Change)     ", hjust=1.15),
#                 left = text_grob("\u2013log10(P-Value)",hjust=0.35, rot = 90))

annotate_figure(figure1,
                          top = text_grob("June vs. April                                 November vs. June                              April vs. November                                  ", face="bold"),
                          left = text_grob("           Male                         Female                    Both Sexes", rot = 90, face="bold"))

outfile1.999<-"/Users/sunshine2.0/Desktop/WORK/output/grid.arrange.final.legend_"
outfile2.999<-".png"
outfile.999<-paste(outfile1.999,padj_cutoff,outfile2.999, sep="")

ggsave(outfile.999, 
       dpi=300, 
       device="png",
       width=29, 
       height=14,
       units="cm")
                
# _______________________________________________________END___________________________________________________________


