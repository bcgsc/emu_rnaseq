library(tximport)
library(tidyverse)
library(rhdf5)
library(DESeq2)

## TO RUN DIFFERENTIAL EXPRESSION ANALYSIS (SEXUAL COMPARISONS)

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
ddsAll <- DESeqDataSetFromTximport(salmon_quant, samples, ~gender)

# Adjusted p-value threshold default, can be modified (default alpha=0.1)
padj_cutoff<-0.05
padj_cutoff1<-0.01

# Log fold change threshold default, can be modified
lfc <-0
lfc_cutoff <-1

## PLOT VOLCANO PLOTS For each condition vs control





# # JUNE VS APRIL
# 
# # Subset the relevant columns from the full dataset and choose the control condition so log2 fold changes are calculated as treatment over control
# dds <- ddsAll[ , (ddsAll$month == "A") ]
# dds$month <- droplevels( dds$month )
# dds$month <- relevel(dds$month, ref = "F")
# # check if subsetted correctly
# #as.data.frame( colData(dds) )
# 
# dds <- DESeq(dds)
# res <- results(dds, alpha=padj_cutoff)
# 
# 
# # number of adjusted p-values less than padj_cutoff
# sum(res$padj < padj_cutoff, na.rm=TRUE)
# # summary (default lfcThreshold = 0)
# summary(res)
# sum(res$baseMean ==0)
# 
# # for each condition vs control, plot a volcano plot
# # convert to dataframe
# 
# df <- as.data.frame(res)
# 
# 
# # make summary table similar to summary(res)
# total<- nrow(df)-sum(df$baseMean ==0, na.rm=TRUE)
# lowCountsSubset<-df[ ((!is.na(df$pvalue))&(is.na(df$padj))),]
# dfsummary1 <- data.frame("LFCup" = c(round(sum(df$log2FoldChange>0& df$padj<padj_cutoff, na.rm=TRUE), digits=0),round(sum(df$log2FoldChange>0& df$padj<padj_cutoff, na.rm=TRUE)/total*100, digits = 1)), 
#                       "LFCdown" = c(round(sum(df$log2FoldChange<0& df$padj<padj_cutoff, na.rm=TRUE), digits = 0),round(sum(df$log2FoldChange<0& df$padj<padj_cutoff, na.rm=TRUE)/total*100, digits = 1)), 
#                       "outliers"= c(round(sum(!is.na(df$log2FoldChange)&is.na(df$pvalue)), digits = 0),round(sum(!is.na(df$log2FoldChange)&is.na(df$pvalue))/round(sum(is.na(df$padj)), digits = 0)*100, digits = 1)),
#                       "lowCounts"= c(round(sum(!is.na(df$pvalue)&is.na(df$padj)), digits = 0),round(sum(!is.na(df$pvalue)&is.na(df$padj))/round(sum(is.na(df$padj)), digits = 0)*100, digits = 1)),
#                       "baseMean"= c(round(sum(df$baseMean ==0, na.rm=TRUE), digits = 0),round(sum(df$baseMean ==0, na.rm=TRUE)/round(sum(is.na(df$padj)), digits = 0)*100, digits = 1)),
#                       "total"= c(round(total, digits = 0), round(sum((sum(df$log2FoldChange>0& df$padj<padj_cutoff, na.rm=TRUE)/total*100),(sum(df$log2FoldChange<0& df$padj<padj_cutoff, na.rm=TRUE)/total*100),(sum(!is.na(df$log2FoldChange)&is.na(df$pvalue))/round(sum(is.na(df$padj)), digits = 0)*100),(sum(!is.na(df$pvalue)&is.na(df$padj))/round(sum(is.na(df$padj)), digits = 0)*100),(sum(df$baseMean ==0, na.rm=TRUE)/round(sum(is.na(df$padj)), digits = 0)*100), digits = 1))),
#                       "treatment" = "MF_JvsA",
#                       "maxlowCounts"=c(round(max(lowCountsSubset$baseMean), digits=0),NA))
#                     
# # make a data table of df statistics
# dfstats1<-data.frame("#NA"=c(round(sum(is.na(df$padj)), digits = 0),NA),
#                      "#sig"=c(sum(df$padj < padj_cutoff, na.rm=TRUE),round(sum(df$padj < padj_cutoff, na.rm=TRUE)/(nrow(df)-sum(is.na(df$padj)))*100, digits = 1)),
#                      "#nonsig"= c(sum(df$padj > padj_cutoff | df$padj == padj_cutoff, na.rm=TRUE),round(sum(df$padj > padj_cutoff | df$padj == padj_cutoff, na.rm=TRUE)/(nrow(df)-sum(is.na(df$padj)))*100, digits = 1)),
#                      "total"= c(nrow(df), NA),
#                      "treatment" = "MF_JvsA")
# 
# 
# # add -log(pvalue) column, and padj category 
# df_log <- mutate(df, logp=-log10(pvalue), padjcat=padj,
#                  lfccat=log2FoldChange)
# 
# df_log <- transform(df_log, padjcat=case_when(padj < padj_cutoff & padj > padj_cutoff1 ~ paste(padj_cutoff1,"< padj <", padj_cutoff),
#                                               padj >= padj_cutoff ~ paste("padj >=", padj_cutoff),
#                                               padj < padj_cutoff1 ~ paste("padj <", padj_cutoff1)))
# 
# df_log<- transform(df_log, lfccat= case_when(lfccat > lfc_cutoff ~ paste("lfc > +", lfc_cutoff), 
#                                              lfccat < -lfc_cutoff ~ paste("lfc < -", lfc_cutoff),
#                                              lfccat <= lfc_cutoff & lfccat >= (-lfc_cutoff) ~ paste("-",lfc_cutoff," =< lfc =< +", lfc_cutoff)))
# 
# df_log <- mutate(df_log, cat=paste(padjcat,";",lfccat))
# df_log$cat[is.na(df_log$padjcat)] <- NA
# df_log$cat[df_log$lfccat == paste("-",lfc_cutoff," =< lfc =< +", lfc_cutoff)] <- NA
# df_log$cat[df_log$padjcat == paste("padj >=", padj_cutoff)] <- NA
# 
# upcat<-df_log %>%
#   filter(cat == paste(padj_cutoff1,"< padj <", padj_cutoff,";","lfc > +", lfc_cutoff) ) %>%
#   nrow()
# 
# upperupcat<-df_log %>%
#   filter(cat == paste("padj <", padj_cutoff1,";","lfc > +", lfc_cutoff) ) %>%
#   nrow()
# 
# 
# downcat<- df_log %>%
#   filter(cat == paste(padj_cutoff1,"< padj <", padj_cutoff,";","lfc < -", lfc_cutoff) ) %>%
#   nrow()
# 
# upperdowncat<- df_log %>%
#   filter(cat == paste("padj <", padj_cutoff1,";","lfc < -", lfc_cutoff) ) %>%
#   nrow()
# 
# # df_log %>%
# #   filter(padj > 0.04 & padj< 0.06)
# 
# #subset significant genes into list of gene names only
# df_log <- mutate(df_log, 'gene'=rownames(df))
# DEsubset1 <- df_log[  !is.na(df_log$cat) , ]
# DEgenelist1 <- DEsubset1$gene
# 
# 
# DEsubset.01 <- df_log %>%
#   filter(cat== paste("padj <", padj_cutoff1,";","lfc < -", lfc_cutoff) | cat== paste("padj <", padj_cutoff1,";","lfc > +", lfc_cutoff) )
# DEgenelist1.01 <- DEsubset.01$gene
# 
# 
# DEsubset.up01 <- df_log %>%
#   filter(cat== paste("padj <", padj_cutoff1,";","lfc > +", lfc_cutoff) )
# DEgenelist1.up01 <- DEsubset.up01$gene
# length(DEgenelist1.up01)
# 
# DEsubset.down05 <- df_log %>%
#   filter(cat== paste("padj <", padj_cutoff,";","lfc < -", lfc_cutoff) )
# DEgenelist1.down01 <- DEsubset.down01$gene
# length(DEgenelist1.down01)
# 
# DEsubset.up05 <-  df_log[  !is.na(df_log$cat) , ]
# DEsubset.up05 <- DEsubset.up05 %>%
#   filter(lfccat== paste("lfc > +", lfc_cutoff) )
# DEgenelist1.up05 <- DEsubset.up05$gene
# length(DEgenelist1.up05)
# 
# DEsubset.down05 <-  df_log[  !is.na(df_log$cat) , ]
# DEsubset.down05 <- DEsubset.down05 %>%
#   filter(lfccat== paste("lfc < -", lfc_cutoff) )
# DEgenelist1.down05 <- DEsubset.down05$gene
# length(DEgenelist1.down05)
# 
# 
# # Default font sizes
# axisfont <- 14
# labelfont <- 5
# allfont <- 18
# 
# # Default for this type of plot
# legendtitle<-"adjusted\np-value"
# 
# # reassign rownames
# rownames(df_log) <- rownames(df)
# 
# # convert padjcat to factor
# df_log$padjcat <- as.factor(df_log$padjcat)
# 
# df_log %>%
#   filter(padj<0.02) %>%
#   filter(padj>0.005)
# 
# 
# # Plot discrete
# p1<-ggplot(df_log, aes(x=log2FoldChange, y=logp, color=cat)) + 
#   geom_point(na.rm = TRUE,alpha=0.5) +
#   geom_vline(xintercept=c(-1,1),linetype="dotted") + 
#   geom_hline(yintercept=c(-log10(padj_cutoff/10),-log10(padj_cutoff1/36)) ,linetype="dotted") + 
#   theme(text = element_text(size=allfont),
#         axis.text = element_text(size=axisfont)) + 
#   ggtitle("June vs. April") +
#   ylab("\u2013log10(P-Value)") + 
#   theme_bw()+
#   xlab("log2(FoldChange)") +
#   labs(color='Adj P-Value' ) +
#   xlim(-8,8)
#  # xlimit() +
#  # ylimit()
# ylim<-ggplot_build(p1)$layout$panel_scales_y[[1]]$range$range[2]
# 
# p1<-p1  + geom_text(x=7.2, y=ylim-3*(ylim/5), label= upcat, color='olivedrab4') 
# p1<-p1 + geom_text(x=-7.2, y=ylim-3*(ylim/5), label= downcat, color='red') 
# p1<-p1 + geom_text(x=7.2, y=ylim-ylim/4, label= upperupcat, color='purple') 
# p1<-p1 + geom_text(x=-7.2, y=ylim-ylim/4, label= upperdowncat, color='darkcyan') 
# 
# outfile1<-"/Users/sunshine2.0/Desktop/WORK/output/volcanoJvsA_"
# outfile2<-".png"
# outfile<-paste(outfile1,padj_cutoff,outfile2, sep="")
# 
# ggsave(outfile, 
#        dpi=300, 
#        device="png",
#        width=23, 
#        height=14,
#        units="cm")
# 
# 
# # Plot continuous
# # ggplot(df_log, aes(x=log2FoldChange, y=logp, color=padj)) +
# #   # geom_point(na.rm = TRUE, alpha=0.3) +
# #   geom_vline(xintercept=lfc, linetype="dotted") +
# #   theme(text = element_text(size=allfont),
# #         axis.text = element_text(size=axisfont)) +
# #   ggtitle("June vs. April") +
# #   ylab("\u2013log10(P-Value)") +
# #   xlab("log2(FoldChange)")  +
# #   labs(color='Adj P-Value' )
# # 
# # outfile1 <- "/home/krwright/Downloads/volcano_continuousJvsA.png"
# # 
# # ggsave(outfile1,
# #        dpi=300,
# #        device="png",
# #        width=23,
# #        height=14,
# #        units="cm")
# 
# # NOVEMBER VS JUNE 
# 
# # Subset the relevant columns from the full dataset and choose the control condition so log2 fold changes are calculated as treatment over control
# dds <- ddsAll[ , (ddsAll$month == "N")|(ddsAll$month == "J") ]
# dds$month <- droplevels( dds$month )
# dds$month <- relevel(dds$month, ref = "J")
# # check if subsetted correctly
# #as.data.frame( colData(dds) )
# 
# dds <- DESeq(dds)
# res <- results(dds, alpha=padj_cutoff)
# 
# # number of adjusted p-values less than 0.1
# sum(res$padj < padj_cutoff, na.rm=TRUE)
# # summary (default lfcThreshold = 0)
# summary(res)
# sum(res$baseMean ==0)
# 
# # for each condition vs control, plot a volcano plot
# # convert to dataframe
# 
# df <- as.data.frame(res)
# 
# # make summary table similar to summary(res)
# total<- nrow(df)-sum(df$baseMean ==0, na.rm=TRUE)
# lowCountsSubset<-df[ ((!is.na(df$pvalue))&(is.na(df$padj))),]
# dfsummary2 <-  data.frame("LFCup" = c(round(sum(df$log2FoldChange>0& df$padj<padj_cutoff, na.rm=TRUE), digits=0),round(sum(df$log2FoldChange>0& df$padj<padj_cutoff, na.rm=TRUE)/total*100, digits = 1)), 
#                           "LFCdown" = c(round(sum(df$log2FoldChange<0& df$padj<padj_cutoff, na.rm=TRUE), digits = 0),round(sum(df$log2FoldChange<0& df$padj<padj_cutoff, na.rm=TRUE)/total*100, digits = 1)), 
#                           "outliers"= c(round(sum(!is.na(df$log2FoldChange)&is.na(df$pvalue)), digits = 0),round(sum(!is.na(df$log2FoldChange)&is.na(df$pvalue))/round(sum(is.na(df$padj)), digits = 0)*100, digits = 1)),
#                           "lowCounts"= c(round(sum(!is.na(df$pvalue)&is.na(df$padj)), digits = 0),round(sum(!is.na(df$pvalue)&is.na(df$padj))/round(sum(is.na(df$padj)), digits = 0)*100, digits = 1)),
#                           "baseMean"= c(round(sum(df$baseMean ==0, na.rm=TRUE), digits = 0),round(sum(df$baseMean ==0, na.rm=TRUE)/round(sum(is.na(df$padj)), digits = 0)*100, digits = 1)),
#                           "total"= c(round(total, digits = 0), round(sum((sum(df$log2FoldChange>0& df$padj<padj_cutoff, na.rm=TRUE)/total*100),(sum(df$log2FoldChange<0& df$padj<padj_cutoff, na.rm=TRUE)/total*100),(sum(!is.na(df$log2FoldChange)&is.na(df$pvalue))/round(sum(is.na(df$padj)), digits = 0)*100),(sum(!is.na(df$pvalue)&is.na(df$padj))/round(sum(is.na(df$padj)), digits = 0)*100),(sum(df$baseMean ==0, na.rm=TRUE)/round(sum(is.na(df$padj)), digits = 0)*100), digits = 1))),
#                           "treatment" = "MF_NvsJ",
#                           "maxlowCounts"=c(round(max(lowCountsSubset$baseMean), digits=0),NA))
#    
# # make a data table of df statistics
# dfstats2<-data.frame("#NA"=c(round(sum(is.na(df$padj)), digits = 0),NA),
#                      "#sig"=c(sum(df$padj < padj_cutoff, na.rm=TRUE),round(sum(df$padj < padj_cutoff, na.rm=TRUE)/(nrow(df)-sum(is.na(df$padj)))*100, digits = 1)),
#                      "#nonsig"= c(sum(df$padj > padj_cutoff | df$padj == padj_cutoff, na.rm=TRUE),round(sum(df$padj > padj_cutoff | df$padj == padj_cutoff, na.rm=TRUE)/(nrow(df)-sum(is.na(df$padj)))*100, digits = 1)),
#                      "total"= c(nrow(df), NA),
#                      "treatment" = "MF_NvsJ")
# 
# # add -log(pvalue) column, and padj category 
# df_log <- mutate(df, logp=-log10(pvalue), padjcat=padj,
#                  lfccat=log2FoldChange)
# 
# df_log <- transform(df_log, padjcat=case_when(padj < padj_cutoff & padj > padj_cutoff1 ~ paste(padj_cutoff1,"< padj <", padj_cutoff),
#                                               padj >= padj_cutoff ~ paste("padj >=", padj_cutoff),
#                                               padj < padj_cutoff1 ~ paste("padj <", padj_cutoff1)))
# 
# df_log<- transform(df_log, lfccat= case_when(lfccat > lfc_cutoff ~ paste("lfc > +", lfc_cutoff), 
#                                              lfccat < -lfc_cutoff ~ paste("lfc < -", lfc_cutoff),
#                                              lfccat <= lfc_cutoff & lfccat >= (-lfc_cutoff) ~ paste("-",lfc_cutoff," =< lfc =< +", lfc_cutoff)))
# 
# df_log <- mutate(df_log, cat=paste(padjcat,";",lfccat))
# df_log$cat[is.na(df_log$padjcat)] <- NA
# df_log$cat[df_log$lfccat == paste("-",lfc_cutoff," =< lfc =< +", lfc_cutoff)] <- NA
# df_log$cat[df_log$padjcat == paste("padj >=", padj_cutoff)] <- NA
# 
# upcat<-df_log %>%
#   filter(cat == paste(padj_cutoff1,"< padj <", padj_cutoff,";","lfc > +", lfc_cutoff) ) %>%
#   nrow()
# 
# upperupcat<-df_log %>%
#   filter(cat == paste("padj <", padj_cutoff1,";","lfc > +", lfc_cutoff) ) %>%
#   nrow()
# 
# 
# downcat<- df_log %>%
#   filter(cat == paste(padj_cutoff1,"< padj <", padj_cutoff,";","lfc < -", lfc_cutoff) ) %>%
#   nrow()
# 
# upperdowncat<- df_log %>%
#   filter(cat == paste("padj <", padj_cutoff1,";","lfc < -", lfc_cutoff) ) %>%
#   nrow()
# 
# #subset significant genes into list of gene names only
# df_log <- mutate(df_log, 'gene'=rownames(df))
# DEsubset2 <- df_log[  !is.na(df_log$cat) , ]
# DEgenelist2 <- DEsubset2$gene
# 
# DEsubset.01 <- df_log %>%
#   filter(cat== paste("padj <", padj_cutoff1,";","lfc < -", lfc_cutoff) | cat== paste("padj <", padj_cutoff1,";","lfc > +", lfc_cutoff) )
# DEgenelist2.01 <- DEsubset.01$gene
# 
# 
# DEsubset.up01 <- df_log %>%
#   filter(cat== paste("padj <", padj_cutoff1,";","lfc > +", lfc_cutoff) )
# DEgenelist2.up01 <- DEsubset.up01$gene
# 
# 
# DEsubset.down01 <- df_log %>%
#   filter(cat== paste("padj <", padj_cutoff1,";","lfc < -", lfc_cutoff) )
# DEgenelist2.down01 <- DEsubset.down01$gene
# 
# DEsubset.up05 <-  df_log[  !is.na(df_log$cat) , ]
# DEsubset.up05 <- DEsubset.up05 %>%
#   filter(lfccat== paste("lfc > +", lfc_cutoff) )
# DEgenelist2.up05 <- DEsubset.up05$gene
# 
# DEsubset.down05 <-  df_log[  !is.na(df_log$cat) , ]
# DEsubset.down05 <- DEsubset.down05 %>%
#   filter(lfccat== paste("lfc < -", lfc_cutoff) )
# DEgenelist2.down05 <- DEsubset.down05$gene
# 
# # Default font sizes
# axisfont <- 14
# labelfont <- 5
# allfont <- 18
# 
# # Default for this type of plot
# legendtitle<-"adjusted\np-value"
# 
# # reassign rownames
# rownames(df_log) <- rownames(df)
# 
# # convert padjcat to factor
# df_log$padjcat <- as.factor(df_log$padjcat)
# 
# 
# 
# # Plot discrete
# p2<-ggplot(df_log, aes(x=log2FoldChange, y=logp, color=cat)) + 
#   geom_point(na.rm = TRUE,alpha=0.5) +
#   geom_vline(xintercept=c(-1,1),linetype="dotted") + 
#   geom_hline(yintercept=c(2.25,3.34) ,linetype="dotted") + 
#   theme(text = element_text(size=allfont),
#         axis.text = element_text(size=axisfont)) + 
#   ggtitle("Nov vs. June") +
#   ylab("\u2013log10(P-Value)") + 
#   theme_bw()+
#   xlab("log2(FoldChange)") +
#   labs(color='Adj P-Value' ) +
#   xlim(-8,8)
# # xlimit() +
# # ylimit()
# ylim<-ggplot_build(p2)$layout$panel_scales_y[[1]]$range$range[2]
# 
# p2<-p2  + geom_text(x=7.2, y=ylim-3*(ylim/5), label= upcat, color='olivedrab4') 
# p2<-p2 + geom_text(x=-7.2, y=ylim-3*(ylim/5), label= downcat, color='red') 
# p2<-p2 + geom_text(x=7.2, y=ylim-ylim/4, label= upperupcat, color='purple') 
# p2<-p2 + geom_text(x=-7.2, y=ylim-ylim/4, label= upperdowncat, color='darkcyan') 
# 
# outfile1<-"/Users/sunshine2.0/Desktop/WORK/output/volcanoNvsJ_"
# outfile2<-".png"
# outfile<-paste(outfile1,padj_cutoff,outfile2, sep="")
# 
# ggsave(outfile, 
#        dpi=300, 
#        device="png",
#        width=23, 
#        height=14,
#        units="cm")
# 
# # 
# # # Plot continuous
# # ggplot(df_log, aes(x=log2FoldChange, y=logp, color=padj)) + 
# #   geom_point(na.rm = TRUE, alpha=0.3) + 
# #   geom_vline(xintercept=lfc, linetype="dotted") + 
# #   theme(text = element_text(size=allfont),
# #         axis.text = element_text(size=axisfont)) + 
# #   ggtitle("November vs. June") +
# #   ylab("\u2013log10(P-Value)") + 
# #   xlab("log2(FoldChange)")  +
# #   labs(color='Adj P-Value' ) +
# #   xlim(-8,8)
# # 
# # outfile1 <- "/Users/sunshine2.0/Desktop/WORK/output/volcano_continuousNvsJ.png"  
# # 
# # ggsave(outfile1, 
# #        dpi=300, 
# #        device="png", 
# #        width=23, 
# #        height=14, 
# #        units="cm")
# 
# # APRIL VS NOVEMBER
# 
# # Subset the relevant columns from the full dataset and choose the control condition so log2 fold changes are calculated as treatment over control
# dds <- ddsAll[ , (ddsAll$month == "N")|(ddsAll$month == "A") ]
# dds$month <- droplevels( dds$month )
# dds$month <- relevel(dds$month, ref = "N")
# # check if subsetted correctly
# #as.data.frame( colData(dds) )
# 
# dds <- DESeq(dds)
# res <- results(dds, alpha=padj_cutoff)
# 
# # number of adjusted p-values less than 0.1
# sum(res$padj < padj_cutoff, na.rm=TRUE)
# # summary (default lfcThreshold = 0)
# summary(res)
# sum(res$baseMean ==0)
# 
# # for each condition vs control, plot a volcano plot
# # convert to dataframe
# 
# df <- as.data.frame(res)
# 
# # make summary table similar to summary(res)
# total<- nrow(df)-sum(df$baseMean ==0, na.rm=TRUE)
# lowCountsSubset<-df[ ((!is.na(df$pvalue))&(is.na(df$padj))),]
# dfsummary3 <-  data.frame("LFCup" = c(round(sum(df$log2FoldChange>0& df$padj<padj_cutoff, na.rm=TRUE), digits=0),round(sum(df$log2FoldChange>0& df$padj<padj_cutoff, na.rm=TRUE)/total*100, digits = 1)), 
#                           "LFCdown" = c(round(sum(df$log2FoldChange<0& df$padj<padj_cutoff, na.rm=TRUE), digits = 0),round(sum(df$log2FoldChange<0& df$padj<padj_cutoff, na.rm=TRUE)/total*100, digits = 1)), 
#                           "outliers"= c(round(sum(!is.na(df$log2FoldChange)&is.na(df$pvalue)), digits = 0),round(sum(!is.na(df$log2FoldChange)&is.na(df$pvalue))/round(sum(is.na(df$padj)), digits = 0)*100, digits = 1)),
#                           "lowCounts"= c(round(sum(!is.na(df$pvalue)&is.na(df$padj)), digits = 0),round(sum(!is.na(df$pvalue)&is.na(df$padj))/round(sum(is.na(df$padj)), digits = 0)*100, digits = 1)),
#                           "baseMean"= c(round(sum(df$baseMean ==0, na.rm=TRUE), digits = 0),round(sum(df$baseMean ==0, na.rm=TRUE)/round(sum(is.na(df$padj)), digits = 0)*100, digits = 1)),
#                           "total"= c(round(total, digits = 0), round(sum((sum(df$log2FoldChange>0& df$padj<padj_cutoff, na.rm=TRUE)/total*100),(sum(df$log2FoldChange<0& df$padj<padj_cutoff, na.rm=TRUE)/total*100),(sum(!is.na(df$log2FoldChange)&is.na(df$pvalue))/round(sum(is.na(df$padj)), digits = 0)*100),(sum(!is.na(df$pvalue)&is.na(df$padj))/round(sum(is.na(df$padj)), digits = 0)*100),(sum(df$baseMean ==0, na.rm=TRUE)/round(sum(is.na(df$padj)), digits = 0)*100), digits = 1))),
#                           "treatment" = "MF_NvsA",
#                           "maxlowCounts"=c(round(max(lowCountsSubset$baseMean), digits=0),NA))
# 
# # make a data table of df statistics
# dfstats3<-data.frame("#NA"=c(round(sum(is.na(df$padj)), digits = 0),NA),
#                      "#sig"=c(sum(df$padj < padj_cutoff, na.rm=TRUE),round(sum(df$padj < padj_cutoff, na.rm=TRUE)/(nrow(df)-sum(is.na(df$padj)))*100, digits = 1)),
#                      "#nonsig"= c(sum(df$padj > padj_cutoff | df$padj == padj_cutoff, na.rm=TRUE),round(sum(df$padj > padj_cutoff | df$padj == padj_cutoff, na.rm=TRUE)/(nrow(df)-sum(is.na(df$padj)))*100, digits = 1)),
#                      "total"= c(nrow(df), NA),
#                      "treatment" = "MF_NvsA")
# 
# # add -log(pvalue) column, and padj category 
# df_log <- mutate(df, logp=-log10(pvalue), padjcat=padj,
#                  lfccat=log2FoldChange)
# 
# df_log <- transform(df_log, padjcat=case_when(padj < padj_cutoff & padj > padj_cutoff1 ~ paste(padj_cutoff1,"< padj <", padj_cutoff),
#                                               padj >= padj_cutoff ~ paste("padj >=", padj_cutoff),
#                                               padj < padj_cutoff1 ~ paste("padj <", padj_cutoff1)))
# 
# df_log<- transform(df_log, lfccat= case_when(lfccat > lfc_cutoff ~ paste("lfc > +", lfc_cutoff), 
#                                              lfccat < -lfc_cutoff ~ paste("lfc < -", lfc_cutoff),
#                                              lfccat <= lfc_cutoff & lfccat >= (-lfc_cutoff) ~ paste("-",lfc_cutoff," =< lfc =< +", lfc_cutoff)))
# 
# df_log <- mutate(df_log, cat=paste(padjcat,";",lfccat))
# df_log$cat[is.na(df_log$padjcat)] <- NA
# df_log$cat[df_log$lfccat == paste("-",lfc_cutoff," =< lfc =< +", lfc_cutoff)] <- NA
# df_log$cat[df_log$padjcat == paste("padj >=", padj_cutoff)] <- NA
# 
# upcat<-df_log %>%
#   filter(cat == paste(padj_cutoff1,"< padj <", padj_cutoff,";","lfc > +", lfc_cutoff) ) %>%
#   nrow()
# 
# upperupcat<-df_log %>%
#   filter(cat == paste("padj <", padj_cutoff1,";","lfc > +", lfc_cutoff) ) %>%
#   nrow()
# 
# 
# downcat<- df_log %>%
#   filter(cat == paste(padj_cutoff1,"< padj <", padj_cutoff,";","lfc < -", lfc_cutoff) ) %>%
#   nrow()
# 
# upperdowncat<- df_log %>%
#   filter(cat == paste("padj <", padj_cutoff1,";","lfc < -", lfc_cutoff) ) %>%
#   nrow()
# 
# 
# #subset significant genes into list of gene names only
# df_log <- mutate(df_log, 'gene'=rownames(df))
# DEsubset3 <- df_log[  !is.na(df_log$cat) , ]
# DEgenelist3 <- DEsubset3$gene
# 
# DEsubset.01 <- df_log %>%
#   filter(cat== paste("padj <", padj_cutoff1,";","lfc < -", lfc_cutoff) | cat== paste("padj <", padj_cutoff1,";","lfc > +", lfc_cutoff) )
# DEgenelist3.01 <- DEsubset.01$gene
# 
# 
# DEsubset.up01 <- df_log %>%
#   filter(cat== paste("padj <", padj_cutoff1,";","lfc > +", lfc_cutoff) )
# DEgenelist3.up01 <- DEsubset.up01$gene
# 
# DEsubset.down01 <- df_log %>%
#   filter(cat== paste("padj <", padj_cutoff1,";","lfc < -", lfc_cutoff) )
# DEgenelist3.down01 <- DEsubset.down01$gene
# 
# DEsubset.up05 <-  df_log[  !is.na(df_log$cat) , ]
# DEsubset.up05 <- DEsubset.up05 %>%
#   filter(lfccat== paste("lfc > +", lfc_cutoff) )
# DEgenelist3.up05 <- DEsubset.up05$gene
# 
# DEsubset.down05 <-  df_log[  !is.na(df_log$cat) , ]
# DEsubset.down05 <- DEsubset.down05 %>%
#   filter(lfccat== paste("lfc < -", lfc_cutoff) )
# DEgenelist3.down05 <- DEsubset.down05$gene
# 
# # Default font sizes
# axisfont <- 14
# labelfont <- 5
# allfont <- 18
# 
# # Default for this type of plot
# legendtitle<-"adjusted\np-value"
# 
# # reassign rownames
# rownames(df_log) <- rownames(df)
# 
# # convert padjcat to factor
# df_log$padjcat <- as.factor(df_log$padjcat)
# 
# # Plot discrete
# p3<-ggplot(df_log, aes(x=log2FoldChange, y=logp, color=cat)) + 
#   geom_point(na.rm = TRUE,alpha=0.5) +
#   geom_vline(xintercept=c(-1,1),linetype="dotted") + 
#   geom_hline(yintercept=c(2.1,3.37) ,linetype="dotted") + 
#   theme(text = element_text(size=allfont),
#         axis.text = element_text(size=axisfont)) + 
#   ggtitle("April vs. Nov") +
#   ylab("\u2013log10(P-Value)") + 
#   theme_bw()+
#   xlab("log2(FoldChange)") +
#   labs(color='Adj P-Value' ) +
#   xlim(-8,8)
# # xlimit() +
# # ylimit()
# ylim<-ggplot_build(p3)$layout$panel_scales_y[[1]]$range$range[2]
# 
# p3<-p3  + geom_text(x=7.2, y=ylim-3*(ylim/5), label= upcat, color='olivedrab4') 
# p3<-p3 + geom_text(x=-7.2, y=ylim-3*(ylim/5), label= downcat, color='red') 
# p3<-p3 + geom_text(x=7.2, y=ylim-ylim/4, label= upperupcat, color='purple') 
# p3<-p3 + geom_text(x=-7.2, y=ylim-ylim/4, label= upperdowncat, color='darkcyan') 
# 
# outfile1<-"/Users/sunshine2.0/Desktop/WORK/output/volcanoAvaN_"
# outfile2<-".png"
# outfile<-paste(outfile1,padj_cutoff,outfile2, sep="")
# 
# ggsave(outfile, 
#        dpi=300, 
#        device="png",
#        width=23, 
#        height=14,
#        units="cm")
# 
# 
# # # Plot continuous
# # ggplot(df_log, aes(x=log2FoldChange, y=logp, color=padj)) + 
# #   geom_point(na.rm = TRUE, alpha=0.3) + 
# #   geom_vline(xintercept=lfc, linetype="dotted") + 
# #   theme(text = element_text(size=allfont),
# #         axis.text = element_text(size=axisfont)) + 
# #   ggtitle("November vs. April") +
# #   ylab("\u2013log10(P-Value)") + 
# #   xlab("log2(FoldChange)")  +
# #   labs(color='Adj P-Value' )
# # 
# # outfile1 <- "/home/krwright/Downloads/volcano_continuousNvsA.png"  
# # 
# # ggsave(outfile1, 
# #        dpi=300, 
# #        device="png", 
# #        width=23, 
# #        height=14, 
# #        units="cm")
# 

#_________________________________________ FEMALE SUBSETS_________________________________________________________

# APRIL

# Subset the relevant columns from the full dataset and choose the control condition so log2 fold changes are calculated as treatment over control
dds <- ddsAll[ , (ddsAll$month == "A")]
dds$gender <- droplevels( dds$gender )
dds$gender <- relevel(dds$gender, ref = "F")
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
DEsubset4 <- df_log[  !is.na(df_log$cat) , ]
DEgenelist4 <- DEsubset4$gene #all significant genes meeting padj cutoff <0.05

DEsubset.01 <- df_log %>%
  filter(cat== paste("P-Value <", padj_cutoff1,";","LFC < -", lfc_cutoff) | cat== paste("P-Value <", padj_cutoff1,";","LFC > +", lfc_cutoff) )
DEgenelist4.01 <- DEsubset.01$gene #all significant genes meeting padj cutoff <0.01

DEsubset.up01 <- df_log %>%
  filter(cat== paste("P-Value <", padj_cutoff1,";","LFC > +", lfc_cutoff) )
DEgenelist4.up01 <- DEsubset.up01$gene #all up-regulated genes meeting padj cutoff <0.01

DEsubset.down01 <- df_log %>%
  filter(cat== paste("P-Value <", padj_cutoff1,";","LFC < -", lfc_cutoff) )
DEgenelist4.down01 <- DEsubset.down01$gene #all down-regulated genes meeting padj cutoff <0.01

DEsubset.up05 <-  df_log[  !is.na(df_log$cat) , ]
DEsubset.up05 <- DEsubset.up05 %>%
  filter(lfccat== paste("LFC > +", lfc_cutoff) )
DEgenelist4.up05 <- DEsubset.up05$gene #all up-regulated genes meeting padj cutoff <0.05

DEsubset.down05 <-  df_log[  !is.na(df_log$cat) , ]
DEsubset.down05 <- DEsubset.down05 %>%
  filter(lfccat== paste("LFC < -", lfc_cutoff) )
DEgenelist4.down05 <- DEsubset.down05$gene #all down-regulated genes meeting padj cutoff <0.05

DE_padj1<-df_log%>% # Obtain list of DE0.01 genes with padj
  dplyr::select(c(gene, padj)) 

# Volcano plot font sizes
axisfont <- 35
labelfont <- 26
allfont <- 39

# reassign rownames
rownames(df_log) <- rownames(df)

# convert padjcat to factor
df_log$padjcat <- as.factor(df_log$padjcat)

# Replace NA with Not Significant or Substantial
df_log<-df_log %>%
  mutate(cat = if_else(is.na(cat), "Not Significant or Substantial", cat))
df_log$cat<-factor(df_log$cat,
                   levels=c("P-Value < 0.01 ; LFC < - 1" ,"P-Value < 0.01 ; LFC > + 1", "0.01 < P-Value < 0.05 ; LFC < - 1" , "0.01 < P-Value < 0.05 ; LFC > + 1", "Not Significant or Substantial"))
cols <- c("0.01 < P-Value < 0.05 ; LFC < - 1" = "#00ace6", "0.01 < P-Value < 0.05 ; LFC > + 1" = "#ff66cc", "P-Value < 0.01 ; LFC < - 1" = "#0000ff", "P-Value < 0.01 ; LFC > + 1" = "#ff0000", "Not Significant or Substantial"="#818181")

# Plot discrete volcano plot
p4<-ggplot(df_log, aes(x=log2FoldChange, y=logp, color=cat)) + 
  geom_point(na.rm = TRUE,alpha=0.5) +
  geom_vline(xintercept=c(-1,1),linetype="dotted") + 
  geom_hline(yintercept=c(4.46,5.2) ,linetype="dotted") + 
  theme(text = element_text(size=allfont),
  axis.text = element_text(size=axisfont)) + 
  ggtitle("April: Male vs Female") +
  ylab(expression(paste("\u2013", log[10], "(P-Value)"))) + 
  theme_bw()+
  xlab(expression(paste(log[2],"(FoldChange)"))) +
  labs(color='Genes' ) +
  xlim(-8,8) +
  scale_colour_manual(values = cols, drop=FALSE, name = "Genes", labels = c("P-Value < 0.01 ; LFC < - 1", "P-Value < 0.01 ; LFC > + 1", "0.01 \u2264 P-Value < 0.05 ; LFC < - 1", "0.01 \u2264 P-Value < 0.05 ; LFC > + 1",  "Not Significant or Substantial"))+
  guides(color = guide_legend(override.aes = list(size = 3)))

ylim<-ggplot_build(p4)$layout$panel_scales_y[[1]]$range$range[2]

p4<-p4  + geom_text(x=7.25, y=ylim-1.3*(ylim/8.25), label= upcat, color='#ff66cc') 
p4<-p4 + geom_text(x=-7.25, y=ylim-1.3*(ylim/8.25), label= downcat, color='#00ace6') 
p4<-p4 + geom_text(x=7.25, y=ylim-.3*(ylim/8.25), label= upperupcat, color='#ff0000') 
p4<-p4 + geom_text(x=-7.25, y=ylim-.3*(ylim/8.25), label= upperdowncat, color='#0000ff') 

outfile1<-"/Users/sunshine2.0/Desktop/WORK/output/volcanoApril_MvsF"
outfile2<-".png"
outfile<-paste(outfile1,padj_cutoff,outfile2, sep="")

ggsave(outfile, 
       dpi=300, 
       device="png",
       width=23, 
       height=14,
       units="cm")

# # Plot continuous volacno plot
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

# JUNE 

# Subset the relevant columns from the full dataset and choose the control condition so log2 fold changes are calculated as treatment over control
dds <- ddsAll[ , (ddsAll$month == "J")]
dds$gender <- droplevels( dds$gender )
dds$gender <- relevel(dds$gender, ref = "F")
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
DEgenelist5 <- DEsubset5$gene #all singificant genes meeting padj cutoff <0.05

DEsubset.01 <- df_log %>%
  filter(cat== paste("padj <", padj_cutoff1,";","lfc < -", lfc_cutoff) | cat== paste("padj <", padj_cutoff1,";","lfc > +", lfc_cutoff) )
DEgenelist5.01 <- DEsubset.01$gene #all significant genes meeting padj cutoff <0.01

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

DE_padj2<-df_log%>% # Obtain list of DE0.01 genes with padj
  dplyr::select(c(gene, padj)) 

# Volcano plot font sizes
axisfont <- 35
labelfont <- 26
allfont <- 39

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
  geom_hline(yintercept=c(4.9, 6.1) ,linetype="dotted") + 
  theme(text = element_text(size=allfont),
        axis.text = element_text(size=axisfont)) + 
  ggtitle("June: Male vs Female") +
  ylab(expression(paste("\u2013", log[10], "(P-Value)"))) + 
  theme_bw()+
  xlab(expression(paste(log[2],"(FoldChange)"))) +
  labs(color='Adj P-Value' ) +
  xlim(-8,8)+
  scale_colour_manual(values = cols, drop=FALSE, name = "Genes", labels = c("P-Value < 0.01 ; LFC < - 1", "P-Value < 0.01 ; LFC > + 1", "0.01 \u2264 P-Value < 0.05 ; LFC < - 1", "0.01 \u2264 P-Value < 0.05 ; LFC > + 1",  "Not Significant or Substantial"))
  
ylim<-ggplot_build(p5)$layout$panel_scales_y[[1]]$range$range[2]

p5<-p5  + geom_text(x=7.25, y=ylim-1.3*(ylim/8.25), label= upcat, color='#ff66cc') 
p5<-p5 + geom_text(x=-7.25, y=ylim-1.3*(ylim/8.25), label= downcat, color='#00ace6') 
p5<-p5 + geom_text(x=7.25, y=ylim-.3*(ylim/8.25), label= upperupcat, color='#ff0000') 
p5<-p5 + geom_text(x=-7.25, y=ylim-.3*(ylim/8.25), label= upperdowncat, color='#0000ff') 

outfile1<-"/Users/sunshine2.0/Desktop/WORK/output/volcanoJune_MvsF"
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

#  NOVEMBER

# Subset the relevant columns from the full dataset and choose the control condition so log2 fold changes are calculated as treatment over control
dds <- ddsAll[ , (ddsAll$month == "N")]
dds$gender <- droplevels( dds$gender )
dds$gender <- relevel(dds$gender, ref = "F")
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

DE_padj3<-df_log%>% # Obtain list of DE0.01 genes with padj
  dplyr::select(c(gene, padj)) 

# Volcano plot font sizes
axisfont <- 35
labelfont <- 26
allfont <- 39

# Default for this type of plot
legendtitle<-"adjusted\np-value"

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
  geom_hline(yintercept=c(3.58,4.65) ,linetype="dotted") + 
  theme(text = element_text(size=allfont),
        axis.text = element_text(size=axisfont)) + 
  ggtitle("Nov: Male vs Female") +
  ylab(expression(paste("\u2013", log[10], "(P-Value)"))) + 
  theme_bw()+
  xlab(expression(paste(log[2],"(FoldChange)"))) +
  labs(color='Adj P-Value' ) +
  xlim(-8,8)+
  scale_colour_manual(values = cols, drop=FALSE, name = "Genes", labels = c("P-Value < 0.01 ; LFC < - 1", "P-Value < 0.01 ; LFC > + 1", "0.01 \u2264 P-Value < 0.05 ; LFC < - 1", "0.01 \u2264 P-Value < 0.05 ; LFC > + 1",  "Not Significant or Substantial"))
  

ylim<-ggplot_build(p6)$layout$panel_scales_y[[1]]$range$range[2]

p6<-p6  + geom_text(x=7.25, y=ylim-1.3*(ylim/8.25), label= upcat, color='#ff66cc') 
p6<-p6 + geom_text(x=-7.25, y=ylim-1.3*(ylim/8.25), label= downcat, color='#00ace6') 
p6<-p6 + geom_text(x=7.25, y=ylim-.3*(ylim/8.25), label= upperupcat, color='#ff0000') 
p6<-p6 + geom_text(x=-7.25, y=ylim-.3*(ylim/8.25), label= upperdowncat, color='#0000ff') 

outfile1<-"/Users/sunshine2.0/Desktop/WORK/output/volcanoNov_MvsF_"
outfile2<-".png"
outfile<-paste(outfile1,padj_cutoff,outfile2, sep="")

ggsave(outfile, 
       dpi=300, 
       device="png",
       width=23, 
       height=14,
       units="cm")

# 
# # Plot continuous volcano plot
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


# 
# #_________________________________________ MALE SUBSETS_________________________________________________________
# 
# # JUNE VS APRIL
# 
# # Subset the relevant columns from the full dataset and choose the control condition so log2 fold changes are calculated as treatment over control
# ddsMF <- ddsAll[ , (ddsAll$month == "A")|(ddsAll$month == "J") ]
# dds <- ddsMF[ , (ddsMF$gender == "M") ]
# dds$month <- droplevels( dds$month )
# dds$month <- relevel(dds$month, ref = "A")
# # check if subsetted correctly
# #as.data.frame( colData(dds) )
# 
# dds <- DESeq(dds)
# res <- results(dds, alpha=padj_cutoff)
# 
# # number of adjusted p-values less than 0.1
# sum(res$padj < padj_cutoff, na.rm=TRUE)
# # summary (default lfcThreshold = 0)
# summary(res)
# sum(res$baseMean ==0)
# 
# # for each condition vs control, plot a volcano plot
# # convert to dataframe
# 
# df <- as.data.frame(res)
# 
# # make summary table similar to summary(res)
# total<- nrow(df)-sum(df$baseMean ==0, na.rm=TRUE)
# lowCountsSubset<-df[ ((!is.na(df$pvalue))&(is.na(df$padj))),]
# dfsummary7 <-  data.frame("LFCup" = c(round(sum(df$log2FoldChange>0& df$padj<padj_cutoff, na.rm=TRUE), digits=0),round(sum(df$log2FoldChange>0& df$padj<padj_cutoff, na.rm=TRUE)/total*100, digits = 1)), 
#                           "LFCdown" = c(round(sum(df$log2FoldChange<0& df$padj<padj_cutoff, na.rm=TRUE), digits = 0),round(sum(df$log2FoldChange<0& df$padj<padj_cutoff, na.rm=TRUE)/total*100, digits = 1)), 
#                           "outliers"= c(round(sum(!is.na(df$log2FoldChange)&is.na(df$pvalue)), digits = 0),round(sum(!is.na(df$log2FoldChange)&is.na(df$pvalue))/round(sum(is.na(df$padj)), digits = 0)*100, digits = 1)),
#                           "lowCounts"= c(round(sum(!is.na(df$pvalue)&is.na(df$padj)), digits = 0),round(sum(!is.na(df$pvalue)&is.na(df$padj))/round(sum(is.na(df$padj)), digits = 0)*100, digits = 1)),
#                           "baseMean"= c(round(sum(df$baseMean ==0, na.rm=TRUE), digits = 0),round(sum(df$baseMean ==0, na.rm=TRUE)/round(sum(is.na(df$padj)), digits = 0)*100, digits = 1)),
#                           "total"= c(round(total, digits = 0), round(sum((sum(df$log2FoldChange>0& df$padj<padj_cutoff, na.rm=TRUE)/total*100),(sum(df$log2FoldChange<0& df$padj<padj_cutoff, na.rm=TRUE)/total*100),(sum(!is.na(df$log2FoldChange)&is.na(df$pvalue))/round(sum(is.na(df$padj)), digits = 0)*100),(sum(!is.na(df$pvalue)&is.na(df$padj))/round(sum(is.na(df$padj)), digits = 0)*100),(sum(df$baseMean ==0, na.rm=TRUE)/round(sum(is.na(df$padj)), digits = 0)*100), digits = 1))),
#                           "treatment" = "M_JvsA",
#                           "maxlowCounts"=c(round(max(lowCountsSubset$baseMean), digits=0),NA))
# 
# # make a data table of df statistics
# dfstats7<-data.frame("#NA"=c(round(sum(is.na(df$padj)), digits = 0),NA),
#                      "#sig"=c(sum(df$padj < padj_cutoff, na.rm=TRUE),round(sum(df$padj < padj_cutoff, na.rm=TRUE)/(nrow(df)-sum(is.na(df$padj)))*100, digits = 1)),
#                      "#nonsig"= c(sum(df$padj > padj_cutoff | df$padj == padj_cutoff, na.rm=TRUE),round(sum(df$padj > padj_cutoff | df$padj == padj_cutoff, na.rm=TRUE)/(nrow(df)-sum(is.na(df$padj)))*100, digits = 1)),
#                      "total"= c(nrow(df), NA),
#                      "treatment" = "M_JvsA")
# 
# # add -log(pvalue) column, and padj category 
# df_log <- mutate(df, logp=-log10(pvalue), padjcat=padj,
#                  lfccat=log2FoldChange)
# 
# df_log <- transform(df_log, padjcat=case_when(padj < padj_cutoff & padj > padj_cutoff1 ~ paste(padj_cutoff1,"< padj <", padj_cutoff),
#                                               padj >= padj_cutoff ~ paste("padj >=", padj_cutoff),
#                                               padj < padj_cutoff1 ~ paste("padj <", padj_cutoff1)))
# 
# df_log<- transform(df_log, lfccat= case_when(lfccat > lfc_cutoff ~ paste("lfc > +", lfc_cutoff), 
#                                              lfccat < -lfc_cutoff ~ paste("lfc < -", lfc_cutoff),
#                                              lfccat <= lfc_cutoff & lfccat >= (-lfc_cutoff) ~ paste("-",lfc_cutoff," =< lfc =< +", lfc_cutoff)))
# 
# df_log <- mutate(df_log, cat=paste(padjcat,";",lfccat))
# df_log$cat[is.na(df_log$padjcat)] <- NA
# df_log$cat[df_log$lfccat == paste("-",lfc_cutoff," =< lfc =< +", lfc_cutoff)] <- NA
# df_log$cat[df_log$padjcat == paste("padj >=", padj_cutoff)] <- NA
# 
# upcat<-df_log %>%
#   filter(cat == paste(padj_cutoff1,"< padj <", padj_cutoff,";","lfc > +", lfc_cutoff) ) %>%
#   nrow()
# 
# upperupcat<-df_log %>%
#   filter(cat == paste("padj <", padj_cutoff1,";","lfc > +", lfc_cutoff) ) %>%
#   nrow()
# 
# 
# downcat<- df_log %>%
#   filter(cat == paste(padj_cutoff1,"< padj <", padj_cutoff,";","lfc < -", lfc_cutoff) ) %>%
#   nrow()
# 
# upperdowncat<- df_log %>%
#   filter(cat == paste("padj <", padj_cutoff1,";","lfc < -", lfc_cutoff) ) %>%
#   nrow()
# 
# 
# #subset significant genes into list of gene names only
# df_log <- mutate(df_log, 'gene'=rownames(df))
# DEsubset7 <- df_log[  !is.na(df_log$cat) , ]
# DEgenelist7 <- DEsubset7$gene
# 
# DEsubset.01 <- df_log %>%
#   filter(cat== paste("padj <", padj_cutoff1,";","lfc < -", lfc_cutoff) | cat== paste("padj <", padj_cutoff1,";","lfc > +", lfc_cutoff) )
# DEgenelist7.01 <- DEsubset.01$gene
# 
# 
# DEsubset.up01 <- df_log %>%
#   filter(cat== paste("padj <", padj_cutoff1,";","lfc > +", lfc_cutoff) )
# DEgenelist7.up01 <- DEsubset.up01$gene
# 
# DEsubset.down01 <- df_log %>%
#   filter(cat== paste("padj <", padj_cutoff1,";","lfc < -", lfc_cutoff) )
# DEgenelist7.down01 <- DEsubset.down01$gene
# 
# DEsubset.up05 <-  df_log[  !is.na(df_log$cat) , ]
# DEsubset.up05 <- DEsubset.up05 %>%
#   filter(lfccat== paste("lfc > +", lfc_cutoff) )
# DEgenelist7.up05 <- DEsubset.up05$gene
# 
# DEsubset.down05 <-  df_log[  !is.na(df_log$cat) , ]
# DEsubset.down05 <- DEsubset.down05 %>%
#   filter(lfccat== paste("lfc < -", lfc_cutoff) )
# DEgenelist7.down05 <- DEsubset.down05$gene
# 
# # Default font sizes
# axisfont <- 14
# labelfont <- 5
# allfont <- 18
# 
# # Default for this type of plot
# legendtitle<-"adjusted\np-value"
# 
# # reassign rownames
# rownames(df_log) <- rownames(df)
# 
# # convert padjcat to factor
# df_log$padjcat <- as.factor(df_log$padjcat)
# 
# 
# 
# # Plot discrete
# p7<-ggplot(df_log, aes(x=log2FoldChange, y=logp, color=cat)) + 
#   geom_point(na.rm = TRUE,alpha=0.5) +
#   geom_vline(xintercept=c(-1,1),linetype="dotted") + 
#   geom_hline(yintercept=c(2.43,3.68) ,linetype="dotted") + 
#   theme(text = element_text(size=allfont),
#         axis.text = element_text(size=axisfont)) + 
#   ggtitle("Male: June vs. April") +
#   ylab("\u2013log10(P-Value)") + 
#   theme_bw()+
#   xlab("log2(FoldChange)") +
#   labs(color='Adj P-Value' ) +
#   xlim(-8,8)
# # xlimit() +
# # ylimit()
# ylim<-ggplot_build(p7)$layout$panel_scales_y[[1]]$range$range[2]
# 
# p7<-p7  + geom_text(x=7.2, y=ylim-3*(ylim/5), label= upcat, color='olivedrab4') 
# p7<-p7 + geom_text(x=-7.2, y=ylim-3*(ylim/5), label= downcat, color='red') 
# p7<-p7 + geom_text(x=7.2, y=ylim-ylim/4, label= upperupcat, color='purple') 
# p7<-p7 + geom_text(x=-7.2, y=ylim-ylim/4, label= upperdowncat, color='darkcyan') 
# 
#   outfile1<-"/Users/sunshine2.0/Desktop/WORK/output/volcanoMaleJvsA_"
#   outfile2<-".png"
#   outfile<-paste(outfile1,padj_cutoff,outfile2, sep="")
# 
# ggsave(outfile, 
#        dpi=300, 
#        device="png",
#        width=23, 
#        height=14,
#        units="cm")
# 
# 
# # Plot continuous
# # ggplot(df_log, aes(x=log2FoldChange, y=logp, color=padj)) + 
# #   geom_point(na.rm = TRUE, alpha=0.3) + 
# #   geom_vline(xintercept=lfc, linetype="dotted") + 
# #   theme(text = element_text(size=allfont),
# #         axis.text = element_text(size=axisfont)) + 
# #   ggtitle("Male: June vs. April") +
# #   ylab("\u2013log10(P-Value)") + 
# #   xlab("log2(FoldChange)")  +
# #   labs(color='Adj P-Value' )
# # 
# # outfile1 <- "/home/krwright/Downloads/volcano_continuousMaleJvsA.png"  
# # 
# # ggsave(outfile1, 
# #        dpi=300, 
# #        device="png", 
# #        width=23, 
# #        height=14, 
# #        units="cm")
# 
# # NOVEMBER VS JUNE 
# 
# # Subset the relevant columns from the full dataset and choose the control condition so log2 fold changes are calculated as treatment over control
# ddsMF <- ddsAll[ , (ddsAll$month == "N")|(ddsAll$month == "J") ]
# dds <- ddsMF[ , (ddsMF$gender == "M") ]
# dds$month <- droplevels( dds$month )
# dds$month <- relevel(dds$month, ref = "J")
# # check if subsetted correctly
# #as.data.frame( colData(dds) )
# 
# dds <- DESeq(dds)
# res <- results(dds, alpha=padj_cutoff)
# 
# 
# # number of adjusted p-values less than 0.1
# sum(res$padj < padj_cutoff, na.rm=TRUE)
# # summary (default lfcThreshold = 0)
# summary(res)
# sum(res$baseMean ==0)
# 
# # for each condition vs control, plot a volcano plot
# # convert to dataframe
# 
# df <- as.data.frame(res)
# 
# # make summary table similar to summary(res)
# total<- nrow(df)-sum(df$baseMean ==0, na.rm=TRUE)
# lowCountsSubset<-df[ ((!is.na(df$pvalue))&(is.na(df$padj))),]
# dfsummary8 <-  data.frame("LFCup" = c(round(sum(df$log2FoldChange>0& df$padj<padj_cutoff, na.rm=TRUE), digits=0),round(sum(df$log2FoldChange>0& df$padj<padj_cutoff, na.rm=TRUE)/total*100, digits = 1)), 
#                           "LFCdown" = c(round(sum(df$log2FoldChange<0& df$padj<padj_cutoff, na.rm=TRUE), digits = 0),round(sum(df$log2FoldChange<0& df$padj<padj_cutoff, na.rm=TRUE)/total*100, digits = 1)), 
#                           "outliers"= c(round(sum(!is.na(df$log2FoldChange)&is.na(df$pvalue)), digits = 0),round(sum(!is.na(df$log2FoldChange)&is.na(df$pvalue))/round(sum(is.na(df$padj)), digits = 0)*100, digits = 1)),
#                           "lowCounts"= c(round(sum(!is.na(df$pvalue)&is.na(df$padj)), digits = 0),round(sum(!is.na(df$pvalue)&is.na(df$padj))/round(sum(is.na(df$padj)), digits = 0)*100, digits = 1)),
#                           "baseMean"= c(round(sum(df$baseMean ==0, na.rm=TRUE), digits = 0),round(sum(df$baseMean ==0, na.rm=TRUE)/round(sum(is.na(df$padj)), digits = 0)*100, digits = 1)),
#                           "total"= c(round(total, digits = 0), round(sum((sum(df$log2FoldChange>0& df$padj<padj_cutoff, na.rm=TRUE)/total*100),(sum(df$log2FoldChange<0& df$padj<padj_cutoff, na.rm=TRUE)/total*100),(sum(!is.na(df$log2FoldChange)&is.na(df$pvalue))/round(sum(is.na(df$padj)), digits = 0)*100),(sum(!is.na(df$pvalue)&is.na(df$padj))/round(sum(is.na(df$padj)), digits = 0)*100),(sum(df$baseMean ==0, na.rm=TRUE)/round(sum(is.na(df$padj)), digits = 0)*100), digits = 1))),
#                           "treatment" = "M_NvsJ",
#                           "maxlowCounts"=c(round(max(lowCountsSubset$baseMean), digits=0),NA))
# 
# # make a data table of df statistics
# dfstats8<-data.frame("#NA"=c(round(sum(is.na(df$padj)), digits = 0),NA),
#                      "#sig"=c(sum(df$padj < padj_cutoff, na.rm=TRUE),round(sum(df$padj < padj_cutoff, na.rm=TRUE)/(nrow(df)-sum(is.na(df$padj)))*100, digits = 1)),
#                      "#nonsig"= c(sum(df$padj > padj_cutoff | df$padj == padj_cutoff, na.rm=TRUE),round(sum(df$padj > padj_cutoff | df$padj == padj_cutoff, na.rm=TRUE)/(nrow(df)-sum(is.na(df$padj)))*100, digits = 1)),
#                      "total"= c(nrow(df), NA),
#                      "treatment" = "M_NvsJ")
# 
# # add -log(pvalue) column, and padj category 
# df_log <- mutate(df, logp=-log10(pvalue), padjcat=padj,
#                  lfccat=log2FoldChange)
# 
# df_log <- transform(df_log, padjcat=case_when(padj < padj_cutoff & padj > padj_cutoff1 ~ paste(padj_cutoff1,"< padj <", padj_cutoff),
#                                               padj >= padj_cutoff ~ paste("padj >=", padj_cutoff),
#                                               padj < padj_cutoff1 ~ paste("padj <", padj_cutoff1)))
# 
# df_log<- transform(df_log, lfccat= case_when(lfccat > lfc_cutoff ~ paste("lfc > +", lfc_cutoff), 
#                                              lfccat < -lfc_cutoff ~ paste("lfc < -", lfc_cutoff),
#                                              lfccat <= lfc_cutoff & lfccat >= (-lfc_cutoff) ~ paste("-",lfc_cutoff," =< lfc =< +", lfc_cutoff)))
# 
# df_log <- mutate(df_log, cat=paste(padjcat,";",lfccat))
# df_log$cat[is.na(df_log$padjcat)] <- NA
# df_log$cat[df_log$lfccat == paste("-",lfc_cutoff," =< lfc =< +", lfc_cutoff)] <- NA
# df_log$cat[df_log$padjcat == paste("padj >=", padj_cutoff)] <- NA
# 
# upcat<-df_log %>%
#   filter(cat == paste(padj_cutoff1,"< padj <", padj_cutoff,";","lfc > +", lfc_cutoff) ) %>%
#   nrow()
# 
# upperupcat<-df_log %>%
#   filter(cat == paste("padj <", padj_cutoff1,";","lfc > +", lfc_cutoff) ) %>%
#   nrow()
# 
# 
# downcat<- df_log %>%
#   filter(cat == paste(padj_cutoff1,"< padj <", padj_cutoff,";","lfc < -", lfc_cutoff) ) %>%
#   nrow()
# 
# upperdowncat<- df_log %>%
#   filter(cat == paste("padj <", padj_cutoff1,";","lfc < -", lfc_cutoff) ) %>%
#   nrow()
# 
# #subset significant genes into list of gene names only
# df_log <- mutate(df_log, 'gene'=rownames(df))
# DEsubset8 <- df_log[  !is.na(df_log$cat) , ]
# DEgenelist8 <- DEsubset8$gene
# 
# DEsubset.01 <- df_log %>%
#   filter(cat== paste("padj <", padj_cutoff1,";","lfc < -", lfc_cutoff) | cat== paste("padj <", padj_cutoff1,";","lfc > +", lfc_cutoff) )
# DEgenelist8.01 <- DEsubset.01$gene
# 
# 
# DEsubset.up01 <- df_log %>%
#   filter(cat== paste("padj <", padj_cutoff1,";","lfc > +", lfc_cutoff) )
# DEgenelist8.up01 <- DEsubset.up01$gene
# 
# DEsubset.down01 <- df_log %>%
#   filter(cat== paste("padj <", padj_cutoff1,";","lfc < -", lfc_cutoff) )
# DEgenelist8.down01 <- DEsubset.down01$gene
# 
# DEsubset.up05 <-  df_log[  !is.na(df_log$cat) , ]
# DEsubset.up05 <- DEsubset.up05 %>%
#   filter(lfccat== paste("lfc > +", lfc_cutoff) )
# DEgenelist8.up05 <- DEsubset.up05$gene
# 
# DEsubset.down05 <-  df_log[  !is.na(df_log$cat) , ]
# DEsubset.down05 <- DEsubset.down05 %>%
#   filter(lfccat== paste("lfc < -", lfc_cutoff) )
# DEgenelist8.down05 <- DEsubset.down05$gene
# 
# # Default font sizes
# axisfont <- 14
# labelfont <- 5
# allfont <- 18
# 
# # Default for this type of plot
# legendtitle<-"adjusted\np-value"
# 
# # reassign rownames
# rownames(df_log) <- rownames(df)
# 
# # convert padjcat to factor
# df_log$padjcat <- as.factor(df_log$padjcat)
# 
# 
# # Plot discrete
# p8<-ggplot(df_log, aes(x=log2FoldChange, y=logp, color=cat)) + 
#   geom_point(na.rm = TRUE,alpha=0.5) +
#   geom_vline(xintercept=c(-1,1),linetype="dotted") + 
#   geom_hline(yintercept=c(2.58,3.78) ,linetype="dotted") + 
#   theme(text = element_text(size=allfont),
#         axis.text = element_text(size=axisfont)) + 
#   ggtitle("Male: Nov vs. June") +
#   ylab("\u2013log10(P-Value)") + 
#   theme_bw()+
#   xlab("log2(FoldChange)") +
#   labs(color='Adj P-Value' ) +
#   xlim(-8,8)
# # xlimit() +
# # ylimit()
# ylim<-ggplot_build(p8)$layout$panel_scales_y[[1]]$range$range[2]
# 
# p8<-p8  + geom_text(x=7.2, y=ylim-3*(ylim/5), label= upcat, color='olivedrab4') 
# p8<-p8 + geom_text(x=-7.2, y=ylim-3*(ylim/5), label= downcat, color='red') 
# p8<-p8 + geom_text(x=7.2, y=ylim-ylim/4, label= upperupcat, color='purple') 
# p8<-p8 + geom_text(x=-7.2, y=ylim-ylim/4, label= upperdowncat, color='darkcyan') 
# 
# outfile1<-"/Users/sunshine2.0/Desktop/WORK/output/volcanoMaleNvsJ_"
# outfile2<-".png"
# outfile<-paste(outfile1,padj_cutoff,outfile2, sep="")
# 
# ggsave(outfile, 
#        dpi=300, 
#        device="png",
#        width=23, 
#        height=14,
#        units="cm")
# 
# 
# # # Plot continuous
# # ggplot(df_log, aes(x=log2FoldChange, y=logp, color=padj)) + 
# #   geom_point(na.rm = TRUE, alpha=0.3) + 
# #   geom_vline(xintercept=lfc, linetype="dotted") + 
# #   theme(text = element_text(size=allfont),
# #         axis.text = element_text(size=axisfont)) + 
# #   ggtitle("Male: November vs. June") +
# #   ylab("\u2013log10(P-Value)") + 
# #   xlab("log2(FoldChange)")  +
# #   labs(color='Adj P-Value' ) +
# #   xlim(-8,8)
# # 
# # outfile1 <- "/home/krwright/Downloads/volcano_continuousMaleNvsJ.png"  
# # 
# # ggsave(outfile1, 
# #        dpi=300, 
# #        device="png", 
# #        width=23, 
# #        height=14, 
# #        units="cm")
# 
# # APRIL vs. NOVEMBER 
# 
# # Subset the relevant columns from the full dataset and choose the control condition so log2 fold changes are calculated as treatment over control
# ddsMF <- ddsAll[ , (ddsAll$month == "N")|(ddsAll$month == "A") ]
# dds <- ddsMF[ , (ddsMF$gender == "M") ]
# dds$month <- droplevels( dds$month )
# dds$month <- relevel(dds$month, ref = "N")
# # check if subsetted correctly
# #as.data.frame( colData(dds) )
# 
# dds <- DESeq(dds)
# res <- results(dds, alpha=padj_cutoff)
# 
# 
# # number of adjusted p-values less than 0.1
# sum(res$padj < padj_cutoff, na.rm=TRUE)
# # summary (default lfcThreshold = 0)
# summary(res)
# sum(res$baseMean ==0)
# 
# # for each condition vs control, plot a volcano plot
# # convert to dataframe
# 
# df <- as.data.frame(res)
# 
# # make summary table similar to summary(res)
# total<- nrow(df)-sum(df$baseMean ==0, na.rm=TRUE)
# lowCountsSubset<-df[ ((!is.na(df$pvalue))&(is.na(df$padj))),]
# dfsummary9 <-  data.frame("LFCup" = c(round(sum(df$log2FoldChange>0& df$padj<padj_cutoff, na.rm=TRUE), digits=0),round(sum(df$log2FoldChange>0& df$padj<padj_cutoff, na.rm=TRUE)/total*100, digits = 1)), 
#                           "LFCdown" = c(round(sum(df$log2FoldChange<0& df$padj<padj_cutoff, na.rm=TRUE), digits = 0),round(sum(df$log2FoldChange<0& df$padj<padj_cutoff, na.rm=TRUE)/total*100, digits = 1)), 
#                           "outliers"= c(round(sum(!is.na(df$log2FoldChange)&is.na(df$pvalue)), digits = 0),round(sum(!is.na(df$log2FoldChange)&is.na(df$pvalue))/round(sum(is.na(df$padj)), digits = 0)*100, digits = 1)),
#                           "lowCounts"= c(round(sum(!is.na(df$pvalue)&is.na(df$padj)), digits = 0),round(sum(!is.na(df$pvalue)&is.na(df$padj))/round(sum(is.na(df$padj)), digits = 0)*100, digits = 1)),
#                           "baseMean"= c(round(sum(df$baseMean ==0, na.rm=TRUE), digits = 0),round(sum(df$baseMean ==0, na.rm=TRUE)/round(sum(is.na(df$padj)), digits = 0)*100, digits = 1)),
#                           "total"= c(round(total, digits = 0), round(sum((sum(df$log2FoldChange>0& df$padj<padj_cutoff, na.rm=TRUE)/total*100),(sum(df$log2FoldChange<0& df$padj<padj_cutoff, na.rm=TRUE)/total*100),(sum(!is.na(df$log2FoldChange)&is.na(df$pvalue))/round(sum(is.na(df$padj)), digits = 0)*100),(sum(!is.na(df$pvalue)&is.na(df$padj))/round(sum(is.na(df$padj)), digits = 0)*100),(sum(df$baseMean ==0, na.rm=TRUE)/round(sum(is.na(df$padj)), digits = 0)*100), digits = 1))),
#                           "treatment" = "M_NvsA",
#                           "maxlowCounts"=c(round(max(lowCountsSubset$baseMean), digits=0),NA))
# 
# # make a data table of df statistics
# dfstats9<-data.frame("#NA"=c(round(sum(is.na(df$padj)), digits = 0),NA),
#                      "#sig"=c(sum(df$padj < padj_cutoff, na.rm=TRUE),round(sum(df$padj < padj_cutoff, na.rm=TRUE)/(nrow(df)-sum(is.na(df$padj)))*100, digits = 1)),
#                      "#nonsig"= c(sum(df$padj > padj_cutoff | df$padj == padj_cutoff, na.rm=TRUE),round(sum(df$padj > padj_cutoff | df$padj == padj_cutoff, na.rm=TRUE)/(nrow(df)-sum(is.na(df$padj)))*100, digits = 1)),
#                      "total"= c(nrow(df), NA),
#                      "treatment" = "M_NvsA")
# 
# # add -log(pvalue) column, and padj category 
# df_log <- mutate(df, logp=-log10(pvalue), padjcat=padj,
#                  lfccat=log2FoldChange)
# 
# df_log <- transform(df_log, padjcat=case_when(padj < padj_cutoff & padj > padj_cutoff1 ~ paste(padj_cutoff1,"< padj <", padj_cutoff),
#                                               padj >= padj_cutoff ~ paste("padj >=", padj_cutoff),
#                                               padj < padj_cutoff1 ~ paste("padj <", padj_cutoff1)))
# 
# df_log<- transform(df_log, lfccat= case_when(lfccat > lfc_cutoff ~ paste("lfc > +", lfc_cutoff), 
#                                              lfccat < -lfc_cutoff ~ paste("lfc < -", lfc_cutoff),
#                                              lfccat <= lfc_cutoff & lfccat >= (-lfc_cutoff) ~ paste("-",lfc_cutoff," =< lfc =< +", lfc_cutoff)))
# 
# df_log <- mutate(df_log, cat=paste(padjcat,";",lfccat))
# df_log$cat[is.na(df_log$padjcat)] <- NA
# df_log$cat[df_log$lfccat == paste("-",lfc_cutoff," =< lfc =< +", lfc_cutoff)] <- NA
# df_log$cat[df_log$padjcat == paste("padj >=", padj_cutoff)] <- NA
# 
# upcat<-df_log %>%
#   filter(cat == paste(padj_cutoff1,"< padj <", padj_cutoff,";","lfc > +", lfc_cutoff) ) %>%
#   nrow()
# 
# upperupcat<-df_log %>%
#   filter(cat == paste("padj <", padj_cutoff1,";","lfc > +", lfc_cutoff) ) %>%
#   nrow()
# 
# 
# downcat<- df_log %>%
#   filter(cat == paste(padj_cutoff1,"< padj <", padj_cutoff,";","lfc < -", lfc_cutoff) ) %>%
#   nrow()
# 
# upperdowncat<- df_log %>%
#   filter(cat == paste("padj <", padj_cutoff1,";","lfc < -", lfc_cutoff) ) %>%
#   nrow()
# 
# 
# #subset significant genes into list of gene names only
# df_log <- mutate(df_log, 'gene'=rownames(df))
# DEsubset9 <- df_log[  !is.na(df_log$cat) , ]
# DEgenelist9 <- DEsubset9$gene
# 
# DEsubset.01 <- df_log %>%
#   filter(cat== paste("padj <", padj_cutoff1,";","lfc < -", lfc_cutoff) | cat== paste("padj <", padj_cutoff1,";","lfc > +", lfc_cutoff) )
# DEgenelist9.01 <- DEsubset.01$gene
# 
# 
# DEsubset.up01 <- df_log %>%
#   filter(cat== paste("padj <", padj_cutoff1,";","lfc > +", lfc_cutoff) )
# DEgenelist9.up01 <- DEsubset.up01$gene
# 
# DEsubset.down01 <- df_log %>%
#   filter(cat== paste("padj <", padj_cutoff1,";","lfc < -", lfc_cutoff) )
# DEgenelist9.down01 <- DEsubset.down01$gene
# 
# DEsubset.up05 <-  df_log[  !is.na(df_log$cat) , ]
# DEsubset.up05 <- DEsubset.up05 %>%
#   filter(lfccat== paste("lfc > +", lfc_cutoff) )
# DEgenelist9.up05 <- DEsubset.up05$gene
# 
# DEsubset.down05 <-  df_log[  !is.na(df_log$cat) , ]
# DEsubset.down05 <- DEsubset.down05 %>%
#   filter(lfccat== paste("lfc < -", lfc_cutoff) )
# DEgenelist9.down05 <- DEsubset.down05$gene
# 
# # Default font sizes
# axisfont <- 14
# labelfont <- 5
# allfont <- 18
# 
# # Default for this type of plot
# legendtitle<-"adjusted\np-value"
# 
# # reassign rownames
# rownames(df_log) <- rownames(df)
# 
# # convert padjcat to factor
# df_log$padjcat <- as.factor(df_log$padjcat)
# 
# 
# 
# # Plot discrete
# p9<-ggplot(df_log, aes(x=log2FoldChange, y=logp, color=cat)) + 
#   geom_point(na.rm = TRUE,alpha=0.5) +
#   geom_vline(xintercept=c(-1,1),linetype="dotted") + 
#   geom_hline(yintercept=c(2.72, 3.9) ,linetype="dotted") + 
#   theme(text = element_text(size=allfont),
#         axis.text = element_text(size=axisfont)) + 
#   ggtitle("Male: April vs. Nov") +
#   ylab("\u2013log10(P-Value)") + 
#   theme_bw()+
#   xlab("log2(FoldChange)") +
#   labs(color='Adj P-Value' ) +
#   xlim(-8,8)
# # xlimit() +
# # ylimit()
# ylim<-ggplot_build(p9)$layout$panel_scales_y[[1]]$range$range[2]
# 
# p9<-p9  + geom_text(x=7.2, y=ylim-3*(ylim/5), label= upcat, color='olivedrab4') 
# p9<-p9 + geom_text(x=-7.2, y=ylim-3*(ylim/5), label= downcat, color='red') 
# p9<-p9 + geom_text(x=7.2, y=ylim-ylim/4, label= upperupcat, color='purple') 
# p9<-p9 + geom_text(x=-7.2, y=ylim-ylim/4, label= upperdowncat, color='darkcyan') 
# 
# outfile1<-"/Users/sunshine2.0/Desktop/WORK/output/volcanoMaleAvaN_"
# outfile2<-".png"
# outfile<-paste(outfile1,padj_cutoff,outfile2, sep="")
# 
# ggsave(outfile, 
#        dpi=300, 
#        device="png",
#        width=23, 
#        height=14,
#        units="cm")


# # Plot continuous
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



## ARRANGE ALL VOLCANO PLOTS INTO ONE FIGURE

library("ggpubr")

figure <-ggarrange(p4 + rremove("xylab")+ labs(title = NULL), p5 + rremove("xylab")+ labs(title = NULL), p6 + rremove("xylab")+ labs(title = NULL),
          common.legend = TRUE,
          legend = "right",
          widths = rep(5), 
          hjust = -0.5,
          vjust = 1,
          ncol = 3, nrow = 1,
          align = "hv") 

figure1 <-annotate_figure(figure,
                          bottom = text_grob(expression(paste(log[2],"(Fold Change)   ")), hjust=1.15),
                          left = text_grob(expression(paste("\u2013", log[10], "(P-Value)")),hjust=0.35, rot = 90))
annotate_figure(figure1,
                          top = text_grob("     April                                               June                                            November                                       ", face="bold"),
                          )  

outfile1.999<-"/Users/sunshine2.0/Desktop/WORK/output/grid.arrange.sexuallyDE_"
outfile2.999<-".png"
outfile.999<-paste(outfile1.999,padj_cutoff,outfile2.999, sep="")

ggsave(outfile.999, 
       dpi=300, 
       device="png",
       width=26, 
       height=7,
       units="cm")
                
# EXPORT LIST OF SEXUALLY DE GENES-----
require('biomaRt')

# Compile list of up and down regulated genes
ensembl_gene_id<-c(DEgenelist4.up01, DEgenelist4.down01, DEgenelist5.up01, DEgenelist5.down01, DEgenelist6.up01, DEgenelist6.down01)
regulated<-c(rep("up", length(DEgenelist4.up01)), rep("down", length(DEgenelist4.down01)), rep("up", length(DEgenelist5.up01)), rep("down", length(DEgenelist5.down01)), rep("up", length(DEgenelist6.up01)), rep("down", length(DEgenelist6.down01)))
month_comparison<-c(rep("April", length(DEgenelist4.up01)+ length(DEgenelist4.down01)), rep("June", length(DEgenelist5.up01)+ length(DEgenelist5.down01)), rep("Nov", length(DEgenelist6.up01)+ length(DEgenelist6.down01)))
DE_genes<-data.frame(ensembl_gene_id, regulated, month_comparison, stringsAsFactors = FALSE)

# Remove .n from gene ID
DE_genes$ensembl_gene_id = substr(DE_genes$ensembl_gene_id,1,nchar(DE_genes$ensembl_gene_id)-2) 
mart <- useEnsembl(biomart = "ensembl", dataset = "dnovaehollandiae_gene_ensembl", version = "102")

annotLookup_DE <-getBM(
  mart = mart,
  attributes = c(
    'ensembl_gene_id',
    'external_gene_name'),
  filter = 'ensembl_gene_id',
  values = DE_genes$ensembl_gene_id,
  uniqueRows = TRUE)

DE_genes<- full_join(DE_genes, annotLookup_DE, key="ensembl_gene_id")

write.csv(DE_genes, "/Users/sunshine2.0/Desktop/WORK/output/2nd_DE_genes.csv")

# _______________________________________________________END___________________________________________________________

