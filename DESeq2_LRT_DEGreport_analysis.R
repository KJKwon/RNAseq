library(DESeq2)
library(DEGreport)
tbl = read.table('Rat_Brain_15M_6M_6W_count_Left+Right_GeneID.txt', header = TRUE, sep = '\t', row.names = 1)
#coldata <- data.frame(Age = factor(c(rep('Control', times = 3), rep('1mM', times = 3), 
                                     #rep('2mM', times= 3),rep('3mM', times= 3))))
#rownames(coldata) = colnames(tbl)
#coldata$Age = factor(coldata$Age, levels = c('Control','1mM','2mM','3mM'))
coldata <- data.frame(Age = factor(c(rep('M15',times = 4), rep('M6',times = 9),rep('W6', times = 5))))
rownames(coldata) = colnames(tbl)
#coldata$Age = factor(coldata$Age, levels = c('M15','M6','W6'))
dds <- DESeqDataSetFromMatrix(countData = tbl, colData = coldata, design = ~Age)
sigLRT_rlog <- rlog(dds)
sigLRT_rlog_mat <- assay(sigLRT_rlog)
dds <- DESeq(dds, test = "LRT", reduced = ~1)
res <- results(dds, name = 'Age_M6_vs_M15')
res.clean <- res[!is.na(res$padj),]
sigLRT_genes <- rownames(res.clean[res.clean$padj < 0.05 & abs(res.clean$log2FoldChange) >= 1,])
cluster_rlog_mat <- sigLRT_rlog_mat[rownames(sigLRT_rlog_mat) %in% sigLRT_genes,]
clusters <- degPatterns(cluster_rlog_mat, metadata = coldata, time = "Age", col = NULL)
tmp.candi <- clusters$df$genes[clusters$df$cluster == 1]
#dds <- DESeq(dds)
resLFC.M15vsM6 <- lfcShrink(dds, coef = "Age_3mM_vs_Control", type = "apeglm")
write.table(resLFC.M15vsM6, 'SH-SY5Y_IronChallenge_3mMvsCtrl_NextSeq_DESeq2_lfcShrink_apeglm_GeneName_output.txt',
            quote = FALSE, sep = '\t')
dds$Age = relevel(dds$Age, ref = 'W6')
dds <- nbinomWaldTest(dds)
resLFC.M6vsW6 <- lfcShrink(dds, coef = "Age_M6_vs_W6", type = "apeglm")
write.table(resLFC.M6vsW6, '6Month_vs_6Week.fastq.trimmed.RAT_ens98_longest_cDNA.BWA.SAMBestCount.DESeq2_lfcShrink_apeglm_GeneID.txt',
            quote = FALSE, sep = '\t')

tbl.total = merge(as.data.frame(resLFC.M15vsM6), as.data.frame(resLFC.M6vsW6), by = 0)
tbl.total = tbl.total[,-c(2,4,5,7,9,10)]
colnames(tbl.total) = c('GeneID','lfc.15Mvs6M','FDR.15Mvs6M','lfc.6Mvs6W','FDR.6Mvs6W')
tbl.total = tbl.total[!is.na(tbl.total$FDR.15Mvs6M) & !is.na(tbl.total$FDR.6Mvs6W),]
tbl.total.clean = tbl.total[tbl.total$lfc.15Mvs6M >= 1 & tbl.total$FDR.15Mvs6M <0.05,]
write.table(tmp.candi,'Rat_Brain_15M_6M_6W_count_fastq.trimmed.RAT_ens98_longest_cDNA.BWA.DEseq2_LRT_DEGreport_cluster_15moUp_Left+Right_geneID.txt',
             ,quote = FALSE, row.names = FALSE)
