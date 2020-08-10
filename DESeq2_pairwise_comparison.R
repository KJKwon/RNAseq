library(DESeq2)
tbl.count = read.table('Rat_Brain_15M_6M_6W_count_GeneName.txt', header = TRUE, row.names = 1)
tbl.rpkm = read.table('Rat_Brain_15M_6M_6W_rpkm_GeneName.txt', header = TRUE, row.names = 1)
tbl.count = tbl.count[,-c(1:16)]
tbl.rpkm = tbl.rpkm[,-c(1:16)]
coldata = data.frame(Age = rep(c('6M','6W'), each = 10))
rownames(coldata) = colnames(tbl.count)
#coldata$Age = factor(coldata$Age, levels = c('Control', '2mM', '10mM'))
coldata$Age = relevel(coldata$Age, ref = '6M')
dds = DESeqDataSetFromMatrix(countData = tbl.count, colData = coldata, design = ~Age)
#Take same standard as edgeR DEG analysis
keep = apply(tbl.rpkm, 1, function(x) sum(x >= 1) >= 10 )
dds = dds[keep,]
dds = DESeq(dds)
res = results(dds, contrast = c("Age","Control","2mM"))
reslfc = lfcShrink(dds, coef = "Age_2mM_vs_10mM", type ='apeglm')
write.table(as.data.frame(reslfc), file = "SH-SY5Y_IronTreat_2mMvs10mM.3pairs.DESeq2lfcShrink_output.txt",quote = FALSE, sep ='\t')
vsd = vst(dds, blind = FALSE)
norm.count = counts(dds, normalized = TRUE)
write.table(assay(vsd), file = "SH-SY5Y_IronTreat.3pairs.DESeq2normalized_GSEA_input_output.txt", quote = FALSE, sep ='\t')
tmp_cor = dist(1-cor(tbl.rpkm, method = 'spearman'))
tmp_clust = hclust(tmp_cor, method = 'average')
tmp_plot = as.dendrogram(tmp_clust)
plot(tmp_plot)
