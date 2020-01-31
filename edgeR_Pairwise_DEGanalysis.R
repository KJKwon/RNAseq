library(edgeR)
tbl_count = read.table('15Mo34_vs_6Mo56.fastq.trimmed.RAT_ens98_dna_rm.STAR.ReadsPerGene.csv',sep = '\t', header = TRUE,row.names = 1)
tmp_cor = dist(1-cor(as.matrix(tbl_clean_cpm), method = 'spearman'))
tmp_clust <- hclust(tmp_cor, method="average")
tmp_clust = as.dendrogram(tmp_clust)
plot(tmp_clust, main = 'Rat_15Mo34_vs_Rat_6Mo56 CPM [Spearman]', cex = 1.5, cex.main = 1.5)
groups = c('15Month','15Month','15Month','15Month','6Month','6Month','6Month','6Month')
y_tn <- DGEList(counts=tbl_count, group=groups)
keep = filterByExpr(y_tn)
y_tn = y_tn[keep, , keep.lib.sizes = FALSE]
y_tn <- calcNormFactors(y_tn)
Age = factor(c('15Mo','15Mo','15Mo','15Mo','6Mo','6Mo','6Mo','6Mo'))
y_tn = estimateCommonDisp(y_tn)
y_tn = estimateTagwiseDisp(y_tn)
write.table(topTags(exactTest(y_tn, pair=c('6Month','15Month')), n = Inf),'15Mo34_vs_6Mo56.fastq.trimmed.RAT_ens98_dna_rm.STAR.ReadsPerGene.edgeR.txt', sep = '\t',
            quote = FALSE)
write.table(cpm(y_tn$counts), '15Mo34Left_vs_6Mo56.fastq.trimmed.RAT_ens98_dna_rm.STAR.ReadsPerGene.newcpm.txt', row.names = TRUE, col.names = TRUE, sep ='\t',
            quote = FALSE)
