library(edgeR)
tbl_count = read.table('15Mo34_vs_6Mo56.fastq.trimmed.RAT_ens98_dna_rm.STAR.ReadsPerGene.csv',sep = '\t', header = TRUE,
                       row.names = 1)
tmp_cor = dist(1-cor(as.matrix(cpm(tbl_clean)), method = 'spearman'))
tmp_clust <- hclust(tmp_cor, method="average")
tmp_clust = as.dendrogram(tmp_clust)
plot(tmp_clust, main = 'Rat_15Mo34_vs_Rat_6Mo56 CPM [Spearman]', cex = 1.5, cex.main = 1.5)

keep <- rowSums(cpm(tbl_count) > 1) >= 4
tbl_clean_count = tbl_count[keep,]
tbl_clean_cpm = cpm(tbl_count)[keep,]
groups = c('15Month','15Month','15Month','15Month','6Month','6Month','6Month','6Month')
y_tn <- DGEList(counts=tbl_clean_count, group=groups)
y_tn <- calcNormFactors(y_tn)
Age = factor(c('15Mo','15Mo','15Mo','15Mo','6Mo','6Mo','6Mo','6Mo'))
y_tn = estimateCommonDisp(y_tn)
y_tn = estimateTagwiseDisp(y_tn)
#design = model.matrix(~Age)
#rownames(design) = colnames(y_tn)
#y_tn <- estimateDisp(y_tn, design, robust = TRUE)
#fit = glmFit(y_tn,design)
#6Month = Up, 15Month = down
#lrt = glmLRT(fit, coef = 2, contrast = c("6Mo","15Mo"))
write.table(topTags(exactTest(y_tn, pair=c('6Month','15Month')), n = Inf),'15Mo34_vs_6Mo56.fastq.trimmed.RAT_ens98_dna_rm.STAR.ReadsPerGene.edgeR.txt', sep = '\t',
            quote = FALSE)
