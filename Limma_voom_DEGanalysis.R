library(edgeR)
library(limma)
tbl = read.table('Rat_Brain_15M_6M_6W_count_GeneID.txt', header = TRUE, sep = '\t', row.names = 1)
tbl_rpkm = read.table('Rat_Brain_15M_6M_6W_rpkm_GeneID.txt', header = TRUE, sep = '\t', row.names = 1)
tbl_male = tbl[,1:16]
tbl_female = tbl[,17:36]
tbl_rpkm_male = tbl_rpkm[,1:16]
tbl_rpkm_female =tbl_rpkm[,17:36]

select = apply(tbl_rpkm_female, 1, function(x) sum(x >= 1) >= 10 )
tbl_female.clean = tbl_female[select,]
tbl_rpkm_female.clean = tbl_rpkm_female[select,]

Age.group = factor(c(rep('6M', times = 10), rep('6W', times = 10)), levels = c('6W','6M'))
y = DGEList(tbl_female.clean, group = Age.group)
y = calcNormFactors(y)
design = model.matrix(~Age.group)
vwts = voomWithQualityWeights(y, design = design)
vfit2 = lmFit(vwts)
vfit2 = eBayes(vfit2)
write.table(topTable(vfit2,coef = 2 , n = Inf),'6Month_vs_6Week.fastq.trimmed.RAT_ens98_longest_cDNA.BWA.SAMBestCount.VoomWeight.limma.txt', sep = '\t',quote = FALSE)
DEG.total = topTable(vfit2 , n = Inf)
DEG.clean = DEG.total[abs(DEG.total$logFC) > 1 & DEG.total$adj.P.Val < 0.05, ]
write.table(DEG.clean,'6Month_vs_6Week.fastq.trimmed.RAT_ens98_longest_cDNA.BWA.SAMBestCount.VoomWeight.limma_clean.txt', sep = '\t',quote = FALSE)

tbl = read.table('SH-SY5Y_IronTreat_24h_MGISeq_count_GeneID.txt', header = TRUE, sep = '\t', row.names = 1)
tbl_rpkm = read.table('SH-SY5Y_IronTreat_24h_MGISeq_rpkm_GeneID.txt', header = TRUE, sep = '\t', row.names = 1)
select = apply(tbl_rpkm[,c(1,2,3,10,11,12)], 1, function(x) sum(x >= 1) >= 3)
tbl_rpkm.clean = tbl_rpkm[select,]
tbl.clean = tbl[select,]
tbl.clean = tbl.clean[,c(1,2,3,10,11,12)]
tbl_rpkm.clean = tbl_rpkm.clean[,c(1,2,3,10,11,12)] 
#pair.select = apply(tbl_rpkm.clean, 1, function(x) sum(x >= 1) >= 3)
tbl.ready = tbl.clean
Conc.group = factor(c('10mM','10mM','10mM','Ctrl','Ctrl','Ctrl'), levels = c('Ctrl','10mM'))
y = DGEList(tbl.ready ,group = Conc.group)
y = calcNormFactors(y)
design = model.matrix(~Conc.group)
vwts = voomWithQualityWeights(y, design = design)
vfit2 = lmFit(vwts)
vfit2 = eBayes(vfit2)
write.table(topTable(vfit2,coef = 2 , n = Inf),'SH-SY5Y_IronTreat.MGISeq.VoomWeight.limma.txt', sep = '\t',quote = FALSE)
DEG.total = topTable(vfit2, coef =2, n = Inf)
DEG.clean = DEG.total[abs(DEG.total$logFC) > 1 & DEG.total$adj.P.Val < 0.05, ]
write.table(DEG.clean,'SH-SY5Y_IronTreat.MGISeq.VoomWeight.limma_clean.txt', sep = '\t',quote = FALSE)
