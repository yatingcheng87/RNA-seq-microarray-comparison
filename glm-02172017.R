dat <-read.csv(file ="stool-count.T.csv", header=TRUE, row.names=1)
dat <-read.csv(file ="stool-count.T.csv", header=TRUE)
rownames(dat) <-dat[,1]
dat<-dat[,-1]
keep = rowSums(cpm(dat) > 3) >= 3
counts = dat[keep,]
key = read.csv(file ="stool.csv", header=TRUE, row.names=1)
sel = grepl("MT-.*", rownames(counts)) + grepl("ERCC-.*", rownames(counts)) + grepl("mt-.*", rownames(counts))
counts = counts[!sel,]

d = DGEList(counts=counts, group = groups)

design.mat <- model.matrix(~ 0 + d$samples$group, data =d$samples)
colnames(design.mat) <- levels(d$samples$group)
colnames(design.mat) <- levels(d$samples$group)

d = calcNormFactors(d)
d = estimateGLMCommonDisp(d, design.mat)
d = estimateGLMTrendedDisp(d, design.mat)
d = estimateGLMTagwiseDisp(d, design.mat)
fit = glmFit(d, design.mat)
lrt12 <- glmLRT(fit, contrast=c(-1,0,1,0))
topTags(lrt12, n=10)
d