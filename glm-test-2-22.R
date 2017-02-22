
# Get counts file from analysis/fname/fname.T.csv
bname = basename(args[1]) 
fname = paste(bname,"-count.T.csv",sep='')
dat = read.csv(file.path("analysis",bname,fname), header = TRUE, row.names=1)
keep = rowSums(cpm(dat) > 3) >= 3
counts = dat[keep,]

sel = grepl("MT-.*", rownames(counts)) + grepl("ERCC-.*", rownames(counts)) + grepl("mt-.*", rownames(counts))

counts = counts[!sel,]


head(counts)
factors <- c("B", "A-tissue", "A", "B-tissue", "A", "B")
groups <-factors



y = DGEList(counts=counts, group = groups)
design <- model.matrix(~ 0 + y$samples$group, data =y$samples)
y = calcNormFactors(y)
y = estimateGLMCommonDisp(y, design)
y = estimateGLMTrendedDisp(y, design)
y = estimateGLMTagwiseDisp(y, design)
fit = glmFit(y, design)

scaled.counts = data.frame(mapply(`*`, counts, y$samples$lib.size *
                                    y$samples$norm.factors/mean(y$samples$lib.size)))
rownames(scaled.counts) = rownames(counts)
dfs = split.data.frame(t(scaled.counts), groups)
dfss = sapply(dfs, colMeans)
scaled.counts = data.frame(mapply(`*`, counts, y$samples$lib.size *
                                    y$samples$norm.factors/mean(y$samples$lib.size)))
rownames(scaled.counts) = rownames(counts)
dfs = split.data.frame(t(scaled.counts), groups)
dfss = sapply(dfs, colMeans)
lrt = glmLRT(fit, contrast=c(-1,0,1,0)) 
topTags(lrt, n =10)


ot1 = topTags(lrt,n=nrow(counts),sort.by="PValue")$table



ot1 = merge(ot1, dfss, by=0)
ot1 = ot1[order(ot1$FDR),] # Sort by ascending FDR
write.csv(ot1,"ot2.csv",row.names=FALSE)