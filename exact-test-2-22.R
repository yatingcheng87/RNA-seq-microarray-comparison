args = c("stool","stool.csv")

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

y = DGEList(counts=counts, group=factors)
y = calcNormFactors(y)
y = estimateCommonDisp(y)
y = estimateTagwiseDisp(y)


scaled.counts = data.frame(mapply(`*`, counts, y$samples$lib.size *
                                    y$samples$norm.factors/mean(y$samples$lib.size)))
rownames(scaled.counts) = rownames(counts)
dfs = split.data.frame(t(scaled.counts), groups)
dfss = sapply(dfs, colMeans)
lrt = exactTest(y, pair = c(1,3))
ot1 = topTags(lrt,n=nrow(counts),sort.by="PValue")$table



ot1 = merge(ot1, dfss, by=0)
ot1 = ot1[order(ot1$FDR),] # Sort by ascending FDR
write.csv(ot1,"ot1.csv",row.names=FALSE)
