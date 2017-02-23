dat1 = read.csv(file = "stool-count.T.csv", header = TRUE)

names=dat1[,1]
rownames(dat1)=make.names(names,unique=TRUE)
dat1 <- dat1[,-1] #this deletes the first column 

