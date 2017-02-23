 countdata <- read.table("ahr.txt", sep="\t", header=T) 
 names=countdata[,1]

 rownames(countdata)=make.names(names,unique=TRUE)
  countdata <- countdata[,-1] #this deletes the first column 
 dim(countdata) 
[1] 58689     6
>  head(countdata) 
           SiCon_DMSO_1 SiCon_DMSO_2 SiCon_DMSO_3 SiAhr_DMSO_1 SiAhr_DMSO_2
X5S_rRNA              2            1            0            1            3
X5_8S_rRNA            0            0            0            0            0
X7SK                  0            0            0            0            0
A1BG                  0            5            1            2            3
A1BG.AS1             38           68           37           66           38
A1CF                  0            4            0            1            0
           SiAhr_DMSO_3
X5S_rRNA              2
X5_8S_rRNA            0
X7SK                  0
A1BG                  3
A1BG.AS1             49
A1CF                  0
 metadata <- read.table("ahrgroup.txt", sep="\t", header=T) 
  D <- DGEList( counts=countdata, group=metadata$group) 
> D
An object of class "DGEList"
$counts
           SiCon_DMSO_1 SiCon_DMSO_2 SiCon_DMSO_3 SiAhr_DMSO_1 SiAhr_DMSO_2
X5S_rRNA              2            1            0            1            3
X5_8S_rRNA            0            0            0            0            0
X7SK                  0            0            0            0            0
A1BG                  0            5            1            2            3
A1BG.AS1             38           68           37           66           38
           SiAhr_DMSO_3
X5S_rRNA              2
X5_8S_rRNA            0
X7SK                  0
A1BG                  3
A1BG.AS1             49
58684 more rows ...

$samples
             group lib.size norm.factors
SiCon_DMSO_1   ahr 24530584            1
SiCon_DMSO_2   ahr 44535047            1
SiCon_DMSO_3   ahr 14950708            1
SiAhr_DMSO_1  DMSO 45033625            1
SiAhr_DMSO_2  DMSO 21110994            1
SiAhr_DMSO_3  DMSO 20089613            1

> sum(rowSums(D$counts)>0) 
[1] 29627
> Dnorm <- calcNormFactors( D ) 
> Dnorm$samples 
             group lib.size norm.factors
SiCon_DMSO_1   ahr 24530584    0.9775597
SiCon_DMSO_2   ahr 44535047    0.9622241
SiCon_DMSO_3   ahr 14950708    0.9872173
SiAhr_DMSO_1  DMSO 45033625    1.0074778
SiAhr_DMSO_2  DMSO 21110994    1.0456790
SiAhr_DMSO_3  DMSO 20089613    1.0221953
>  lin.norm.factors <- (sum(D$samples$lib.size)/nrow(D$samples))/ D$samples$lib.size 
>  lin.norm.factors
[1] 1.1567232 0.6371408 1.8979098 0.6300869 1.3440909 1.4124262
>  Dc <- estimateCommonDisp( Dnorm ) 
> Dc$common.dispersion
[1] 0.008540497
> Dc.test <- exactTest( Dc ) 
> topTags(Dc.test)
Comparison of groups:  DMSO-ahr 
               logFC   logCPM       PValue          FDR
SCG2        2.381025 2.650193 4.394851e-50 2.579294e-45
TMEM255A   -1.691290 3.122411 2.954272e-33 8.669165e-29
CHGB        1.665798 3.186863 1.226693e-31 2.399779e-27
ATP6V1A     1.312144 6.936145 3.183383e-31 4.670739e-27
FOXQ1      -1.305898 4.769411 5.213729e-28 6.119771e-24
GDF15       1.576898 2.979769 6.648803e-28 6.503527e-24
HIST1H2BD   1.586562 2.664667 1.404858e-24 1.177853e-20
SLC7A10     1.509196 2.805245 1.627210e-24 1.193741e-20
CSGALNACT2  1.180856 5.566456 2.004863e-24 1.307371e-20
TMEM2       1.324758 3.506094 5.094798e-23 2.990086e-19
>  Dc.table <- topTags( Dc.test, n=nrow(Dc$counts) )$table 
> Dt <- estimateTagwiseDisp( Dnorm, prior.df=10, grid.length=500 ) 
Running estimateCommonDisp() on DGEList object before proceeding with estimateTagwiseDisp().
>  summary(Dt$tagwise.dispersion) 
     Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
0.0001334 0.0001616 0.0055990 0.0121600 0.0118900 0.5466000 
>  Dt.test <- exactTest( Dt) 
>  Dt.table <- topTags( Dt.test, n=nrow(Dt$counts) )$table
> write.table(Dt.table, file="top_table_edgeR_Plus_vs_Minus.txt", sep="\t") 
> Dt <- estimateTagwiseDisp( Dnorm, prior.df=13, grid.length=500 ) 
Running estimateCommonDisp() on DGEList object before proceeding with estimateTagwiseDisp().
>  summary(Dt$tagwise.dispersion) 
     Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
0.0001334 0.0001616 0.0062310 0.0118700 0.0122600 0.5466000 
>  Dt.test <- exactTest( Dt) 
> Dt.table <- topTags( Dt.test, n=nrow(Dt$counts) )$table 
>  write.table(Dt.table, file="top_table_edgeR_Plus_vs_Minus1.txt", sep="\t") 
> 
