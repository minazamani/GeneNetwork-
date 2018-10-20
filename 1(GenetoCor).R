#autor Mina Zamani...First Partial of My thesis Code about Collective bahavior of Gene Expression in Cancerous and Normal cells
# last Update 12 oct 2018 

#libraries
library(Hmisc)



#Load data from https://github.com/mskcc/RNAseqDB/blob/master/data/normalized/brca-rsem-fpkm-tcga-t.txt.gz
#the normalized gene expression levels (FPKM). This set of data files was not only quantile normalized, but also was corrected for batch effects (using tool ComBat).
#https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5903355/

#a <-  read.delim("brca-rsem-fpkm-tcga-t.txt.gz")
#b <-  read.delim("brca-rsem-fpkm-tcga.txt.gz")
a <-  read.delim("kirc-rsem-fpkm-tcga-t.txt.gz")
b <-  read.delim("kirc-rsem-fpkm-tcga.txt.gz")

rownames(a) <- a$Hugo_Symbol
rownames(b) <- b$Hugo_Symbol

Cancera <- as.matrix(a[,-1:-2])
Normalb <- as.matrix(b[,-1:-2])

#Extract Significant Genes using Statistical Tests (Pvalue and Fold Change)

m = as.matrix(cbind(Cancera[, 1:dim(Cancera)[2]], Normalb[, 1:dim(Normalb)[2]]))# Combine Tow set of Data

p.value.all.genes = as.matrix(apply(m, 1, function(x) { t.test(x[1:dim(Cancera)[2]], x[dim(Cancera)[2]+1:(dim(Cancera)[2]+dim(Normalb)[2]+1)]) $p.value } ))
SignificantGenesTtest <-subset(p.value.all.genes , p.value.all.genes[,1] <0.01)
RowGeneTtest <- as.matrix(rownames(SignificantGenesTtest))
SigRowTtest <-(match(RowGeneTtest , rownames(b)))

Cancera <- as.matrix(Cancera[SigRowTtest, ])
Normalb <- as.matrix(Normalb[SigRowTtest, ])

Cancermedian <- as.matrix(apply(Cancera, 1, median))
Noramalmedian <-as.matrix(apply(Normalb, 1, median))
fold.changes <- Cancermedian/Noramalmedian #*****************************************
fold.changes <- na.omit(fold.changes)
#Or use this function >> foldChanges(y, labels, matrixALL, matrixNoL, middle=mean.middle, log=TRUE)


SignificantGenesFC <-subset(fold.changes , fold.changes[,1] > 4)
SignificantGenesFC <- as.matrix(SignificantGenesFC[!is.infinite(SignificantGenesFC),])

x <- as.matrix(rownames(SignificantGenesFC))
SigRow <-  as.matrix(match(x, rownames(Normalb)))
Cancera <- as.matrix((Cancera[SigRow, ]))
Normalb <- as.matrix((Normalb[SigRow, ]))
Cancera <- na.omit(Cancera)
Normalb <- na.omit(Normalb)


#Calculate Correlation

CorCancer <- rcorr(as.matrix(t(Cancera)))
CorCancer <- CorCancer$r
CorNormal <- rcorr(as.matrix(t(Normalb)))
CorNormal <- as.matrix(CorNormal$r) 

write.table(CorCancer ,"CorCancerKirc.txt",sep = "\t")
write.table(CorNormal ,"CorNormalKirc.txt",sep = "\t")


