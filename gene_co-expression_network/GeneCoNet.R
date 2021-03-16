#!/usr/bin/env Rscript

DESCRIPTION="
##
## DESCRIPTION
##

Construct metabolite co-expression network.

## Example:
Rscript GeneCoNet.R TPM_matrix diffExprGenes timepoint

# TPM_matrix (tab delimited; TPM values for each gene across samples of interest):
Name	T1	T1	T2	T2	T3	T3
gene1	51.7	54.8	68.4	164.4	150.3	121.3
gene2	0.0	0.0	1.1	129.9	155.4	101.6
gene3	1000.3	1025.9	1236.8	2005.6	6038.7	3384.1
..
..

# diffExprGenes (list of differentially expressed genes):
gene1
gene2
gene3
..
..

# timepoint
Time point from TPM_matrix to extract (e.g. T2)


# Output files:
timepoint.corData.txt - network edges

timepoint.modules.txt - network modules

"

args = commandArgs(trailingOnly=TRUE)
#args <- c("DAM_Pos_T5_ATAC_vs_HTHC", "T5", "ATAC", "HTHC")

# Test if there is three arguments: if not, return an error
if (length(args)!=3) {
        cat(DESCRIPTION)
        stop("Three arguments must be supplied!", call.=FALSE)
}


## Main code
expression.file=args[1]
diffExprGenes=args[2]
timepoint=args[3]

#Load libraries
suppressMessages(suppressWarnings(library(DGCA, quietly = TRUE)))
suppressMessages(suppressWarnings(library(WGCNA, quietly = TRUE)))
suppressMessages(suppressWarnings(library(matrixStats, quietly = TRUE)))
#Parameters
maxPval=0.05
#Read data
tpm.data=as.matrix(read.table(expression.file,h=T,row.names=1,sep="\t"))
deg.data=read.table(diffExprGenes,h=T,row.names=1,sep="\t")
deg.ids=rownames(deg.data)

#Remove non numerical data and reformat the matrix
indx <- grepl(timepoint, colnames(tpm.data))
TP=which(indx==T)
tpm.data=tpm.data[deg.ids,TP]
C=colnames(tpm.data)
R=rownames(tpm.data)
dims=dim(tpm.data)
tpm.data=as.numeric(tpm.data)
dim(tpm.data)=dims
colnames(tpm.data)=C
rownames(tpm.data)=R

##Start computing co-expression networks##
t.data=t(tpm.data)
cor.data=matCorr(t.data,corrType="pearson")
pairs.cor.data=data.frame(n1=rownames(cor.data)[row(cor.data)],n2=colnames(cor.data)[col(cor.data)],cor=c(cor.data))
cor.data.row=rownames(cor.data)
cor.data.col=colnames(cor.data)
nsample.data=matNSamp(t.data)
cor.data.pval=matCorSig(cor.data,nsample.data)
colnames(cor.data.pval)=cor.data.col
rownames(cor.data.pval)=cor.data.row
pairsPval=data.frame(n1=rownames(cor.data.pval)[row(cor.data.pval)],n2=colnames(cor.data.pval)[col(cor.data.pval)],pval=c(cor.data.pval))
cor.data.pval.vec=as.vector(cor.data.pval)
cor.data.adjPval=adjustPVals(cor.data.pval.vec,adjust="BH")
cor.data.adjPval=as.numeric(format.pval(cor.data.adjPval,digits=2,nsmall=3))
dim(cor.data.adjPval)=dim(cor.data.pval)
colnames(cor.data.adjPval)=cor.data.col
rownames(cor.data.adjPval)=cor.data.row
pairsAdjPval=data.frame(n1=rownames(cor.data.adjPval)[row(cor.data.adjPval)],n2=colnames(cor.data.adjPval)[col(cor.data.adjPval)],adjPval=c(cor.data.adjPval))
#
cor.data.val=cbind(pairs.cor.data,pval=pairsPval$pval,adjPval=pairsAdjPval$adjPval)
cor.data.val.final=cor.data.val[complete.cases(cor.data.val),]
cor.data.val.final.filtered=cor.data.val.final[cor.data.val.final$adjPval <= maxPval,]
#
output=paste(timepoint,".corData.txt",sep = "")
write.table(cor.data.val.final.filtered,file=output,quote=FALSE,sep="\t",row.names=FALSE)
#
gene_tree<-hclust(as.dist(1-cor.data),method="average")
module_labels <- cutreeDynamicTree(dendro=gene_tree, minModuleSize=10,deepSplit=TRUE)
m <- as.data.frame(module_labels)
rownames(m)=rownames(cor.data)
output=paste(timepoint,".modules.txt",sep = "")
write.table(m,file=output,quote=FALSE,sep="\t",row.names=T,col.names=F)

#END
