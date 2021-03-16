#!/usr/bin/env Rscript

DESCRIPTION="
##
## DESCRIPTION
##

Construct metabolite co-expression network.

## Example:
Rscript GeneCoNet.R expression_matrix samplesInfo_file timepoint treatment1 treatment2

# expression_matrix (tab delimited; header line required, sample IDs should match the 'Sample.ID' columns in samplesInfo_file):
groupId,goodPeakCount,medMz,medRt,maxQuality,compoundId,ATAC_2wk_1058,ATAC_2wk_1455,ATAC_2wk_1499
1,5,156.001175,1.182,0.738764,156.001175@1.181877,145133.81,89625.4,94933.16
2,3,133.014267,16.027,0.756112,133.014267@16.027122,57118.36,14854.14,22026.89
3,41,132.867844,6.946,0.718885,132.867844@6.945603,128486.34,182182.98,176913.56
4,3,201.037231,6.186,0.784493,201.037231@6.186069,0,0,0
..
..

# samplesInfo_file (comma delimited; 'Sample.ID', 'Weight', 'Time' and 'Treatment' headers required):
Sample.ID,Time,TP,Weight,Species,Treatment
ATAC_12wk_1103,12wk,T11,140,Pacuta,ATAC
ATAC_12wk_2306,12wk,T11,150,Pacuta,ATAC
ATAC_12wk_1777,12wk,T11,130,Pacuta,ATAC
ATAC_24hrs_1059,24hrs,T5,240,Pacuta,ATHC
ATAC_24hrs_1757,24hrs,T5,190,Pacuta,ATHC
ATAC_24hrs_1563,24hrs,T5,300,Pacuta,ATHC
ATAC_2wk_1047,2wk,T7,190,Pacuta,HTHC
..
..

# timepoint
Time point from sampleInfo_file to run test on (e.g. T5)

# treatment1 & treatment2
The two 'Treatment' IDs from sampleInfo_file to use for differential accumulation testing (e.g. ATAC ATHC)


# Output files:
metabolites_file.cor.txt - network edges

metabolites_file.modules.txt - network modules

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
