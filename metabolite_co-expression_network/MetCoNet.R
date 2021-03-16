#!/usr/bin/env Rscript

DESCRIPTION="
##
## DESCRIPTION
##

Construct metabolite co-expression network.

## Example:
Rscript MetCoNet.R metabolites_file samplesInfo_file timepoint treatment1 treatment2

# metabolites_file (comma delimited; header line required, sample IDs should match the 'Sample.ID' columns in samplesInfo_file):
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
if (length(args)!=4) {
        cat(DESCRIPTION)
        stop("Four arguments must be supplied!", call.=FALSE)
}


## Main code
input=args[1]
samplefile=args[2]
timepoint=args[3]
Treatment1=args[4]
Treatment2=args[5]

#Load libraries
suppressWarnings(suppressMessages(library(DGCA, quietly = TRUE)))
suppressWarnings(suppressMessages(library(WGCNA, quietly = TRUE)))
suppressWarnings(suppressMessages(library(matrixStats, quietly = TRUE)))
#Parameters
maxPval=0.05
#Read data
metabolites=read.table(input,h=T,row.names=1,sep="\t",check.names=F)
samplesInfo=read.table(samplefile,h=T,row.names=1,sep=",")

A=samplesInfo[which(samplesInfo$TP==timepoint & samplesInfo$Treatment==Treatment1),]
H=samplesInfo[which(samplesInfo$TP==timepoint & samplesInfo$Treatment==Treatment2),]
A.norm=sweep(metabolites[,rownames(A)],2,as.numeric(A$Weight),FUN='/')
H.norm=sweep(metabolites[,rownames(H)],2,as.numeric(H$Weight),FUN='/')
data=cbind(A.norm,H.norm)
dataA=t(data)
corAMB=matCorr(dataA,corrType="pearson")
pairscorAMB=data.frame(n1=rownames(corAMB)[row(corAMB)],n2=colnames(corAMB)[col(corAMB)],cor=c(corAMB))
corAMBrow=rownames(corAMB)
corAMBcol=colnames(corAMB)
nsampleA=matNSamp(dataA)
corAMBpval=matCorSig(corAMB,nsampleA)
colnames(corAMBpval)=corAMBcol
rownames(corAMBpval)=corAMBrow
pairsPval=data.frame(n1=rownames(corAMBpval)[row(corAMBpval)],n2=colnames(corAMBpval)[col(corAMBpval)],pval=c(corAMBpval))
corAMBpvalVec=as.vector(corAMBpval)
corAMBAdjPval=adjustPVals(corAMBpvalVec,adjust="BH")
corAMBAdjPval=as.numeric(format.pval(corAMBAdjPval,digits=2,nsmall=3))
dim(corAMBAdjPval)=dim(corAMBpval)
colnames(corAMBAdjPval)=corAMBcol
rownames(corAMBAdjPval)=corAMBrow
pairsAdjPval=data.frame(n1=rownames(corAMBAdjPval)[row(corAMBAdjPval)],n2=colnames(corAMBAdjPval)[col(corAMBAdjPval)],adjPval=c(corAMBAdjPval))
corAMBval=cbind(pairscorAMB,pval=pairsPval$pval,adjPval=pairsAdjPval$adjPval)
corAMBvalFinal=corAMBval[complete.cases(corAMBval),]
corAMBvalFinalFiltered=corAMBvalFinal[corAMBvalFinal$adjPval <= maxPval,]
output=paste(input,".cor.txt",sep = "")
write.table(corAMBvalFinalFiltered,file=output,quote=FALSE,sep="\t",row.names=FALSE)
tree<-hclust(as.dist(1-corAMB),method="average")
module_labels <- cutreeDynamicTree(dendro=tree, minModuleSize=10,deepSplit=TRUE)
m <- as.data.frame(module_labels)
rownames(m)=rownames(corAMB)
output=paste(input,".modules.txt",sep = "")
write.table(m,file=output,quote=FALSE,sep="\t",row.names=T,col.names=F)
