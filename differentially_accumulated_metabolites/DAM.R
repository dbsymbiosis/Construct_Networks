#!/usr/bin/env Rscript

DESCRIPTION="
##
## DESCRIPTION
##
R script to identify metabolites that are differentially accumulated at a specific timepoint between two conditions.
## Example:
Rscript DAM.R metabolites_file samplesInfo_file output_file timepoint treatment1 treatment2
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
# output_file
..
..
# timepoint
Time point from sampleInfo_file to run test on (e.g. T5)
# treatment1 & treatment2
The two 'Treatment' IDs from sampleInfo_file to use for differential accumulation testing (e.g. ATAC ATHC)
"

args = commandArgs(trailingOnly=TRUE)
#args <- c("MC_Neg_ATAC_vs_HTAC.csv.filtered.csv", "Mcap_Sample_Info.csv", "DAM_Neg_T5_ATAC_vs_HTAC", "T5", "ATAC", "HTAC")

# Test if there is six arguments: if not, return an error
if (length(args)!=2) {
  cat(DESCRIPTION)
  stop("Six arguments must be supplied!", call.=FALSE)
}


metabolites=read.table(args[1],h=T,row.names=1,sep=",",check.names=F)
samplesInfo=read.table(args[2],h=T,row.names=1,sep=",")
outputfile=args[3]
timepoint=args[4]
Treatment1=args[5]
Treatment2=args[6]

A=samplesInfo[which(samplesInfo$Time==timepoint & samplesInfo$Treatment==Treatment1),]
H=samplesInfo[which(samplesInfo$Time==timepoint & samplesInfo$Treatment==Treatment2),]

A.norm=sweep(metabolites[,rownames(A)],2,as.numeric(A$Weight),FUN='/')
H.norm=sweep(metabolites[,rownames(H)],2,as.numeric(H$Weight),FUN='/')

a=rownames(data.frame(which(rowMeans(A.norm==0)<0.10)))
h=rownames(data.frame(which(rowMeans(H.norm==0)<0.10)))
tokeep=intersect(a,h)

A.norm=A.norm[tokeep,]
H.norm=H.norm[tokeep,]

# Compute t.test
pval=NULL

for(i in 1:length(tokeep))
{
  test=t.test(A.norm[i,],H.norm[i,])
  pval=c(pval,test$p.value)
}

pvalAdj=p.adjust(pval,method = "BH", n = length(pval))

# Compute fold change
A.log2 = log2(A.norm)
A.log2[A.log2=="-Inf"] <-0
control = apply(A.log2,1,mean)

H.log2 = log2(H.norm)
H.log2[H.log2=="-Inf"] <-0
test = apply(H.log2,1,mean)


foldchange = test - control


# Compute VIP score
suppressWarnings(suppressMessages(library(mixOmics)))

tmp=metabolites[tokeep,c(colnames(A.norm),colnames(H.norm))]

s=c(rep(Treatment1,length(colnames(A.norm))),rep(Treatment2,length(colnames(H.norm))))

suppressWarnings(suppressMessages(p<-plsda(t(tmp),s)))

plotIndiv(p,ind.names=s,ellipse=T,legend=T)

vp=vip(p)

# Create table with all data
stats=cbind("FC"= foldchange,"VIP"=vp[,1],"PVALUE"=pval,"ADJPVALUE"=pvalAdj)
tmp=cbind(stats,tmp)

# Filter the data and keep only significant metabolites -> FoldChange >= 2 and VIP score >= 1
sig=which(abs(stats[,1]) >= 2 & stats[,2] >=1)

# Get siginificant metabolites data
final=metabolites[sig,1:5]

# Save table of significant metabolites or DAM
final=cbind(final,tmp[sig,])
write.table(data.frame("molID"=rownames(final),final,check.names=F),file=outputfile,quote=F,sep="\t",row.names=F)

#END