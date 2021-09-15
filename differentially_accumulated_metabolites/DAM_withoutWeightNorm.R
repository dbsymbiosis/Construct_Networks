#!/usr/bin/env Rscript

DESCRIPTION="
##
## DESCRIPTION
##

R script to identify metabolites that are differentially accumulated between two conditions.

## Example:
Rscript DAM.R metabolites_file samplesInfo_file output_file treatment1 treatment2

# metabolites_file (comma delimited; first row needs to be the column names; sample IDs in first row should match those in the 'Sample.ID' columns in the samplesInfo_file):
groupId,goodPeakCount,medMz,medRt,maxQuality,compoundId,ATAC_2wk_1058,ATAC_2wk_1455,ATAC_2wk_1499
1,5,156.001175,1.182,0.738764,156.001175@1.181877,145133.81,89625.4,94933.16
2,3,133.014267,16.027,0.756112,133.014267@16.027122,57118.36,14854.14,22026.89
3,41,132.867844,6.946,0.718885,132.867844@6.945603,128486.34,182182.98,176913.56
4,3,201.037231,6.186,0.784493,201.037231@6.186069,0,0,0
..
..

# samplesInfo_file (tab delimited; 'Sample.ID', and 'Treatment' column names are required):
Sample.ID	Treatment
ATAC_12wk_1103	ATAC_12wk
ATAC_12wk_2306	ATAC_12wk
ATAC_12wk_1777	ATAC_12wk
ATAC_24hrs_1059	ATAC_24hrs
ATAC_24hrs_1757	ATAC_24hrs
ATAC_24hrs_1563	ATAC_24hrs
ATAC_2wk_1047	ATAC_2wk
..
..


# output_file
..
..

# treatment1 & treatment2
The two 'Treatment' IDs from sampleInfo_file to use for differential accumulation testing (e.g. ATAC_12wk HTHC_12wk)

"

args = commandArgs(trailingOnly=TRUE)
#args <- c("MC_Neg_ATAC_vs_HTAC.csv.filtered.csv", "Mcap_Sample_Info.csv", "DAM_Neg_T5_ATAC_vs_HTAC", "ATAC", "HTAC")

# Test if there is five arguments: if not, return an error
if (length(args)!=5) {
        cat(DESCRIPTION)
        stop("Five arguments must be supplied!", call.=FALSE)
}

# Parse command line arguments
metabolites=read.table(args[1],h=T,row.names=1,sep=",",check.names=F)
samplesInfo=read.table(args[2],h=T,sep="\t")
outputfile=args[3]
Treatment1=args[5]
Treatment2=args[6]

# Extract samples for each treatment
T1.samples=samplesInfo[which(samplesInfo$Treatment==Treatment1),]
T2.samples=samplesInfo[which(samplesInfo$Treatment==Treatment2),]

T1.metabolites=metabolites[,T1.samples$Sample.ID]
T2.metabolites=metabolites[,T2.samples$Sample.ID]

# Keep rows which have at least 1 sample that is > 0
t1=rownames(data.frame(which(rowSums(T1.metabolites>0)>0)))
t2=rownames(data.frame(which(rowSums(T2.metabolites>0)>0)))
tokeep=union(t1,t2)

T1.metabolites=T1.metabolites[tokeep,]
T2.metabolites=T2.metabolites[tokeep,]

# Compute t.test
pval=NULL

for(i in 1:length(tokeep))
{
  test=t.test(T1.metabolites[i,],T2.metabolites[i,])
  pval=c(pval,test$p.value)
}

pvalAdj=p.adjust(pval,method = "BH", n = length(pval))

# Compute fold change
T1.log2 = log2(T1.metabolites)
T1.log2[T1.log2=="-Inf"] <-0
control = apply(T1.log2,1,mean)

T2.log2 = log2(T2.metabolites)
T2.log2[T2.log2=="-Inf"] <-0
test = apply(T2.log2,1,mean)


foldchange = test - control


# Compute VIP score
suppressWarnings(suppressMessages(library(mixOmics)))

tmp=metabolites[tokeep,c(colnames(T1.metabolites),colnames(T2.metabolites))]

s=c(rep(Treatment1,length(colnames(T1.metabolites))),rep(Treatment2,length(colnames(T2.metabolites))))

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
