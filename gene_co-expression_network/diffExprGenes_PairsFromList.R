library(DESeq2)

args = commandArgs(trailingOnly=TRUE)
# test if there is at least one argument: if not, return an error
if (length(args)!=3) {
  stop("Three arguments must be supplied (matrix, sample files, pairs2run).n", call.=FALSE)
}

matrix.file <- args[1]
sample.file <- args[2]
pairs.to.run <- args[3]
# matrix.file <- "test_data.matrix"
# sample.file <- "test_samples.txt"
# pairs.to.run <- "test_data.pairs2run"

## Load matrix file
data = read.table(matrix.file, header=T, com='', row.names=1, check.names=F, sep='\t')
data = as.matrix(data)
data <- round(data)
mode(data) <- "integer"

## Load samples file [sample \t replicate]
samples_data = read.table(sample.file, header=F, check.names=F, fill=T)
samples_data = samples_data[samples_data[,2] != '',]
colnames(samples_data) = c('sample_name', 'replicate_name')
sample_types = as.character(unique(samples_data[,1]))
rep_names = as.character(samples_data[,2])

## Load pairs to run file [sample_ID1 \t sample_ID2]
pairs.to.run_data = read.table(pairs.to.run, header=F, check.names=F, fill=T, stringsAsFactors = FALSE)
colnames(pairs.to.run_data) = c('ID1', 'ID2')

## Get just the columns 
data = data[, colnames(data) %in% rep_names, drop=F ] # Get just data columns in samples file
nsamples = length(sample_types)

## Get list of sample<->replicate {S1:[S1R1, S1R2, S1R3], S2:[.....]}
sample_type_list = list()
for (i in 1:nsamples) {
  samples_want = samples_data[samples_data[,1]==sample_types[i], 2]
  sample_type_list[[sample_types[i]]] = as.vector(samples_want)
}

## Get list of sample order [S1, S1, S1, S2, S2, S2, S3, ....]
sample_factoring = colnames(data)
for (i in 1:nsamples) {
  sample_type = sample_types[i]
  replicates_want = sample_type_list[[sample_type]]
  sample_factoring[ colnames(data) %in% replicates_want ] = sample_type
}


## Check if all samples are in *.matrix file.
if (all(rep_names %in% colnames(data))) {
  if (length(rep_names) == length(colnames(data))) {
    print("All columns in *.matrix file are in samples.txt")
  } else {
    warning("WARNING: **NOT** all columns in *.matrix file are in samples.txt")
  }
} else {
  stop("ERROR: Not all samples in samples.txt are in the *.matrix file")
}

## Check if all pairs2run are in *.matrix file.
if (all(unique(c(as.vector(pairs.to.run_data$ID1), as.vector(pairs.to.run_data$ID2))) %in% unique(samples_data$sample_name))) {
  if (length(unique(c(as.vector(pairs.to.run_data$ID1), as.vector(pairs.to.run_data$ID2)))) == length(unique(samples_data$sample_name))) {
    print("All samples in samples.txt are in pairs2run.txt")
  } else {
    warning("WARNING: **NOT** all columns in samples.txt are in pairs2run.txt")
  }
} else {
  stop("ERROR: Not all samples in pairs2run.txt are in the samples.txt")
}


## For each pairwise comparison. 
matrix.file.results <- paste(matrix.file,"_DiffExprResults.txt",sep='')
write(paste(c("seqName", "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj", "controlSample", "treatmentSample"), collapse='\t'), matrix.file.results, append = FALSE)
diff.expressed.genes <- c()
for (i in 1:nrow(pairs.to.run_data)) {
  print (paste("Paires to compare", pairs.to.run_data[i,1], pairs.to.run_data[i,2]))
  
  # Get columns which match each sample we want to test
  d1 <- data[,sample_factoring %in% pairs.to.run_data[i,1]]
  d2 <- data[,sample_factoring %in% pairs.to.run_data[i,2]]
  
  ## 
  ## DESeq2 analysis
  ## 
  
  # Object describing conditions of each column - control (c, c, c) or treatment (t, t, t)
  conditions <- c(rep("c", length(colnames(d1))), rep("t", length(colnames(d2))))
  coldata <- data.frame(row.names=c(colnames(d1), colnames(d2)), conditions)
  
  # "dds" object for DESeq2 - contains counts and condition information
  dds <- DESeqDataSetFromMatrix(countData=cbind(d1, d2), colData=coldata, design=~conditions)
  # Run DESeq2 - returns normalized expression level in logarithmic base 2 scale
  results <- results(DESeq(dds))
  
  # Add names of samples being compared to table.
  results$controlSample <- pairs.to.run_data[i,1]
  results$treatmentSample <- pairs.to.run_data[i,2]
  
  # Filter results - fold change is log2FC >1 with adjusted p-values<0.05
  res <- na.exclude(as.data.frame(results))
  filter <- res[(abs(res$log2FoldChange)>1 & res$padj<0.05),]
  
  print (paste("Number of Diff Exppressed genes:", length(rownames(filter))))
  write.table(filter, matrix.file.results, append = TRUE, quote = FALSE, col.names = FALSE, sep = '\t')
  diff.expressed.genes <- c(diff.expressed.genes, rownames(filter))
}

# Write names of genes that are differentially expressed between any two paire-wise comparisons. 
write(unique(diff.expressed.genes), paste(matrix.file,"_DiffExprGeneNames.txt",sep=''))

