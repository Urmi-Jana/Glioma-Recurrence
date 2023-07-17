library(DESeq2)
# library(TCGAbiolinks)
# BiocManager::install("TCGAbiolinks")
# library(tidyverse)

co <- read.csv("rec.csv", header=TRUE) #read in csv files

# 1. select the stranded_first column 2. filter out 'lncRNA'
cr <- co[, 1:5 ] 
df1 <- co[co$gene_type=='lncRNA' , c(1, 2, 3, 5)]
df1<- df1[, c(1, 4)]

# same for non-recurrence
non_rec <- read.csv("non-rec.csv", header=TRUE)
df2 <- non_rec[non_rec$gene_type=='lncRNA',]
df2 <- df2[, c(1,4) ]

# merge values to make a count matrix
df <- merge(df1, df2, by='gene_id')

rec2 <- read.csv("rec2.csv", header=TRUE)
df3 <- rec2[rec2$gene_type=='lncRNA',]
df3 <- df3[, c(1,4) ]


df <- merge(df, df3, by="gene_id")

# change rownames and remove redundant column
rownames(df) <- df$gene_id
df <- df[, -c(1)]



# change column names
colnames(df) <- c( "1", "2", "3")
# print( head(df))

# metadata

mr <- data.frame(
  sample_name = c("1", "2", "3"),
  status = c("rec", "non-rec", "rec")
)

rownames(mr) <- mr$sample_name
colnames(df) <- mr$sample_name
print(head(mr))

# DESeq2 object
dds <- DESeqDataSetFromMatrix(
  countData=df,
  colData=mr,
  design=~status)


# filter to remove genes with 0 counts
dds <- dds[rowSums(counts(dds)) >= 5]

# factor levels to tell deseq about factor to compare against
dds$status = relevel(dds$status, ref='non-rec')

# run seq function

dds <- DESeq(dds)

res <- results(dds,  pAdjustMethod = "none", alpha = "0.5")

print(res)
print(summary(res))