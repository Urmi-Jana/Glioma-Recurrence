library(DESeq2)
library(TCGAbiolinks)

clinical_data <- GDCquery_clinic("CPTAC-3", type = "clinical")
data <- clinical_data[clinical_data$primary_site=='Brain',]
# print(data[1:10, 1:5])
which(colnames(data) %in% c("progression_or_recurrence"))

data2 <- data[, c(2,79)]
data2 <- data2[rowSums(is.na(data2)) == 0,]
data2 <- data2[!data2$progression_or_recurrence=='not reported',]

query <- GDCquery(
    project = "CPTAC-3",
    data.category = "Transcriptome Profiling",
    experimental.strategy = "RNA-Seq",
    workflow.type = "STAR - Counts",
    data.type = "Gene Expression Quantification",
    access = "open",
)

output_query <- (getResults(query))
all(data2$submitter_id %in% output_query$cases.submitter_id)

gliomas_query_data<- output_query[output_query$cases.submitter_id %in% data2$submitter_id,]

samples <- gliomas_query_data$cases

gliomas_query <- GDCquery(
    project = "CPTAC-3",
    data.category = "Transcriptome Profiling",
    experimental.strategy = "RNA-Seq",
    workflow.type = "STAR - Counts",
    data.type = "Gene Expression Quantification",
    access = "open",
    barcode = samples,
)

# GDCdownload(gliomas_query)
gliomas_data <- GDCprepare(gliomas_query, summarizedExperiment = TRUE)

# count matrix
gliomas_data_assay <- assay(gliomas_data, 'unstranded')

gene_metadata <- as.data.frame(rowData(gliomas_data))
coldata <- as.data.frame(colData(gliomas_data))
# coldata <- coldata$progression_or_recurrence

# write.table(coldata, file="coldata.csv", sep = ",", row.names=FALSE)

all(rownames(coldata) %in% colnames(gliomas_data_assay))

dds <- DESeqDataSetFromMatrix(
  countData=gliomas_data_assay,
  colData=coldata,
  design=~progression_or_recurrence)

# remove genes with less than 10 reads across all samples
keep <- rowSums(counts(dds) >= 5)
dds <- dds[keep,]

# set factor level
dds$progression_or_recurrence <- relevel(dds$progression_or_recurrence, ref = 'no')


# run deseq function
dds <- DESeq(dds)
res <- results(dds, pAdjustMethod = "none", alpha = "0.3")

summary(res)

head(res)
# ordered_res <- res[order(res$padj),]
# head(ordered_res)

results <- list(coldata = coldata, gliomas_data_assay = gliomas_data_assay, 
        gene_metadata = gene_metadata, res = res, dds = dds)
saveRDS (file = 'results.RDS', results)

data <- readRDS('results.RDS')
print(data$dds)
