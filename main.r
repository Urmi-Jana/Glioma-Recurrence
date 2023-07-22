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

GDCdownload(gliomas_query)
gliomas_data <- GDCprepare(query, summarizedExperiment = TRUE)

data2$cases <- samples
data2 <- data2[, c(3, 2)]
colnames(df) <- c('','id', 'recurrence')
# count matrix
gliomas_data_assay <- assay(gliomas_data, 'unstranded')
write.table(gliomas_data_assay, file="counts.csv", sep = ",", row.names=FALSE)


