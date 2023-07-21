library(DESeq2)
library(TCGAbiolinks)

clinical_data <- GDCquery_clinic("CPTAC-3")
data <- clinical_data[clinical_data$primary_site=='Brain',]
# print(data[1:10, 1:5])
any(colnames(data) %in% c("progression_or_recurrence"))
# print(colnames(data))
empty_columns <- colSums(is.na(data)) == 0
data2 <- data[, !empty_columns]

# data2 <- data2[, c('submitter_id','progression_or_recurrence')]
# data2 <- data2[rowSums(is.na(data2)) == 0,]

query <- GDCquery(
    project = "CPTAC-3",
    data.category = "Transcriptome Profiling",
    experimental.strategy = "RNA-Seq",
    workflow.type = "STAR - Counts",
    data.type = "Gene Expression Quantification",
    access = "open",
)

output_query <- (getResults(query))

data <- GDCprepare(gliomas_query, summarizedExperiment = TRUE)
# output <- apply(output_query,2,as.character)
# # print(output)
# print(length(output))

# samples <- output_query$cases[1:20]

gliomas_query <- GDCquery(
    project = "CPTAC-3",
    data.category = "Transcriptome Profiling",
    experimental.strategy = "RNA-Seq",
    workflow.type = "STAR - Counts",
    data.type = "Gene Expression Quantification",
    access = "open",
    barcode = samples,
)

output <- getResults(gliomas_query)
length(samples)

# GDCdownload(gliomas_query)
gliomas_data <- GDCprepare(gliomas_query, summarizedExperiment = TRUE)

df <- gliomas_data[gliomas_data$primary_diagnosis=='Glioblastoma']

# count matrix
gliomas_data_assay <- assay(gliomas_data, 'unstranded')
# write.table(gliomas_data_assay, file="assay.csv", sep = ",", row.names=FALSE)

all("C3L-02544" %in% data$submitter_id)
