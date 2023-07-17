library(DESeq2)
library(TCGAbiolinks)

clinical_data <- GDCquery_clinic("CPTAC-3")
data <- clinical_data[clinical_data$primary_site=='Brain',]
# print(data[1:10, 1:5])
any(colnames(data) %in% c("progression_or_recurrence"))
# print(colnames(data))
empty_columns <- colSums(is.na(data)) == 0
data2 <- data[, !empty_columns]

data2 <- data2[, c('submitter_id','progression_or_recurrence')]
data2 <- data2[rowSums(is.na(data2)) == 0,]

query <- GDCquery(
    project = "CPTAC-3",
    data.category = "Transcriptome Profiling",
    experimental.strategy = "RNA-Seq",
    workflow.type = "STAR - Counts",
    data.type = "Gene Expression Quantification",
    access = "open",
    # barcode = data2[data2$submitter_id],
)


output_query <- (getResults(query))
# output <- apply(output_query,2,as.character)
# # print(output)
# print(length(output))

samples <- output_query$cases[1:300]

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
# write.table(output, file="res.csv", sep = ",", row.names=FALSE)
length(output)
