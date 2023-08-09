# Glioma-Recurrence
Predicting the recurrence of cancerous gliomas using RNA-sequences

## Abstract: 
Recurrence of tumors in the brain can be predicted by analysing the RNA-seq of the tumors.

#### Dataset used: 
CPTAC-3 from the TCGA portal. Data is filtered to include only the cases with primary site as the Brain.

#### Steps planned:
1. Extract the raw counts for the tumors using GDC query.
2. Construct the count matrix and column data for differential gene expression with DESeq2.
3. Extract the up- and down-regulated genes.
4. The up- and down-regulated genes are considered to be in the positive (recurrent) class (as their behaviour increases the probability of recurrence. The rest of the genes are considered as the negative (non-recurrenct) class.
5. Use a ML or DL model to predict the possibility that the test case (containing the RNA seq counts of all the genes as in the training cases) will have a recurrent tumor.

#### Paper that is being referenced: 
Wang X, Han L, Zhou L, Wang L, Zhang LM. Prediction of candidate RNA signatures for recurrent ovarian cancer prognosis by the construction of an integrated competing endogenous RNA network. Oncol Rep. 2018 Nov;40(5):2659-2673. doi: 10.3892/or.2018.6707. Epub 2018 Sep 13. PMID: 30226545; PMCID: PMC6151886.
