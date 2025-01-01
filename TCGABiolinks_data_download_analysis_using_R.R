
install.packages("BiocManager")
BiocManager::install("TCGAbiolinks")


library(TCGAbiolinks)
library(tidyverse)

install.packages("BiocManager")
BiocManager::install("maftools")

library(maftools)
library(pheatmap)
library(SummarizedExperiment)

# get list of projects
gdcprojects <- getGDCprojects()
getProjectSummary('TCGA-BRCA')

# building a query
query_TCGA <- GDCquery(project = 'TCGA-BRCA', data.category = 'Transcriptome Profiling')
output_query_TCGA <- getResults(query_TCGA)


# build a query to retrieve gene expression data ------------
query_TCGA <- GDCquery(project = 'TCGA-BRCA',
                       data.category = 'Transcriptome Profiling',
                       experimental.strategy = 'RNA-Seq',
                       workflow.type = 'STAR - Counts',
                       access = 'open',
                       barcode = c('TCGA-5L-AAT0-01A-12R-A41B-07', 'TCGA-A2-A04U-01A-11R-A115-07','TCGA-AN-A04A-01A-21R-A034-07'))

getResults(query_TCGA)

# download data - GDCdownload
GDCdownload(query_TCGA)

# prepare data
tcga_brca_data <- GDCprepare(query_TCGA, summarizedExperiment = TRUE)
brca_matrix <- assay(tcga_brca_data, 'fpkm_unstrand')


# build a query to retrieve DNA methylation data --------------

getProjectSummary('TCGA-LUAD')

query_LUAD <- GDCquery(project = 'TCGA-LUAD', data.category = 'DNA Methylation')
output_query_LUAD <- getResults(query_LUAD)

query_LUAD_methyl <- GDCquery(project = 'TCGA-LUAD',
                         data.category = 'DNA Methylation',
                         platform = 'Illumina Human Methylation 27',
                         access = 'open',
                         data.type = 'Methylation Beta Value',
                         barcode = c('TCGA-38-4630-01A-01D-1205-05', 'TCGA-55-1594-11A-01D-0945-05', 'TCGA-17-Z004-01A-01D-0752-05'))

output_query_methyl <- getResults(query_LUAD_methyl )

# download data - GDCdownload
GDCdownload(query_LUAD_methyl)


# plot probes showing differences in beta values between samples

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
BiocManager::install("sesame")

library(sesame)

 ## prepare data
tcga_luad_meth_data <- GDCprepare(query_LUAD_methyl, summarizedExperiment = TRUE)
assay(tcga_luad_meth_data)

yes

# removing rows with NA values in all 3 samples
na.omit(luad_matrix)
  
# top 10 plot probes showing differences in beta values (methylation values) between samples                                            
analyse_methy <- tcga_luad_meth_data %>% 
  assay %>% 
  rowVars() %>% 
  order(decreasing = TRUE) %>% 
  head(10)

# heatmap showing the visualization of top 10 probes having high variance amongst samples
pheatmap(assay(tcga_luad_meth_data)[analyse_methy,])



# download and visualize mutation data from TCGA ----------------------
query_mutation <- GDCquery(project = 'TCGA-BRCA',
                           data.category = 'Simple Nucleotide Variation',
                           access = 'open',
                           barcode = c('TCGA-AR-A24U-01A-11D-A167-09,TCGA-AR-A24U-10A-01D-A167-09',
                                       'TCGA-AQ-A04J-01A-02W-A050-09,TCGA-AQ-A04J-10A-01W-A055-09',
                                       'TCGA-BH-A0HK-01A-11W-A071-09,TCGA-BH-A0HK-10A-01W-A071-09'))

output_query_mutation <- getResults(query_mutation)

GDCdownload(query_mutation)

maf <- GDCprepare(query_mutation, summarizedExperiment = TRUE)

# maftools utils to read and create dashboard for analyzing mutations
maftools.input <- read.maf(maf)

plotmafSummary(maf = maftools.input,
               addStat = 'median',
               rmOutlier = TRUE,
               dashboard = TRUE)


# oncoprint
oncoplot(maf = maftools.input,
         top = 10,
         removeNonMutated = TRUE)
