# prepare matrices for Data Fusion

# RNA data

# normalise RNA-seq matrix function
# make RNA function to generate normalised count matrix
library(DESeq2)

RNAdeseq <- function(data.table.RNA = "data/RNAbyGENESKCM_RSEM.txt",
                     results_file = "results/vst_normalized_counts.txt", label.samples = cbind(c(rep("x", 40), rep("y",40)))){
  SKCM.RNAseq <- read.table(data.table.RNA,
                           sep = "\t", header =T, row.names = 1)

  corrected.names <- sapply(1:length(colnames(SKCM.RNAseq)), function(x){ gsub("\\.", "-",
                                                                              colnames(SKCM.RNAseq)[x])})
  colnames(SKCM.RNAseq) <- corrected.names
  
  if(length(which(SKCM.RNAseq =="NA", arr.ind = T)) != 0){
    cat(paste("RNA Matrix has the following NA values, and they will be converted to zero to avoid errors:\n"))
    print(which(SKCM.RNAseq=="NA", arr.ind = T))
    cts[which(SKCM.RNAseq=="NA", arr.ind = T)] <- 0
  }
  
  
  cts <- matrix(data= round(as.numeric(unlist(SKCM.RNAseq)), digits = 0), ncol=ncol(SKCM.RNAseq), nrow=nrow(SKCM.RNAseq))
  colnames(cts) <- colnames(SKCM.RNAseq)
  rownames(cts) <- rownames(SKCM.RNAseq)
  
  # start deseq part 
  coldata <- label.samples
  rownames(coldata) <- colnames(SKCM.RNAseq)
  coldata <- data.frame(coldata)
  colnames(coldata) <- c("condition")
  
  
  # remove seq with  less than 100 reads
  cts <- cts[which(rowSums(cts) > 100), ]
  
  # remove seq with more than 1000000 reads
  cts <- cts[which(rowSums(cts) < 1000000), ]
  
  # remove genes with more than 20% of zeros
  numb.zero <- sapply(1:nrow(cts), function(x){length(which(cts[x, ] == 0))})
  
  cts <- cts[which(numb.zero < round(ncol(cts)*0.2, digits=0)), ]
  
  
  dds <- DESeqDataSetFromMatrix(countData = cts,
                                colData = coldata,
                                design = ~ condition)
  
  
  vsd <- vst(dds, blind=T)
  
  dat_norm <- assay(vsd)
  # https://hbctraining.github.io/DGE_workshop_salmon_online/lessons/02_DGE_count_normalization.html
  
  write.table(dat_norm, file= results_file, sep="\t", quote=F, col.names=T, row.names = T)
}


# take rnaseq data
# it can be downloaded from http://firebrowse.org/?cohort=SKCM&download_dialog=true
# and take RSEM counts to produce the RNAbyGENESKCM_SELECTED_RSEM.txt
# file RNAbyGENESKCM_SELECTED_RSEM.txt should have this format
#  TCGA-BF-A1PU-01	TCGA-BF-A1PV-01
# A1BG 2222 3333




library("RTCGAToolbox")


# this is supplementary table 1D, it can be downloaded from  https://www.cell.com/cell/fulltext/S0092-8674(15)00634-0?_returnURL=https%3A%2F%2Flinkinghub.elsevier.com%2Fretrieve%2Fpii%2FS0092867415006340%3Fshowall%3Dtrue
# check Table S1 xlx file (table D, Patient-centric table) 

clinical.data <- read.table("data/S1D_table_SKCM_clinical_Molecular_data.txt", header = T, sep = "\t")

SKCMData <- getFirehoseData(dataset="SKCM", runDate="20160128",
                            forceDownload=TRUE, RNASeq2Gene = T)
RNA.data <- biocExtract(SKCMData, "RNASeq2Gene")
RNA.tab <- getData(RNA.data, type = "assays")@"data"[[1]]

# write.table(RNA.tab, "data/RNAbyGENESKCM_RSEM.txt", sep = "\t", quote = F)

name.RNA <- colnames(RNA.tab)
# correct RNA names
name.RNA <- colnames(RNA.tab)
name.RNA1 <- sapply(1:length(name.RNA), function(x){paste(strsplit(name.RNA[x], "-")[[1]][1:3], collapse = "-")})
name.RNA2 <- sapply(1:length(name.RNA), function(x){paste(strsplit(name.RNA[x], "-")[[1]][4], collapse = "-")})
name.RNA2 <- sapply(1:length(name.RNA), function(x){paste(strsplit(name.RNA2[x], "*")[[1]][1:2], collapse = "")})
name.RNA <-  sapply(1:length(name.RNA), function(x){paste(name.RNA1[x], name.RNA2[x], sep = "-")})
availableRNA <- match(clinical.data[,1], name.RNA)
availableRNA <- availableRNA[which(is.na(availableRNA)==F)]
name.RNA <- name.RNA[availableRNA[which(is.na(availableRNA)==F)]]
RNA.tab <- RNA.tab[, availableRNA]
colnames(RNA.tab) <- name.RNA

write.table(RNA.tab, "data/RNAbyGENESKCM_SELECTED_RSEM.txt", sep = "\t", quote = F)

# file RNAbyGENESKCM_SELECTED_RSEM.txt should have this format
#  TCGA-BF-A1PU-01	TCGA-BF-A1PV-01
# A1BG 2222 3333




RNAdeseq(data.table.RNA = "data/RNAbyGENESKCM_SELECTED_RSEM.txt", results_file = "results/vsd_normalised_SKCM.txt",
         label.samples = cbind(c(rep("x", 166), rep("y", 166))))



############################


# add methylation


# download methylation file from http://firebrowse.org/?cohort=SKCM&download_dialog=true

# cat SKCM.methylation__humanmethylation450__jhu_usc_edu__Level_3__within_bioassay_data_set_function__data.data.txt  | awk -F '\t' '$4!="X"{print $0}' > filter1_methyl.txt
# cat filter1_methyl.txt  | awk -F '\t' '$4!="Y"{print $0}' > filter2_methyl.txt

#BiocManager::install("minfi")
#BiocManager::install("IlluminaHumanMethylation450kanno.ilmn12.hg19")

#install.packages("data.table")

# https://www.bioconductor.org/help/course-materials/2015/BioC2015/methylation450k.html#genetic-variants-and-cell-type-composition
# use minfi to read file and remove SNP probes

library("minfi")
library("data.table")
library("IlluminaHumanMethylation450kanno.ilmn12.hg19")
methylation.tab <- readTCGA("data/filter2_methyl.txt", sep = "\t", keyName = "Composite Element REF", Betaname = "Beta_value",
                            pData = NULL, array = "IlluminaHumanMethylation450k",
                            mergeManifest = FALSE,
                            showProgress = TRUE)

methylation.tab <- addSnpInfo(methylation.tab)
# remove all probes mapping on SNPs
methylation.tab <- dropLociWithSnps(methylation.tab, snps=c("Probe","SBE","CpG"), maf=0)

# produce matrix to be analysed with consensus cluster Plus
beta.val <- methylation.tab@assays@data$Beta

#add probe name
rownames(beta.val) <- methylation.tab@rowRanges@ranges@NAMES

# detect and remove rows with NA
NA.rows <- sapply(1:nrow(beta.val), function(x){any(is.na(beta.val[x, ]))})

beta.val <- beta.val[which(NA.rows==F), ]


# remove not used samples
clinical.data <- read.table("data/S1D_table_SKCM_clinical_Molecular_data.txt", header = T, sep = "\t")
pat.names1 <- sapply(1:length(colnames(beta.val)), function(x){paste(strsplit(colnames(beta.val)[x], "-")[[1]][c(1,2,3)], collapse = "-")})
pat.names2 <- sapply(1:length(colnames(beta.val)), function(x){gsub("([A-Z]).*$", "", strsplit(colnames(beta.val)[x], "-")[[1]][4])})
colnames(beta.val) <-  paste(pat.names1, pat.names2, sep = "-")
meth.table <- beta.val[, match(clinical.data[, 1], colnames(beta.val))]
all(colnames(meth.table)==clinical.data[,1])
# take 1% most variable probes
ngen <- round(0.01*nrow(meth.table), digits=0)
mads <- apply(meth.table, 1, mad)
var.meth <- meth.table[rev(order(mads))[1:ngen],]


rm(methylation.tab)

write.table(var.meth, "results/Methylation_SKCM_TCGA_01_part.txt", row.names = T, col.names = T, quote = F, sep = "\t")


######################################

# Code to prepare data for supplementary figure 1



# download mutation data
library("TCGAbiolinks")


query.maf <- GDCquery(project = "TCGA-SKCM", 
                           data.category = "Simple Nucleotide Variation", 
                           data.type = "Masked Somatic Mutation",
                           workflow.type = "Aliquot Ensemble Somatic Variant Merging and Masking",
                           access = "open",
                           legacy = FALSE)


GDCdownload(query.maf)


# then all downloaded files have to be copy and pasted in one directory
# the package should have downloaded the files inside the GDCdata/TCGA-SKCM folder  
# then we want to extract all files and merge in one
# gunzip -c * > SKCM.maf
# we remove all lines with "#"
# grep -v "#" SKCM.maf  > SKCM_1.maf
# we want to take maf files header to later remove the header of all maf files from the file
# head -8 SKCM.maf  > head_SKCM.txt
# grep -v "Hugo_Symbol" SKCM_1.maf  > SKCM_2.maf
# take only used columns
# awk 'BEGIN { FS = "\t" } ; {print $1,$2,$3, $4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$16,$51,$76,$77,$122}' SKCM_2.maf > SKCM_3.maf
# tail -1 head_SKCM.txt | awk '{print $1,$2,$3, $4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$16,$51,$76,$77,$122}' > last_line.txt
# grep "#" head_SKCM.txt > forst_part.vcf
# cat forst_part.vcf last_line.txt SKCM_3.maf > reduced_SKCM_hg38.maf
# cat reduced_SKCM_hg38.maf  | tr " " "\t" > SKCM_hg38.maf

