# make DF on SKCM dataset

# filtering functions

# get patient annotation 
patients.table <- read.table("./data/S1D_table_SKCM_clinical_Molecular_data.txt", header = T, sep = "\t")

# add RNA data and take 1500 most variable genes

RNA.data <- read.table("./results/vsd_normalised_SKCM.txt", header = T, sep = "\t", row.names = 1)
colnames(RNA.data) <- gsub("\\.", "-", colnames(RNA.data))

mads <- apply(RNA.data,1,mad)
RNA.red <- RNA.data[rev(order(mads))[1:1500],]




# take id of patients without RNA, miRNA, CNA
noRNA <- patients.table[which(is.na(match(patients.table[,1], colnames(RNA.red)))==T), 1]

patients.list <- setdiff(patients.table[,1], unique(noRNA))


meth.red <- read.table("./results/Methylation_SKCM_TCGA_01_part.txt", sep = "\t", header = T)

colnames(meth.red) <- gsub("\\.","-", colnames(meth.red))

meth.red <- meth.red[, match(patients.list, colnames(meth.red))]

RNA.red <- RNA.red[, match(patients.list, colnames(RNA.red))]

patients.table <- patients.table[match(patients.list, patients.table[, 1]), ]

if(F){
  # check all labels match
  all(colnames(RNA.red)==colnames(meth.red))
  all(colnames(RNA.red)==patients.table[,1])
}

# make matrices to be used later for heatmaps
write.table(RNA.red, "results/RNA_reduced_SKCM.txt", col.names = T, row.names = T, sep = "\t", quote = F)
write.table(meth.red, "results/Meth_reduced_SKCM.txt", col.names = T, row.names = T, sep = "\t", quote = F)



if(T){
  # write table for jSVD
  dir.create("results/dataSKCMjsvd_2D")
  dir.create("results/dataSKCMjsvd_2D/preproc")
  write.table(t(RNA.red), "results/dataSKCMjsvd_2D/preproc/RNA_prep.txt", col.names = F, row.names = F, sep = "\t")
  write.table(t(meth.red), "results/dataSKCMjsvd_2D/preproc/Met_prep.txt", col.names = F, row.names = F, sep = "\t")
  write.table(patients.table, "results/dataSKCMjsvd_2D/patients_SKCM_clinical.txt", sep = "\t", quote = F)
}


##################
# table 1 ########
##################


# execute analysis on all domains with Spectrum and silhuette plot

Matrix.score <- matrix(data=NA, ncol=10, nrow = 7)


rownames(Matrix.score) <- c("RNA-seq", "Methylation", "Spectrum", "SNF", "NEMO", "jSVD", "MUT_CLASS")

# clinical/molecular/cluster class to check
rev.class <- c(3:11)

colnames(Matrix.score) <- c(colnames(patients.table)[rev.class], "Silhouette")

# compute corrected rand index
library("pdfCluster")

# silhouette score
library("cluster")


library(Spectrum)
# RNA-seq
Spectrum.clus <- Spectrum(data = data.frame(RNA.red), fontsize=8, dotsize=2)
# relevant patient classes
Matrix.score[1,(1:length(rev.class)) ] <- sapply((1:length(rev.class)), function(j){round(adj.rand.index(patients.table[,rev.class[j]], Spectrum.clus$assignments), digits = 3)})
sil.score <- summary(silhouette(Spectrum.clus$assignments, dist(Spectrum.clus$similarity_matrix)), FUN = mean)
Matrix.score[1, ncol(Matrix.score)] <- paste(round(sil.score$clus.avg.widths, digits=3), collapse = ",")


# Methylation
Spectrum.clus <- Spectrum(data = meth.red, fontsize=8, dotsize=2)
Matrix.score[2,(1:length(rev.class)) ] <- sapply((1:length(rev.class)), function(j){round(adj.rand.index(patients.table[,rev.class[j]], Spectrum.clus$assignments), digits = 3)})
sil.score <- summary(silhouette(Spectrum.clus$assignments, dist(Spectrum.clus$similarity_matrix)), FUN = mean)
Matrix.score[2, ncol(Matrix.score)] <- paste(round(sil.score$clus.avg.widths, digits=3), collapse = ",")


# make DF 
SKCM.data <- vector(mode = "list", length = 2)
SKCM.data[[1]] <- data.frame(RNA.red)
SKCM.data[[2]] <- data.frame(meth.red)



Spectrum.clus <- Spectrum(data = SKCM.data, fontsize=8, dotsize=2)
Matrix.score[3,(1:length(rev.class)) ] <- sapply((1:length(rev.class)), function(j){round(adj.rand.index(patients.table[,rev.class[j]], Spectrum.clus$assignments), digits = 3)})
sil.score <- summary(silhouette(Spectrum.clus$assignments, dist(Spectrum.clus$similarity_matrix)), FUN = mean)
Matrix.score[3, ncol(Matrix.score)] <- paste(round(sil.score$clus.avg.widths, digits=3), collapse = ",")




# NEMO

library(SNFtool) # this package requires heatmap.plus that is no more available in recent versions of R
library(NEMO)

# compute SNF SNFtools
## First, set all the parameters:
K = 20;		# number of neighbors, usually (10~30)
alpha = 0.5;  	# hyperparameter, usually (0.3~0.8)
T = 20; 	# Number of Iterations, usually (10~20)

Dist1 = dist2(as.matrix(t(SKCM.data[[1]])),as.matrix(t(SKCM.data[[1]])))
Dist2 = dist2(as.matrix(t(SKCM.data[[2]])),as.matrix(t(SKCM.data[[2]])))

W1 = affinityMatrix(Dist1, K, alpha)
W2 = affinityMatrix(Dist2, K, alpha)

W = SNF(list(W1,W2), K, T)


C = 3 
SNF.class = spectralClustering(W, C); 	# the final subtypes information
Matrix.score[4,(1:length(rev.class)) ] <- sapply((1:length(rev.class)), function(j){round(adj.rand.index(patients.table[,rev.class[j]], SNF.class), digits = 3)})
sil.score <- summary(silhouette(SNF.class, dist(W)), FUN = mean)
Matrix.score[4, ncol(Matrix.score)] <- paste(round(sil.score$clus.avg.widths, digits=3), collapse = ",")




nemo.clus <- nemo.clustering(SKCM.data)
Matrix.score[5,(1:length(rev.class)) ] <- sapply((1:length(rev.class)), function(j){round(adj.rand.index(patients.table[,rev.class[j]], nemo.clus), digits = 3)})
affinity.graph <- nemo.affinity.graph(SKCM.data, k = length(unique(nemo.clus)))
sil.score <- summary(silhouette(nemo.clus, dist(affinity.graph)), FUN = mean)
Matrix.score[5, ncol(Matrix.score)] <- paste(round(sil.score$clus.avg.widths, digits=3), collapse = ",")




# JSVD

library("ConsensusClusterPlus")
system("python3 JSVD_Pfeffer_et_al_2019.py 4")
JSVD.data <- read.table("results/U_JSVD_Pfeffer_et_al_2019.txt", header = F, sep = "\t")


dmeth <- matrix(data = unlist(JSVD.data), nrow=nrow(JSVD.data), ncol = ncol(JSVD.data))
resultsCCP <- ConsensusClusterPlus(t(dmeth),maxK=4,reps=10000,pItem=0.8,pFeature=1,title="example3",
                                   distance="euclidean",clusterAlg="km", plot = "png")
# check you have same classes
class.jsvd <- resultsCCP[[3]]$consensusClass

# make table of classes



write.table(cbind(patients.table[,1], class.jsvd), "./results/jsvd_classes.txt", sep = "\t", quote = F, col.names = F, row.names = F)

Matrix.score[6,(1:length(rev.class)) ] <- sapply((1:length(rev.class)), function(j){round(adj.rand.index(patients.table[, rev.class[j]], class.jsvd), digits = 3)})
sil.score <- summary(silhouette(class.jsvd, dist(JSVD.data)), FUN = mean)
Matrix.score[6, ncol(Matrix.score)] <- paste(round(sil.score$clus.avg.widths, digits=3), collapse = ",")

Matrix.score[7, ] <- c(sapply((1:length(rev.class)), function(j){round(adj.rand.index(patients.table[, rev.class[j]], patients.table[, 3]), digits = 3)}), NA)

write.table(Matrix.score[, c(1,5,6,ncol(Matrix.score))], file="results/SKCM_assessment_matrix_jSVD_SKCM_2domain.txt", sep ="\t", quote = F)

library(ComplexHeatmap)

MUT <- patients.table$MUTATIONSUBTYPES
RNA <- patients.table$RNASEQ.CLUSTER_CONSENHIER
Met <- patients.table$MethTypes.201408

feat.mat <- cbind(MUT, RNA, Met, class.jsvd)

feat.mat <- feat.mat[order(-class.jsvd), ]

rownames(feat.mat) <- c(1:nrow(feat.mat))

feat.df <- cbind(feat.mat[,1], feat.mat[,2], feat.mat[,3], feat.mat[,4])

feat.df <- data.frame(feat.df)

colnames(feat.df) <- c("MUT", "RNA", "Met", "class.jSVD")

# make legend heatmap
ha <- HeatmapAnnotation(df=feat.df, gp = gpar(col = "white"),show_annotation_name = F, simple_anno_size = unit(0.3, "cm"),
                        annotation_label = c("Mutation", "RNA", "Methylation", "jSVD"),show_legend = F,
                        col = list(MUT = c("-" = "white", "Triple_WT" = "green", "RAS_Hotspot_Mutants" = "red", "NF1_Any_Mutants" = "blue",
                                           "BRAF_Hotspot_Mutants" = "purple"),
                                   RNA = c("-" = "white","immune" = "purple", "keratin" = "grey", "MITF-low" = "brown"),
                                   Met = c("CpG island-methylated" = "orange", "hyper-methylated" = "yellow", "normal-like" = "grey", "hypo-methylated" = "blue"),
                                   class.jSVD = c("1" = "red", "2" = "blue","3" = "orange", "4" = "brown")))
# legend for draw
lgd1 = Legend(at = c("-", "Triple_WT", "RAS_Hotspot_Mutants", "NF1_Any_Mutants", "BRAF_Hotspot_Mutants"),  legend_gp = gpar(fill=c("white", "green", "red", "blue", "purple")), title = "Mutation")
lgd2 = Legend(at = c("immune", "keratin", "MITF-low"),  legend_gp = gpar(fill=c("purple", "grey", "brown")), title = "RNA")
lgd3 = Legend(at = c("CpG island-methylated", "hyper-methylated", "normal-like", "hypo-methylated"),  legend_gp = gpar(fill=c("orange", "yellow", "grey", "blue")), title = "Methylation")
lgd4 = Legend(at = c("1", "2", "3"),  legend_gp = gpar(fill=c("red", "blue", "orange")), title = "jSVD class")

pd = packLegend(lgd1, lgd2, lgd3, lgd4, max_height = unit(10, "cm"), 
                column_gap = unit(1, "cm"))


##################
# figure 1 #######
##################

# make heatmap function
# make survival analysis
library("survminer")
library("survival")
library("ggplot2")
library("ComplexHeatmap")
library("samr")


MakeGeneHeatmap <- function(gene.matrix = "data/UM_vst_normalized_counts.txt", patient.class = class.jsvd,
                            results.signif.genes = "./results/sign.genes.txt",
                            results.filtered.matrix = "./results/filt_matrix.txt",
                            number.genes = NA,
                            fold.change.par = 2,
                            make.heatmap = T,
                            results.heatmap = "./results/heatmap.png",
                            Methyl.plot = F,
                            first.row = 2,
                            return.plot.heat = F,
                            annotation.add = T,
                            heat.legend = ha,
                            col.rampSel = c(-2, 0, 2)){
  TGCA_dat <- strsplit(scan(file = gene.matrix, what = character(), sep = "\n", quiet = T) , "\t")
  
  # remove first 2 rows and first column
  # NB sapply automatically transpose matrix (it was made for symmetric matrices)
  G_mat <- sapply(first.row:length(TGCA_dat), function(x) TGCA_dat[[x]][2:length(TGCA_dat[[first.row]])])
  row.nam <- sapply(first.row:length(TGCA_dat), function(x) {TGCA_dat[[x]][1]})
  G_mat <- t(G_mat)
  num_col <- ncol(G_mat)
  num_row <- nrow(G_mat)
  
  # remove NA values that would cause problems later
  if(length(which(G_mat =="NA", arr.ind = T)) != 0){
    cat(paste("Matrix", gene.matrix, "has the following NA values, and they will be converted to zero to avoid errors:\n"))
    print(which(G_mat=="NA", arr.ind = T))
    G_mat[which(G_mat=="NA", arr.ind = T)] <- 0
  }
  
  G_mat <- matrix(as.numeric(G_mat), nrow = num_row, ncol = num_col)
  
  
  
  rownames(G_mat) <- row.nam
  
  # convert in correct format patient.class vector
  #patient.class <- convertBinaryClass(patient.class)
  # censoring status 1 censored 2 dead 
  #patient.class <- as.numeric(patient.class)+1
  
  data2 <- list(x = G_mat, y = patient.class, logged2 = TRUE, genenames = row.names(G_mat), geneid = row.names(G_mat))
  
  ##Calculate all significant genes:
  if(length(unique(patient.class))==2){
    samr.obj <- samr(data2, resp.type = "Two class unpaired", nperms = 100)
  } else {
    samr.obj <- samr(data2, resp.type = "Multiclass", nperms = 100)
    
  }
  
  # min.foldchange -> consider genes with at least 2 fold change
  delta.table <- samr.compute.delta.table(samr.obj, min.foldchange = fold.change.par , nvals = 50)
  
  # find first value of delta with median FDR = 0
  delta.value <- delta.table[which(delta.table[, 5] == 0)[1], 1]
  if(is.na(delta.value)){
    delta.value <- 0
  }   
  
  
  if(is.na(number.genes) !=T ){
    # find delta which considers around 300 genes 
    delta.Value <- delta.table[which(delta.table[, 4] <= number.genes)[1], 1]
    if(is.na(delta.Value)){ # if there is no delta according to the selected threshold
      cat(paste0("no delta values for that threshold of genes, using delta = ", delta.value, "instead\n"))
      cat("More Info:\n")
      cat(delta.table[which(delta.table[, 5] == 0)[1], ])
    } else {
      # find highest delta.value
      delta.value <- max(delta.value, delta.Value)
    }
  }
  
  cat(paste(delta.table[which(delta.table[, 1] == delta.value)[1], 4], "genes selected\n"))
  
  
  siggenes.table <- samr.compute.siggenes.table(samr.obj, del = delta.value, data2, delta.table, all.genes = FALSE)
  #Note: the "min.foldchange=0.1" means that the fold change for the two groups should be >0.1.
  
  a <- siggenes.table$genes.up; # all up regulated genes
  b <- siggenes.table$genes.lo; # all down regulated genes
  c_vect <- rbind(a,b)
  if(is.na(results.signif.genes)!=T){
    write.table(c_vect, results.signif.genes, sep = "\t")
  }
  sign.genes <- c_vect[, 2] # take significant genes
  # take rows with significant genes
  sign.row <- sapply(1:length(sign.genes), function(x) {match(sign.genes[x], rownames(G_mat))})
  sign_mat <- G_mat[sign.row, ]
  # order sign_mat according to patient class and transpose
  colnames(sign_mat) <- patient.class
  sign_mat <- sign_mat[, order(-as.numeric(patient.class))]
  if(is.na(results.filtered.matrix)!=T){
    write.table(sign_mat, results.filtered.matrix, sep = "\t")
  }
  #print(sign_mat)
  # call function that makes heatmap
  if(make.heatmap == T){
    if(Methyl.plot){
      my_palette <- c("blue","white","yellow") #colorRampPalette(c("blue","white","yellow"))(n = 50)  
    } else {
      my_palette <- c("blue","white","red") # colorRampPalette(c("blue","white","red"))(n = 50)
    }
    #data.name <- "data/TCGA_expr_SAM_genes_SVD_R.txt"
    #dataset = read.delim(data.name, header=TRUE, sep="\t")
    #results.name <- "results/SVD_prov.png"
    #data.only = data.matrix(dataset[,2:ncol(dataset)])
    
    
    
    hr <- hclust(as.dist(1-cor(t(sign_mat), method="pearson")), method="average")
    #hc <- hclust(as.dist(1-cor(data.only, method="pearson")), method="average")
    
    distance <- dist(sign_mat)
    cluster <- hclust(distance, method = "average")
    
    ordered.mat <- sign_mat[order(hr$order), ]
    ordered.mat <- t(scale(t(ordered.mat)))
    
    
    library(circlize)
    col_fun = colorRamp2(col.rampSel, my_palette)
    
    if(annotation.add){
      p <- Heatmap(ordered.mat,cluster_columns = F,show_row_dend = F,show_column_dend = F,
                   col = col_fun,show_row_names = F,show_column_names = F,show_heatmap_legend = F, top_annotation = heat.legend)
    }else{
      p <- Heatmap(ordered.mat,cluster_columns = F,show_row_dend = F,show_column_dend = F,
                   col = col_fun,show_row_names = F,show_column_names = F,show_heatmap_legend = F)
    }
    png(results.heatmap, width = 7, height = 7, units = 'in', res = 300)
    print(p)
    dev.off()
    
    
    
  }
  if(return.plot.heat){
    #library(cowplot)
    library(grid)
    if(annotation.add){
      Heatmap(ordered.mat,cluster_columns = F,show_row_dend = F,show_column_dend = F,
              col = col_fun,show_row_names = F,show_column_names = F,show_heatmap_legend = F, top_annotation = heat.legend)
    }else{
      Heatmap(ordered.mat,cluster_columns = F,show_row_dend = F,show_column_dend = F,
              col = col_fun,show_row_names = F,show_column_names = F,show_heatmap_legend = F)
    }
  }
  
}
# BiocManager::install('samr')
library(samr)

Rna.heat <- MakeGeneHeatmap(gene.matrix= "results/RNA_reduced_SKCM.txt",results.heatmap = "./results/RNA_reduced_SKCM_MP_weiht_01.png", results.filtered.matrix = "./results/RNA_SKCM_SIGN_SAM_MP_weight.txt",
                            patient.class = class.jsvd, return.plot.heat = T, annotation.add = T, heat.legend = ha, fold.change.par = 2, number.genes = 100,
                            col.rampSel = c(-2, 0, 2), results.signif.genes = "results/RNA_seq_sam_genes.txt")

Met.heat <- MakeGeneHeatmap(gene.matrix = "results/Meth_reduced_SKCM.txt", patient.class = class.jsvd,
                            fold.change.par = 2, Methyl.plot = T, return.plot.heat = T, annotation.add = F,results.filtered.matrix = "./results/Methyl_SKCM_SIGN_SAM_MP_weight.txt",
                            results.heatmap = "./results/Met_jSVD_heatmapSKCM_MP_weight1.png", first.row = 2, number.genes = 100, col.rampSel = c(-1.5, 0, 1.5)) #, number.genes = 300)

png("./results/joinedHeatmap_SKCM_JSVD_full.png", width = 8, height = 5, units = 'in', res = 300)
ht_list = Rna.heat %v% Met.heat
draw(ht_list, annotation_legend_list = pd, annotation_legend_side = "right", align_annotation_legend = "heatmap_top")
dev.off()

##############################################


##################
# figure 4 #######
##################



# make survival cluster

patient.table <- read.table("./results/dataSKCMjsvd_2D/patients_SKCM_clinical.txt", header = T, sep = "\t")
class.jsvd <- read.table("results/jsvd_classes.txt")
class.jsvd_full <- class.jsvd

# remove patients without FU
prob.id <- setdiff(which(patient.table$CURATED_DAYS_TO_DEATH_OR_LAST_FU!="-"), grep("a", patient.table$CURATED_DAYS_TO_DEATH_OR_LAST_FU))
patient.table <- patient.table[prob.id, ]
class.jsvd <- class.jsvd[prob.id, ]

# remove patients with not available state
class.jsvd <- class.jsvd[grep("Not", patient.table$CURATED_MELANOMA_SPECIFIC_VITAL_STATUS..0...ALIVE.OR.CENSORED..1...DEAD.OF.MELANOMA., invert = T), ]
patient.table <- patient.table[grep("Not", patient.table$CURATED_MELANOMA_SPECIFIC_VITAL_STATUS..0...ALIVE.OR.CENSORED..1...DEAD.OF.MELANOMA., invert = T), ]



# set as dead
patient.table$CURATED_MELANOMA_SPECIFIC_VITAL_STATUS..0...ALIVE.OR.CENSORED..1...DEAD.OF.MELANOMA.[which(patient.table$CURATED_MELANOMA_SPECIFIC_VITAL_STATUS..0...ALIVE.OR.CENSORED..1...DEAD.OF.MELANOMA.=="1")] <- 2
# set as censored
patient.table$CURATED_MELANOMA_SPECIFIC_VITAL_STATUS..0...ALIVE.OR.CENSORED..1...DEAD.OF.MELANOMA.[which(patient.table$CURATED_MELANOMA_SPECIFIC_VITAL_STATUS..0...ALIVE.OR.CENSORED..1...DEAD.OF.MELANOMA.=="0")] <- 1




surv.table <- data.frame(cbind(patient.table$CURATED_MELANOMA_SPECIFIC_VITAL_STATUS..0...ALIVE.OR.CENSORED..1...DEAD.OF.MELANOMA.,
                               patient.table$CURATED_DAYS_TO_DEATH_OR_LAST_FU))

colnames(surv.table) <- c("MD","Time")
# survival times
futimes <-  round(as.numeric(patient.table$CURATED_DAYS_TO_DEATH_OR_LAST_FU)/30, digits = 0) 
patient.table$CURATED_DAYS_TO_DEATH_OR_LAST_FU <- round(as.numeric(patient.table$CURATED_DAYS_TO_DEATH_OR_LAST_FU)/30, digits = 0) 
# censoring status 1 censored 2 dead 
fustat <- patient.table$CURATED_MELANOMA_SPECIFIC_VITAL_STATUS..0...ALIVE.OR.CENSORED..1...DEAD.OF.MELANOMA.


library(survival)
library(survminer)

makeVectSelect<-function(selected.group=3){
  patient.class <- rep(0, length(fustat))
  patient.class[which(class.jsvd[, 2]==selected.group)] <- 1
  return(patient.class)
}

class.1 <- makeVectSelect(1)
class.2 <- makeVectSelect(2)
class.3 <- makeVectSelect()

surv.table <- data.frame(cbind(patient.table$MD, as.numeric(patient.table$Time), class.1, class.2, class.3))

colnames(surv.table)[c(1,2)] <- c("status", "time")
class1 <- survfit( Surv(time=futimes, event= as.numeric(fustat)) ~ class.1, data = surv.table)
class2 <- survfit( Surv(time=futimes, event= as.numeric(fustat)) ~ class.2, data = surv.table)
class3 <- survfit( Surv(time=futimes, event= as.numeric(fustat)) ~ class.3, data = surv.table)

fit <- list(class1, class2, class3)
p <- ggsurvplot_combine(fit, surv.table,pval = TRUE, palette =c("000000000","red","NA",
                                                                "blue","0000000","orange"), risk.table = FALSE,
                        legend.labs =  c(" ", paste0("class 1 (", sum(class.1),"/",length(which(class.jsvd_full[,2]==1)) ,"),"),"",
                                         paste0("class 2 (", sum(class.2),"/",length(which(class.jsvd_full[,2]==2)), ")"),"   ",
                                         paste0("class 3 (", sum(class.3), "/",length(which(class.jsvd_full[,2]==3)), ")")))

result.surv <- "results/single_groups_survival_jSVD.png"
png(result.surv, width = 7, height = 7, units = 'in', res = 300)
par(mar = c(3, 3, 3, 3))
print(p)
dev.off()

surv_pvalue(fit)

# check survival of 1 vs 3
survdiff(Surv(time=futimes[which(class.jsvd[,2]!=2)], event= as.numeric(fustat[which(class.jsvd[,2]!=2)])) ~ class.1[which(class.jsvd[,2]!=2)], data = surv.table[which(class.jsvd[,2]!=2), ])

# check survival 1 vs 2
survdiff(Surv(time=futimes[which(class.jsvd[,2]!=3)], event= as.numeric(fustat[which(class.jsvd[,2]!=3)])) ~ class.1[which(class.jsvd[,2]!=3)], data = surv.table[which(class.jsvd[,2]!=3), ])
# check survival 3 vs 2
survdiff(Surv(time=futimes[which(class.jsvd[,2]!=1)], event= as.numeric(fustat[which(class.jsvd[,2]!=1)])) ~ class.2[which(class.jsvd[,2]!=1)], data = surv.table[which(class.jsvd[,2]!=1), ])

#######################################################


##################
# figure 2 #######
##################

# check enriched pathway 


library(clusterProfiler)
library(org.Hs.eg.db)



RNA_SAM <- read.table("results/RNA_seq_sam_genes.txt", sep ="\t", header  =T)
# convert from entrezID to ensembl and symbol ID
gene.df <- bitr(RNA_SAM[,2], fromType = "SYMBOL",
                toType = c("ENSEMBL"),
                OrgDb = org.Hs.eg.db)

RNA_all <- read.table("results/RNA_reduced_SKCM.txt", sep ="\t", header  =T)

# convert from entrezID to ensembl and symbol ID
gene.dfsymbol <- bitr(rownames(RNA_all), fromType = "SYMBOL",
                      toType = c("ENSEMBL"),
                      OrgDb = org.Hs.eg.db)


ego <- enrichGO(gene = gene.df[,2],#significant genes in our analysis
                universe= gene.dfsymbol[,2],# all genes that were considered in the analysis
                keyType = "ENSEMBL",
                OrgDb = org.Hs.eg.db,
                ont = "BP",
                pAdjustMethod = "fdr",
                qvalueCutoff = 0.05,
                readable = TRUE)
cluster_summary <- data.frame(ego)
write.table(cluster_summary, "results/RNA_sam_significant_GO_BP.txt", sep  ="\t", quote = F, col.names = T, row.names = F)

dotplot(ego, showCategory=10)

#########################################################################

##################
# figure 3 #######
##################



# check methylation on significant genes

methylation.data <- read.table("data/gdac.broadinstitute.org_SKCM.Methylation_Preprocess.Level_3.2016012800.0.0/SKCM.meth.by_mean.data.txt", sep = "\t", header = T)

# load full jsvd class 

class.jsvd <- read.table("results/jsvd_classes.txt")

# take only SAM genes

meth.reduced <- methylation.data[match(RNA_SAM[, 2], methylation.data[, 1]), ]

# remove NA genes

meth.reduced <- meth.reduced[which(is.na(meth.reduced[,1])==F), ]

colnames(meth.reduced) <- gsub("\\.", "-", colnames(meth.reduced))

# take only same patients

meth.reduced <- meth.reduced[, c(1, match(class.jsvd[,1], colnames(meth.reduced)))]


# make plot of 1 vs 2 

RNA_all <- read.table("results/RNA_reduced_SKCM.txt", sep ="\t", header  =T)

RNA_selected <- RNA_all[match(meth.reduced[,1], rownames(RNA_all)), ]


rownames(meth.reduced) <- meth.reduced[,1]

meth.reduced <- meth.reduced[, (2:ncol(meth.reduced))]
# for all genes compute logFC

logFC_RNA21 <- sapply(1:nrow(RNA_selected), function(x){log(mean(as.numeric(unlist(RNA_selected[x, which(class.jsvd[,2]==2)])))/mean(unlist(as.numeric(RNA_selected[x, which(class.jsvd[,2]!=2)]))))})
logFC_Meth21 <- sapply(1:nrow(meth.reduced), function(x){log(mean(as.numeric(unlist(meth.reduced[x, which(class.jsvd[,2]==2)])))/mean(as.numeric(unlist(meth.reduced[x, which(class.jsvd[,2]!=2)]))))})

results.mat <- cbind(rownames(RNA_selected), logFC_RNA21, logFC_Meth21)

png("results/scatter_plot_methylation_expressionII.png", width = 8, height = 5, units = 'in', res = 300)

color.vect <- rep("azure4", length(logFC_RNA21))
color.vect[which(logFC_Meth21<0)] <- "red"

plot(logFC_Meth21, logFC_RNA21, xlab= expression(paste("Methylation ", log(mu[2]/mu[13]))), ylab =  expression(paste("RNA-seq ", log(mu[2]/mu[13]))),col = color.vect, pch = 16)

abline(v=0)

legend( "topright", 
        legend=c(as.character(length(which(logFC_Meth21<0))),
                 as.character(length(which(logFC_Meth21>0)))),
        col=c("red","azure4"), bg="transparent",
        pch=16, cex = 2, box.lty = 0)
dev.off()

#######################################################################


#################################
# supplementary figure 1  #######
#################################




# make three files with only class mutations

all.variants <- read.table("data/SKCM_hg38.maf", sep ="\t", header =T, quote = "")



#make variant files of class1

class1.maf <- sapply(1:length(which(class.jsvd[,2]==1)), function(x){grep(class.jsvd[which(class.jsvd[,2]==1)[x], 1], all.variants[,15])})
class1.maf <- all.variants[unlist(class1.maf),]

write.table(class1.maf, "results/class1_maf.maf", sep = "\t",row.names = F, quote = F)


class2.maf <- sapply(1:length(which(class.jsvd[,2]==2)), function(x){grep(class.jsvd[which(class.jsvd[,2]==2)[x], 1], all.variants[,15])})
class2.maf <- all.variants[unlist(class2.maf),]
write.table(class2.maf, "results/class2_maf.maf", sep = "\t",row.names = F, quote = F)


class3.maf <- sapply(1:length(which(class.jsvd[,2]==3)), function(x){grep(class.jsvd[which(class.jsvd[,2]==3)[x], 1], all.variants[,15])})
class3.maf <- all.variants[unlist(class3.maf),]
write.table(class3.maf, "results/class3_maf.maf", sep = "\t", row.names = F, quote = F)



# add variants domain and find most important SNPS for classification
library(maftools)


laml = read.maf(maf = "results/class1_maf.maf", isTCGA = T) #, vc_nonSyn = var.cat, isTCGA = F) #, clinicalData = laml.clin)



png("results/class1_SKCM.png", width = 10, height = 5, units = 'in', res = 300)
oncoplot(maf = laml,genes = c("BRAF", "NRAS","NF1", "CDKN2A")) # , top = 10)
dev.off()




laml = read.maf(maf = "results/class2_maf.maf", isTCGA = T) #, vc_nonSyn = var.cat, isTCGA = F) #, clinicalData = laml.clin)


png("results/class2_SKCM.png", width = 10, height = 5, units = 'in', res = 300)
oncoplot(maf = laml,genes = c("BRAF", "NRAS","NF1", "CDKN2A")) # , top = 10)
dev.off()


laml = read.maf(maf = "results/class3_maf.maf", isTCGA = T) #, vc_nonSyn = var.cat, isTCGA = F) #, clinicalData = laml.clin)


png("results/class3_SKCM.png", width = 10, height = 5, units = 'in', res = 300)
oncoplot(maf = laml,genes = c("BRAF", "NRAS","NF1", "CDKN2A")) # , top = 10)
dev.off()



