## load in packages
# setwd('D:/Fan_Files/ICA/CMML 5.23/')
.libPaths("D:/Programs/R/R-4.4.1/library")
library(Matrix)
library(Seurat)
library(scds)
library(SingleCellExperiment)
library(PRROC)
library(reticulate)
library(dplyr)
library(pbapply)
set.seed(123)

## calculate doublet score（ pbmc-ch, hm-12k）
# location of target dataset
locs <- c('./02_data/real_data/pbmc-ch.rds', './02_data/real_data/hm-12k.rds')
test <- readRDS('./02_data/real_data/pbmc-ch.rds')
summary(test[[2]])
sum(test[[2]]=='doublet')
2545/15272
730/12820



score.list.cxds <- list()
score.list.bcds <- list()
score.list.hybrid <- list()


for(loc in locs){
  data <- readRDS(loc)
  count <- data[[1]]; dim(count)
  label <- data[[2]]; table(label)
  label <- ifelse(label == 'doublet', 1, 0); table(label)
  doublet.rate <- sum(label==1) / length(label); doublet.rate
  sce <- SingleCellExperiment(assays = list(counts = count))
  # cxds, bcds, hybrid
  scdstime <- system.time({
    sce <- cxds_bcds_hybrid(sce)
  })
  CD <- colData(sce)
  score.cxds <- CD$cxds_score
  score.bcds <- CD$bcds_score
  score.hybrid <- CD$hybrid_score
  # save scores
  score.list.cxds <- append(score.list.cxds, list(score.cxds))
  score.list.bcds <- append(score.list.bcds, list(score.bcds))
  score.list.hybrid <- append(score.list.hybrid, list(score.hybrid))
}
# save results
saveRDS(score.list.cxds, './04_output/paper_result/cxds_real_score_pbmc_hm6k.rds')
saveRDS(score.list.bcds, './04_output/paper_result/bcds_real_score_pbmc_hm6k.rds')
saveRDS(score.list.hybrid, './04_output/paper_result/hybrid_real_score_pbmc_hm6k.rds')


## -----------------------------------------------------------------------------
## 性能评估：在 10%, 20%, 40% 识别率下的精确度、召回率、真阴性率


## assessment
methods <- list(
  cxds = score.list.cxds,
  bcds = score.list.bcds,
  hybrid = score.list.hybrid
)

rs <- c(0.1, 0.2, 0.4)
results_3 <- data.frame()
results_4 <- data.frame()
results_5 <- data.frame()

for(method_name in names(methods)){
  score.list <- methods[[method_name]]
  
  for(i in seq_along(locs)){
    data <- readRDS(locs[i])
    dataset_name <- gsub(".*/(.*)\\.rds", "\\1", locs[i])
    label <- ifelse(data[[2]] == 'doublet', 1, 0)
    score <- score.list[[i]]
    
    for(r in rs){
      # calculate threshold
      d <- floor(length(label) * r)
      thresh <- sort(score, decreasing = TRUE)[d]
      pred <- score > thresh
      
      # 计算指标
      tp <- sum(pred[label == 1])
      fp <- sum(pred[label == 0])
      fn <- sum(!pred[label == 1])
      tn <- sum(!pred[label == 0])
      
      precision <- tp / (tp + fp)
      recall <- tp / (tp + fn)
      tnr <- tn / (tn + fp)
      
      # 构建结果行
      results_3 <- rbind(results_3, data.frame(
        method = method_name,
        dataset = dataset_name,
        detection_rate = r,
        precision = round(precision, 4),
        recall = round(recall, 4),
        tnr = round(tnr, 4)
      ))
    }
  }
}

# 保存评估结果为CSV
write.csv(results_4, "./04_output/paper_result/03_performance_metrics.csv", row.names = FALSE)

###-----------------------------------------------------------------------------
## scrublet

# read python module
scr <- import('scrublet')
scipy.io <- import('scipy.io')
np <- import('numpy')
os <- import('os')

# list to save doublet scores
score.list.scrublet <- list()

# loop over each dataset
for(loc in locs){
  data <- readRDS(loc)
  count <- data[[1]]; dim(count)
  label <- data[[2]]; table(label)
  label <- ifelse(label == 'doublet', 1, 0); table(label)
  doublet.rate <- sum(label==1) / length(label); doublet.rate
  # scrublet
  scrub_time <- system.time({
    result <- scr$Scrublet(counts_matrix = t(count), expected_doublet_rate = doublet.rate, random_state = 10L)
    results <- result$scrub_doublets(min_counts=2, 
                                     min_cells=3, 
                                     min_gene_variability_pctl=85, 
                                     n_prin_comps=30L)
  })
  score <- as.vector(results[[1]])
  score.list.scrublet <- append(score.list.scrublet, list(score))
}
# save results, change the location accordingly
saveRDS(score.list.scrublet, '04_output/paper_result/scrublet_real_score.rds')

results_4 <- results_3
# loop over identification rates
for(r in rs){
  print('====================')
  print(r)
  precisions <- c()
  recalls <- c()
  tnrs <- c()
  result <- matrix(data = 0, nrow = length(locs), ncol=3)
  for(i in 1:length(locs)){
    print(locs[i])
    data <- readRDS(locs[i])
    dataset_name <- gsub(".*/(.*)\\.rds", "\\1", locs[i])
    # obtain the doublet labels
    label <- data[[2]]; table(label)
    label <- ifelse(label == 'doublet', 1, 0); table(label)
    # calculate threshold based on identification rate
    score <- score.list.scrublet[[i]]
    d <- floor(length(label) * r); d
    thresh <- sort(score, decreasing = T)[d]; thresh
    # predict doublet based on threshold
    pred <- score > thresh; table(pred)
    # result
    tp <- sum(pred[which(label==1)]==1); tp
    fp <- sum(pred[which(label==0)]==1); fp
    fn <- sum(pred[which(label==1)]==0); fn
    tn <- sum(pred[which(label==0)]==0); tn
    
    precision <- tp/(tp + fp); precision
    recall <- tp/(tp + fn); recall
    tnr <- tn/(tn + fp); tnr
    
    precisions[i] <- precision
    recalls[i] <- recall
    tnrs[i] <- tnr
    
    results_4 <- rbind(results_4, data.frame(
      method = 'Scrublet',
      dataset = dataset_name,
      detection_rate = r,
      precision = round(precision, 4),
      recall = round(recall, 4),
      tnr = round(tnr, 4)
    ))
    }
}

##------------------------------------------------------------------------------
## doubletFinder
score.list <- list()

# loop over each dataset
system.time({
  for(loc in locs){
    # loc <- locs[i]
    print(loc)
    data <- readRDS(loc)
    count <- data[[1]]; dim(count)
    label <- data[[2]]; table(label)
    label <- ifelse(label == 'doublet', 1, 0); table(label)
    doublet.rate <- sum(label==1) / length(label); doublet.rate
    
    # doubletfinder
    dblf_time <- system.time({
      ## Pre-process Seurat object (standard) --------------------------------------------------------------------------------------
      doublet.seurat <- CreateSeuratObject(counts = count, project = "doublet", min.cells = 1); doublet.seurat
      doublet.seurat <- NormalizeData(doublet.seurat)
      doublet.seurat <- ScaleData(doublet.seurat)
      doublet.seurat <- FindVariableFeatures(doublet.seurat, selection.method = "vst", nfeatures = 2000)
      doublet.seurat <- RunPCA(doublet.seurat)
      
      ## pK Identification (no ground-truth) ---------------------------------------------------------------------------------------
      sweep.res.doublet <- paramSweep(doublet.seurat, PCs = 1:10)
      sweep.stats.doublet <- summarizeSweep(sweep.res.doublet, GT = FALSE)
      bcmvn.doublet <- find.pK(sweep.stats.doublet)
      pK <- bcmvn.doublet$pK[which.max(bcmvn.doublet$BCmetric)]; pK <- as.numeric(levels(pK))[pK]; pK
      doublet.seurat <- doubletFinder(doublet.seurat, PCs = 1:10, pN = 0.25, pK = pK, nExp = sum(label==1))
      attribute <- paste('pANN', 0.25, pK, sum(label==1), sep = '_'); attribute
    })
    score <- doublet.seurat@meta.data[[attribute]]
    score.list <- append(score.list, list(score))
  }
})
# save results, change the location accordingly
saveRDS(score.list, '04_output/paper_result/doubletfinder_real_score.rds')
results_5 <- results_4
# loop over identification rates
for(r in rs){
  print('====================')
  print(r)
  precisions <- c()
  recalls <- c()
  tnrs <- c()
  result <- matrix(data = 0, nrow = length(locs), ncol=3)
  for(i in 1:length(locs)){
    print(locs[i])
    data <- readRDS(locs[i])
    dataset_name <- gsub(".*/(.*)\\.rds", "\\1", locs[i])
    
    # obtain the doublet labels
    label <- data[[2]]; table(label)
    label <- ifelse(label == 'doublet', 1, 0); table(label)
    # calculate threshold based on identification rate
    score <- score.list[[i]]
    d <- floor(length(label) * r); d
    thresh <- sort(score, decreasing = T)[d]; thresh
    # predict doublet based on threshold
    pred <- score > thresh; table(pred)
    # result
    tp <- sum(pred[which(label==1)]==1); tp
    fp <- sum(pred[which(label==0)]==1); fp
    fn <- sum(pred[which(label==1)]==0); fn
    tn <- sum(pred[which(label==0)]==0); tn
    
    precision <- tp/(tp + fp); precision
    recall <- tp/(tp + fn); recall
    tnr <- tn/(tn + fp); tnr
    
    precisions[i] <- precision
    recalls[i] <- recall
    tnrs[i] <- tnr
    results_5 <- rbind(results_5, data.frame(
      method = 'DoubletFinder',
      dataset = dataset_name,
      detection_rate = r,
      precision = round(precision, 4),
      recall = round(recall, 4),
      tnr = round(tnr, 4)
    ))
    }
}


