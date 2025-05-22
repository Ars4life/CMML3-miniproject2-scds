library(Matrix)
library(Seurat)
library(scds)
library(SingleCellExperiment)
library(PRROC)

# simulation data with different doublet rates


sim.data <- readRDS("./02_data/sim_data/sim_rate.rds")
rt1 <- c(0.02,0.04,0.08,0.1,0.2, 0.4)
result_4 <- data.frame()
score.list.cxds <- list()
score.list.bcds <- list()
score.list.hybrid <- list()

for (r in rt1) {
  data <- sim.data[[as.character(r)]]
  if (is.null(data)) {
    warning(paste("No data found for rate:", rate))
    next
  }
  data.matrix <- data[[1]][[1]]; dim(data.matrix)
  data.type <- data[[2]][[1]]; table(data.type)
  data.type <- ifelse(data.type=='doublet', 1, 0); table(data.type)
  index.doublet <- which(data.type==1)
  data.clean <- data.matrix[,-index.doublet]; dim(data.clean)
  sce <- SingleCellExperiment(assays = list(counts = data.matrix))
  sce <- cxds_bcds_hybrid(sce)
  CD <- colData(sce)
  score.cxds <- CD$cxds_score
  score.bcds <- CD$bcds_score
  score.hybrid <- CD$hybrid_score
  # save scores
  methods <- list(
  cxds = score.cxds,
  bcds = score.bcds,
  hybrid = score.hybrid)
  score.list.cxds <- append(score.list.cxds, list(score.cxds))
  score.list.bcds <- append(score.list.bcds, list(score.bcds))
  score.list.hybrid <- append(score.list.hybrid, list(score.hybrid))
  for(method_name in names(methods)){
    label <- ifelse(data[[2]][[1]] == 'doublet', 1, 0)
    score <- methods[method_name][[1]]
    d <- floor(length(label) * r)
    thresh <- sort(score, decreasing = TRUE)[d]
    pred <- score > thresh
    # 计算指标
    tp <- sum(pred[label == 1]);tp
    fp <- sum(pred[label == 0]);fp
    fn <- sum(!pred[label == 1]);fn
    tn <- sum(!pred[label == 0]);tn
    print(method_name)
    precision <- tp / (tp + fp)
    recall <- tp / (tp + fn)
    tnr <- tn / (tn + fp)
    
    result_4 <- rbind(result_4, data.frame(
      method = method_name,
      detection_rate = r,
      precision = round(precision, 4),
      recall = round(recall, 4),
      tnr = round(tnr, 4)
    ))

  }
}
write.csv(result_4,'./04_output/sim_rate2.csv')

for (r in rt1) {
  data <- sim.data[[as.character(r)]]
  if (is.null(data)) {
    warning(paste("No data found for rate:", rate))
    next
  }
  data.matrix <- data[[1]][[1]]; dim(data.matrix)
  data.type <- data[[2]][[1]]; table(data.type)
  data.type <- ifelse(data.type=='doublet', 1, 0); table(data.type)
  index.doublet <- which(data.type==1)
  data.clean <- data.matrix[,-index.doublet]; dim(data.clean)
  sce <- SingleCellExperiment(assays = list(counts = data.matrix))
  sce <- cxds_bcds_hybrid(sce)
  CD <- colData(sce)
  score.cxds <- CD$cxds_score
  score.bcds <- CD$bcds_score
  score.hybrid <- CD$hybrid_score
  # save scores
  methods <- list(
    cxds = score.cxds,
    bcds = score.bcds,
    hybrid = score.hybrid)
  score.list.cxds <- append(score.list.cxds, list(score.cxds))
  score.list.bcds <- append(score.list.bcds, list(score.bcds))
  score.list.hybrid <- append(score.list.hybrid, list(score.hybrid))
  for(method_name in names(methods)){
    label <- ifelse(data[[2]][[1]] == 'doublet', 1, 0)
    score <- methods[method_name][[1]]
    d <- floor(length(label) * r)
    thresh <- sort(score, decreasing = TRUE)[d]
    pred <- score > thresh
    # 计算指标
    tp <- sum(pred[label == 1]);tp
    fp <- sum(pred[label == 0]);fp
    fn <- sum(!pred[label == 1]);fn
    tn <- sum(!pred[label == 0]);tn
    print(method_name)
    precision <- tp / (tp + fp)
    recall <- tp / (tp + fn)
    tnr <- tn / (tn + fp)
    
    result_4 <- rbind(result_4, data.frame(
      method = method_name,
      detection_rate = r,
      precision = round(precision, 4),
      recall = round(recall, 4),
      tnr = round(tnr, 4)
    ))
    
  }
}