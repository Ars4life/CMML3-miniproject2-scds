# setwd('D:/Fan_Files/ICA/CMML 5.23/')
library(Seurat)
library(Matrix)

# 数据集路径
locs <- c(
  './02_data/real_data/pbmc-ch.rds',
  './02_data/real_data/hm-12k.rds',
  './02_data/sim_data/sim_DE.rds',
  './02_data/sim_data/sim_type.rds'
)

# 创建输出目录（如果不存在）
output_dir <- "./02_data/03_processed_data"

if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# 预处理函数
preprocess_dataset <- function(loc) {
  tryCatch({
    # 读取数据
    raw_data <- readRDS(loc)
    count_matrix <- raw_data[[1]]  # 假设数据存储结构为list(counts, labels)
    
    # 记录开始时间
    start_time <- Sys.time()
    
    # 创建Seurat对象
    seu <- CreateSeuratObject(
      counts = count_matrix,
      project = basename(loc),
      min.cells = 1,
      min.features = 200
    )
    
    # 数据标准化
    seu <- NormalizeData(seu, normalization.method = "LogNormalize")
    
    # 特征缩放
    seu <- ScaleData(seu, features = rownames(seu))
    
    # 查找高变基因
    seu <- FindVariableFeatures(
      seu,
      selection.method = "vst",
      nfeatures = 2000,
      verbose = FALSE
    )
    
    # PCA降维
    seu <- RunPCA(seu, npcs = 20, verbose = FALSE)
    
    # 计算处理耗时
    processing_time <- round(as.numeric(difftime(Sys.time(), start_time, units = "mins")), 1)
    
    # 保存处理后的数据
    output_path <- file.path(output_dir, paste0("processed_", basename(loc)))
    saveRDS(seu, output_path)
    
    # 返回处理信息
    return(list(
      dataset = basename(loc),
      status = "Success",
      time_mins = processing_time,
      n_cells = ncol(seu),
      n_features = nrow(seu)
    ))
  }, error = function(e) {
    return(list(
      dataset = basename(loc),
      status = paste("Error:", e$message),
      time_mins = NA,
      n_cells = NA,
      n_features = NA
    ))
  })
}

# 执行批量处理
log_df <- data.frame()
for (i in seq_along(locs)) {
  cat("\nProcessing", i, "/", length(locs), ":", locs[i], "\n")
  
  # 执行预处理
  log_entry <- preprocess_dataset(locs[i])
  
  # 记录日志
  log_df <- rbind(log_df, as.data.frame(log_entry))
  
  # 打印进度
  cat("-> Status:", log_entry$status, "\n")
  if (log_entry$status == "Success") {
    cat("Processed", log_entry$n_cells, "cells with", 
        log_entry$n_features, "features\n")
    cat("Time elapsed:", log_entry$time_mins, "minutes\n")
  }
}


# 最终状态报告
success_rate <- mean(log_df$status == "Success")
cat("\nProcessing completed!\nSuccess rate:", 
    scales::percent(success_rate), "\n")
