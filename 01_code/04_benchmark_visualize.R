## setwd('D:/Fan_Files/ICA/CMML 5.23/')
.libPaths("D:/Programs/R/R-4.4.1/library")
library(ggplot2)
library(dplyr)
library(gridExtra)


# 读取数据
data <- read.csv("./04_output/paper_result/03_performance_metrics.csv")

# 将detection_rate转换为分类变量
data$detection_rate <- as.factor(data$detection_rate)

# 分别绘制两个数据集的图
plot_for_dataset <- function(dataset_name, metric = "precision") {
  df <- data %>% 
    filter(dataset == dataset_name) %>% 
    select(method, detection_rate, !!sym(metric))  # 动态选择指标
  
  ggplot(df, aes(x = method, y = !!sym(metric), fill = detection_rate)) +
    geom_col(position = position_dodge(width = 0.8), width = 0.7) +
    scale_fill_brewer(palette = "Set2") +
    labs(
      title = paste0("Dataset: ", dataset_name),
      x = "Method",
      y = toupper(metric),
      fill = "Detection Rate"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      axis.text.x = element_text(angle = 45, hjust = 1)
    )
}

# 生成两个数据集的precision图
p1 <- plot_for_dataset("pbmc-ch", "precision")
p2 <- plot_for_dataset("hm-12k", "precision")

# 并排显示（需要gridExtra包）
grid.arrange(p1, p2, ncol = 2)



library(ggplot2)
library(dplyr)

# 读取数据
data <- read.csv("./04_output/DE_result.csv") %>% select(-1)  # 移除第一列序号列


# 绘制分组柱状图
ggplot(data, aes(x = Method, y = Value, fill = Metric)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +  # 分组柱状图
  geom_text(
    aes(label = round(Value, 3)), 
    position = position_dodge(width = 0.8), 
    vjust = -0.3, 
    size = 3.5, 
    color = "black"
  ) +  # 添加数值标签
  scale_fill_manual(
    values = c("Precision" = "#1f77b4", "Recall" = "#ff7f0e", "TNR" = "#2ca02c")
  ) +  # 自定义颜色
  labs(
    title = "Method Performance Comparison Across Metrics",
    x = "Method",
    y = "Value",
    fill = "Metric"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),  # 调整X轴标签角度
    legend.position = "top",  # 图例置于顶部
    panel.grid.major.x = element_blank()  # 隐藏X轴主要网格线
  ) +
  scale_y_continuous(
    limits = c(0, 1.05), 
    breaks = seq(0, 1, 0.2), 
    expand = c(0, 0)  # 避免Y轴空白
  )



library(ggplot2)
library(tidyr)
library(dplyr)

# 读取数据并整理格式
data <- read.csv("./04_output/sim_rate2.csv") %>% 
  select(-1) %>%  # 移除第一列序号
  pivot_longer(
    cols = c(precision, recall, tnr),
    names_to = "metric",
    values_to = "value"
  )

# 统一规范列名（可选）
colnames(data) <- c("method", "detection_rate", "metric", "value")

# 生成分面折线图
ggplot(data, aes(x = detection_rate, y = value, color = method)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  facet_wrap(~ metric, ncol = 3, scales = "free_y") +
  scale_color_manual(values = c("#1f77b4", "#ff7f0e", "#2ca02c")) +  # 蓝/橙/绿
  labs(
    title = "Method Performance Across Doublet Rates",
    x = "Detection Rate",
    y = "Value",
    color = "Method"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    strip.text = element_text(face = "bold"),  # 分面标题加粗
    legend.position = "bottom"
  ) +
  scale_x_continuous(breaks = unique(data$detection_rate))  # 显示所有检测率刻度




# 将detection_rate转换为因子（强制等间距）
data <- data %>% 
  mutate(
    detection_rate_factor = factor(
      detection_rate,
      levels = unique(sort(detection_rate)),  # 保持原始顺序
      labels = as.character(unique(sort(detection_rate)))  # 显示原始数值标签
    )
  )

# 生成分面折线图
ggplot(data, aes(x = detection_rate_factor, y = value, color = method, group = method)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  facet_wrap(~ metric, ncol = 3, scales = "free_y") +
  scale_color_manual(values = c("#1f77b4", "#ff7f0e", "#2ca02c")) +
  labs(
    title = "Method Performance Across Doublet Rates (Uniform X-axis)",
    x = "Detection Rate",
    y = "Value",
    color = "Method"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    strip.text = element_text(face = "bold"),
    legend.position = "bottom",
    axis.text.x = element_text(angle = 45, hjust = 1)  # 调整标签角度防止重叠
  )
