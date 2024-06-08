library(TwoSampleMR)

bmi_file <- system.file("extdata", "bmi.txt", package = "TwoSampleMR")
bmi_exp_dat <- read_exposure_data(bmi_file)


a <- mtcars

# 创建示例数据集
data <- data.frame(
  x = c(1, 2, 3, 4, 5), 
  y = c(2, 4, 6, 8, 10)
)


model <- lm(y ~ x, data = data)
summary(model)

coef_x <- model$coefficients[2]  # x变量的回归系数
se_x <- model$coefficients[2,2]  # x变量的标准误差



library(tidyverse)

# 生成模拟数据
protein_levels <- c(15, 16, 17, 18, 19)
weight_means <- c(8.3, 8.66, 8.92, 9.06, 9.21)
weight_sd <- 0.5

df <- data.frame(
  protein_level = rep(protein_levels, each = 50),
  weight = c(
    rnorm(50, weight_means[1], weight_sd),
    rnorm(50, weight_means[2], weight_sd),
    rnorm(50, weight_means[3], weight_sd),
    rnorm(50, weight_means[4], weight_sd),
    rnorm(50, weight_means[5], weight_sd)
  )
)

# 绘制图表
ggplot(df, aes(x = protein_level, y = weight, fill = as.factor(protein_level))) +
  geom_boxplot() +
  labs(x = "Protein Level (%)", y = "Weight (kg)") +
  theme_minimal()

# 绘制图表
ggplot(df, aes(x = protein_level, y = weight, fill = as.factor(protein_level))) +
  geom_boxplot() +
  geom_line(data = data.frame(protein_level = protein_levels, weight_mean = weight_means),
            aes(x = protein_level, y = weight_mean),
            color = "red",
            size = 1.5) +
  labs(x = "Protein Level (%)", y = "Weight (kg)") +
  theme_minimal()

# 计算每个蛋白质含量下的平均体重
weight_summary <- df %>%
  group_by(protein_level) %>%
  summarize(weight_mean = mean(weight))

# 绘制图表
ggplot(df, aes(x = protein_level, y = weight, fill = as.factor(protein_level))) +
  geom_boxplot() +
  geom_line(data = weight_summary,
            aes(x = protein_level, y = weight_mean, group = 1),
            color = "red",
            size = 1.5) +
  labs(x = "Protein Level (%)", y = "Weight (kg)") +
  theme_minimal()
