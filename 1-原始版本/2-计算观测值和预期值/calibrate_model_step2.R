#' @title 线性模型校准与相关性显著性检验
#' 
#' @description
#' 本脚本用于对观测到的中性变异水平与突变似然得分之间的关系进行线性模型拟合，并保存模型系数与截距，供后续计算预期值使用。同时可选地对相关性进行显著性检验，比较不同变量与观测值之间的相关性强度。
#' 
#' @param prefix 字符串，文件名前缀，用于区分不同数据集或分析流程。
#' @param do_sig_test 字符串，是否进行相关性显著性检验（"yes" 或 "no"），默认为 "yes"。
#' 
#' @details
#' 1. 读取观测值与突变得分的输入文件（分为参考区和OriB-OriH区）。
#' 2. 对不同突变组（G>A_and_T>C 和 other）分别拟合线性模型（obs_max_het ~ sum_likelihood），提取系数与截距并保存。
#' 3. 若选择进行显著性检验，则：
#'    - 对各组数据拟合多元线性模型（obs_max_het ~ sum_likelihood + length），输出模型摘要。
#'    - 使用 cocor 包对相关系数进行依赖性检验，比较 sum_likelihood 与 length 对 obs_max_het 的相关性显著性。
#'    - 检验结果保存至指定输出文件。
#' 
#' @return 无返回值，结果写入指定输出文件。
#' 
#' @examples
#' calibrate_step2(prefix = "example_", do_sig_test = "yes")
#' 
#' @author GitHub Copilot
# 保存线性模型拟合结果，供后续计算预期值使用。

args <- commandArgs(trailingOnly = TRUE)

calibrate_step2 <- function(prefix = NA, do_sig_test = NA){
  # 设置默认值
  if(is.na(prefix)) {
    prefix = ""
  } 
  if(is.na(do_sig_test)) {
    do_sig_test = "yes"
  } 
  
  file <- read.delim(file = sprintf('output/calibration/%sloci_obs_vs_scores.txt', prefix), header = TRUE, sep = "\t")
  ori_file <- read.delim(file = sprintf('output/calibration/%sloci_obs_vs_scores_ori.txt', prefix), header = TRUE, sep = "\t")
  
  model_fits <- rbind(c("region", "mutation_group", "item", "value"),
                      c("ref_exc_ori", "G>A_and_T>C", "coefficient",  coef(lm(obs_max_het ~ sum_likelihood, data = file[file$mutation_group == "G>A_and_T>C",]))["sum_likelihood"]),
                      c("ref_exc_ori", "G>A_and_T>C", "intercept", coef(lm(obs_max_het ~ sum_likelihood, data = file[file$mutation_group == "G>A_and_T>C",]))["(Intercept)"]),
                      c("ref_exc_ori", "other", "coefficient", coef(lm(obs_max_het ~ sum_likelihood, data = file[file$mutation_group == "other",]))["sum_likelihood"]),
                      c("ref_exc_ori", "other", "intercept", coef(lm(obs_max_het ~ sum_likelihood, data = file[file$mutation_group == "other",]))["(Intercept)"]),
                      c("ori", "G>A_and_T>C", "coefficient", coef(lm(obs_max_het ~ sum_likelihood, data = ori_file[ori_file$mutation_group == "G>A_and_T>C",]))["sum_likelihood"]),
                      c("ori", "G>A_and_T>C", "intercept", coef(lm(obs_max_het ~ sum_likelihood, data = ori_file[ori_file$mutation_group == "G>A_and_T>C",]))["(Intercept)"]),
                      c("ori", "other", "coefficient", coef(lm(obs_max_het ~ sum_likelihood, data = ori_file[ori_file$mutation_group == "other",]))["sum_likelihood"]),
                      c("ori", "other", "intercept", coef(lm(obs_max_het ~ sum_likelihood, data = ori_file[ori_file$mutation_group == "other",]))["(Intercept)"]))
  
  write.table(model_fits, file = sprintf("output/calibration/%slinear_model_fits.txt", prefix), sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
  
  
  if (do_sig_test == "yes"){
    
    # 初始化输出文件
    sink(sprintf("output/calibration/%scorrelation_significance_testing.txt", prefix), append = FALSE)
    
    # 检验哪个变量与观测到的中性变异相关性更强
    print("G>A and T>C mutation group (for reference excluding OriB-OriH)")
    print(summary(lm(obs_max_het ~ sum_likelihood + length, data = file[file$mutation_group == "G>A_and_T>C",])))
    print("All other mutation group (for reference excluding OriB-OriH)")
    print(summary(lm(obs_max_het ~ sum_likelihood + length, data = file[file$mutation_group == "other",])))
    print("G>A and T>C mutation group for OriB-OriH region")
    print(summary(lm(obs_max_het ~ sum_likelihood + length, data = ori_file[ori_file$mutation_group == "G>A_and_T>C",])))
    print("All other mutation group for OriB-OriH region")
    print(summary(lm(obs_max_het ~ sum_likelihood + length, data = ori_file[ori_file$mutation_group == "other",])))
    
    # 借助 cocor 包展示：观测到的中性变异水平与突变似然得分的相关性比与位点长度的相关性更显著。
    # 这是针对相关变量存在共享自变量（即相关系数依赖）的两组相关性进行单侧假设检验。
    # 使用 cor.test 提取用于检验的相关系数。
    
    library(cocor)
    
    print("cocor package analysis")
    print("G>A and T>C mutation group (for reference excluding OriB-OriH)")
    print(cocor.dep.groups.overlap(cor.test(file[file$mutation_group == "G>A_and_T>C", c("obs_max_het")], file[file$mutation_group == "G>A_and_T>C", c("sum_likelihood")], method = "pearson")$estimate, 
                             cor.test(file[file$mutation_group == "G>A_and_T>C", c("obs_max_het")], file[file$mutation_group == "G>A_and_T>C", c("length")], method = "pearson")$estimate, 
                             cor.test(file[file$mutation_group == "G>A_and_T>C", c("sum_likelihood")], file[file$mutation_group == "G>A_and_T>C", c("length")], method = "pearson")$estimate, 
                             n = 39, alternative = "greater", alpha = 0.05, conf.level = 0.95, null.value = 0))
    
    print("All other mutation group (for reference excluding OriB-OriH)")
    print(cocor.dep.groups.overlap(cor.test(file[file$mutation_group == "other", c("obs_max_het")], file[file$mutation_group == "other", c("sum_likelihood")], method = "pearson")$estimate, 
                             cor.test(file[file$mutation_group == "other", c("obs_max_het")], file[file$mutation_group == "other", c("length")], method = "pearson")$estimate, 
                             cor.test(file[file$mutation_group == "other", c("sum_likelihood")], file[file$mutation_group == "other", c("length")], method = "pearson")$estimate, 
                             n = 39, alternative = "greater", alpha = 0.05, conf.level = 0.95, null.value = 0))
    
    print("G>A and T>C mutation group for OriB-OriH region")
    print(cocor.dep.groups.overlap(cor.test(ori_file[ori_file$mutation_group == "G>A_and_T>C", c("obs_max_het")], ori_file[ori_file$mutation_group == "G>A_and_T>C", c("sum_likelihood")], method = "pearson")$estimate, 
                             cor.test(ori_file[ori_file$mutation_group == "G>A_and_T>C", c("obs_max_het")], ori_file[ori_file$mutation_group == "G>A_and_T>C", c("length")], method = "pearson")$estimate, 
                             cor.test(ori_file[ori_file$mutation_group == "G>A_and_T>C", c("sum_likelihood")], ori_file[ori_file$mutation_group == "G>A_and_T>C", c("length")], method = "pearson")$estimate, 
                             n = 8, alternative = "greater", alpha = 0.05, conf.level = 0.95, null.value = 0))
    
    print("All other mutation group for OriB-OriH region")
    print(cocor.dep.groups.overlap(cor.test(ori_file[ori_file$mutation_group == "other", c("obs_max_het")], ori_file[ori_file$mutation_group == "other", c("sum_likelihood")], method = "pearson")$estimate, 
                             cor.test(ori_file[ori_file$mutation_group == "other", c("obs_max_het")], ori_file[ori_file$mutation_group == "other", c("length")], method = "pearson")$estimate, 
                             cor.test(ori_file[ori_file$mutation_group == "other", c("sum_likelihood")], ori_file[ori_file$mutation_group == "other", c("length")], method = "pearson")$estimate, 
                             n = 8, alternative = "greater", alpha = 0.05, conf.level = 0.95, null.value = 0))
    sink()
  }
}

calibrate_step2(prefix = args[1], do_sig_test = args[2])
