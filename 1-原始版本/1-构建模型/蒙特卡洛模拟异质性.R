library(doParallel)
library(dplyr)
library(parallel)
library(stats)
library(tidyr)

setwd("/mnt/f/OneDrive/文档（科研）/脚本/Download/16-Constraint-mtDNA-model/2-整理后/1-构建模型")
dir.create("../output/蒙特卡洛模拟/")

# 用于跨世代（卵母细胞发生）模拟杂合度的函数
# 需要用户提供一系列输入参数

# 这段 R 代码是一个**线粒体突变与瓶颈效应的模拟函数**。
# 它不是分析 gnomAD 数据的主流程，而是一个**遗传动力学模型**，用来模拟母系传递中 mtDNA 异质性随世代的变化。
# 简要说明如下：
# ### 一、函数目的
# `mtDNA_simulate()` 通过蒙特卡洛法模拟：
# * 受精卵中的 mtDNA 在经历细胞分裂、瓶颈、突变、回突变、再扩增等阶段时，
# * **异质性（heteroplasmy）** 如何变化。
# 它输出每一代的 heteroplasmy 水平，用于研究：
# * 为什么有些突变能在几代后固定；
# * 为什么有些突变被瓶颈随机消失；
# * 不同突变率、瓶颈大小等参数如何影响这些过程。

# ### 二、主要参数与循环层级
# 外层循环：
# * `back_mutation_rates`：回突变概率比例。
# * `num_mt_bottleneck`：瓶颈阶段的 mtDNA 拷贝数（例如 8 – 512）。
# * `mutation_rates`：突变率。
# * `starting_heteroplasmy`：初始异质性。
# * `num_generations`：模拟的代数。
# 每个组合都会独立运行一次完整的模拟。
# ### 三、关键生物阶段建模
# 1. **瓶颈阶段（原始生殖细胞形成）**
#    模拟 mtDNA 在连续细胞分裂中的随机分配（`rhyper` 超几何分布）。
#    杂合度随拷贝随机分布而漂移，形成“遗传漂变”。
# 2. **突变与回突变**
#    * `rbinom` 用于随机生成正向突变（参考碱基→突变碱基）。
#    * `back_mut` 模拟突变碱基回到参考型。
#    * 突变率由 `mutation_rates` 控制，回突率为 `back_mutation_rates * mutation_rate`。
# 3. **细胞扩增与分裂**
#    每个生殖细胞经历 18 次有丝分裂，每次复制与分裂都会影响突变比例。
# 4. **克隆扩增（卵母细胞形成）**
#    在瓶颈后，mtDNA 拷贝回到 2^19 数量级。
#    再经历若干复制与回突，得到下一代卵母细胞的最终异质性。
# ### 四、输出结果
# 每一代都会写入结果文件：
# ```
# individual   heteroplasmy   generation   mutation_rate   bottleneck_size   back_mutation_rate   starting_heteroplasmy
# ```
# 可用于：
# * 画出异质性在多代的漂移轨迹；
# * 比较不同参数对异质性维持或丢失的影响。


mtDNA_simulate <- function(p){
  for(z in back_mutation_rates){
    for(i in num_mt_bottleneck){
      for(j in mutation_rates){
        for(start_het in starting_heteroplasmy){
          for(gen in 1:num_generations){
            
            if(gen == 1){
              heteroplasmy = start_het
            }
            print(paste("generation:", gen, "heteroplasmy:", heteroplasmy))
            
            # 瓶颈阶段：从受精卵形成原始生殖细胞（PGC）
            if(heteroplasmy > 0){
              # 步骤 1 - 瓶颈，仅在发生突变时才执行，否则杂合度将保持为 0
              # 在未发生 mtDNA 复制的情况下，每次细胞分裂 mtDNA 拷贝数减半
              # 需要模拟的分裂次数取决于最终瓶颈大小 i
              num_divisions = ifelse(i == 8, 16, 
                                     ifelse(i == 16, 15, 
                                            ifelse(i == 32, 14, 
                                                   ifelse(i == 64, 13, 
                                                          ifelse(i == 128, 12, 
                                                                 ifelse(i==256, 11, 
                                                                        ifelse(i==512, 10, 
                                                                               print("error"))))))))
              # 设定初始合子中的 mtDNA 拷贝数
              k = num_mt_zygote
              
              # 模拟 mtDNA 拷贝在子细胞之间的随机分配
              for(division in 1:num_divisions){
                nn = 1  # 观测次数
                m = k * heteroplasmy  # 替换等位基因数量
                n = k - m  # 参考等位基因数量
                k = k / 2  # 分配给子细胞的 mtDNA 拷贝数
                
                # rhyper 从超几何分布中生成随机数
                # 细胞分裂后的杂合度随机分配，rhyper 返回被选中的替换等位基因数
                heteroplasmy = rhyper(nn, m, n, k) / k
                
                print(paste("heteroplasmy after bottleneck division ", division, "is", heteroplasmy))
              }
              
              if(k != i){
                stop("bottleneck size after cell divisions does not equal bottleneck")
              }
            } 
            print(paste("heteroplasmy after bottleneck", heteroplasmy))
            
            # 模拟 mtDNA 复制过程中的突变
            g = 16568  # 基因组大小，不含 m.3107 处的间隔
            back_size = i * heteroplasmy  # 瓶颈后可能发生回突变的替换等位基因数量
            size = i - back_size  # 瓶颈后可能发生突变的参考等位基因数量
            
            # PGC 扩增形成初级卵母细胞：单个 PGC 将经历 18 次有丝分裂，在此过程中可能发生突变
            for(cycle in 1:18){
              
              # 复制阶段
              if(back_size == 0){  # 杂合度 = 0
                prob = g * j  # mtDNA 中发生突变的概率（碱基数 * 突变率）
                back_prob = 0
              } else {
                prob = j  # 发生突变的概率，仅针对该位点，即剩余参考等位基因再次发生相同突变
                back_prob = z * j  # 回突变的概率，z 为“正向”突变概率 j 的倍数
              }
              
              # rbinom 从二项分布中生成随机数
              # 正向突变
              new_mut = rbinom(1, size, prob)  # 新增突变等位基因，1 表示试验次数
              # 回突变
              back_mut = rbinom(1, back_size, back_prob)  # 新增参考等位基因，1 表示试验次数
              # 为下一次复制准备输入，size 为参考等位基因数，back_size 为替换等位基因数
              size = (2 * size) - new_mut + back_mut  # 复制后的参考等位基因数
              back_size = (2 * back_size) + new_mut - back_mut  # 复制后的突变等位基因数
              
              if((size + back_size) != (i*2)){
                stop("number copies after proliferation does not equal doubling of bottleneck")
              }
              
              # 分裂阶段
              nn = 1  # 观测次数
              m = back_size  # 替换等位基因数
              n = size  # 参考等位基因数
              k = i  # 分配给子细胞的 mtDNA 拷贝数
              
              # rhyper 从超几何分布中生成随机数
              # 细胞分裂后随机分配杂合度，rhyper 返回被选中的替换等位基因数
              # 其中 size 为参考等位基因数，back_size 为替换等位基因数
              back_size = rhyper(nn, m, n, k)  # 更新后的替换等位基因数
              size = i - back_size
              
              print(paste("heteroplasmy after replication + division for ", cycle, back_size / i))
              
              if((size + back_size) != i){
                stop("number copies after proliferation does not equal bottleneck")
              }
            }
            
            # PGC 形成初级卵母细胞后的杂合度
            heteroplasmy = back_size / i
            print(paste("heteroplasmy after proliferation", heteroplasmy))
            
            # 复制回 2^19 拷贝的次数取决于瓶颈大小 i，以形成成熟卵母细胞
            num_replications = ifelse(i == 8, 16, 
                                      ifelse(i == 16, 15, 
                                             ifelse(i == 32, 14, 
                                                    ifelse(i == 64, 13, 
                                                           ifelse(i == 128, 12, 
                                                                  ifelse(i==256, 11, 
                                                                         ifelse(i==512, 10, 
                                                                                print("error"))))))))
            
            for(replication in 1:num_replications){
              # 设定突变率
              if(back_size == 0){  # 杂合度 = 0
                prob = g * j  # mtDNA 中发生突变的概率（碱基数 * 突变率）
                back_prob = 0
              } else {
                prob = j  # 发生突变的概率，仅针对该位点，即剩余参考等位基因再次发生相同突变
                back_prob = z * j  # 回突变的概率，z 为“正向”突变概率 j 的倍数
              }
              
              # rbinom 从二项分布中生成随机数
              # 正向突变
              new_mut = rbinom(1, size, prob)  # 新增突变等位基因，1 表示试验次数
              # 回突变
              back_mut = rbinom(1, back_size, back_prob)  # 新增参考等位基因，1 表示试验次数
              # 为下一次复制准备输入，size 为参考等位基因数，back_size 为替换等位基因数
              size = (2 * size) - new_mut + back_mut  # 复制后的参考等位基因数
              back_size = (2 * back_size) + new_mut - back_mut  # 复制后的突变等位基因数
              
              print(paste("heteroplasmy after replication for ", replication, back_size / (size + back_size)))
            }
            
            if((size + back_size) != num_mt_zygote){
              stop("size after clonal amplification does not equal 2^19")
            }
            
            # 成熟卵母细胞的最终杂合度，用于下一代模拟
            heteroplasmy = back_size / num_mt_zygote  # 替换等位基因数量除以总 mtDNA 拷贝数
            print(paste("heteroplasmy after clonal amplification", heteroplasmy))
            
            # 写入结果文件
            # 参考的表头为 "individual\theteroplasmy\tgeneration\tmutation_rate\tbottleneck_size\tback_mutation_rate\tstarting_heteroplasmy\n"
            line = sprintf("%s\t%s\t%s\t%s\t%s\t%s\t%s\n", p, heteroplasmy, gen, j, i, z, start_het)
            cat(line, file = sim_file_path, append = TRUE, sep ="\n")
          }
        }
      }
    }
  }
}

# 设置并行计算环境
numCores <- detectCores()
registerDoParallel(numCores)
options(digits=22)  # 让 R 显示尽可能多的小数位


### 首先复现 Colgnahi 等人在图 2B 中报道的成熟卵母细胞数据

# 设定参数
mutation_rates = 10^-8  # 每个碱基的突变率
population_size = 10000  # 模拟人群的个体数
num_generations = 1  # 需要模拟的世代数
num_mt_bottleneck = list(8, 16, 32, 64, 128)  # 瓶颈阶段的 mtDNA 拷贝数
num_mt_zygote = 2^19  # 合子中的 mtDNA 拷贝数
back_mutation_rates = 1  # 回突变概率的倍数系数，z = 1 表示正反向突变概率相同
starting_heteroplasmy = 0.1  # 第 0 代的初始杂合度

# 输出结果文件
sim_file_path <-  "../output/蒙特卡洛模拟/Colnaghi等人复制实验数据文件.txt"
header = "individual\theteroplasmy\tgeneration\tmutation_rate\tbottleneck_size\tback_mutation_rate\tstarting_heteroplasmy\n"
cat(header, file = sim_file_path, append = FALSE, sep = "\t")

# 并行调用 mtDNA_simulate，对每个个体进行模拟
mclapply(1:population_size, FUN = mtDNA_simulate, mc.cores = numCores)  # 使用 mclapply 并行化



### 接下来在多组突变率范围内运行模拟

# 设定参数
mutation_rates = list(10^-7, 7.5 * 10^-8, 5 * 10^-8, 2.5 * 10^-8, 10^-8, 7.5 * 10^-9, 5 * 10^-9, 2.5 * 10^-9, 10^-9)
population_size = 10000
num_generations = 5
num_mt_bottleneck = 128
num_mt_zygote = 2^19
back_mutation_rates = 1
starting_heteroplasmy = 0 

# 输出结果文件
sim_file_path <-  "../output/蒙特卡洛模拟/蒙特卡洛模拟_结果.txt"
header = "individual\theteroplasmy\tgeneration\tmutation_rate\tbottleneck_size\tback_mutation_rate\tstarting_heteroplasmy\n"
cat(header, file = sim_file_path, append = FALSE, sep = "\t")

# 并行调用 mtDNA_simulate，对每个个体进行模拟
mclapply(1:population_size, FUN = mtDNA_simulate, mc.cores = numCores) 


# 使用自助法评估模拟数据中的最大杂合度 

# 自助法函数
sim_bootstrap <- function(mutation_rate, df) {  
  # 初始化数据框
  output <- data.frame(matrix(vector(), ncol = 2))
  colnames(output) <-c("mutation_rate", "max_heteroplasmy")
  output[1,] <- c(NA, NA)
  
  # 针对每个突变率执行 1000 次自助抽样，每次从总群体（10,000）中抽取 100 个样本
  for(k in 1:1000){
    j <- trimws(mutation_rate) 
    table <- df %>% filter(mutation_rate == !!j) %>%
      sample_n(size = 100, replace = T) # 从 10,000 个个体中抽样 100 个
    max_het <- max(table$heteroplasmy)
    row <- as.data.frame(t(c(j, max_het)))
    colnames(row) <- c("mutation_rate", "max_heteroplasmy")
    output <- rbind(output, row)
  }
  
  # 移除初始化的 NA 行
  output <- output[!is.na(output$mutation_rate),]
} 

results <- read.table("../output/蒙特卡洛模拟/蒙特卡洛模拟_结果.txt", header = TRUE)
gen_number = 5

output <- mclapply(mutation_rates, FUN = sim_bootstrap, df = results[results$generation == gen_number,], mc.cores = numCores) 
output <- do.call("rbind", output)

write.table(output, file = "../output/蒙特卡洛模拟/采样最大异质性.txt", row.names = FALSE, sep = "\t")

