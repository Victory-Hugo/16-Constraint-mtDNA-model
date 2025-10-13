# Save linear model fits, for use in calculating expected.

args <- commandArgs(trailingOnly = TRUE)

calibrate_step2 <- function(prefix = NA, do_sig_test = NA){
  # set defaults
  if(is.na(prefix)) {
    prefix = ""
  } 
  if(is.na(do_sig_test)) {
    do_sig_test = "yes"
  } 
  
  file <- read.delim(file = sprintf('output_files/calibration/%sloci_obs_vs_scores.txt', prefix), header = TRUE, sep = "\t")
  ori_file <- read.delim(file = sprintf('output_files/calibration/%sloci_obs_vs_scores_ori.txt', prefix), header = TRUE, sep = "\t")
  
  model_fits <- rbind(c("region", "mutation_group", "item", "value"),
                      c("ref_exc_ori", "G>A_and_T>C", "coefficient",  coef(lm(obs_max_het ~ sum_likelihood, data = file[file$mutation_group == "G>A_and_T>C",]))["sum_likelihood"]),
                      c("ref_exc_ori", "G>A_and_T>C", "intercept", coef(lm(obs_max_het ~ sum_likelihood, data = file[file$mutation_group == "G>A_and_T>C",]))["(Intercept)"]),
                      c("ref_exc_ori", "other", "coefficient", coef(lm(obs_max_het ~ sum_likelihood, data = file[file$mutation_group == "other",]))["sum_likelihood"]),
                      c("ref_exc_ori", "other", "intercept", coef(lm(obs_max_het ~ sum_likelihood, data = file[file$mutation_group == "other",]))["(Intercept)"]),
                      c("ori", "G>A_and_T>C", "coefficient", coef(lm(obs_max_het ~ sum_likelihood, data = ori_file[ori_file$mutation_group == "G>A_and_T>C",]))["sum_likelihood"]),
                      c("ori", "G>A_and_T>C", "intercept", coef(lm(obs_max_het ~ sum_likelihood, data = ori_file[ori_file$mutation_group == "G>A_and_T>C",]))["(Intercept)"]),
                      c("ori", "other", "coefficient", coef(lm(obs_max_het ~ sum_likelihood, data = ori_file[ori_file$mutation_group == "other",]))["sum_likelihood"]),
                      c("ori", "other", "intercept", coef(lm(obs_max_het ~ sum_likelihood, data = ori_file[ori_file$mutation_group == "other",]))["(Intercept)"]))
  
  write.table(model_fits, file = sprintf("output_files/calibration/%slinear_model_fits.txt", prefix), sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
  
  
  if (do_sig_test == "yes"){
    
    # inititate output file
    sink(sprintf("output_files/calibration/%scorrelation_significance_testing.txt", prefix), append = FALSE)
    
    # Test which variable is more significantly associated with observed neutral variation
    print("G>A and T>C mutation group (for reference excluding OriB-OriH)")
    print(summary(lm(obs_max_het ~ sum_likelihood + length, data = file[file$mutation_group == "G>A_and_T>C",])))
    print("All other mutation group (for reference excluding OriB-OriH)")
    print(summary(lm(obs_max_het ~ sum_likelihood + length, data = file[file$mutation_group == "other",])))
    print("G>A and T>C mutation group for OriB-OriH region")
    print(summary(lm(obs_max_het ~ sum_likelihood + length, data = ori_file[ori_file$mutation_group == "G>A_and_T>C",])))
    print("All other mutation group for OriB-OriH region")
    print(summary(lm(obs_max_het ~ sum_likelihood + length, data = ori_file[ori_file$mutation_group == "other",])))
    
    # Show that the observed level of neutral variation is more significantly correlated with mutation likelihood scores than locus length, using the cocor package. 
    # Testing two correlations which have dependent groups (ie they have one variable in common), and a one-sided hypothesis
    # cor.test is used to extract the correlation coefficient for testing
    
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
