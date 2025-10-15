library(doParallel)
library(dplyr)
library(parallel)
library(stats)
library(tidyr)


# simulation to show correlation between mutation rates and maximum heteroplasmy

dir.create("../output_files/simulation/")

# function to simulate heteroplasmy across generations (oogenesis)
# requires input parameters, which are user supplied
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
            
            # bottleneck, formation of a primordial germ cell (PGC) from the zygote
            if(heteroplasmy > 0){
              # step 1 - bottleneck, only apply if a mutation has occured - as otherwise heteroplasmy level remains to be 0
              # mtDNA copy number is halved each cell division, without mtDNA replication
              # the number of cell divisions to model depends on the final bottleneck size (i)
              num_divisions = ifelse(i == 8, 16, 
                                     ifelse(i == 16, 15, 
                                            ifelse(i == 32, 14, 
                                                   ifelse(i == 64, 13, 
                                                          ifelse(i == 128, 12, 
                                                                 ifelse(i==256, 11, 
                                                                        ifelse(i==512, 10, 
                                                                               print("error"))))))))
              # set the number of mtDNA copies in the starting zygote
              k = num_mt_zygote
              
              # model the random segregation of mtDNA copies into daughter cells
              for(division in 1:num_divisions){
                nn = 1  # the number of observations
                m = k * heteroplasmy  # number of alternate alleles
                n = k - m  # number of reference alleles
                k = k / 2  # number of mtDNA copies to partion to daughter cell
                
                # rhyper simulates random number from the hypergeometric distribution
                # random heteroplasmy assignment after cell division, rhyper is number of alternate alleles selected
                heteroplasmy = rhyper(nn, m, n, k) / k
                
                print(paste("heteroplasmy after bottleneck division ", division, "is", heteroplasmy))
              }
              
              if(k != i){
                stop("bottleneck size after cell divisions does not equal bottleneck")
              }
            } 
            print(paste("heteroplasmy after bottleneck", heteroplasmy))
            
            # model mutation during mtDNA replication
            g = 16568  # genome size, excluding the spacer at m.3107
            back_size = i * heteroplasmy  # number alternate alleles after bottleneck that could back mutate
            size = i - back_size  # number reference alleles after bottleneck that could mutate
            
            # PGC proliferation to form a primary oocyte - a single PGC will undergo 18 cell divisions by mitosis, during which mutation can happen
            for(cycle in 1:18){
              
              # REPLICATION
              if(back_size == 0){  # heteroplasmy = 0
                prob = g * j  # probability of mutation in the mtDNA (number of base pairs * mutation rate)
                back_prob = 0
              } else {
                prob = j  # probability of mutation, specific to the position with the mutation, ie remaining reference alleles having the same mutation occur
                back_prob = z * j  # probability of back mutation, where z is a multiplier of the 'forward' probability j
              }
              
              # rbinom simulate random number from the binomial distribution
              # forward mutation
              new_mut = rbinom(1, size, prob)  # new mutant alleles, note the 1 represents the number of observations
              # back mutation
              back_mut = rbinom(1, back_size, back_prob)  # new reference alleles, note the 1 represents the number of observations
              # inputs for next replication, size is the number of reference alleles, and back_size is the number of alternate alleles
              size = (2 * size) - new_mut + back_mut  # number reference alleles after replication
              back_size = (2 * back_size) + new_mut - back_mut  # number mutant alleles after replication
              
              if((size + back_size) != (i*2)){
                stop("number copies after proliferation does not equal doubling of bottleneck")
              }
              
              # DIVISION
              nn = 1  # the number of observations
              m = back_size  # number of alternate alleles
              n = size  # number of reference alleles
              k = i  # number of mtDNA copies to partion to daughter cell 
              
              # rhyper simulates random number from the hypergeometric distribution
              # random heteroplasmy assignment after cell division, rhyper is number of alternate alleles selected
              # size is the number of reference alleles, and back_size is the number of alternate alleles
              back_size = rhyper(nn, m, n, k)  # new number of alternate alleles 
              size = i - back_size
              
              print(paste("heteroplasmy after replication + division for ", cycle, back_size / i))
              
              if((size + back_size) != i){
                stop("number copies after proliferation does not equal bottleneck")
              }
            }
            
            # heteroplasmy after formation of a primary oocyte from a PGC
            heteroplasmy = back_size / i
            print(paste("heteroplasmy after proliferation", heteroplasmy))
            
            # per replication event back to 2^19, mature oocyte
            # the number of replication events needed depends on the bottleneck size (i)
            num_replications = ifelse(i == 8, 16, 
                                      ifelse(i == 16, 15, 
                                             ifelse(i == 32, 14, 
                                                    ifelse(i == 64, 13, 
                                                           ifelse(i == 128, 12, 
                                                                  ifelse(i==256, 11, 
                                                                         ifelse(i==512, 10, 
                                                                                print("error"))))))))
            
            for(replication in 1:num_replications){
              # set mutation rate
              if(back_size == 0){  # heteroplasmy = 0
                prob = g * j  # probability of mutation in the mtDNA (number of base pairs * mutation rate)
                back_prob = 0
              } else {
                prob = j  # probability of mutation - specific to the position with the mutation, ie remaining reference alleles having the same mutation occur
                back_prob = z * j  # probability of back mutation, where z is a multiplier of the 'forward' probability j
              }
              
              # rbinom simulate random number from the binomial distribution
              # forward mutation
              new_mut = rbinom(1, size, prob)  # new mutant alleles, note the 1 represents the number of observations
              # back mutation
              back_mut = rbinom(1, back_size, back_prob)  # new reference alleles, note the 1 represents the number of observations
              # inputs for next replication, size is the number of reference alleles, and back_size is the number of alternate alleles
              size = (2 * size) - new_mut + back_mut  # number reference alleles after replication
              back_size = (2 * back_size) + new_mut - back_mut  # number mutant alleles after replication
              
              print(paste("heteroplasmy after replication for ", replication, back_size / (size + back_size)))
            }
            
            if((size + back_size) != num_mt_zygote){
              stop("size after clonal amplification does not equal 2^19")
            }
            
            # final heteroplasmy in the mature oocyte, that will be used for the next generation
            heteroplasmy = back_size / num_mt_zygote  # number alternate alleles divided by total mtDNA copies
            print(paste("heteroplasmy after clonal amplification", heteroplasmy))
            
            # append to file
            # header for reference is "individual\theteroplasmy\tgeneration\tmutation_rate\tbottleneck_size\tback_mutation_rate\tstarting_heteroplasmy\n"
            line = sprintf("%s\t%s\t%s\t%s\t%s\t%s\t%s\n", p, heteroplasmy, gen, j, i, z, start_het)
            cat(line, file = sim_file_path, append = TRUE, sep ="\n")
          }
        }
      }
    }
  }
}

# set up for parallelization
numCores <- detectCores()
registerDoParallel(numCores)
options(digits=22)  # so R displays maximum decimal places


### first replicate data from Colgnahi et al, Figure 2B for mature oocytes

# establish parameters
mutation_rates = 10^-8  # mutation rate per bp
population_size = 10000  # number of individuals in population
num_generations = 1  # number of generations to model across
num_mt_bottleneck = list(8, 16, 32, 64, 128)  # number of mtDNA copies in the bottleneck
num_mt_zygote = 2^19  # number of mtDNA copies in the zygote
back_mutation_rates = 1  # a multiplier to apply to the 'forward' mutation rate, ie z = 1 applies the same back and forward mutation rates
starting_heteroplasmy = 0.1  # heteroplasmy level in generation 0

# file to write results to
sim_file_path <-  "../output_files/simulation/colnaghi_etal_replication.txt"
header = "individual\theteroplasmy\tgeneration\tmutation_rate\tbottleneck_size\tback_mutation_rate\tstarting_heteroplasmy\n"
cat(header, file = sim_file_path, append = FALSE, sep = "\t")

# apply the function mtDNA_simulate to each element of population size in parallel, number define by numCores
mclapply(1:population_size, FUN = mtDNA_simulate, mc.cores = numCores)  # mclapply to parallize



### now run our simulation for range of mutation rates

# establish parameters
mutation_rates = list(10^-7, 7.5 * 10^-8, 5 * 10^-8, 2.5 * 10^-8, 10^-8, 7.5 * 10^-9, 5 * 10^-9, 2.5 * 10^-9, 10^-9)
population_size = 10000
num_generations = 5
num_mt_bottleneck = 128
num_mt_zygote = 2^19
back_mutation_rates = 1
starting_heteroplasmy = 0 

# file to write results to
sim_file_path <-  "../output_files/simulation/simulation_results.txt"
header = "individual\theteroplasmy\tgeneration\tmutation_rate\tbottleneck_size\tback_mutation_rate\tstarting_heteroplasmy\n"
cat(header, file = sim_file_path, append = FALSE, sep = "\t")

# apply the function mtDNA_simulate to each element of population size in parallel, number define by numCores
mclapply(1:population_size, FUN = mtDNA_simulate, mc.cores = numCores) 


# using bootstrap to assess maximum heteroplasmy in simulated data 

# bootstrap function
sim_bootstrap <- function(mutation_rate, df) {  
  # initialize the dataframe
  output <- data.frame(matrix(vector(), ncol = 2))
  colnames(output) <-c("mutation_rate", "max_heteroplasmy")
  output[1,] <- c(NA, NA)
  
  # for each mutation rate, do 1000 bootstrap replicates sampling 100 from the total population (10,000)
  for(k in 1:1000){
    j <- trimws(mutation_rate) 
    table <- df %>% filter(mutation_rate == !!j) %>%
      sample_n(size = 100, replace = T) # sampling 100 from total 10,000
    max_het <- max(table$heteroplasmy)
    row <- as.data.frame(t(c(j, max_het)))
    colnames(row) <- c("mutation_rate", "max_heteroplasmy")
    output <- rbind(output, row)
  }
  
  # remove the initialized NA row
  output <- output[!is.na(output$mutation_rate),]
} 

results <- read.table("../output_files/simulation/simulation_results.txt", header = TRUE)
gen_number = 5

output <- mclapply(mutation_rates, FUN = sim_bootstrap, df = results[results$generation == gen_number,], mc.cores = numCores) 
output <- do.call("rbind", output)

write.table(output, file = "../output_files/simulation/sampling_max_heteroplasmy.txt", row.names = FALSE, sep = "\t")

