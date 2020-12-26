library(dplyr)
source("./functions/cv_analysis.R")

results <- data.frame()

for(fold in 1:45){
  for(tile in 1:20){
    filename <- paste0("./results/cv_results_full_", fold, "_", tile, ".Rdata")
    load(filename)
    results <- rbind(results,
                     summarize_cv_results(cv_results, get_reps_ids))
    print(paste("done", fold, tile))
  }
  
}

outFileName <- "summarized_results.Rdata"
save(results, file=outFileName)


