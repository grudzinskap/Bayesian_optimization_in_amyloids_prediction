#real labels
pep424_tab <- read.table("./benchmark/pep424_evaluation.txt", sep = "\t")
levels(pep424_tab[[3]]) <- unlist(strsplit(levels(pep424_tab[[3]]), " "))

pep424_prots <- lapply(as.character(pep424_tab[[2]]) %>% strsplit(split = ""), function(i)
  i[i != " "])

real_labels <- (as.character(pep424_tab[[3]]) == "+")[lengths(pep424_prots) > 5]


#fold amyloid preds
foldAmyloid_data <- readLines("./benchmark/FoldAmyloid_pred.txt")
foldAmyloid_preds <- grepl("1", 
                           foldAmyloid_data[grep("Options: Scale = Expected number of contacts 8A", 
                                                 foldAmyloid_data) + 1])

#PASTA 2.0 preds
file_names <- list.files("./benchmark/pasta2_preds/")[grepl("fasta.seq.best_pairings_list.dat", list.files("./benchmark/pasta2_preds/"))]

pasta2_prot_id <- as.numeric(sapply(strsplit(file_names, ".", fixed = TRUE), function(i) substr(i[1], 5, nchar(i[1]))))

# pasta2_preds <- sapply(file_names, function(i) {
#   all_lines <- readLines(paste0("./benchmark/pasta2_preds/", i))
#   as.numeric(strsplit(strsplit(all_lines, "PASTA energy ")[[1]][2], "*[ ]")[[1]][1])
# }) > -5

pasta2_energy <- sapply(file_names, function(i) {
  all_lines <- readLines(paste0("./benchmark/pasta2_preds/", i))
  as.numeric(strsplit(strsplit(all_lines, "PASTA energy ")[[1]][2], "*[ ]")[[1]][1])
})

pasta2_preds <- 1 - (pasta2_energy - min(pasta2_energy))/max(pasta2_energy - min(pasta2_energy))

# pasta 2 cutoff: 0.6226317
1 - (-5 - min(pasta2_energy))/max(pasta2_energy - min(pasta2_energy))

#appnn preds

library(appnn)
appnn_preds <- read.fasta("./benchmark/pep424_better_names.fasta", seqtype = "AA") %>% 
  sapply(paste0, collapse = "") %>% 
  appnn %>% 
  sapply(function(i) i[["overall"]])


write.csv(data.frame(real_labels = real_labels,
                     PASTA2 = pasta2_preds[order(pasta2_prot_id)],
                     FoldAmyloid = foldAmyloid_preds,
                     appnn = appnn_preds), "./results/benchmark_otherpreds.csv",
          row.names = FALSE)
