source("./functions/choose_properties.R")
source("./functions/create_encodings.R")
source("./functions/encode_amyloids.R")
source("./functions/cv.R")
source("./functions/cv_analysis.R")
source("./functions/make_classifier.R")


require(seqinr)
require(dplyr)
require(pbapply)
require(biogram)
require(ranger)
require(hmeasure)
require(pbapply)

load("./data/aa_groups.RData")
aa_groups <- string2list(aa_groups)

raw_seqs_list <- c(read.fasta("./data/amyloid_pos_benchmark.fasta", seqtype = "AA"),
                   read.fasta("./data/amyloid_neg_benchmark.fasta", seqtype = "AA"))
#sequences longer than 5 aa and shorter than 26 aa
purified_seqs_id <- lengths(raw_seqs_list) > 5 & lengths(raw_seqs_list) < 16
seqs_list <- raw_seqs_list[purified_seqs_id]

seqs_m <- tolower(t(sapply(seqs_list, function(i)
  c(i, rep(NA, max(lengths(seqs_list)) - length(i))))))

ets <- c(rep(1, length(read.fasta("./data/amyloid_pos_benchmark.fasta", seqtype = "AA"))),
         rep(0, length(read.fasta("./data/amyloid_neg_benchmark.fasta", seqtype = "AA"))))
ets <- ets[purified_seqs_id]

seq_lengths <- unname(lengths(seqs_list))

test_dat <- read.fasta("./benchmark/pep424_better_names.fasta", seqtype = "AA")
test_dat_m <- tolower(t(sapply(test_dat, function(i)
  c(i, rep(NA, max(lengths(test_dat)) - length(i))))))

pasted_seqs <- sapply(seqs_list, paste0, collapse = "")
pasted_test <- sapply(test_dat, paste0, collapse = "")

# not in the test set
nit <- sapply(pasted_seqs, function(i) sum(grepl(i, pasted_test, fixed = TRUE))) == 0

# time around 40 [s]
# tmp <- make_classifier(seqs_m, ets, seq_lengths, 6, aa_groups[9], test_dat_m)
# system.time(make_classifier(seqs_m, ets, seq_lengths, 6, aa_groups[9], test_dat_m))

# Get list of classifiers

full_aa <- tolower(a()[-1]) %>% as.list
names(full_aa) <- 1L:20
full_aa <- list(full_aa)

for (seed_ in c(1, 5, 1995, 268855, 123)) {
  print(seed_)
  class_list <- pblapply(c(6, 10, 15), function(single_length)
    data.frame(class14592 = make_classifier(seqs_m, ets, seq_lengths, single_length, 
                                            aa_groups[14592], test_dat_m, seed_)
    )
  )
  
  learn_lengths <- c(6, 10, 15)
  dat <- cbind(read.csv("./results/benchmark_otherpreds.csv")['real_labels'],
               do.call(cbind, lapply(1L:3, function(i) {
                 single_preds_df <- class_list[[i]]
                 colnames(single_preds_df) <- paste0(colnames(single_preds_df), "_", learn_lengths[i])
                 single_preds_df
               })))
  
  bench_res <- HMeasure(dat[[1]], dat[-1], threshold = c(0.6226317, rep(0.5, ncol(dat) - 2)))[["metrics"]] %>%
    mutate(MCC = calc_mcc(TP, TN, FP, FN), classifier = rownames(.)) %>%
    select(classifier, AUC, MCC, Sens, Spec) %>%
    group_by(classifier) 
  
  csv_name <- paste("./results/benchmark_amyloGram_PG_seed_", seed_, ".csv", sep = "")
  write.csv(bench_res, csv_name, row.names = FALSE)
}
# mtry = 13, 16, 17

