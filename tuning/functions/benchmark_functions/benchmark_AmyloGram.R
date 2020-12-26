set.seed(1)

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

 class_list <- pblapply(c(6, 10, 15), function(single_length)
   data.frame(class14592 = make_classifier(seqs_m, ets, seq_lengths, single_length, aa_groups[14592], test_dat_m),
              class14596 = make_classifier(seqs_m, ets, seq_lengths, single_length, aa_groups[14596], test_dat_m),
              class14533 = make_classifier(seqs_m, ets, seq_lengths, single_length, aa_groups[14533], test_dat_m),
              class18297 = make_classifier(seqs_m, ets, seq_lengths, single_length, aa_groups[18297], test_dat_m),
              class6704 = make_classifier(seqs_m, ets, seq_lengths, single_length, aa_groups[6704], test_dat_m),
              raw_aa = make_classifier(seqs_m, ets, seq_lengths, single_length, full_aa[1], test_dat_m)
   )
 )
# save(class_list, file = "./results/class_list_best5.RData")


raw_benchmark <- pblapply(c(6, 10, 15), function(single_length)
  data.frame(raw_aa = make_classifier(seqs_m, ets, seq_lengths, single_length, full_aa[1], test_dat_m)
   )
 )

# save(raw_benchmark, file = "./results/class_raw_benchmark.RData")
learn_lengths <- c(6, 10, 15)
# load("./results/class_list_best5.RData")
# load("./results/class_raw_benchmark.RData")

dat <- cbind(read.csv("./results/benchmark_otherpreds.csv"),
             do.call(cbind, lapply(1L:3, function(i) {
               single_preds_df <- class_list[[i]]
               colnames(single_preds_df) <- paste0(colnames(single_preds_df), "_", learn_lengths[i])
               single_preds_df
             })),
             do.call(cbind, lapply(1L:3, function(i) {
               single_preds_df <- raw_benchmark[[i]]
               colnames(single_preds_df) <- paste0(colnames(single_preds_df), "_", learn_lengths[i])
               single_preds_df
             })))

#write.csv(dat, "./results/benchmark_allpreds_raw.csv", row.names = FALSE)

get_last <- function(x, what) {
  splitted <- strsplit(as.character(x), what)[[1]]
  splitted[length(splitted)]
}

bench_res <- HMeasure(dat[[1]], dat[-1], threshold = c(0.6226317, rep(0.5, ncol(dat) - 2)))[["metrics"]] %>%
  mutate(MCC = calc_mcc(TP, TN, FP, FN), classifier = rownames(.)) %>%
  select(classifier, AUC, MCC, Sens, Spec) %>%
  group_by(classifier) %>%
  mutate(pos = get_last(classifier, what = "_"),
         nice_name = strsplit(as.character(classifier), "_")[[1]][1]) %>%
  mutate(nice_name = get_last(nice_name, "class")) %>%
  mutate(nice_name = ifelse(nice_name == "raw", "full alphabet", nice_name)) %>%
  ungroup %>% 
  mutate(pos = as.numeric(pos),
         nice_name = factor(nice_name, 
                            levels = c("PASTA2", "FoldAmyloid", "appnn", "14592", "14596", "14533", "18297", 
                                       "full alphabet", "16548"))) %>%
  arrange(nice_name) %>%
  slice(1L:18)

write.csv(bench_res, "./results/benchmark_allpreds_PG_seed1.csv", row.names = FALSE)
