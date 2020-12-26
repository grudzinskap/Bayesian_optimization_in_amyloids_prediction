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

jackknife_res <- pblapply(1L:nrow(seqs_m), function(seq_id) {
  class1 <- make_classifier(seqs_m[-seq_id, ], 
                            ets[-seq_id], seq_lengths[-seq_id], 
                            10, aa_groups[14592], seqs_m[seq_id, , drop = FALSE])
  c(pred = class1, et = ets[seq_id])
})

save(jackknife_res, file = "./results/jackknife_res.RData")


class_full <- make_classifier(seqs_m, ets, seq_lengths, 10, aa_groups[14592], seqs_m)

load("./results/jackknife_res.RData")

jackknife_res
jackknife_df <- data.frame(cbind(do.call(rbind, jackknife_res), normal = class_full))

HMeasure(jackknife_df[, "et"], jackknife_df[, -2])[["metrics"]] %>% 
  data.frame() %>% 
  mutate(MCC = calc_mcc(as.numeric(TP), TN, FP, FN)) %>% 
  select(AUC, MCC)
