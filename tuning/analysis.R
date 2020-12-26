# code below creates and cross-validates classifiers for all reduced amino acid alphabets

source("./functions/choose_properties.R")
source("./functions/create_encodings.R")
source("./functions/encode_amyloids.R")
source("./functions/cv.R")
source("./functions/cv_analysis.R")

require(seqinr)
require(dplyr)

# require(biogram, lib.loc = "~/R/library")
# require(cvTools, lib.loc = "~/R/library")
# require(ranger, lib.loc = "~/R/library")
# require(hmeasure, lib.loc = "~/R/library")
# 
# 
# args <- commandArgs(trailingOnly = TRUE)
# idx <- as.numeric(args[1])
# idx2 <- as.numeric(args[2])

require(pbapply)
require(biogram)
require(cvTools)
require(ranger)
require(hmeasure)
require(e1071)

if(Sys.info()["nodename"] %in% c("phobos", "michal-XPS14"))
  pathway <- "/home/michal/Dropbox/amyloid_cv_results/"

if(Sys.info()["nodename"] == "MICHALKOMP" )
  pathway <- "C:/Users/Michal/Dropbox/amyloid_cv_results/"

ftraits <- choose_properties()

# 131071 combinations for each length of reduced alphabet: 524284 total

#time-consuming step commented out. Instead just load the results.
#aa_groups <- create_encodings(ftraits)
#save(aa_groups, file = "aa_groups.RData")
load("./data/aa_groups.RData")
aa_groups <- string2list(aa_groups)

raw_seqs_list <- c(read.fasta("./data/amyloid_pos_full.fasta",seqtype = "AA"),
               read.fasta("./data/amyloid_neg_full.fasta",seqtype = "AA"))
#sequences longer than 5 aa and shorter than 26 aa
purified_seqs_id <- lengths(raw_seqs_list) > 5 & lengths(raw_seqs_list) < 26
seqs_list <- raw_seqs_list[purified_seqs_id]

#save(purified_seqs_id, seqs_list, file = "./results/purified_seqs.RData")

seqs_m <- tolower(t(sapply(seqs_list, function(i)
  c(i, rep(NA, max(lengths(seqs_list)) - length(i))))))


ids <- split(1:length(aa_groups), rep(1L:ntiles, each=length(aa_groups)/ntiles+1))
#for single aa_group around 85 sec
id <- "ID9711K6"
extracted_ngrams <- extract_ngrams(seqs_m, aa_groups[id])

#extracted_ngrams <- extract_ngrams(seqs_m, aa_groups[ids[[idx2]]])

for (i in names(aa_groups)){
  first <- all(unlist(aa_groups[[i]]['1']) == unlist(list('g')))
  second <- all(unlist(aa_groups[[i]]['2']) == unlist(list('k', 'p', 'r')))
  third <- all(unlist(aa_groups[[i]]['3']) == unlist(list('i', 'l', 'v')))
  if (all(c(first, second, third))) {
    print(i)
    print(aa_groups[[i]])
    print(' ')
  }
}


#for single aa_group around 85 sec
#extracted_ngrams <- extract_ngrams(seqs_m, aa_groups)
#add to extracted_ngrams raw aa n-grams

extracted_ngrams[["ID9711K6"]][["dimnames"]]

#create folds for cv
ets <- c(rep(1, length(read.fasta("./data/amyloid_pos_full.fasta",seqtype = "AA"))),
         rep(0, length(read.fasta("./data/amyloid_neg_full.fasta",seqtype = "AA"))))
ets <- ets[purified_seqs_id]
seq_lengths <- unname(lengths(seqs_list))

all_folds <- create_all_folds(ets, seq_lengths)

cv_results <- do_cv(all_folds[idx], extracted_ngrams, create_hv(seqs_m))
filename <- paste0("./results/cv_results_full_", idx, "_", idx2, ".Rdata")
print(filename)
save(cv_results, file = filename)
print("success")

cv_results <- do_cv(all_folds, extracted_ngrams, create_hv(seqs_m))
cv_summary <- summarize_cv_results(cv_results, reps_ids = get_reps_ids())

#cv_summary <- summarize_cv_results(cv_results, get_reps_ids)

# negative control - full amino acid alphabet

# full_aa <- tolower(a()[-1]) %>% as.list
# names(full_aa) <- 1L:20
# full_aa <- list(full_aa)
# extracted_ngrams_full <- extract_ngrams(seqs_m, full_aa)
# cv_results_full <- do_cv(all_folds, extracted_ngrams_full, create_hv(seqs_m))
# save(cv_results_full, file = "./results/cv_results_full.RData")

# alternative alphabets from literature

# stand_enc <- c(string2list(gsub("-", "_", tolower("ADEGHKNPQRST-C-FY-ILMV-W"))),
#   string2list(gsub("-", "_", tolower("AG-C-DEKNPQRST-FILMVWY-H"))))
# 
# extracted_ngrams_stand <- extract_ngrams(seqs_m, stand_enc)
# cv_results_stand <- do_cv(all_folds, extracted_ngrams_stand, create_hv(seqs_m))
# save(cv_results_stand, file = "./results/cv_results_stand.RData")

# compute bootstrap significance for benchmark results
source("functions/benchmark_significance.R")

# compute pairwise identity between proteins
source("functions/compute_pid.R")
