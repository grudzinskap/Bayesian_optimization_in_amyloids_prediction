# library(Biostrings)

remove_dash <- function(x)
  x[x != "-"]

seq2bio <- function(x)
  Biostrings::AAString(paste0(x, collapse = ""))

source("./functions/choose_properties.R")
source("./functions/create_encodings.R")
source("./functions/encode_amyloids.R")
source("./functions/cv.R")
source("./functions/cv_analysis.R")
source("./functions/make_classifier.R")

compare_seq <- function(shorter, longer) {
  windows <- seq2ngrams(as.vector(longer), n = length(shorter), u = a()[-1]) %>% 
    decode_ngrams %>% 
    strsplit(split = "")
  
  sapply(windows, function(single_window) {
    unname(mean(shorter == single_window))
  }) %>% max
}

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

# all_prots <- c(read.fasta("data/amyloid_neg_benchmark.fasta", seqtype = "AA"), 
#                read.fasta("data/amyloid_pos_benchmark.fasta", seqtype = "AA")) %>% 
#   lapply(toupper) %>% 
#   lapply(remove_dash) %>% 
#   lapply(seq2bio)
# 
# pep424 <- read.fasta("benchmark/pep424_better_names.fasta", seqtype = "AA") %>% 
#   lapply(seq2bio)
# 
# pep424_pid <- lapply(pep424, function(i) 
#   pbsapply(all_prots, function(j) 
#     Biostrings::pid(Biostrings::pairwiseAlignment(i, j, type="global"), type = "PID1")
#   )
# )
# 
# save(pep424_pid, file = "/results/pep424_pid.RData")

# alns <- pbsapply(combn(1L:length(seqs_list), 2, simplify = FALSE), function(i) {
#   seq1 <- seq2bio(seqs_list[[i[[1]]]])
#   seq2 <- seq2bio(seqs_list[[i[[2]]]])
#   aln <- Biostrings::pairwiseAlignment(seq1, seq2, type="global", substitutionMatrix = NULL)
#   slot(aln, "score")
# })
# 
# save(alns, file = "./results/alns.RData")

set.seed(5)
for(i in 1L:15) {
  perm_alns <- pbsapply(combn(1L:length(seqs_list), 2, simplify = FALSE), function(i) {
    seq1 <- seq2bio(sample(seqs_list[[i[[1]]]]))
    seq2 <- seq2bio(sample(seqs_list[[i[[2]]]]))
    aln <- Biostrings::pairwiseAlignment(seq1, seq2, type="global", substitutionMatrix = NULL)
    slot(aln, "score")
  })
  cat(c(perm_alns, "\n"), sep = ";", file = "./results/perm_alns.txt", append = TRUE)
}
