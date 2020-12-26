source("./functions/choose_properties.R")
source("./functions/create_encodings.R")
source("./functions/encode_amyloids.R")


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

ets <- c(rep(1, length(read.fasta("./data/amyloid_pos_full.fasta",seqtype = "AA"))),
         rep(0, length(read.fasta("./data/amyloid_neg_full.fasta",seqtype = "AA"))))
ets <- ets[purified_seqs_id]

load("./results/alns.RData")

perm_alns <- do.call(cbind, unlist(lapply(c("./results/perm_alns5.txt",
                        "./results/perm_alns.txt"), function(i)
  lapply(strsplit(readLines(i), ";"), as.numeric)), recursive = FALSE))

p_scores <- unlist(lapply(1L:length(alns), function(i)
  mean(perm_alns[i, ] >= alns[i])))

table(p_scores)
mean(p_scores > 0.8)

# what do we align? 
# 0:non-amyloid/non-amyloid
# 1: amyloid/non-amyloid
# 2: amyloid/amyloid, 

aln_type <- unlist(lapply(combn(1L:length(seqs_list), 2, simplify = FALSE), 
              function(single_combn) sum(ets[single_combn])))

p_score_df <- data.frame(p_score = p_scores, type = aln_type)

library(ggplot2)

mutate(p_score_df, p_score = p_score == 0) %>% 
  group_by(type) %>% 
  summarise(mean(p_score))

ggplot(p_score_df, aes(x = p_score, fill = factor(type))) +
  geom_density() +
  facet_wrap(~ type)

