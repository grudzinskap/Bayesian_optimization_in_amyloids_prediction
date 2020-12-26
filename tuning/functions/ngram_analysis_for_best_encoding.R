source("./functions/choose_properties.R")
source("./functions/create_encodings.R")
source("./functions/encode_amyloids.R")
source("./functions/cv.R")
source("./functions/cv_analysis.R")

require(seqinr)
require(dplyr)
require(pbapply)
require(biogram)
require(cvTools)
require(ranger)
require(hmeasure)
require(e1071)
library(reshape2)

if(Sys.info()["nodename"] %in% c("phobos", "michal-XPS14"))
  pathway <- "/home/michal/Dropbox/amyloid_cv_results/"

if(Sys.info()["nodename"] == "MICHALKOMP" )
  pathway <- "C:/Users/Michal/Dropbox/amyloid_cv_results/"

if(Sys.info()["nodename"] == "tobit" )
  pathway <- "~/Dropbox/amyloid_cv_results/"

amyloids <- read.csv("results/all_summaries.csv")
load("./data/aa_groups.RData")
aa_groups <- string2list(aa_groups)

#we select best encodings, by their ranks of AUC, for different lengths of sequences
best_positions <- amyloids %>% 
  select(len_range, enc_adj, AUC_mean) %>%
  group_by(len_range) %>%
  mutate(position = rank(AUC_mean)) %>%
  ungroup %>%
  group_by(enc_adj) %>%
  summarise(cum_rank = sum(position)) %>%
  arrange(desc(cum_rank)) %>%
  slice(1L:5)

best_enc <- best_positions[["enc_adj"]][1]

########### properties - best encodings ############

load("./results/enc_dupes.RData")
#remember, indices in enc_adj and are equal to indices in enc_dupes + 2 and trait_tab  + 2

best_dupes <- lapply(enc_dupes[["aa_duplicates"]][best_enc], function(single_enc)
  single_enc %>%
    sapply(function(i) substr(i, 3, nchar(i))) %>%
    strsplit("K") %>%
    lapply(as.numeric) %>%
    do.call(rbind, .) %>%
    data.frame %>%
    select(X1) %>%
    unlist(use.names = FALSE) - 2)

frequencies <- data.frame(trait_tab[unlist(best_dupes, use.names = FALSE), ] %>%
                            apply(1, na.omit) %>% unlist %>% table)
colnames(frequencies) <- c("X", "Freq")
frequencies$X <- as.numeric(as.character(frequencies$X))

trait_props_all <- prop_MK %>%
  inner_join(frequencies) %>%
  arrange(desc(Freq)) %>%
  select(name, Freq) %>%
  mutate(prop = Freq/length(unlist(best_dupes, use.names = FALSE))) %>%
  droplevels()

trait_props_all[["name"]] <- factor(trait_props_all[["name"]],
                                    levels = rev(as.character(trait_props_all[["name"]])))

levels(trait_props_all[["name"]]) <- sapply(strwrap(levels(trait_props_all[["name"]]), 30, simplify = FALSE),
                                            function(i) paste0(i, collapse = "\n"))

ggplot(trait_props_all, aes(x = name, y = prop)) +
  geom_bar(stat = "identity") +
  coord_flip()

####### ngrams for the best encoding ############

partition_enc_id <- filter(amyloids, enc_adj %in% best_enc) %>%
  select(partition, enc, enc_adj) %>%
  filter(!duplicated(.)) %>%
  slice(1) %>%
  select(partition, enc) %>%
  unlist


ngrams_best_enc <- lapply(1L:45, function(single_replicate_id) {
  load(paste0(pathway, "results/cv_results_full_", single_replicate_id, "_", 
              partition_enc_id[["partition"]],".Rdata"))
  
  sapply(cv_results[[partition_enc_id[["enc"]]]][[1]], function(single_fold)
    single_fold[[2]]) 
}) %>% 
  unlist %>% 
  table(dnn = "ngram") %>%
  data.frame %>%
  filter(Freq == 225) %>%
  select(ngram) %>%
  unlist %>%
  as.character


load("./data/aa_groups.RData")
aa_groups <- string2list(aa_groups)

raw_seqs_list <- c(read.fasta("./data/amyloid_pos_full.fasta",seqtype = "AA"),
                   read.fasta("./data/amyloid_neg_full.fasta",seqtype = "AA"))
#sequences longer than 5 aa and shorter than 26 aa
purified_seqs_id <- lengths(raw_seqs_list) > 5 & lengths(raw_seqs_list) < 26

ets <- c(rep(1, length(read.fasta("./data/amyloid_pos_full.fasta",seqtype = "AA"))),
         rep(0, length(read.fasta("./data/amyloid_neg_full.fasta",seqtype = "AA"))))
ets <- ets[purified_seqs_id]

seqs_list <- raw_seqs_list[purified_seqs_id]

seqs_m <- tolower(t(sapply(seqs_list, function(i)
  c(i, rep(NA, max(lengths(seqs_list)) - length(i))))))


extracted_ngrams <- extract_ngrams(seqs_m, aa_groups[[best_enc]])

seqs_deg <- t(apply(seqs_m, 1, function(seq) degenerate(seq, aa_groups[[best_enc]])))

positive_ngrams_count_best_encoding <- count_specified(seqs_deg[which(ets==1),], ngrams_best_enc) 
negative_ngrams_count_best_encoding <- count_specified(seqs_deg[which(ets==0),], ngrams_best_enc) 

dane=rbind(
  data.frame(v = positive_ngrams_count_best_encoding$v, 
             j = positive_ngrams_count_best_encoding$j,
             type = "pos", 
             nseq = 397),
  data.frame(v = negative_ngrams_count_best_encoding$v, 
             j = negative_ngrams_count_best_encoding$j,
             type = "neg", 
             nseq = 1033)
)

save(dane, file = "./results/imp_ngrams_counts.RData")

ngram_count_df <- dane
# lengths of n-grams
imp_ngram_len <- nchar(decode_ngrams(ngrams_best_enc))

# protein lengths
pos_prot_len <- apply(seqs_deg[which(ets==1),], 1, function(i) sum(!is.na(i)))
neg_prot_len <- apply(seqs_deg[which(ets==0),], 1, function(i) sum(!is.na(i))) 

total_ngram_number <- c(sapply(imp_ngram_len, function(i)
  pos_prot_len - i + 1) %>% 
    colSums,
  sapply(imp_ngram_len, function(i)
    neg_prot_len - i + 1) %>% 
    colSums
)

ngram_freq_final_tab <- dane %>% 
  group_by(type,j) %>%
  summarise(count_raw=sum(v>0)) %>% 
  ungroup %>% 
  mutate(ngram_number = total_ngram_number,
         count = count_raw/ngram_number) %>%
  dcast(j ~ type) %>%
  mutate(name=ngrams_best_enc,
         decoded_name=decode_ngrams(name),
         diff_freq = pos - neg) %>%
  arrange(diff_freq*pos) %>%
  inner_join(relative_frequency_p_vals) %>%
  select(decoded_name, diff_freq, pval, pos, neg)


# n-gram patterns from literature -----------------------------------------
# "Sequence determinants of amyloid fibril formation"
# http://www.pnas.org/content/101/1/87.full

# write output file


amyl_raw <- list(X1 = character(0), 
                 X2 = character(0), 
                 X3 = c("v", "l", "s", "c", "w", "f", "n", "q", "e"), 
                 X4 = c("i", "l", "t", "y", "w", "f", "n", "e"), 
                 X5 = c("f", "i", "y"), 
                 X6 = character(0))


namyl_raw <- list(X1 = "p", 
                  X2 = c("p", "k", "r", "h", "w"), 
                  X3 = character(0), 
                  X4 = character(0), 
                  X5 = character(0), 
                  X6 = c("p", "k", "r", "h"))

degenerate_pattern <- function(x)
  lapply(x, function(tmp_seq) {
    for (i in 1L:length(aa_groups[[best_enc]])) {
      tmp_seq[tmp_seq %in% aa_groups[[best_enc]][[i]]] <- names(aa_groups[[best_enc]])[i]
    }
    
    res <- sort(unique(tmp_seq))
    if(length(res) == 0) {
      "|"
    } else {
      res
    }
  })

amyl_deg <- as.matrix(expand.grid(degenerate_pattern(amyl_raw)))
namyl_deg <- as.matrix(expand.grid(degenerate_pattern(namyl_raw)))

remove_redundant_space <- function(ngrams) 
  unlist(lapply(ngrams, function(x) {
    if(grepl(pattern = "\\_?\\|\\_?", x)) {
      sub("_*$", "", sub("^_*", "", gsub("|", "", x, fixed = TRUE)))
    } else {
      x
    }
  }))

get_patterns <- function(ns, ds, seqs) {
  lapply(1L:length(ns), function(i) {
    seq2ngrams(seqs, ns[[i]], c(1L:6, "|"), ds[[i]]) %>% 
      decode_ngrams %>% 
      unique %>% 
      remove_redundant_space %>% 
      unique
  })
}

ns = c(1, rep(2, 4), rep(3, 3))
ds = list(0, 0, 1, 2, 3, c(0, 0), c(0, 1), c(1, 0))

amyl_ngrams <- get_patterns(ns, ds, amyl_deg) %>% unlist %>% unique
namyl_ngrams <- get_patterns(ns, ds, namyl_deg) %>% unlist %>% unique

ngram_freq_final_tab[["association"]] <- ((ngram_freq[["decoded_name"]] %in% amyl_ngrams) +
                                            (ngram_freq[["decoded_name"]] %in% namyl_ngrams) * -1) %>% 
  factor(labels = c("Non-amyloidogenic", "Not found", "Amyloidogenic"))

write.csv(ngram_freq_final_tab, file = "./results/ngram_freq.csv", row.names = FALSE)
