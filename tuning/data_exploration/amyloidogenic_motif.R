#incorporated into ngram_analysis_for_best_encoding.R

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

ngram_freq <- read.csv("./results/ngram_freq.csv")

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

((ngram_freq[["decoded_name"]] %in% amyl_ngrams) +
  (ngram_freq[["decoded_name"]] %in% namyl_ngrams) * -1) %>% 
  factor()
