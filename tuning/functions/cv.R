#' Create vector of hexamer names
#'
#' Creates a vector of hexamer indices required for cross-validation.
#' @param seq matrix of sequences
#' 
#' @return Vector of peptide names.
create_hv <- function(seq) 
  unlist(lapply(1L:nrow(seq), function(i) {
    res <- do.call(rbind, strsplit(decode_ngrams(seq2ngrams(seq[i, ][!is.na(seq[i, ])], 6, a()[-1])), ""))
    rep(i, nrow(res))
  }))


#' Create folds for cv
#' 
#' Returns indices of folds for hexagrams used in cross-validation.
#' 
#' @param ets a vector of etiquettes (0 or 1).
#' @param seq_lengths a vector of lengths of sequences (must be provided, since we do not use
#' exact sequences and calculating length from n-gram is inconvenient)).
#' 
#' @return List of 90 lists. Each element contains four data frames containing indices for
#' amyloids in training set, nonamyloids in training set, amyloids in test set and 
#' nonamyloidsin test set.

create_all_folds <- function(ets, seq_lengths) {
  seq_label <- cut(seq_lengths, breaks = c(5, 6, 10, 15, 25))
  
  #set the seed to create the same folds
  set.seed(1)
  lapply(1L:3, function(constant_pos) {
    lapply(1L:3, function(constant_neg) {

      lapply(1L:5, function(dummy) {
        splitted_seqs <- lapply(levels(seq_label), function(single_label) {
          pos_seqs <- which(ets == 1 & seq_label == single_label)
          neg_seqs <- which(ets == 0 & seq_label == single_label)
          
          fold_list <- lapply(list(pos = pos_seqs, neg = neg_seqs), function(single_n) {
            folded <- cvFolds(length(single_n), K = 5)
            data.frame(id = single_n[folded[["subsets"]]], which = folded[["which"]])
          })
        })
        
        list(pos_train = do.call(rbind, lapply(splitted_seqs[1L:constant_pos], 
                                               function(single_split) single_split[["pos"]])),
             pos_test = do.call(rbind, lapply(splitted_seqs[-(1L:constant_pos)], 
                                              function(single_split) single_split[["pos"]])),
             neg_train = do.call(rbind, lapply(splitted_seqs[1L:constant_pos], 
                                               function(single_split) single_split[["neg"]])),
             neg_test = do.call(rbind, lapply(splitted_seqs[-(1L:constant_pos)], 
                                              function(single_split) single_split[["neg"]])))
      })
    })
  }) %>% unlist(recursive = FALSE) %>% 
    unlist(recursive = FALSE)
}

# for future reference, positive and negative sequence length in cv
# lapply(c(6, 10, 15), function(constant_pos) {
#   lapply(c(6, 10, 15), function(constant_neg) {
#     lapply(1L:10, function(dummy) {
#       c(pos = constant_pos, neg = constant_neg)
#     })
#   })
# }) %>% unlist %>% matrix(ncol = 2, byrow = TRUE)


#' Perform cross-validation
#'
#' Performs cross-validation using the n-grams extracted from amyloids and non-amyloids. 
#' @param all_folds a list of folds as given by \code{\link{create_all_folds}}.
#' @param extracted_ngrams a list of n-grams as given by \code{\link{extract_ngrams}}.
#' @param hv a vector of hexamers indices as given by \code{\link{create_hv}}.
#' 
#' @note The classifier of choice is the \code{\link{[ranger]ranger}} implementation of 
#' random forests algorithm.
#' 
#' @return List of length equal to the number of repetitions of cross-validation.
#' Each element of the list is a list with the number of elements equal to the number 
#' of encodings. Each element of this list contains a list with two elements:
#' a) data frame with the prediction, real value of target and id of the protein,
#' b) list of important n-grams as found by \code{\link[biogram]{test_features}}.

do_cv <- function(all_folds, extracted_ngrams, hv)
  lapply(extracted_ngrams, function(encoded_group)
    lapply(all_folds, function(fold_list) 
      do_single_cv(fold_list, encoded_group, hv)
    )
  )


do_single_cv <- function(fold_list, encoded_group, hv) {
  lapply(1L:5, function(fold) {
    #training data
    train_pos <- encoded_group[hv %in% fold_list[[1]][fold_list[[1]][, "which"] != fold, "id"], ]
    train_neg <- encoded_group[hv %in% fold_list[[3]][fold_list[[3]][, "which"] != fold, "id"], ]
    
    test_bis <- test_features(c(rep(1, nrow(train_pos)), rep(0, nrow(train_neg))),
                              rbind(train_pos, train_neg), adjust = NULL)
    imp_bigrams <- cut(test_bis, breaks = c(0, 0.05, 1))[[1]]
    
    train_data <- data.frame(as.matrix(rbind(train_pos, train_neg)[, imp_bigrams]),
                             tar = as.factor(c(rep(1, nrow(train_pos)), rep(0, nrow(train_neg)))))
    model_cv <- ranger(tar ~ ., train_data, write.forest = TRUE, probability = TRUE)
    
    #test data
    test_pos <- rbind(encoded_group[hv %in% fold_list[[1]][fold_list[[1]][, "which"] == fold, "id"], ],
                      encoded_group[hv %in% fold_list[[2]][fold_list[[2]][, "which"] == fold, "id"], ])
    test_neg <- rbind(encoded_group[hv %in% fold_list[[3]][fold_list[[3]][, "which"] == fold, "id"], ], 
                      encoded_group[hv %in% fold_list[[4]][fold_list[[4]][, "which"] == fold, "id"], ])
    
    #number of n-grams from protein
    ngram_prots_pos <- c(as.vector(table(hv[hv %in% fold_list[[1]][fold_list[[1]][, "which"] == fold, "id"]])),
                         as.vector(table(hv[hv %in% fold_list[[2]][fold_list[[2]][, "which"] == fold, "id"]])))
    ngram_prots_neg <- c(as.vector(table(hv[hv %in% fold_list[[3]][fold_list[[3]][, "which"] == fold, "id"]])),
                         as.vector(table(hv[hv %in% fold_list[[4]][fold_list[[4]][, "which"] == fold, "id"]])))
    
    preds <- cbind(predict(model_cv, data.frame(as.matrix(rbind(test_pos, test_neg)[, imp_bigrams])))[["predictions"]][, 2], 
                   # probability of being amyloid
                   c(rep(1, nrow(test_pos)), rep(0, nrow(test_neg))),
                   c(unlist(lapply(1L:length(ngram_prots_pos), function(prot_id)
                     rep(prot_id, ngram_prots_pos[prot_id]))), 
                     unlist(lapply(1L:length(ngram_prots_neg), function(prot_id)
                       rep(prot_id, ngram_prots_neg[prot_id])))))
    
    # is randomForest consistent with ranger?
    # model_rf <- randomForest(tar ~ ., train_data)
    # predict(model_rf, data.frame(rbind(test_pos, test_neg)[, imp_bigrams])) == 
    #   predict(model_cv, data.frame(rbind(test_pos, test_neg)[, imp_bigrams]))[["predictions"]]
    
    preds_df <- preds %>% 
      data.frame %>% 
      rename(prob = X1, tar = X2, prot = X3) %>% 
      mutate(prot = paste0(tar, "_", prot)) %>%
      group_by(prot) %>%
      # assumption - peptide is amyloid if at least one hexagram has prob > 0.5, 
      # so we take maximum probabilities for all hexagrams belonging to the peptide
      summarise(prob = max(prob), tar = unique(tar), len = 5 + length(prot)) %>%
      mutate(len_range = cut(len, include.lowest = TRUE, breaks = c(5, 6, 10, 15, 25)))
    
    perf_measures <- lapply(levels(preds_df[["len_range"]]), function(single_range) {
      dat <- preds_df[preds_df[["len_range"]] == single_range, ]
      tryCatch(HMeasure(dat[["tar"]], dat[["prob"]])[["metrics"]], 
               warning = function(w) 
                 list(HMeasure(dat[["tar"]], dat[["prob"]])[["metrics"]],w))
    })
    
    #check if scores were switched
    reverted <- !sapply(perf_measures, is.data.frame)

    aggregated_measures <- data.frame(len_range = levels(preds_df[["len_range"]]), 
               n = as.vector(table(preds_df[["len_range"]])),
               reverted = reverted,
               do.call(rbind, lapply(1L:length(reverted), function(single_range_id) {
                 if(reverted[single_range_id]) {
                   perf_measures[[single_range_id]][[1]]
                 } else {
                   perf_measures[[single_range_id]]
                 }
               }))
    )
    
    list(aggregated_measures, imp_bigrams)
  })
}
