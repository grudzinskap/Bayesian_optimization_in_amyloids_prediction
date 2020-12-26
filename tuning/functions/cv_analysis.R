# Functions here should not be used, because they were written with the different data structure in mind. 
# They are preserved only for a potential future use.

#' Explain replicates of cross-validation
#'
#' Since the integration with WCSS, the exact structure of cross-validation results 
#' might be obstruse, so this function returns additional information about the cross-validation 
#' replicates.
#' 
#' @return List of maximum lengths of sequences in training set.

get_reps_ids <- function() 
  lapply(c(6, 10, 15), function(constant_pos) {
    lapply(c(6, 10, 15), function(constant_neg) {
      lapply(1L:5, function(dummy) {
        c(pos = constant_pos, neg = constant_neg)
      })
    })
  }) %>% unlist %>% 
    matrix(ncol = 2, byrow = TRUE) %>%
    data.frame %>%
    rename(pos = X1, neg = X2)




#' Summarize results of cross-validation
#'
#' Computes mean performance measure of cross-validation results.
#'
#' @param cv_results results of cross-valdiation as per \code{\link{do_cv}} function.
#' @param reps_ids sizes of positive and negative data sets as provided by
#' \code{\link{get_reps_ids}}
#' 
#' @return Data frame containing mean and standard deviation of AUC, MCC, 
#' sensitivity and specificity.
summarize_cv_results <- function(cv_results, reps_ids) {
  res <- do.call(rbind, lapply(1L:length(cv_results), function(single_enc_id) {
    all_reps_res <- lapply(cv_results[[single_enc_id]], function(single_rep) {
      do.call(rbind, lapply(single_rep, function(single_fold_res) {
        single_fold_res[[1]]
      })) %>% mutate(MCC = calc_mcc(TP, TN, FP, FN)) %>%
        select(len_range, AUC, MCC, Sens, Spec)
    }) %>% do.call(rbind, .) %>% 
      cbind(., get_reps_ids()) %>% 
      group_by(len_range, pos, neg) %>%
      summarize_each(funs(mean, sd)) %>%
      cbind(., enc = single_enc_id)
  }))
  
  # dplyr changes order of len_ranger, need to get it back
  res[["len_range"]] <- factor(res[["len_range"]], levels(res[["len_range"]])[c(4, 3, 1, 2)])
  res[["pos"]] <- as.factor(res[["pos"]])
  res[["neg"]] <- as.factor(res[["neg"]])
  res
}


#helper function calculating MCC
calc_mcc <- function(TP, TN, FP, FN)
  (TP*TN - FP*FN)/sqrt((TP + FP)*(TP + FN)*(TN + FP)*(TN + FN))

# levels(cv_summary[["len_range"]]) <- paste0("Test peptide length: ", levels(cv_summary[["len_range"]]))
# library(ggplot2)
# ggplot(cv_summary, aes(x = pos, y = neg, fill = mauc, label = round(mauc, 4))) +
#   geom_tile() +
#   geom_text(color = "red") +
#   facet_wrap(~ len_range) +
#   scale_x_discrete("Positive training set") +
#   scale_y_discrete("Positive negative set")



#' Get important n-grams from cross-validation
#'
#' Extracts important n-grams (n-grams that occured in all folds of every replicate of
#' cross-validation for a single combination of positive and negative training data 
#' set sizes).
#' @param cv_results results of cross-valdiation as per \code{\link{do_cv}} function.
#' @param reps_ids sizes of positive and negative data sets as provided by
#' \code{\link{get_reps_ids}}
#' 
#' @return List of data frames containing important n-grams for a single encoding.
summarize_cv_ngrams <- function(cv_results, reps_ids) {
  lapply(1L:length(cv_results), function(single_enc_id) {
    all_reps_res <- lapply(cv_results[[single_enc_id]], function(single_rep) {
      unlist(lapply(single_rep, function(single_fold) single_fold[[2]]))
    })
    
    # ids of replicates
    reps_ids_factor <- factor(apply(reps_ids(), 1, paste, collapse = "_"))
    
    # concatenated n-grams from all replicates, still splitted by training sets
    conc_res <- lapply(levels(reps_ids_factor), function(rep_group) {
      all_reps_res[reps_ids_factor == rep_group] %>%
        unlist %>%
        table(dnn = "ngram") %>%
        as.data.frame(responseName = "count", stringsAsFactors = FALSE) %>%
        filter(count == 25) %>%
        select(ngram) %>%
        unlist(use.names = FALSE)
    })
    
    do.call(rbind, lapply(1L:length(conc_res), function(group_id) {
      train_sizes <- as.numeric(strsplit(levels(reps_ids_factor)[group_id], "_")[[1]])
      data.frame(pos = train_sizes[1], neg = train_sizes[2], ngram = conc_res[[group_id]],
                 stringsAsFactors = FALSE)
    }))
  })
}

liberal_mean <- function(x)
  mean(x, na.rm = TRUE)

liberal_sd <- function(x)
  sd(x, na.rm = TRUE)

# imp_ngrams <- summarize_cv_ngrams(cv_results_full, get_reps_ids)
# 
# tmp <- imp_ngrams[[1]][["ngram"]]
# 
# gap_ngrams(tmp)

