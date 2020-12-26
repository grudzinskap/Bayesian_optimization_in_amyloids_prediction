#aa_group must a list of length 1
make_classifier <- function(dat, ets, seq_lengths, max_len, aa_group, test_dat, seed_) {
  fdat <- extract_ngrams(dat[seq_lengths <= max_len, ], aa_group)[[1]]
  fets_raw <- ets[seq_lengths <= max_len]
  flens <- seq_lengths[seq_lengths <= max_len] - 5
  fets <- unlist(lapply(1L:length(flens), function(i) rep(fets_raw[i], flens[i])))
  
  test_bis <- test_features(fets, fdat, adjust = NULL)
  imp_bigrams <- cut(test_bis, breaks = c(0, 0.05, 1))[[1]]
  
  train_data <- data.frame(as.matrix(fdat[, imp_bigrams]), tar = factor(fets))
  
  model <- ranger(tar ~ ., train_data, write.forest = TRUE, probability = TRUE, seed = seed_)

  test_ngrams <- extract_ngrams(test_dat, aa_group)[[1]]
  
  test_lengths <- apply(test_dat, 1, function(i) sum(!is.na(i))) - 5
  
  raw_pred <- predict(model, data.frame(as.matrix(test_ngrams)[, imp_bigrams, drop = FALSE]))[["predictions"]]
  if(!is.matrix(raw_pred)) {
    raw_pred <- matrix(raw_pred, ncol = 2)
  }
  
  preds <- cbind(raw_pred[, 2], 
                 unlist(lapply(1L:length(test_lengths), function(i) rep(i, test_lengths[i]))))
  
  preds %>% 
    data.frame %>% 
    rename(prob = X1, prot = X2) %>%
    group_by(prot) %>%
    # assumption - peptide is amyloid if at least one hexagram has prob > 0.5, 
    # so we take maximum probabilities for all hexagrams belonging to the peptide
    summarise(prob = max(prob)) %>%
    select(prob) %>%
    unlist
}

make_classifier_MBO <- function(dat, ets, seq_lengths, max_len, aa_group, test_dat, n_, mtry_, seed_) {
  fdat <- extract_ngrams(dat[seq_lengths <= max_len, ], aa_group)[[1]]
  fets_raw <- ets[seq_lengths <= max_len]
  flens <- seq_lengths[seq_lengths <= max_len] - 5
  fets <- unlist(lapply(1L:length(flens), function(i) rep(fets_raw[i], flens[i])))
  
  test_bis <- test_features(fets, fdat, adjust = NULL)
  imp_bigrams <- cut(test_bis, breaks = c(0, 0.05, 1))[[1]]
  
  train_data <- data.frame(as.matrix(fdat[, imp_bigrams]), tar = factor(fets))
  
  model <- ranger(tar ~ ., train_data, num.trees = n_, mtry = mtry_, seed = seed_,
                  write.forest = TRUE, probability = TRUE)

  test_ngrams <- extract_ngrams(test_dat, aa_group)[[1]]
  
  test_lengths <- apply(test_dat, 1, function(i) sum(!is.na(i))) - 5
  
  raw_pred <- predict(model, data.frame(as.matrix(test_ngrams)[, imp_bigrams, drop = FALSE]))[["predictions"]]
  if(!is.matrix(raw_pred)) {
    raw_pred <- matrix(raw_pred, ncol = 2)
  }
  
  preds <- cbind(raw_pred[, 2], 
                 unlist(lapply(1L:length(test_lengths), function(i) rep(i, test_lengths[i]))))
  
  preds %>% 
    data.frame %>% 
    rename(prob = X1, prot = X2) %>%
    group_by(prot) %>%
    # assumption - peptide is amyloid if at least one hexagram has prob > 0.5, 
    # so we take maximum probabilities for all hexagrams belonging to the peptide
    summarise(prob = max(prob)) %>%
    select(prob) %>%
    unlist
}

#time around 50 [s]
#tmp <- make_classifier(seqs_m, ets, seq_lengths, 6, aa_groups[15608], test_dat_m)
# system.time(make_classifier(seqs_m, ets, seq_lengths, 6, aa_groups[9], test_dat_m))


make_classifier_whole_protein <- function(dat, ets, seq_lengths, max_len, aa_group, test_dat) {
  fdat <- extract_ngrams(dat[seq_lengths <= max_len, ], aa_group)[[1]]
  fets_raw <- ets[seq_lengths <= max_len]
  flens <- seq_lengths[seq_lengths <= max_len] - 5
  fets <- unlist(lapply(1L:length(flens), function(i) rep(fets_raw[i], flens[i])))
  
  test_bis <- test_features(fets, fdat, adjust = NULL)
  imp_bigrams <- cut(test_bis, breaks = c(0, 0.05, 1))[[1]]
  
  train_data <- data.frame(as.matrix(fdat[, imp_bigrams]), tar = factor(fets))
  
  model <- ranger(tar ~ ., train_data, write.forest = TRUE, probability = TRUE)
  
  test_ngrams <- extract_ngrams(test_dat, aa_group)[[1]]
  
  predict(model, data.frame(as.matrix(test_ngrams)[, imp_bigrams]))[["predictions"]][, 2]
}