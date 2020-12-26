# raw data for bootstrap
b_raw <- read.csv("/results/benchmark_allpreds_raw.csv") %>% 
  select(real_labels, PASTA2, FoldAmyloid, appnn, class14592_6, 
         class14592_10, class14592_15, raw_aa_6, raw_aa_10, raw_aa_15)

HMeasure(b_raw[[1]], b_raw[-1], threshold = c(0.6226317, rep(0.5, ncol(b_raw) - 2)))[["metrics"]] %>%
  mutate(MCC = calc_mcc(TP, TN, FP, FN), classifier = rownames(.)) %>%
  select(classifier, AUC, MCC, Sens, Spec)


set.seed(1)

b_dat <- lapply(1L:1000, function(replic_id) {
  b_sample <- b_raw[sample(1L:nrow(b_raw), replace = TRUE), ]
  HMeasure(b_sample[[1]], b_sample[-1], threshold = c(0.6226317, rep(0.5, ncol(b_sample) - 2)))[["metrics"]] %>%
    mutate(MCC = calc_mcc(TP, TN, FP, FN), classifier = rownames(.)) %>%
    select(classifier, AUC, MCC, Sens, Spec) %>%
    data.frame(replic = replic_id, .)
}) %>%
  do.call(rbind, .)

write.csv(b_dat, file = "/results/predictor_boot.csv", row.names = FALSE)
