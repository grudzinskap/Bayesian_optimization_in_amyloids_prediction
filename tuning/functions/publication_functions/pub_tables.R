format_table <- function(x, caption, label, range) {
  bold_max <- function(x) {
    nx <- as.numeric(x)
    x[nx == max(nx)] <- paste0("\\textbf{", x[nx == max(nx)], "}")
    x
  }
  
  pretty_format <- function(x)
    formatC(x, digits = 4, format = "f")
  
  tab_bold <- x[range] %>%
    mutate_each(funs(pretty_format)) %>% 
    mutate_each(funs(bold_max))
  
  rws <- seq(1, nrow(tab_bold) - 1, by = 2)
  col <- rep("\\rowcolor[gray]{0.85}", length(rws))
  
  tab_full <- cbind(x[-range], tab_bold)
  
  res <- print(xtable(tab_full, caption = caption, label = label, align = rep("c", ncol(tab_full) + 1)), 
               include.rownames = FALSE, booktabs = TRUE,
               add.to.row = list(pos = as.list(rws), command = col), print.results = FALSE, 
               caption.placement = "top", sanitize.text.function = identity, 
               sanitize.rownames.function = identity)
  res
}

group2df <- function(group_list, caption = NULL, label = NULL) {
  tab <- data.frame(Groups = sapply(group_list, function(i)
    paste0(toupper(sort(i)), collapse = ", ")))
  tab <- cbind(ID = 1L:nrow(tab), tab)
  rws <- seq(1, nrow(tab) - 1, by = 2)
  col <- rep("\\rowcolor[gray]{0.85}", length(rws))
  colnames(tab) <- c("Subgroup ID", "Amino acids")
  res <- print(xtable(tab, caption = caption, label = label), 
               include.rownames = FALSE, booktabs = TRUE,
               add.to.row = list(pos = as.list(rws), command = col), print.results = FALSE, 
               caption.placement = "top", sanitize.text.function = identity, 
               sanitize.rownames.function = identity)
  res
}

# cross-validation summary table ----------------------------

pub_full_summary_table <- rbind(select(amyloids, AUC_mean, MCC_mean, Sens_mean, Spec_mean, pos, len_range) %>%
                                  group_by(pos, len_range) %>%
                                  summarise_each(funs(mean)) %>%
                                  mutate(classifiers = "All reduced alphabets"),
                                select(amyloids, pos, len_range, AUC_mean, MCC_mean, Sens_mean, Spec_mean, enc_adj) %>%
                                  filter(enc_adj == best_enc) %>%
                                  select(-enc_adj) %>%
                                  mutate(classifiers = "Best reduced alphabet"),
                                select(amyloids, AUC_mean, MCC_mean, Sens_mean, Spec_mean, pos, len_range, enc_adj) %>%
                                  filter(enc_adj %in% c(1L:2)) %>%
                                  select(-enc_adj) %>%
                                  group_by(pos, len_range) %>%
                                  summarise_each(funs(mean)) %>%
                                  mutate(classifiers = "Reduced alphabet from literature"),
                                select(full_alphabet, pos, len_range, AUC_mean, MCC_mean, Sens_mean, Spec_mean) %>%
                                  mutate(classifiers = "Full alphabet")
)

colnames(pub_full_summary_table) <- c("pos", "len_range", "AUC_mean", "MCC_mean", "Sens_mean", "Spec_mean", 
                                      "classifiers")



mutate(pub_full_summary_table, 
       len_range = factor(len_range, levels = levels(len_range)[c(3L:4, 1L:2)])) %>%
  group_by(len_range, classifiers) %>%
  summarise_each(funs(mean)) %>%
  ungroup %>%
  arrange(classifiers, len_range) %>%
  select(-pos) %>%
  rename("Length of peptides in testing set" = len_range, Classifier = classifiers, "Mean AUC" = AUC_mean,
         "Mean MCC" = MCC_mean, "Mean sensitivity" = Sens_mean, "Mean specificity" = Spec_mean) %>%
  select(c(2, 1, 3L:6)) %>%
  format_table("Summarized results of 15 repeats of cross-validation. The results of all reduced 
               alphabets (without the best performing encoding, the best encoding, the full 
               alphabet and the reduced alphabet are presented separately.", "tab:cv_summary", 3L:6) %>%
  cat

# best encoding ----------------------------

cat(group2df(aa_groups[[best_enc]], caption = "The best-performing encoding.", 
             label = "tab:best_enc"))

# benchmark table --------------------------------------------

slice(bench_measures, c(1L:3, 5, 16L:18)) %>%
  mutate(pos = factor(pos, labels = c("6", "6-10", "6-15"))) %>%
  mutate(nice_name = paste0(nice_name, " (", pos, ")")) %>%
  select(nice_name, AUC, MCC, Sens, Spec) %>%
  mutate(nice_name = sub(" (NA)", "", as.character(nice_name), fixed = TRUE)) %>%
  mutate(nice_name = sub("14592 (6-10)", "AmyloGram", as.character(nice_name), fixed = TRUE)) %>%
  rename(Classifier = nice_name, "AUC" = AUC, "MCC" = MCC, 
         "Sensitivity" = Sens, "Specificity" = Spec) %>%
  slice(c(4, 5L:7, 1L:3)) %>%
  format_table("Results of benchmark on \\textit{pep424} data set for AmyloGram, PASTA2, FoldAmyloid and 
               random forest predictor learned on n-grams extracted without any amino acid 
               encoding from the sequences of the length specified in the brackets.", 
               "tab:bench_summary", 2L:5) %>%
  cat

# sizes of training sets, requires analysis.R --------------

test_sizes <- data.frame(len = seq_lengths, et = ets) %>% 
  mutate(len_group = cut(seq_lengths, breaks = c(5, 6, 10, 15, 25)),
         nhexamers = len - 5) %>% 
  group_by(len_group, et) %>% 
  summarise(total_hexamers = sum(nhexamers), total = length(nhexamers)) %>% 
  mutate(set = "Test",
           etnice = ifelse(et == 1, "Amyloid", "Non-amyloid")) %>% 
  ungroup %>% 
  select(set, len_group, etnice, total, total_hexamers) %>% 
  na.omit

filter(test_sizes, len_group != "(15,25]") %>% 
  mutate(set = "Train") %>% 
  group_by(etnice) %>% 
  mutate(total = cumsum(total), total_hexamers = cumsum(total_hexamers),
         len_group = factor(len_group, labels = c("(5,6]", "(5,10]", "(5,15]"))) %>% 
  rbind(., test_sizes) %>% 
  xtable(digits = 0) %>% 
  print(include.rownames = FALSE)

# supplemental table, benchmark statistical significance

rws <- seq(1, nrow(b_res) - 1, by = 2)
col <- rep("\\rowcolor[gray]{0.85}", length(rws))

rename(b_res[1L:5], Measure = measure, Classifier = classifier,
       Mean = m, `Lower bound` = l, `Upper bound` = u) %>% 
  arrange(Measure, desc(Classifier)) %>% 
  xtable(digits = 4, 
         caption = "The mean values of performance measures and their confidence intervals obtained in the bootstrap.",
         label = "tab:meas") %>% 
  print(include.rownames = FALSE, booktabs = TRUE,
        add.to.row = list(pos = as.list(rws), command = col))
