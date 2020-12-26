source("./functions/choose_properties.R")
source("./functions/create_encodings.R")
source("./functions/encode_amyloids.R")
source("./functions/cv.R")
source("./functions/cv_analysis.R")
source("./functions/make_classifier.R")
source("./functions/plot_tools.R")

if(Sys.info()["nodename"] == "tobit" )
  setwd("~/Dropbox/doktorat/moje_prace/amyloid_ngram2/")

library(biogram) #calc_ed
require(seqinr) #for choose_properties
require(xtable)
require(reshape2)
library(dplyr)

# read data ----------------------------------------

amyloids <- read.csv("results/all_summaries.csv") %>% 
  select(-enc, -partition) %>% 
  arrange(enc_adj, pos, len_range) 

load("./data/aa_groups.RData")
aa_groups <- string2list(aa_groups)
aa_groups[1L:2] <- c(aa1 = string2list(gsub("-", "_", tolower("ADEGHKNPQRST-C-FY-ILMV-W"))),
                     aa2 = string2list(gsub("-", "_", tolower("AG-C-DEKNPQRST-FILMVWY-H"))))

#############this part is copied from function choose_properties ################
prop_MK <- read.csv2("./data/AA_index_mk2.csv") %>% filter(!is.na(chosen))

years <- prop_MK %>% select(name) %>% unlist %>% as.character %>% sapply(function(i) 
  strsplit(last(strsplit(i, ", ")[[1]]), ")", fixed = TRUE)[[1]][1]) %>%
  as.numeric

prop_MK <- cbind(prop_MK, years = years) %>% filter(years >= 1980) %>%
  select(X, name) %>% rbind(., 
                            data.frame(X = 545:550, 
                                       name =  c("Values of Wc in proteins from class Beta, cutoff 6 A, separation 5 (Wozniak-Kotulska, 2014)",
                                                 "Values of Wc in proteins from class Beta, cutoff 8 A, separation 5 (Wozniak-Kotulska, 2014)",
                                                 "Values of Wc in proteins from class Beta, cutoff 12 A, separation 5 (Wozniak-Kotulska, 2014)",
                                                 "Values of Wc in proteins from class Beta, cutoff 6 A, separation 15 (Wozniak-Kotulska, 2014)",
                                                 "Values of Wc in proteins from class Beta, cutoff 8 A, separation 15 (Wozniak-Kotulska, 2014)",
                                                 "Values of Wc in proteins from class Beta, cutoff 12 A, separation 15 (Wozniak-Kotulska, 2014)")))
#################################################################################

trait_tab <- read.csv(file = "./results/trait_tab.csv")

# best ecodings ----------------------------------------------

best_positions <- amyloids %>% 
  select(len_range, enc_adj, AUC_mean) %>%
  group_by(len_range) %>%
  mutate(position = rank(AUC_mean)) %>%
  ungroup %>%
  group_by(enc_adj) %>%
  summarise(cum_rank = sum(position)) %>%
  arrange(desc(cum_rank)) %>%
  slice(1)

best_enc <- best_positions[["enc_adj"]]

# amyloids %>% 
#   filter(enc_adj == best_enc) %>% 
#   arrange(AUC_mean) %>% 
#   print(head = FALSE)

load("./results/cv_results_full.RData")

full_alphabet <- lapply(1L:length(cv_results_full[[1]]), function(single_replicate_id) {
  lapply(1L:5, function(single_fold) 
    data.frame(cv_results_full[[1]][[single_replicate_id]][[single_fold]][[1]], replicate = floor((single_replicate_id - 1)/15) + 1)
  )
}) %>% unlist(recursive = FALSE) %>%
  do.call(rbind, .) %>% 
  mutate(MCC = calc_mcc(TP, TN, FP, FN)) %>%
  select(replicate, len_range, AUC, MCC, Sens, Spec) %>%
  group_by(replicate, len_range) %>%
  summarize_each(funs(mean = liberal_mean, sd = liberal_sd), AUC, MCC, Sens, Spec) %>%
  ungroup %>% 
  mutate(pos = c(6, 10, 15)[replicate])

load("./results/cv_results_stand.RData")

stand_alphabet <- lapply(cv_results_stand, function(single_enc)
  lapply(1L:length(single_enc), function(single_replicate_id) {
    lapply(1L:5, function(single_fold) 
      data.frame(single_enc[[single_replicate_id]][[single_fold]][[1]], replicate = floor((single_replicate_id - 1)/15) + 1)
    )
  }) %>% unlist(recursive = FALSE) %>%
    do.call(rbind, .) %>% 
    mutate(MCC = calc_mcc(TP, TN, FP, FN)) %>%
    select(replicate, len_range, AUC, MCC, Sens, Spec) %>%
    group_by(replicate, len_range) %>%
    summarize_each(funs(mean = liberal_mean, sd = liberal_sd), AUC, MCC, Sens, Spec) %>%
    ungroup %>% 
    mutate(pos = c(6, 10, 15)[replicate])) %>% 
  do.call(rbind, .) %>% 
  mutate(enc_adj = c(rep(1, 12), rep(2, 12))) %>% 
  select(replicate, len_range, enc_adj, pos, AUC_mean, MCC_mean, 
         Sens_mean, Spec_mean, AUC_sd, MCC_sd, Sens_sd, Spec_sd) %>% 
  arrange(enc_adj, pos, len_range) 

amyloids[amyloids[["enc_adj"]] %in% c(1, 2), ] <- stand_alphabet

# properites of benc encoding  --------------------------------------

load("./results/enc_dupes.RData")

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
  droplevels() %>%
  filter(Freq == 11) %>%
  select(name) %>%
  unlist() %>%
  as.character()

traits_best_ids <- filter(frequencies, Freq == 11) %>%
  select(X) %>%
  unlist(use.names = FALSE)

best_enc_aa <- do.call(rbind, lapply(1L:length(aa_groups[[best_enc]]), function(i)
  data.frame(id = i, aa = toupper(aa_groups[[best_enc]][[i]]), stringsAsFactors = FALSE)
))

norm_props <- normalize_properties()[traits_best_ids, ] %>% 
  data.frame %>% 
  cbind(gr = rownames(.), .) %>% 
  melt(variable = "aa") %>%
  mutate(aa = as.character(aa))
  
best_enc_props <- inner_join(best_enc_aa, norm_props)

# best n-grams
ngram_freq <- read.csv("./results/ngram_freq.csv")

# distance from the best encoding --------------------------------------

# very time consuming, load from file instead
# distances = lapply(aa_groups, function(enc) calc_ed(aa_groups[[best_enc]], enc))
# distances_si = lapply(aa_groups, function(enc) calc_si(aa_groups[[best_enc]], enc))
# save(distances, distances_si, file="results/distances.Rdata")
load(file = "results/distances.Rdata")

selected_props <- aaprop[frequencies %>% top_n(5) %>% select(X) %>% unlist, ]

props_normalization <- lapply(aa_groups, function(enc) {
  coords_a <- lapply(aa_groups[[best_enc]], function(single_subgroup) rowMeans(selected_props[, single_subgroup, drop = FALSE]))
  coords_b <- lapply(enc, function(single_subgroup) rowMeans(selected_props[, single_subgroup, drop = FALSE]))
  
  norm_factor <- sum(sapply(coords_a, function(single_coords_a) {
    distances <- sapply(coords_b, function(single_coords_b) 
      # vector of distances between groups
      sqrt(sum((single_coords_a - single_coords_b)^2))
    )
    # c(dist = min(distances), id = unname(which.min(distances)))
    min(distances)
  }))
  norm_factor
})

normalized_distances <- data.frame(enc_adj=seq_along(aa_groups),
                                   ed=unlist(distances)*unlist(props_normalization))

#for full alphabet
full_alphabet_groups = lapply(unlist(aa_groups[[1]]), function(x) x)
coords_a <- lapply(aa_groups[[best_enc]], function(single_subgroup) rowMeans(selected_props[, single_subgroup, drop = FALSE]))
coords_b <- lapply(full_alphabet_groups, function(single_subgroup) rowMeans(selected_props[, single_subgroup, drop = FALSE]))

norm_factor <- sum(sapply(coords_a, function(single_coords_a) {
  distances <- sapply(coords_b, function(single_coords_b) 
    # vector of distances between groups
    sqrt(sum((single_coords_a - single_coords_b)^2))
  )
  # c(dist = min(distances), id = unname(which.min(distances)))
  min(distances)
}))
normalized_distances <- rbind(normalized_distances,c(0, 14*norm_factor)) 


ed_dat <- inner_join(amyloids %>% 
                       mutate(et = factor(ifelse(enc_adj %in% best_enc, "best", ifelse(enc_adj %in% 1L:2, "literature", "")))) %>%
                       select(len_range, enc_adj, AUC_mean, et) %>%
                       rbind(select(full_alphabet, AUC_mean, len_range) %>% 
                               mutate(et = "full alphabet", enc_adj = 0)) %>%
                       group_by(enc_adj) %>%
                       summarise(AUC_mean = mean(AUC_mean), et=et[1]), normalized_distances)

ed_dat[["et"]] <- factor(ed_dat[["et"]], labels = c("Encoding", "Best-performing encoding",
                                                    "Encoding\nfrom literature", "Full alphabet"))
cor.test(~ AUC_mean + ed, ed_dat)

# similarity index distances -------------------------------

df_si <- data.frame(enc_adj = c(seq_along(aa_groups), 0), 
                    si = c(unlist(distances_si), calc_si(aa_groups[[best_enc]], full_alphabet_groups))
)

si_dat <- inner_join(amyloids %>% 
                       mutate(et = factor(ifelse(enc_adj %in% best_enc, "best", ifelse(enc_adj %in% 1L:2, "literature", "")))) %>%
                       select(len_range, enc_adj, AUC_mean, et) %>%
                       rbind(select(full_alphabet, AUC_mean, len_range) %>% 
                               mutate(et = "full alphabet", enc_adj = 0)) %>%
                       group_by(enc_adj) %>%
                       summarise(AUC_mean = mean(AUC_mean), et=et[1]), df_si)
si_dat[["et"]] <- factor(si_dat[["et"]], labels = c("Encoding", "Best-performing encoding",
                                                    "Encoding\nfrom literature", "Full alphabet"))

cor.test(~ AUC_mean + si, si_dat)

# Pearson's product-moment correlation
# 
# data:  AUC_mean and si
# t = 80.641, df = 18536, p-value < 2.2e-16
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#  0.4988838 0.5201983
# sample estimates:
#       cor 
# 0.5096192 

# results of benchmark --------------------------------------

bench_measures <- read.csv("./results/benchmark_allpreds.csv")

# benchmark - statistical sigificance ------------------------------------------------

b_dat <- read.csv("./results/predictor_boot.csv")

levels(b_dat[["classifier"]]) <- sub("class14592", "AmyloGram", levels(b_dat[["classifier"]]))
levels(b_dat[["classifier"]]) <- sub("raw_aa", "Full alphabet", levels(b_dat[["classifier"]]))
levels(b_dat[["classifier"]]) <- sub("_", " (", levels(b_dat[["classifier"]]), fixed = TRUE)
levels(b_dat[["classifier"]]) <- sub("6", "6)", levels(b_dat[["classifier"]]), fixed = TRUE)
levels(b_dat[["classifier"]]) <- sub("10", "10)", levels(b_dat[["classifier"]]), fixed = TRUE)
levels(b_dat[["classifier"]]) <- sub("15", "15)", levels(b_dat[["classifier"]]), fixed = TRUE)
levels(b_dat[["classifier"]]) <- sub("PASTA2", "PASTA 2.0", levels(b_dat[["classifier"]]), fixed = TRUE)

b_res <- lapply(colnames(b_dat)[-c(1L:2)], function(i) {
  tmp_dat <- b_dat[, c(colnames(b_dat)[1L:2], i)]
  colnames(tmp_dat)[3] <- "measure"
  
  #n_compar <- nlevels(b_dat[["classifier"]]) * (nlevels(b_dat[["classifier"]]) - 1)/2
  n_compar <- 4
  
  group_by(tmp_dat, classifier) %>%
    summarise(m = mean(measure), 
              # l = quantile(measure, 1 - 0.95^(1/n_compar)),
              # u = quantile(measure, 0.95^(1/n_compar))
              l = mean(measure) - qt(1 - 0.975^(1/n_compar), df = length(measure) - 1)*sd(measure)/sqrt(length(measure)),
              u = mean(measure) + qt(1 - 0.975^(1/n_compar), df = length(measure) - 1)*sd(measure)/sqrt(length(measure))
    ) %>% 
    data.frame(measure = i, .)
}) %>%
  do.call(rbind, .) %>% 
  mutate(classifier = factor(classifier,
                             levels = c("Full alphabet (6)", "Full alphabet (10)", "Full alphabet (15)", 
                                        "AmyloGram (6)", "AmyloGram (10)", "AmyloGram (15)", 
                                        "appnn", "FoldAmyloid", "PASTA 2.0")),
         type = ifelse(grepl("alphabet", classifier), "Full alphabet", 
                       ifelse(grepl("AmyloGram", classifier), "AmyloGram", "Other classifier")))

levels(b_res[["measure"]]) <- c("AUC", "MCC", "Sensitivity", "Specificity")

# pairwise identity ---------------------------

load("./results/pep424_pid.RData")

ets <- c(rep(1, length(read.fasta("./data/amyloid_pos_benchmark.fasta", seqtype = "AA"))),
         rep(0, length(read.fasta("./data/amyloid_neg_benchmark.fasta", seqtype = "AA"))))

pep424_ets <- read.csv("results/benchmark_allpreds_raw.csv")[[1]]

pid100 <- lapply(1L:length(pep424_pid), function(i) {
  single_pep_dat <- cut(pep424_pid[[i]], breaks = c(0, 99, 100), include.lowest = TRUE)
  cbind(pep = i, et = pep424_ets[i], 
        as.data.frame(table(PID = single_pep_dat, ets = ets), responseName = "nprot")[c(2, 4), 2L:3])
}) %>% 
  do.call(rbind, .) %>% 
  inner_join(as.data.frame(table(ets), responseName = "count")) %>% 
  mutate(fprot = nprot/count) %>% 
  group_by(et, ets) %>% 
  summarise(fprot = mean(fprot), nprot = sum(nprot)) %>% 
  ungroup %>% 
  mutate(et = factor(et, labels = c("Non-amyloidogenic (pep424)", "Amyloidogenic (pep424)")),
         ets = factor(ets, labels = c("Non-amyloidogenic (training dataset)", "Amyloidogenic (training dataset)")))

