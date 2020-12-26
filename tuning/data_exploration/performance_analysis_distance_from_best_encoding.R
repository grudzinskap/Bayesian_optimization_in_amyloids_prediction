#' Analysis of the encodings close to the best encoding
#' 

#### loading functions and pacakges ######
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

#### loading data and identifying best encoding #####

amyloids <- read.csv("results/all_summaries.csv")
trait_tab <- read.csv(file = "./results/trait_tab.csv")
load("./results/enc_dupes.RData")
load("./data/aa_groups.RData")
aa_groups <- string2list(aa_groups)

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

##### calculating encoding distances ######

# very time consuming, load from file instead
# distances = lapply(aa_groups, function(enc) calc_ed(aa_groups[[best_enc]], enc))
# save(distances, file="results/distances.Rdata")
load(file = "results/distances.Rdata")

# group_sample <- with(set.seed(0), sample(1:length(aa_groups), 10))
# distances = lapply(aa_groups[group_sample], function(enc) calc_ed(aa_groups[[best_enc]], enc))

#first we get info about properties that contribute to this encoding.
#we take into account all duplicates

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
                            apply(1, na.omit) %>% unlist %>% table) %>% arrange(desc(Freq)) 
colnames(frequencies) <- c("X", "Freq")
frequencies$X <- as.numeric(as.character(frequencies$X))


#now we load info about values for each aminoacid for properties
data(aaprop)

#we need to extract only selected properties
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

#normalized distances
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

source("functions/pub_functions.R")
#we join information about distance with mean AUC
inner_join(amyloids %>% 
             mutate(et = factor(ifelse(enc_adj %in% best_enc, "best", ifelse(enc_adj %in% 1L:2, "literature", "")))) %>%
             select(len_range, enc_adj, AUC_mean, et) %>%
             rbind(select(full_alphabet, AUC_mean, len_range) %>% 
                     mutate(et = "full alphabet", enc_adj = 0)) %>%
             group_by(enc_adj) %>%
             summarise(AUC_mean = mean(AUC_mean), et=et[1]), normalized_distances) -> dat

dat$et <- factor(dat$et, labels = c("Reduced alphabet", "Best performing reduced alphabet",
                      "Reduced alphabet from literature", "Full alphabet"))

#this shows really strong correlation
library(ggplot2)
ggplot(dat, aes(x=ed, y=AUC_mean, color=et)) + 
  geom_point() +
  scale_color_manual("", values = c("grey", "red", "blue", "green")) +
  xlab("Normalized encoding distance") +
  ylab("AUC") +
  my_theme +
  geom_point(data = filter(dat, et != "Reduced alphabet"), 
             aes(x = ed, y = AUC_mean, color = et))
