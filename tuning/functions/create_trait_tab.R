require(seqinr)
require(dplyr)
require(pbapply)
require(biogram)

source("./functions/choose_properties.R")
source("./functions/create_encodings.R")

ftraits <- choose_properties()
trait_combn <- create_traits_combination(ftraits)
load("./data/aa_groups.RData")
aa_groups <- string2list(aa_groups)

idk_tab <- sapply(names(aa_groups), function(i)
  substr(i, 3, nchar(i))) %>%
  strsplit("K") %>%
  lapply(as.numeric) %>%
  do.call(rbind, .) 

head(idk_tab)
colnames(idk_tab) <- c("ID", "K")

idk_tab[-c(1L:2), ] %>%
  data.frame() %>%
  group_by(K) %>%
  summarise(max(ID))

trait_tab <- do.call(rbind, lapply(trait_combn, function(i) {
  cbind(i, matrix(NA, ncol = 17 - ncol(i), nrow = nrow(i)))
}))

write.csv(trait_tab, file = "./results/trait_tab.csv", row.names = FALSE)

# trait_tab[best_enc,] %>%
#   apply(1, na.omit) %>%
#   unlist %>%
#   table

#encodings and duplicates
enc_dupes <- create_encodings(ftraits, TRUE)
save(enc_dupes, file = "./results/enc_dupes.RData")
