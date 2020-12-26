# create feature-encoding tables - on which features are encodings based

source("create_encodings.R")

#table of features-encodings (features on which encodings are based)
fenc_table <- do.call(rbind, lapply(all_traits_combn_list, function(i)
  cbind(i, matrix(NA, nrow = nrow(i), ncol = length(vtraits) - ncol(i)))))
rownames(fenc_table) <- paste0("ID", 1L:nrow(fenc_table))

#names of identical encodings
dupl_encodings <- unlist(lapply(all_aa_groups, function(i) {
  pasted_enc <- sapply(i, function(i) paste0(unlist(i), collapse = ""))
  unique_enc <- unique(pasted_enc)
  lapply(unique_enc, function(j)
    names(pasted_enc[pasted_enc == j]))
}), recursive = FALSE)

#features duplicating encodings
dupl_features <- lapply(dupl_encodings, function(duplicates) {
  enc_id <- strsplit(duplicates, "K") %>% 
    sapply(., first) %>%
    unlist
  fenc_table[enc_id, ] %>% 
    as.vector %>% 
    na.omit %>% 
    as.vector %>% 
    table %>% 
    data.frame
})