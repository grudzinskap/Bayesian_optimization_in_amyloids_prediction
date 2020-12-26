
ngram_status <- select(ngram_freq, decoded_name, diff_freq) %>% 
  mutate(status = ifelse(diff_freq > 0, "Amyloid", "Non-amyloid")) %>% 
  select(decoded_name, status) %>% 
  rename(ngram = decoded_name)

# trans_df <- expand.grid(list(c(1L:6, "_"), c(1L:6, "_")), stringsAsFactors = TRUE) %>% 
#   apply(1, paste0, collapse = "") %>% 
#   sapply(function(single_transition) grepl(single_transition, levels(ngram_freq[["decoded_name"]]))) %>% 
#   melt %>% 
#   mutate(Var1 = factor(Var1, labels = levels(ngram_freq[["decoded_name"]]))) %>% 
#   rename(ngram = Var1, trans = Var2) %>% 
#   mutate(len = nchar(as.character(ngram))) %>% 
#   filter(len > 1) %>%
#   inner_join(ngram_status) %>% 
#   group_by(trans, status) %>% 
#   summarise(freq = mean(value)) %>% 
#   mutate(start = substr(as.character(trans), 0, 1),
#          end = substr(as.character(trans), 2, 2))
# 
# ggplot(trans_df, aes(x = status, y = freq, fill = freq, label = round(freq, 4))) +
#   geom_bar(stat = "identity") +
#   geom_text(color = "red") +
#   facet_grid(start ~ end, labeller = label_both)


pos_df <- expand.grid(list(c(1L:6, "_"), c(1L:6, "_"),  c(1L:6, "_")), stringsAsFactors = TRUE) %>% 
  apply(1, paste0, collapse = "") %>% 
  sapply(function(single_transition) grepl(single_transition, levels(ngram_freq[["decoded_name"]]))) %>% 
  melt %>% 
  mutate(Var1 = factor(Var1, labels = levels(ngram_freq[["decoded_name"]]))) %>% 
  rename(ngram = Var1, trans = Var2) %>% 
  mutate(len = nchar(as.character(ngram))) %>% 
  filter(len > 2) %>%
  inner_join(ngram_status) %>% 
  mutate(p1 = substr(trans, 0, 1), 
         p2 = substr(trans, 2, 2),
         p3 = substr(trans, 3, 3)) 

lapply(list(group_by(pos_df, p1, status) %>% 
       summarise(freq = mean(value)) %>% 
       mutate(trans = as.character(p1)), 
     group_by(pos_df, p1, p2, status) %>% 
       summarise(freq = mean(value)) %>% 
       mutate(trans = paste0(as.character(p1), as.character(p2))),
     group_by(pos_df, p1, p2, p3, status) %>% 
       summarise(freq = mean(value)) %>% 
       mutate(trans = paste0(as.character(p1), as.character(p2), as.character(p3)))
), function(i) select(ungroup(i), trans, status, freq)) %>% 
  do.call(rbind, .)