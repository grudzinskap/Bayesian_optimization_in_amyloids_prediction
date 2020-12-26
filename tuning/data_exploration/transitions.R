get_transition_df <- function(transition_vector, dist)
  sapply(transition_vector, function(single_transition)
    grepl(single_transition, levels(ngram_freq[["decoded_name"]]))
  ) %>% melt %>% 
  mutate(Var1 = factor(Var1, labels = levels(ngram_freq[["decoded_name"]]))) %>% 
  rename(ngram = Var1, trans = Var2) %>% 
  mutate(len = nchar(as.character(ngram))) %>% 
  filter(len > (dist + 1)) %>%
  inner_join(ngram_status) %>% 
  group_by(trans, status) %>% 
  summarise(freq = mean(value)) %>% 
  mutate(dist = dist)

ngram_status <- select(ngram_freq, decoded_name, diff_freq) %>% 
  mutate(status = ifelse(diff_freq > 0, "Amyloid", "Non-amyloid")) %>% 
  select(decoded_name, status) %>% 
  rename(ngram = decoded_name)

transitions0 <- expand.grid(list(c(1L:6, "_"), c(1L:6, "_")), stringsAsFactors = TRUE) %>% 
  apply(1, paste0, collapse = "")

transitions1 <- expand.grid(list(c(1L:6, "_"), ".", c(1L:6, "_")), stringsAsFactors = TRUE) %>% 
  apply(1, paste0, collapse = "")

trans_df0 <- get_transition_df(transitions0, 0) 

trans_df1 <- get_transition_df(transitions1, 1) 

levels(trans_df1[["trans"]]) <- levels(trans_df0[["trans"]])

trans_df <- rbind(trans_df0, trans_df1) %>% 
  group_by(trans) %>% 
  mutate(dist = factor(dist), 
         all_zero = sum(freq) == 0) %>% 
  filter(!all_zero) %>% 
  mutate(start = substr(as.character(trans), 0, 1),
         end = substr(as.character(trans), 2, 2))

ggplot(trans_df, aes(x = dist, y = freq, fill = status)) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_wrap(~ trans)

adj_matrix <- trans_df0 %>% 
  mutate(start = substr(as.character(trans), 0, 1),
         end = substr(as.character(trans), 2, 2)) %>% 
  select(freq, status, start, end) %>% 
  acast(start ~ end ~ status, value.var = "freq")

library(circlize)

chod_colors <- adjustcolor(rainbow(length(adj_matrix[, , 1]), start = 0.1, end = 0.9), alpha.f = 0.5)
dim(chod_colors) = dim(adj_matrix[, , 1]) 

chordDiagram(adj_matrix[, , 1], 
             col = chod_colors)


chordDiagram(adj_matrix[, , 1], 
             grid.col = chod_colors[colSums(adj_matrix[, , 1]) != 0],
             grid.border = chod_colors[colSums(adj_matrix[, , 1]) != 0],
             row.col = chod_colors[colSums(adj_matrix[, , 1]) != 0],
             column.col = chod_colors[colSums(adj_matrix[, , 1]) != 0])
chordDiagram(adj_matrix[, , 2], 
             grid.col = chod_colors[colSums(adj_matrix[, , 2]) != 0], 
             row.col = chod_colors[colSums(adj_matrix[, , 2]) != 0],
             column.col = chod_colors[colSums(adj_matrix[, , 2]) != 0])
