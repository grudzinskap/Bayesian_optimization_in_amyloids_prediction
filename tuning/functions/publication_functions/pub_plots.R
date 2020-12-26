# plot data ----------------------------------------------

amyloids_plot <- select(amyloids, AUC_mean, MCC_mean, Sens_mean, Spec_mean, pos, len_range, enc_adj) %>%
  mutate(et = factor(ifelse(enc_adj %in% best_enc, "best", ifelse(enc_adj %in% 1L:2, "literature", "")))) %>%
  select(AUC_mean, MCC_mean, Sens_mean, Spec_mean, pos, len_range, et, enc_adj) %>%
  rbind(select(full_alphabet, AUC_mean, MCC_mean, Sens_mean, Spec_mean, pos, len_range) %>% 
          mutate(et = "full alphabet", enc_adj = 0)) %>%
  mutate(len_range = factor(len_range, levels = c("[5,6]", "(6,10]", "(10,15]", "(15,25]")),
         pos = factor(pos, labels = paste0("Training peptide length: ", 
                                           c("6", "6-10", "6-15"))),
         et = factor(et, labels = c("Encoding", "Best-performing encoding",
                                    "Standard encoding", "Full alphabet"))) %>%
  mutate(len_range = factor(len_range, 
                            labels = paste0("Test peptide length: ", c("6 ", "7-10", "11-15", "16-25"))),
         et2 = ifelse(enc_adj == 1L, "Standard encoding (Kosiol, et al., 2004)", as.character(et)),
         et2 = ifelse(enc_adj == 2L, "Standard encoding (Melo and Marti-Renom, 2006)", as.character(et2)),
         et2 = factor(et2, levels = c("Encoding", 
                                      "Best-performing encoding", 
                                      "Full alphabet", 
                                      "Standard encoding (Kosiol, et al., 2004)", 
                                      "Standard encoding (Melo and Marti-Renom, 2006)"),
                      labels = c("Encoding", 
                                 "Best-performing encoding", 
                                 "Full alphabet", 
                                 "Standard encoding [21]", 
                                 "Standard encoding [22]")),
         et = et2)

# write.csv(amyloids_plot, file = "./results/amyloid_plot_data.csv")

# Fig 1 all encodings sens/spec  ----------------------------------------

sesp_dat <- amyloids_plot
levels(sesp_dat[["pos"]]) <- c("Training peptide length: 6", "Training peptide length: 6-10", 
                               "Training peptide length: 6-15")

# sesp_plot <- ggplot(sesp_dat, aes(x = Spec_mean, y = Sens_mean, color = et)) +
#   #geom_density_2d(color = "grey", contour = TRUE) +
#   #stat_density2d(aes(fill=..level..), color = "grey", contour = TRUE, geom="polygon") +
#   scale_alpha_continuous(range = c(0.35, 1)) +
#   scale_y_continuous("Mean sensitivity") +
#   scale_x_continuous("Mean specificity") +
#   geom_point(data = sesp_dat, aes(shape = et)) +
#   geom_point(data = filter(sesp_dat, et != "Encoding"), aes(shape = et)) +
#   guides(color = guide_legend(nrow = 2), shape = guide_legend(nrow = 2)) +
#   scale_shape_manual("", values = c(1, 16, 16, 17, 15), drop = FALSE) +
#   scale_color_manual("", values = c("grey", "firebrick1", "green3", "dodgerblue", "dodgerblue"), drop = FALSE) +
#   scale_size_manual("", values = c(1, 1, 1, 1.5, 1.5), drop = FALSE) +
#   facet_grid(pos ~ len_range) +
#   my_theme 


# sesp_plot <- ggplot(sesp_dat, aes(x = Spec_mean, y = Sens_mean, color = et)) +
#   geom_bin2d(bins = 20, color = "black") + 
#   scale_fill_gradient2("Number of encodings", low = "beige", mid = "orange3", 
#                        high = "deeppink1", midpoint = 4500) +
#   scale_y_continuous("Mean sensitivity") +
#   scale_x_continuous("Mean specificity\n") +
#   geom_point(data = droplevels(filter(sesp_dat, et != "Encoding")),
#              aes(x = Spec_mean, y = Sens_mean, color = et2, shape = et2),
#              fill = NA) +
#   guides(color = guide_legend(nrow = 2), shape = guide_legend(nrow = 2), 
#          fill = guide_colorbar(barwidth = unit(10, "line"))) +
#   scale_shape_manual("", values = c(21, 23, 24, 25), drop = FALSE) +
#   scale_color_manual("", values = c("firebrick1", "green3", "dodgerblue", "dodgerblue"), drop = FALSE) +
#   scale_size_manual("", values = c(0.5, 0.5, 0.5, 0.5) + 0.5, drop = FALSE) +
#   facet_grid(pos ~ len_range) +
#   my_theme

# sesp_plot <- ggplot(sesp_dat, aes(x = Spec_mean, y = Sens_mean, color = et)) +
#   geom_bin2d(bins = 20, color = "black") + 
#   scale_fill_gradient2("Number of encodings", low = "beige", mid = "orange3", 
#                        high = "deeppink1", midpoint = 4500) +
#   scale_y_continuous("Mean sensitivity") +
#   scale_x_continuous("Mean specificity\n") +
#   geom_point(data = droplevels(filter(sesp_dat, et != "Encoding")),
#              aes(x = Spec_mean, y = Sens_mean, color = et2, shape = et2),
#              fill = NA) +
#   guides(color = guide_legend(nrow = 2), shape = guide_legend(nrow = 2), 
#          fill = guide_colorbar(barwidth = unit(10, "line"))) +
#   scale_shape_manual("", values = c(21, 23, 24, 25), drop = FALSE) +
#   scale_color_manual("", values = c("firebrick1", "green3", "dodgerblue", "dodgerblue"), drop = FALSE) +
#   scale_size_manual("", values = c(0.5, 0.5, 0.5, 0.5) + 0.5, drop = FALSE) +
#   facet_grid(pos ~ len_range) +
#   my_theme
# 
# 
# #png("./publication/figures/sesp_plot.png", height = 4, width = 6.5, unit = "in", res = 200)
# cairo_ps("./publication/figures/sesp_plot.eps", height = 4.5, width = 6.5)
# # should be eps, but it's too big for overleaf
# print(sesp_plot)
# dev.off()


sesp_plots_names <- unlist(lapply(levels(sesp_dat[["pos"]]), function(single_pos) 
  lapply(levels(sesp_dat[["len_range"]]), function(single_len_range) {
    paste0(single_pos, ", ", tolower(single_len_range))
  })
))

cat(sapply(1L:length(sesp_plots_names), function(i) 
  paste0("\\subsection{", sesp_plots_names[i], "}\n", 
         "\\centerline{\\includegraphics{figures/sesp_plot", i, ".eps}}")
), sep = "\n\n")

sesp_plots <- lapply(levels(sesp_dat[["pos"]]), function(single_pos) 
  lapply(levels(sesp_dat[["len_range"]]), function(single_len_range) {
    subdat <- filter(sesp_dat, pos == single_pos, len_range == single_len_range)
    ggplot(subdat,  aes(x = Spec_mean, y = Sens_mean, color = et)) +
      geom_bin2d(bins = 20, color = "black") + 
      scale_fill_gradient2("Number of encodings", low = "beige", mid = "orange3", 
                           high = "deeppink1", midpoint = 750) +
      scale_y_continuous("Mean sensitivity") +
      scale_x_continuous("Mean specificity\n") +
      geom_point(data = droplevels(filter(subdat, et != "Encoding")),
                 aes(x = Spec_mean, y = Sens_mean, color = et2, shape = et2),
                 fill = NA, stroke = 1, size = 4) +
      guides(color = guide_legend(nrow = 2), shape = guide_legend(nrow = 2), 
             fill = guide_colorbar(barwidth = unit(10, "line"))) +
      scale_shape_manual("", values = c(21, 23, 24, 25), drop = FALSE) +
      scale_color_manual("", values = c("firebrick1", "green3", "dodgerblue", "dodgerblue"), drop = FALSE) +
      scale_size_manual("", values = c(0.5, 0.5, 0.5, 0.5) + 3.5, drop = FALSE) +
      #ggtitle(paste0(single_pos, "\n", single_len_range)) +
      my_theme
  })
) %>% unlist(recursive = FALSE)

for(i in 1L:length(sesp_plots)) {
  #png("./publication/figures/sesp_plot.png", height = 4, width = 6.5, unit = "in", res = 200)
  cairo_ps(paste0("./supplements/figures/sesp_plot", i, ".eps"), height = 7.5, width = 6)
  # should be eps, but it's too big for overleaf
  print(sesp_plots[[i]])
  dev.off()
}

# Fig 2 AUC boxplot  ----------------------------------------

AUC_boxplot <- ggplot(amyloids_plot, aes(x = len_range, y = AUC_mean)) +
  geom_boxplot(outlier.color = "grey", outlier.shape = 1, outlier.size = 1) +
  geom_point(data = filter(amyloids_plot, et2 != "Encoding"), 
             aes(x = len_range, y = AUC_mean, shape = et2, color = et2, size = et2),
             fill = NA) +
  scale_x_discrete("") +
  scale_y_continuous("Mean AUC") +
  guides(color = guide_legend(nrow = 3), shape = guide_legend(nrow = 3)) +
  scale_shape_manual("", values = c(1, 21, 23, 24, 25), drop = FALSE) +
  scale_color_manual("", values = c("grey", "firebrick1", "green3", "dodgerblue", "dodgerblue"), drop = FALSE) +
  #scale_color_manual("", values = c("grey", rep("black", 4)), drop = FALSE) +
  scale_size_manual("", values = c(1, 1.5, 1.5, 1.5, 1.5), drop = FALSE) +
  facet_wrap(~ pos, nrow = 3) +
  my_theme + 
  coord_flip() 

cairo_ps("./publication/figures/AUC_boxplot.eps", height = 3.1, width = 6.5)
#png("./pub_figures/AUC_boxplot.png", height = 648, width = 648)
print(AUC_boxplot)
dev.off()

# Fig 3 MCC boxplot  ----------------------------------------

MCC_boxplot <- ggplot(amyloids_plot, aes(x = len_range, y = MCC_mean)) +
  geom_boxplot(outlier.color = "grey", outlier.shape = 1) +
  geom_point(data = filter(amyloids_plot, et != "Encoding"), 
             aes(x = len_range, y = MCC_mean, color = et, shape = et)) +
  scale_x_discrete("") +
  scale_y_continuous("Mean MCC") +
  scale_shape_manual("", values = c(1, 16, 15, 15), drop = FALSE) +
  scale_color_manual("", values = c("grey", "firebrick1", "dodgerblue", "green3"), drop = FALSE) +
  facet_wrap(~ pos, nrow = 3) +
  my_theme + 
  coord_flip()

# cairo_ps("./publication/figures/MCC_boxplot.eps", height = 9, width = 8)
# #png("./pub_figures/MCC_boxplot.png", height = 648, width = 648)
# MCC_boxplot
# dev.off()

# Fig 4 properties  ----------------------------------------

ggplot(best_enc_props, aes(x = as.factor(id), y = value, label = aa)) +
  geom_quasirandom() +
  geom_text_repel() +
  facet_wrap(~ gr, ncol = 2)



# Fig 5 n-grams  ----------------------------------------

g_legend<-function(a.gplot) {
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  tmp$grobs[[leg]]
}

gr_aa <- group_by(best_enc_aa, id) %>% 
  summarise(gr = paste0("{", paste0(aa, collapse = ""), "}")) 

ngram_freq_plot <- mutate(ngram_freq, decoded_name = gsub("_", "-", decoded_name)) %>%
  mutate(decoded_name = factor(decoded_name, levels = as.character(decoded_name)),
         amyloid = diff_freq > 0) %>%
  melt() %>%
  filter(variable %in% c("pos", "neg")) %>%
  droplevels %>%
  mutate(variable = factor(variable, labels = c("Amyloid", "Non-amyloid")))

for(i in 1L:6)
  levels(ngram_freq_plot[["decoded_name"]]) <- gsub(as.character(i), gr_aa[i, "gr"],
                                                    levels(ngram_freq_plot[["decoded_name"]]))

labels_colors <- c("black", "chartreuse3", "dodgerblue2", "firebrick1", "darkorange", "darkseagreen4", "darkorchid3")

gen_labels <- function(single_gr, x, gr_aa) {
  new_lab <- x
  for (other_gr in gr_aa[["gr"]][-single_gr])
    levels(new_lab) <- gsub(as.character(other_gr), paste0(rep(" ", nchar(other_gr)), collapse = ""), 
                            levels(new_lab), fixed = TRUE)
  levels(new_lab) <- gsub("-", " ", levels(new_lab), fixed = TRUE)
  new_lab
}

# generate list of labels where every element (amino acid group or dash) is present only once

all_labels <- c(eval({
  new_lab <- ngram_freq_plot[["decoded_name"]]
  levels(new_lab) <- gsub("[,A-Z\\{}]", " ", levels(new_lab))
  list(new_lab)
}),
lapply(1L:6, function(i) gen_labels(i, ngram_freq_plot[["decoded_name"]], gr_aa)))

levels(all_labels[[1]]) <- gsub("-", "X", levels(all_labels[[1]]))

# create series of plots where only one element (group of amino acids or dash) is plotted 

ngram_plots <- lapply(1L:7, function(i) {
  p <- ggplot(ngram_freq_plot, aes(x = decoded_name, y = value)) +
    geom_bar(aes(fill = variable), position = "dodge", 
             stat = "identity", color = "black", size = 0.1,
             width = 0.7) +
    geom_point(data = group_by(ngram_freq_plot, decoded_name)  %>% filter(value == max(value)),
               aes(y = value + 0.004, shape = association), size = 2, stroke = 0.2, fill = "white") +
    scale_fill_manual("", values = c("white", "black")) +
    scale_shape_manual("Experimentally validated motif:", breaks = c("Amyloidogenic", "Non-amyloidogenic"), 
                       values = c(21, 16, NA)) +
    scale_y_continuous("Frequency\n") +
    scale_x_discrete("", labels = all_labels[[i]]) + 
    theme(axis.text.y = element_text(size = 6, colour = labels_colors[i], family = "mono", face = "bold"),
          legend.box = "vertical") +
    coord_flip() +
    my_theme
})

amyl_legend <- g_legend(ngram_plots[[1]])

ngrams_plots_final <- lapply(1L:length(ngram_plots), function(i)
  if(i < 7) {
    arrangeGrob(ngram_plots[[i]] + 
                  theme(legend.position="none",
                        axis.title.x = element_text(size=14 + size_mod, vjust = -1, color = "white")),
                rectGrob(x = unit(0.5, "npc"), y = unit(0.5, "npc"), gp = gpar(col = "white")), 
                nrow = 2, heights = c(0.96, 0.04))
  } else {
    arrangeGrob(ngram_plots[[i]] + theme(legend.position="none"),
                amyl_legend, 
                nrow = 2, heights = c(0.96, 0.04))
    
  }
)

# combine plots

cairo_ps("./publication/figures/ngrams.eps", height = 8.1, width = 5.3)
for(i in 1L:7) {
  grid.draw(ngrams_plots_final[[i]])
}
dev.off()


# convert_factor <- function(x)
#   factor(x, levels = c(1L:6, "|", "_"))
# 


# in case we need to get n-grams in a tabular format
#writeLines(as.character(ngram_freq_plot[["decoded_name"]]), "n_gramy_Ania.txt")

# convert_factor <- function(x)
#   factor(x, levels = c(1L:6, "|", "_"))
# 
# block_df <- unique(ngram_freq[["decoded_name"]]) %>% 
#   as.character() %>% 
#   strsplit(split = "") %>% 
#   lapply(function(i) c(rep("", 5 - length(i)), i)) %>% 
#   do.call(rbind, .) %>% 
#   data.frame() %>% 
#   mutate_each(funs(convert_factor)) %>% 
#   mutate(y = factor(1L:nrow(.))) %>% 
#   melt(measure.vars = paste0("X", 1L:5), variable.name = "x", value.name = "block") %>% 
#   mutate(block = factor(block))
# 
# block_plot <- ggplot(droplevels(filter(block_df, block != "")), aes(x = x, y = y, fill = block)) +
#   geom_tile(color = "black") +
#   scale_fill_manual(values = c("lightgrey", "chartreuse3", "dodgerblue2", 
#                                "firebrick1", "darkorange", "darkseagreen4", "cyan3"), 
#                     labels = c("-", gr_aa[["gr"]])) +
#   scale_x_discrete("") +
#   scale_y_discrete("") +
#   #theme_void() +
#   my_theme +
#   theme(axis.text = element_blank(),
#         axis.ticks = element_blank(),
#         panel.grid.major = element_line(color= NA),
#         panel.grid.major = element_line(color= NA),
#         panel.background = element_rect(color = NA))
# 
# ngram_plot <- ggplot(ngram_freq_plot, aes(x = decoded_name, y = value)) +
#   geom_bar(aes(fill = variable), position = "dodge", stat = "identity") +
#   geom_point(data = group_by(ngram_freq_plot, decoded_name)  %>% filter(value == max(value)),
#              aes(y = value + 0.007, shape = association), size = 2) +
#   scale_fill_manual("", values = c("firebrick1", "dodgerblue2")) +
#   scale_shape_manual("Experimentally tested motif:", 
#                      breaks = c("Amyloidogenic", "Non-amyloidogenic"), 
#                      values = c(16, 17, NA)) +
#   scale_y_continuous("Frequency") +
#   scale_x_discrete("") +
#   my_theme +
#   theme(axis.text.y = element_blank()) +
#   coord_flip() 
# 
# block_plot_nol <- ggplot_gtable(ggplot_build(block_plot + theme(legend.position="none")))
# ngram_plot_nol <- ggplot_gtable(ggplot_build(ngram_plot + theme(legend.position="none")))
# 
# max_par = unit.pmax(block_plot_nol$heights[3:4], ngram_plot_nol$heights[3:4])
# block_plot_nol$heights[3:4] <- max_par
# block_plot_nol$heights[3:4] <- max_par
# 
# pdf("tmp.pdf", height = 8.5, width = 3.5)
# grid.draw(arrangeGrob(block_plot_nol,
#                       ngram_plot_nol,
#                       g_legend(block_plot),
#                       g_legend(ngram_plot), nrow = 2, ncol = 2, 
#                       heights=c(0.90, 0.1)))
# dev.off()

# Fig 6 alternative (similarity index)  ----------------------------------------

# careful - check if similarity index is used instead of the encoding distance

si_dat <- si_dat %>% 
  mutate(et2 = ifelse(enc_adj == 1L, "Standard encoding (Kosiol, et al., 2004)", as.character(et)),
         et2 = ifelse(enc_adj == 2L, "Standard encoding (Melo and Marti-Renom, 2006)", as.character(et2)),
         et2 = factor(et2, levels = c("Encoding", 
                                      "Best-performing encoding", 
                                      "Full alphabet", 
                                      "Standard encoding (Kosiol, et al., 2004)", 
                                      "Standard encoding (Melo and Marti-Renom, 2006)"),
                      labels = c("Encoding", 
                                 "Best-performing encoding", 
                                 "Full alphabet", 
                                 "Standard encoding [21]", 
                                 "Standard encoding [22]")),
         et = et2)

write.csv2(si_dat, row.names = FALSE, file = "./results/si_dat.csv")

# si_AUC_plot <- ggplot(si_dat, aes(x=si, y=AUC_mean, color=et, shape = et)) + 
#   geom_point() +
#   scale_color_manual("", values = c("grey", "firebrick1", "dodgerblue", "green3")) +
#   scale_shape_manual("", values = c(1, 16, 15, 15), drop = FALSE) +
#   xlab("Similarity index") +
#   ylab("AUC") +
#   my_theme +
#   geom_point(data = filter(si_dat, et != "Encoding"), 
#              aes(x = si, y = AUC_mean, color = et)) +
#   guides(color = guide_legend(nrow = 2), shape = guide_legend(nrow = 5)) +
#   scale_shape_manual("", values = c(1, 16, 16, 17, 15), drop = FALSE) +
#   scale_color_manual("", values = c("grey", "firebrick1", "green3", "dodgerblue", "dodgerblue"), drop = FALSE) +
#   scale_size_manual("", values = c(1, 1, 1, 1.5, 1.5), drop = FALSE) 

si_AUC_plot <- ggplot(si_dat, aes(x=si, y=AUC_mean)) + 
  geom_bin2d(bins = 25, color = "black") + 
  scale_fill_gradient2("Number of encodings", low = "white", high = "orange",
                       midpoint = 500) +
  xlab("Similarity to the best-performing encoding\n") +
  ylab("AUC") +
  my_theme +
  geom_point(data = droplevels(filter(si_dat, et != "Encoding")),
             aes(x = si, y = AUC_mean, color = et, shape = et),
             fill = NA) +
  guides(color = guide_legend(nrow = 4), shape = guide_legend(nrow = 4), 
         fill = guide_colorbar(barwidth = unit(6, "line"))) +
  scale_shape_manual("", values = c(21, 23, 24, 25), drop = FALSE) +
  scale_color_manual("", values = c("firebrick1", "green3", "dodgerblue", "dodgerblue"), drop = FALSE) +
  #scale_color_manual("", values = rep("black", 4), drop = FALSE) +
  scale_size_manual("", values = c(0.5, 0.5, 0.5, 0.5) + 1.5, drop = FALSE) 

cairo_ps("./publication/figures/ed_AUC.eps", height = 3.5, width = 5)
print(si_AUC_plot)
dev.off()


# save(amyloids_plot, best_enc_props, ngram_freq_plot, ed_dat,
#      file = "./presentation/presentation.RData")

# supplemental figure: pairwise identity -----------------------

cairo_ps("./supplements/figures/pid.eps", height = 5.1, width = 5)
ggplot(pid100, aes(x = et, y = nprot, fill = ets, label = nprot, color = ets)) +
  geom_bar(stat = "identity", position = position_dodge(width = 1), color = "black") +
  geom_text(aes(y = nprot - 10), position = position_dodge(width = 1)) + 
  scale_x_discrete("") +
  scale_y_continuous("Total number of sequences with 100% pairwise identity") +
  scale_fill_manual("", values = c("black", "white")) +
  scale_color_manual(values = c("white", "black")) +
  guides(color = "none") +
  my_theme
dev.off()


# supplemental figure: benchmark significance -----------------------

cairo_ps("./supplements/figures/signif.eps", height = 7.9, width = 6)
ggplot(b_res, aes(y = m, x = classifier, ymin = l, ymax = u, color = type)) +
  #geom_point() +
  geom_errorbar() +
  facet_wrap(~ measure, ncol = 2, scales = "free_y") +
  scale_color_discrete("") +
  my_theme + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  scale_x_discrete("") +
  scale_y_continuous("Value")
dev.off()

# supplemental figure: size and flexibility -----------------------


nice_prop_names <- data.frame(property = c(levels(best_enc_props[["gr"]]),
                                           "DAWD720101", "BIGC670101", "ZIMJ680102"),
                              nice_name = c("Hydrophobicity index\n(Argos et al., 1982)",
                                            "Average flexibility indices\n(Bhaskaran-Ponnuswamy, 1988)",
                                            "Polarizability parameter\n(Charton-Charton, 1982)",
                                            "Thermodynamic beta sheet propensity\n(Kim-Berg, 1993)",
                                            "Size\n(Dawson, 1972)",
                                            "Residue volume\n(Bigelow, 1967)",
                                            "Bulkiness\n(Zimmerman et al., 1968)"))

prop_plot_dat <- lapply(c(levels(best_enc_props[["gr"]]), 
                          "DAWD720101", "BIGC670101", "ZIMJ680102"), function(single_prop) {
  amyl_flex <- read.fasta("./data/amyloid_pos_full.fasta",seqtype = "AA") %>% 
    sapply(function(i) {
      i <- i[i %in% a()]
      mean(aaprop[single_prop, as.vector(tolower(i))])
    })
  
  nonamyl_flex <- read.fasta("./data/amyloid_neg_full.fasta",seqtype = "AA") %>% 
    sapply(function(i) {
      i <- i[i %in% a()]
      mean(aaprop[single_prop, as.vector(tolower(i))])
    })
  
  rbind(data.frame(value = amyl_flex, status = "Amyloidogenic", property = single_prop),
        data.frame(value = nonamyl_flex, status = "Non-amyloidogenic", property = single_prop))
}) %>% 
  do.call(rbind, .) %>% 
  inner_join(nice_prop_names) %>% 
  mutate(size_prop = property %in% c("DAWD720101", "BIGC670101", "ZIMJ680102"),
         size_prop = ifelse(size_prop, "Size-related properties", "Properties chosen by AmyloGram"))

cairo_ps("./supplements/figures/violins.eps", height = 6.9, width = 6)
ggplot(prop_plot_dat, aes(x = nice_name, color = status, y = value)) +
  geom_violin(fill = NA, position = position_dodge(1)) +
  geom_boxplot(width = 0.25, position = position_dodge(1), outlier.shape = NA) +
  facet_wrap(~ size_prop, ncol = 1, scales = "free_y") +
  scale_x_discrete("") +
  scale_y_continuous("Normalized value") +
  coord_flip() +
  scale_color_discrete("") +
  my_theme
dev.off()




cairo_ps("./supplements/figures/flex_bulk.eps", height = 7.9, width = 6)
lapply(c("DAWD720101", "BIGC670101", "ZIMJ680102"), function(i) {
  data.frame(status = filter(prop_plot_dat, property == c("BHAR880101"))[["status"]],
             flex = filter(prop_plot_dat, property == c("BHAR880101"))[["value"]],
             size = filter(prop_plot_dat, property == c(i))[["value"]],
             nice = filter(prop_plot_dat, property == c(i))[["nice_name"]])
}) %>% do.call(rbind, .) %>% 
  ggplot(aes(x = flex, y = size)) +
  stat_density2d(aes(fill=status, alpha=..level..), 
                 color = "black",
                 contour = TRUE, geom="polygon") +
  scale_x_continuous("Average flexibility indices (Bhaskaran-Ponnuswamy, 1988)") +
  scale_y_continuous("Size-related property") +
  scale_fill_discrete("") +
  guides(alpha = FALSE) +
  #facet_wrap(~ nice, ncol = 1, scales = "free_x") +
  facet_wrap(~ nice, ncol = 1) +
  my_theme
dev.off()

two_prop <- aaprop[c("BHAR880101", "DAWD720101", "BIGC670101", "ZIMJ680102"), ] %>% 
  t %>% 
  data.frame %>% 
  mutate(aa = rownames(.))

size_prop_names <- data.frame(variable = c("DAWD720101", "BIGC670101", "ZIMJ680102"),
                              nice_name = c("Size\n(Dawson, 1972)",
                                            "Residue volume\n(Bigelow, 1967)",
                                            "Bulkiness\n(Zimmerman et al., 1968)"))


lapply(c("DAWD720101", "BIGC670101", "ZIMJ680102"), function(i) {
  test_res <- cor.test(two_prop[, 1], two_prop[, i])
  data.frame(p_val = test_res[["p.value"]],
             est = unname(test_res[["estimate"]]),
             variable = i)
}) %>% 
  do.call(rbind, .) %>% 
  inner_join(size_prop_names) %>% 
  select(Property = nice_name, `Pearson's correlation` = est, `p-value` = p_val) %>% 
  xtable
  

cairo_ps("./supplements/figures/size_correlations.eps", height = 7.9, width = 6)
melt(two_prop, id.vars = c("BHAR880101", "aa")) %>% 
  inner_join(size_prop_names) %>% 
  ggplot(aes(x = BHAR880101, y = value, label = toupper(aa))) +
  geom_point() +
  geom_text_repel(size = 5) +
  theme_bw() +
  scale_x_continuous("Average flexibility indices (Bhaskaran-Ponnuswamy, 1988)") +
  scale_y_continuous("Normalized value") +
  geom_smooth(se = FALSE, method = "lm", color = "red") +
  facet_wrap(~ nice_name, nrow = 1)
dev.off()