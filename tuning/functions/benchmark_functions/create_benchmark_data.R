library(seqinr)
pep424_tab <- read.table("./benchmark/pep424_evaluation.txt", sep = "\t")
levels(pep424_tab[[3]]) <- c(0, 0, 1, 1)


pep424_prots <- lapply(as.character(pep424_tab[[2]]) %>% strsplit(split = ""), function(i)
  i[i != " "])

# write pasta files with seqs
write.fasta(pep424_prots[lengths(pep424_prots) > 5], 
            names = as.character(pep424_tab[[1]]), file.out = "pep424.fasta")

write.fasta(pep424_prots[lengths(pep424_prots) > 5], 
            names = paste0("prot", 1L:sum(lengths(pep424_prots) > 5)), file.out = "pep424_better_names.fasta")

#reads data and removes seqs already present in test set
get_train_seqs <- function(file_name, pattern_vector) {
  seq <- read.fasta(file_name, seqtype = "AA")
  
  seq_pasted <- sapply(seq, function(i) paste0(i, collapse = ""))
  
  seq_duplicated <- lapply(seq_pasted, function(i)
    grep(i, pattern_vector)) %>% unlist %>% unique
  
  seq[-seq_duplicated]
}

seq_pos <- get_train_seqs("./data/amyloid_pos_full.fasta", sapply(pep424_prots, paste0, collapse = ""))
seq_neg <- get_train_seqs("./data/amyloid_neg_full.fasta", sapply(pep424_prots, paste0, collapse = ""))

write.fasta(seq_pos, names = names(seq_pos), file.out = "./data/amyloid_pos_benchmark.fasta")
write.fasta(seq_neg, names = names(seq_neg), file.out = "./data/amyloid_neg_benchmark.fasta")