#' Normalize properties
#'
#' Normalizes amino acid properties from AAIndex and adds properties 
#' from Wozniak 2014 (also normalized)
#'
#' @return a 550 x 20 matrix of normalized properties.

normalize_properties <- function() {
  data("aaindex")
  
  aa_props <- sapply(aaindex, function(i) i[["I"]])
  tableA <- read.table("./data/tableA.csv", sep = ";", dec = ".", head = TRUE)
  
  #names for additional properties, all names commented out
#   add_names <- c("Values of Wc in proteins from class Beta, cutoff 6 A, separation 5 (Wozniak-Kotulska, 2014)",
#                  "Values of Wc in proteins from class Beta, cutoff 8 A, separation 5 (Wozniak-Kotulska, 2014)",
#                  "Values of Wc in proteins from class Beta, cutoff 12 A, separation 5 (Wozniak-Kotulska, 2014)",
#                  "Values of Wc in proteins from class Beta, cutoff 6 A, separation 15 (Wozniak-Kotulska, 2014)",
#                  "Values of Wc in proteins from class Beta, cutoff 8 A, separation 15 (Wozniak-Kotulska, 2014)",
#                  "Values of Wc in proteins from class Beta, cutoff 12 A, separation 15 (Wozniak-Kotulska, 2014)")
#   traits_names <- c(sapply(prop_MK[["X"]], function(i) aaindex[[i]][["D"]]), 
#                     add_names)
  
  aa_nprop <- t(apply(cbind(aa_props, tableA[tableA[["X"]] %>% as.character %>% aaa %>% order, 2L:7]), 2, function(i) {
    res <- i - min(i, na.rm = TRUE)
    res/max(res, na.rm = TRUE)
  }))
  
  colnames(aa_nprop) <- a(colnames(aa_nprop))
  
  #properties below are hydrophobicity scales that should be reversed
  aa_nprop[c(252, 519, 543, 544), ] <- 1 - aa_nprop[c(252, 519, 543, 544), ]
  
  aa_nprop
}

#' Choose properties
#'
#' Creates a list of traits for encodings
#'
#' @return a vector of trait indices in the expanded aaindex table

choose_properties <- function() {
  
  aa_nprop <- normalize_properties()
  
  #key for selecting properties - new (younger than 1980 years properties)
  prop_MK <- read.csv2("./data/AA_index_mk2.csv") %>% filter(!is.na(chosen))
  
  years <- prop_MK %>% select(name) %>% unlist %>% as.character %>% sapply(function(i) 
    strsplit(last(strsplit(i, ", ")[[1]]), ")", fixed = TRUE)[[1]][1]) %>%
    as.numeric
  
  prop_MK <- cbind(prop_MK, years = years) %>% filter(years >= 1980)
  
  traits <- c(prop_MK[["X"]], 545:550)
  
  #final traits
  #removed redundant hydrophobicities and interactivities. Left interactivities and hydrophobicities
  #with the highest correlation with other parameters
  traits[-c(6, 13, 18, 19, 20, 21, 22, 24, 25)]
}