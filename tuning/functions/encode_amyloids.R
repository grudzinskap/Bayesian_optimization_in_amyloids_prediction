#' Convert encoding from string to list
#'
#' Converts encodings from string format to list.
#' @param x vector of encodings in the string format.
#' 
#' @return List of encodings in the list format.

string2list <- function(x)
  lapply(x, function(i) {
    pasted_group <- strsplit(i, "_", fixed = TRUE)[[1]]
    res <- strsplit(pasted_group, "")
    names(res) <- 1L:length(res)
    res
  }) 

#' Extract n-grams
#'
#' Extracts n-grams from sequences.
#' @param seq matrix of sequences
#' @param aa_groups list of encodings
#' 
#' @return List of length equal to the length of \code{aa_groups}. Each element 
#' is a list of matrices, where each matrix contains n-grams from a single sequence. 

extract_ngrams <- function(seq, aa_groups) {
  gl <- lapply(1L:nrow(seq), function(i) {
    res <- do.call(rbind, strsplit(decode_ngrams(seq2ngrams(seq[i, ][!is.na(seq[i, ])], 6, a()[-1])), ""))
    cbind(res, id = paste0("P", rep(i, nrow(res))))
  })
  
  #very time consuming
  lapply(aa_groups, function(single_group) 
    as.simple_triplet_matrix(do.call(rbind, lapply(gl, function(single_protein) {
      bitrigrams <- as.matrix(count_multigrams(ns = c(1, rep(2, 4), rep(3, 3)), 
                                               ds = list(0, 0, 1, 2, 3, c(0, 0), c(0, 1), c(1, 0)),
                                               seq = degenerate(single_protein[, -7], single_group),
                                               u = as.character(1L:length(single_group))))
      
      bitrigrams <- bitrigrams > 0
      storage.mode(bitrigrams) <- "integer"
      
      bitrigrams
    })))
  )
}
