#' AmyloGram Graphical User Interface
#'
#' Launches graphical user interface that predicts presence of amyloids.
#'
#' @section Warning : Any ad-blocking software may cause malfunctions.
#' @export AmyloGram_gui

library(seqinr)
library(biogram)

data(AmyloGram_model)

AmyloGram_gui_tuned <- function()
  runApp("inst/AmyloGram")

AmyloGram_gui_tuned()
