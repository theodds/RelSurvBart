tikzprint <- function(gg, filename, ...) {
  filename_tex <- paste0(filename, ".tex")
  filename_pdf <- paste0(filename, ".pdf")
  tikz(file = filename_tex, standAlone = TRUE, ...)
  plot(gg)
  dev.off()
  tools::texi2pdf(filename_tex, clean = TRUE)
  file.rename(filename_pdf, paste0("figures/", filename_pdf))
  file.remove(filename_tex)
  return(TRUE)
}

