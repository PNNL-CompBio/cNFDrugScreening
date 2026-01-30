#00_cNF_helper_code.R
##standard metadata across all cNFs, including colors if possible


library(synapser)
synLogin()
syn <- list(get = synapser::synGet, store = synapser::synStore)
library(readxl)
library(tidyr)
library(dplyr)

meta1 <- readxl::read_xlsx(syn$get('syn65595365')$path) |>
  tidyr::separate(Specimen,into=c('Patient','Tumor'),sep='_',remove = FALSE)|>
  dplyr::select(Specimen,Patient,Tumor,aliquot)|>
  mutate(cohort=1)

meta2 <- readxl::read_xlsx(syn$get('syn69920464')$path,sheet='Sheet2')|>
  tidyr::separate(Specimen,into=c('Patient','Tumor'),sep='_',remove = FALSE)|>
  dplyr::select(Specimen,Patient,Tumor,aliquot='SampleAlias')|>
  mutate(cohort=2)

meta <- rbind(meta1,meta2) %>%
  filter(!(cohort == 1 & aliquot %in% c(2, 5, 6)))


pcols <- c(NF0017='steelblue',NF0021='orange2',NF0019='orchid4',
           NF0022='goldenrod4',NF0018='olivedrab',NF0020='darkred', NF0022='tan',
           NF0023='darkgrey',NF0025='lightblue',NF0026='yellow3',NF0027='magenta3',
           NF0028='lightgreen',NF0031='pink2')



# This function is used in the notebooks.
# All it does is convert a file from pdf to png so it can be easily displayed after Knitting.
pdf_to_png_if_possible <- function(pdf, dpi = 200) {
  if (!file.exists(pdf)) return(NULL)

  png <- sub("\\.pdf$", ".png", pdf, ignore.case = TRUE)

  # Rebuild png if it doesn't exist, or if the pdf is newer
  rebuild <- !file.exists(png) || (file.info(pdf)$mtime > file.info(png)$mtime)

  if (rebuild) {
    if (requireNamespace("pdftools", quietly = TRUE)) {
      out <- pdftools::pdf_convert(pdf, format = "png", dpi = dpi, pages = 1)
      if (length(out) && file.exists(out[1])) {
        if (file.exists(png)) file.remove(png)
        file.rename(out[1], png)
      }
    } else if (requireNamespace("magick", quietly = TRUE)) {
      img <- magick::image_read_pdf(pdf, density = dpi)
      magick::image_write(img[1], path = png, format = "png")
    } else {
      # If no converter is available, just fall back to the PDF.
      return(pdf)
    }
  }

  if (file.exists(png)) png else pdf
}

# This is used to display plots in the Markdown Notebooks
show_plots <- function(paths, dpi = 200) {
  ok <- paths[file.exists(paths)]
  if (!length(ok)) return(NULL)
  show_paths <- vapply(ok, pdf_to_png_if_possible, FUN.VALUE = character(1), dpi = dpi)
  knitr::include_graphics(show_paths)
}
