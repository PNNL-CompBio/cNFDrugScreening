##standard metadata across all cNFs, including colors if possible


library(synapser)
synLogin()
library(readxl)


meta1 <- readxl::read_xlsx(syn$get('syn65595365')$path) |>
  tidyr::separate(Specimen,into=c('Patient','Tumor'),sep='_',remove = FALSE)|>
  dplyr::select(Specimen,Patient,Tumor,aliquot)|>
  mutate(cohort=1)

meta2 <- readxl::read_xlsx(syn$get('syn69920464')$path,sheet='Sheet2')|>
  tidyr::separate(Specimen,into=c('Patient','Tumor'),sep='_',remove = FALSE)|>
  dplyr::select(Specimen,Patient,Tumor,aliquot='SampleAlias')|>
  mutate(cohort=2)

meta <- rbind(meta1,meta2)



pcols <- c(NF0017='steelblue',NF0021='orange2',NF0019='orchid4',
           NF0022='goldenrod4',NF0018='olivedrab',NF0020='darkred', NF0022='tan',
           NF0023='darkgrey',NF0025='lightblue',NF0026='yellow3',NF0027='magenta3',
           NF0028='lightgreen',NF0031='pink2')

