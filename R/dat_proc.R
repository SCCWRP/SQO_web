######
# dat proc

library(tidyverse)


# indicator species lookup table ------------------------------------------

indic_lookup <- read.csv('data/raw/indic_lookup.csv', stringsAsFactors = F)

save(indic_lookup, file = 'data/indic_lookup.RData', compress = 'xz')


# write shiny contaminant widgets to txt ----------------------------------

# write shiny widgets for contaminants to text file for each copy/paste to index

# file with 'default' contamint inputs
cnt_bld <- read.csv('data/raw/contam_build.csv', stringsAsFactors = F)
cnt_bld[is.na(cnt_bld)] <- 'NULL'

sink('data/raw/contamwidg.txt', append = T)
for(i in 1:nrow(cnt_bld)){
  
  grp <- cnt_bld[i, 'ChemGroup', drop = T]
  cnt <- cnt_bld[i, 'Chem', drop = T]
  sed <- cnt_bld[i, 'Sediment.conc', drop = T]
  dis <- cnt_bld[i, 'Dissolved.surface', drop = T]
  por <- cnt_bld[i, 'Porewater', drop = T]
  
  out <- paste0('# ', grp, '\n', 
                'fluidRow(\n\tcolumn(2, h6(\'', grp, '\')),\n',
                '\tcolumn(2, h6(\'', cnt, '\')),\n',
                '\tcolumn(2, numericInput(\'', cnt, '_sed\', NULL, value = ', sed, ', min = 0, width = \'100%\')),\n',
                '\tcolumn(2, numericInput(\'', cnt, '_dis\', NULL, value = ', dis, ', min = 0, width = \'100%\')),\n',
                '\tcolumn(2, numericInput(\'', cnt, '_por\', NULL, value = ', por, ', min = 0, width = \'100%\'))\n',
                '\t),\n\n'
  )
  
  cat(out)
  
}
sink()