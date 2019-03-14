######
# dat proc

library(tidyverse)

##
# indicator species lookup table

indic_lookup <- read.csv('data/raw/indic_lookup.csv', stringsAsFactors = F)

save(indic_lookup, file = 'data/indic_lookup.RData', compress = 'xz')