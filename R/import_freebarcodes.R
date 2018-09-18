# Script to import free-barcodes pregenerated barcodes list and put them in data/ as usable objects

library(tidyverse)

bc_files = list.files('/mnt/bigData2/resources/freebarcodes/freebarcodes-master/barcodes',
                      pattern = 'txt',
                      full.names = TRUE)


read_bc_file_write_rdata = function(bc_file) {
  bc_vec = read_tsv(bc_file,
                    col_names = 'bc') %>%
    pull(bc)

  save(bc_vec,
       file = paste0('~/mpradesigntools/data/',
                     stringr::str_extract(bc_file, 'barcodes[0-9]+-[0-9]'),
                     '.RData'))

}

map(bc_files,
    read_bc_file_write_rdata)
