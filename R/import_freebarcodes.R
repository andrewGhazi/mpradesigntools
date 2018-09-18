# Script to import free-barcodes pregenerated barcodes list and put them in data/ as usable objects

library(tidyverse)
library(parallel)

bc_files = list.files('/mnt/bigData2/resources/freebarcodes/freebarcodes-master/barcodes',
                      pattern = 'txt',
                      full.names = TRUE)


read_bc_file_write_rdata = function(bc_file) {
  set_name = stringr::str_extract(bc_file, 'barcodes[0-9]+-[0-9]')
  assign(set_name,
         read_tsv(bc_file,
                  col_names = 'bc') %>%
           pull(bc))

  save(list = c(set_name),
       file = paste0('~/mpradesigntools/data/',
                     set_name,
                     '.RData'))

}

mclapply(bc_files,
         read_bc_file_write_rdata,
         mc.cores = 10)

#### Make the table for the readme ----

count_bc_file_lines = function(bc_file) {
  system(paste0('wc -l ',
                bc_file),
         intern = TRUE) %>%
    stringr::str_extract('^[0-9]+') %>%
    as.numeric
}

data_frame(barcode_set = map_chr(bc_files,
                                 ~stringr::str_extract(.x, 'barcodes[0-9]+-[0-9]')),
           n_barcodes = map_dbl(bc_files,
                                count_bc_file_lines)) %>%
  add_case(barcode_set = 'twelvemers',
           n_barcodes = 1140292) %>%
  arrange(barcode_set) %>%
  knitr::kable()

#### Add documentation for barcode sets ----

write_data_documentation = function(bc_file){
  set_name = stringr::str_extract(bc_file, 'barcodes[0-9]+-[0-9]')

  write(paste0('#\' ', set_name, '\n', '\'', set_name, '\'\n'),
        file = 'R/data.R',
        append = TRUE)

  return('doneskies')
}

map(bc_files,
    write_data_documentation)
