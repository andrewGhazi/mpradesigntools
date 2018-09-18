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

# mclapply(bc_files,
#          read_bc_file_write_rdata,
#          mc.cores = 10)

#### Make the table for the readme ----

count_bc_file_lines = function(bc_file) {
  system(paste0('wc -l ',
                bc_file),
         intern = TRUE) %>%
    stringr::str_extract('^[0-9]+') %>%
    as.numeric
}

# data_frame(barcode_set = map_chr(bc_files,
#                                  ~stringr::str_extract(.x, 'barcodes[0-9]+-[0-9]')),
#            n_barcodes = map_dbl(bc_files,
#                                 count_bc_file_lines)) %>%
#   add_case(barcode_set = 'twelvemers',
#            n_barcodes = 1140292) %>%
#   arrange(barcode_set) %>%
#   knitr::kable()

#### Add documentation for barcode sets ----

write_data_documentation = function(bc_file){
  set_name = stringr::str_extract(bc_file, 'barcodes[0-9]+-[0-9]')

  write(paste0('#\' ', set_name, '\n', '\'', set_name, '\'\n'),
        file = 'R/data.R',
        append = TRUE)

  return('doneskies')
}

# map(bc_files,
#     write_data_documentation)

#### scan for mpra-intert ----
# let's double check that these new barcodes

matches_all_nucleotides = function(bc){
  # https://stackoverflow.com/questions/469913/regular-expressions-is-there-an-and-operator
  grepl('(?=.*A)(?=.*C)(?=.*G)(?=.*T)', bc,
        perl = TRUE)
}

load_and_count = function(set_name){
  load(list.files('~/mpradesigntools/data',
                  pattern = set_name,
                  full.names = TRUE))
  bc_set = get(set_name)

  sum(unlist(mclapply(bc_set,
               matches_all_nucleotides,
               mc.cores = 12)))
}

# library(Biostrings)
# nucruns = vector(mode = 'character', length = 4) %>% DNAStringSet
# ni = 1
# for (i in 4) {
#   for (j in c('A', 'G', 'T', 'C')) {
#     nucruns[ni] = rep(j, i) %>% paste(collapse = '') %>% DNAStringSet
#     ni = ni + 1
#   }
# }

count_nuc_runs = function(set_name){
  load(list.files('~/mpradesigntools/data',
                  pattern = set_name,
                  full.names = TRUE))

  bc_set = get(set_name)

  Biostrings::vcountPDict(nucruns, DNAStringSet(bc_set)) %>%
    colSums %>%
    {. > 0} %>%
    sum
}

count_mir_seeds = function(set_name){

  human_seeds = readRNAStringSet('~/plateletMPRA/data/mature.fa') %>%
    {.[grepl('Homo sapiens', names(.))]} %>%
    unique %>%
    subseq(2,7) %>%
    DNAStringSet %>%
    as.data.frame() %>%
    rownames_to_column() %>%
    as_tibble()

  load(list.files('~/mpradesigntools/data',
                  pattern = set_name,
                  full.names = TRUE))

  bc_set = get(set_name)

  vwhichPDict(DNAStringSet(unique(human_seeds$x)),
              DNAStringSet(bc_set)) %>% #this takes ~40 minutes. All seeds takes ~1h45m
    map_lgl(~length(.x) != 0) %>%
    sum
}


# set_properties = data_frame(barcode_set = map_chr(bc_files,
#                                  ~stringr::str_extract(.x, 'barcodes[0-9]+-[0-9]')),
#            n_barcodes = map_dbl(bc_files,
#                                 count_bc_file_lines)) %>%
#   add_case(barcode_set = 'twelvemers',
#            n_barcodes = 1140292) %>%
#   mutate(frac_with_all = map_dbl(barcode_set,
#                                  load_and_count) / n_barcodes,
#          frac_with_nuc_runs = unlist(mclapply(barcode_set,
#                                               count_nuc_runs,
#                                               mc.cores = 12)) / n_barcodes,
#          frac_with_mir_seeds = unlist(mclapply(barcode_set,
#                                                count_mir_seeds,
#                                                mc.cores = 15)) / n_barcodes)
#
# save(set_properties,
#      file = '~/plateletMPRA/outputs/set_properties.RData')

# So none of the new sets have nucleotide runs, and the larger ones have all 4
# nucleotides ~85% - 98% of the time

