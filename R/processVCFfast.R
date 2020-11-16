#Function to take an input vcf and generate an output vcf with the specified up/downstream context ranges and # of barcodes per snp
# And now let's do it quickly using dplyr and purrr. processVCF

#' @importFrom stringr str_split
spreadAllelesAcrossRows = function(snp){
  #snp is a row from a vcf data_frame
  # if the 'ALT' column has commas in it, spread those out across otherwise identical rows

  if (!grepl(',', snp$ALT)) {
    return(snp)
  } else {
    altAlleles = snp$ALT %>% stringr::str_split(',', simplify = TRUE)
    res = snp[rep(1, length(altAlleles)),]
    res$ALT = altAlleles %>% as.vector
    return(res)
  }
}

#' @importFrom Biostrings countPattern
#' @importFrom Biostrings DNAString
#' @importFrom Biostrings reverse
countDigSites = function(biostring, enzyme1, enzyme2, enzyme3) {
  #Count the number of times enzyme1, XbaI, and SfiI sites occur in a given biostring

  # enzyme1 = DNAString('GGTACC') #KpnI
  # xba = DNAString('TCTAGA') #XbaI
  # sfi = DNAString('GGCCNNNNNGGCC') #SfiI

  sum(countPattern(enzyme1, biostring, fixed = FALSE),
      countPattern(enzyme2, biostring, fixed = FALSE),
      countPattern(enzyme3, biostring, fixed = FALSE),
      countPattern(enzyme1 %>% reverse, biostring, fixed = FALSE),
      countPattern(enzyme2 %>% reverse, biostring, fixed = FALSE),
      countPattern(enzyme3 %>% reverse, biostring, fixed = FALSE))
}

#' @importFrom Biostrings subseq
#' @importFrom Biostrings toString
generateInsConstruct = function(snpseq, mid, reverseGene, upstreamContextRange, downstreamContextRange){
  #This function generates a mutant construct sesquence based on
  # snpseq - the genomic context
  # mid - the insertion allele
  # reverseGene - a logical indicating whether or not the SNP is for a gene that's reversed
  # upstreamContextRange - the length of context upstream of the SNP
  # downstreamContextRange - the length of context downstream of the SNP
  #   this is necessary for reverse genes on the - strand to accurately specify the start point of the insertion

  # If the insertion is for a gene that's transcribed from the reverse strand, it needs to be the COMPLEMENT of the allele to the LEFT of the position
  # Otherwise it's the mutant allele to the right of the position

  # if (reverseGene) { # Deprecating this check because snpseq is pre-flipped now
  #   toString(c(subseq(snpseq, 1, downstreamContextRange), #The insertion goes to the right of the position given
  #              Biostrings::complement(mid),
  #              subseq(snpseq, downstreamContextRange + 1, length(snpseq))))
  # } else {
  #   toString(c(subseq(snpseq, 1, upstreamContextRange + 1), #The insertion goes to the right of the position given
  #              mid,
  #              subseq(snpseq, upstreamContextRange + 2, length(snpseq))))
  # }
  toString(c(subseq(snpseq, 1, upstreamContextRange + 1), #The insertion goes to the right of the position given
             mid,
             subseq(snpseq, upstreamContextRange + 2, length(snpseq))))
}

#' @importFrom Biostrings subseq
generateDelConstruct = function(snpseq, refwidth, upstreamContextRange) {
  c(subseq(snpseq,
           1,
           upstreamContextRange),
    subseq(snpseq,
           upstreamContextRange + refwidth + 1,
           length(snpseq)))
}

randomly_change_pattern = function(dig_pattern){
  site_to_change = sample(1:nchar(dig_pattern), size = 1)
  old_allele = stringr::str_sub(dig_pattern, site_to_change, site_to_change)
  allele_options = c('A', 'T', 'C', 'G')[!(c('A', 'T', 'C', 'G') %in% old_allele)]

  new_allele = sample(allele_options, size = 1)

  stringr::str_sub(dig_pattern, start = site_to_change, end = site_to_change) = new_allele
  return(dig_pattern)
}

reassign_pattern = function(construct_seq, aberrant_site_loc, replacement) {
  stringr::str_sub(construct_seq,
          start = BiocGenerics::start(aberrant_site_loc),
          end = BiocGenerics::end(aberrant_site_loc)) = replacement
  return(construct_seq)
}

change_pattern = function(ab_pattern,
                          pos_to_change,
                          allele_to_use){
  substr(ab_pattern,
         pos_to_change,
         pos_to_change) = allele_to_use

  return(ab_pattern)
}

#' Randomly correct aberrant digestion sites
#'
#' For a SNP with aberrant digestion sites in the context, randomly change bases
#' in the site across barcodes
#'
#' @param snp a data_frame of one SNP
#' #' @param nper The number of barcoded sequences to be generated per allele per
#'   SNP
#' @param upstreamContextRange the amount of sequence context to acquire upstream of the SNP
#' @param downstreamContextRange the amount of sequence context to acquire downstream of the SNP
#' @param fwprimer a string containing the forward PCR primer to be used
#' @param revprimer a string containing the reverse PCR primer to be used
randomly_fix = function(snp,
                        res_df,
                        dig_patterns,
                        dig_site_locations){

  dig_sites_present = purrr::map_lgl(dig_site_locations, ~length(BiocGenerics::width(.x)) != 0)

  if (sum(dig_sites_present) > 1) {
    res_df = data_frame(ID = snp$ID,
                            CHROM = snp$CHROM,
                            POS = snp$POS,
                            REF = snp$REF,
                            ALT = snp$ALT,
                            result = 'Failed - More than one aberrant digestion site in context')
    return(res_df)
  }

  aberrant_in_context = any(dig_sites_present)

  aberrant_site = dig_site_locations[which(dig_sites_present)]  %>% .[[1]]

  if (aberrant_in_context) {

    aberrant_pattern = dig_patterns[which(dig_sites_present)]

    # This assures that the changes are unique, if possible
    altered_patterns = data_frame(pos_to_change = 1:nchar(aberrant_pattern),
                                  possible_alleles = map(pos_to_change, ~dplyr::setdiff(c('A', 'C', 'G', 'T'),
                                                                                        substr(aberrant_pattern,
                                                                                               .x, .x)))) %>%
      tidyr::unnest_legacy() %>%
      {dplyr::sample_n(.,
                       nrow(res_df) / 2,
                       replace = (nrow(res_df) / 2 > nrow(.)))} %>%
      mutate(altered_pattern = map2_chr(pos_to_change, possible_alleles,
                                        change_pattern,
                                        ab_pattern = aberrant_pattern))

    res_df %<>%
      mutate(aberrant_pattern = aberrant_pattern,
             fixed_pattern = rep(altered_patterns$altered_pattern, times = 2), # give the same altered patterns to both alleles
             fixed_pattern_index = rep(1:(dplyr::n()/2), times = 2),
             constrseq_fixed = map2_chr(constrseq, fixed_pattern,
                                        reassign_pattern,
                                        aberrant_site_loc = aberrant_site))
  } else {
    res_df = data_frame(ID = snp$ID,
                            CHROM = snp$CHROM,
                            POS = snp$POS,
                            REF = snp$REF,
                            ALT = snp$ALT,
                            result = 'Failed - aberrant site generated at interface between sequence elements')
    return(res_df)
  }

  return(res_df)

}

#' process an individual SNP
#'
#' Take one SNP, get its genomic context, concatenate the parts, and check it
#' doesn't contain aberrant digestion sites. If it does, resample the barcodes
#' and try a few more times. If it still does, return a failure stating why.
#'
#' @param snp a data_frame containing the VCF information for one SNP as well as
#'   a barcode pool to sample from.
#' @param nper The number of barcoded sequences to be generated per allele per
#'   SNP
#' @param upstreamContextRange the amount of sequence context to acquire upstream of the SNP
#' @param downstreamContextRange the amount of sequence context to acquire downstream of the SNP
#' @param fwprimer a string containing the forward PCR primer to be used
#' @param revprimer a string containing the reverse PCR primer to be used
#' @param enzyme1 the first restriction enzyme's recognition pattern
#' @param enzyme2 the first restriction enzyme's recognition pattern
#' @param enzyme3 the first restriction enzyme's recognition pattern
#' @return a data_frame of labeled sequences with appropriate information on the changes made
#' @import BSgenome.Hsapiens.UCSC.hg38
#' @importFrom Biostrings reverseComplement
#' @importFrom Biostrings replaceLetterAt
#' @importFrom tibble data_frame
processSnp = function(snp,
                      nper,
                      upstreamContextRange,
                      downstreamContextRange,
                      fwprimer,
                      revprimer,
                      enzyme1,
                      enzyme2,
                      enzyme3,
                      alter_aberrant = FALSE,
                      extra_elements = FALSE,
                      max_construct_size = NULL,
                      flip_RV = TRUE){
  # snp is one row from the expanded vcf, including the reverseGene column which
  # indicates whether or not to use the reverse complement genomic context. It
  # also has a dedicated pool of barcodes to select from

  genome = BSgenome.Hsapiens.UCSC.hg38
  # kpn = 'GGTACC' #KpnI
  # xba = 'TCTAGA' #XbaI
  # sfi = 'GGCCNNNNNGGCC' #SfiI

  isSNV = (snp$REF %in% c('A', 'C', 'G', 'T') && snp$ALT %in% c('A', 'C', 'G', 'T')) #look at me interchanging between SNV and SNP willy-nilly ####
  isINS = snp$REF == '-'
  isDEL = snp$ALT == '-'

  notes = NULL # initialize


  #There are three code blocks below for each of these cases. Use RStudio's code folding to open up only the one of interest.
  # This could probably be made to be 1/3 the length by writing a function that
  # adapts to the type of SNP when generating the sequences, but that might be
  # hard. There are subtle differences in each of the types of variants.

  # these get passed later to randomly_fix
  dig_patterns = c(enzyme1 = enzyme1,
                   enzyme2 = enzyme2,
                   enzyme3 = enzyme3,
                   enzyme1_rev = enzyme1 %>% reverse,
                   enzyme2_rev = enzyme2 %>% reverse,
                   enzyme3_rev = enzyme3 %>% reverse)

  bc_length = nchar(snp$bcPools[[1]][1])
  #### check the construct size, shorten if applicable ----
  if (isINS) {
    # insertions are 1bp longer because the ref is empty
    tot_construct_length = sum(nchar(fwprimer),
                               nchar(revprimer),
                               upstreamContextRange,
                               downstreamContextRange,
                               max(nchar(snp$REF) + 1 , nchar(snp$ALT) + 1),
                               nchar(enzyme1),
                               nchar(enzyme2),
                               bc_length)
  } else if (isDEL) {
    tot_construct_length = sum(nchar(fwprimer),
                               nchar(revprimer),
                               upstreamContextRange,
                               downstreamContextRange,
                               max(-nchar(snp$REF), -nchar(snp$ALT)),
                               nchar(enzyme1),
                               nchar(enzyme2),
                               bc_length)
  } else {
    tot_construct_length = sum(nchar(fwprimer),
                               nchar(revprimer),
                               upstreamContextRange,
                               downstreamContextRange,
                               max(nchar(snp$REF), nchar(snp$ALT)),
                               nchar(enzyme1),
                               nchar(enzyme2),
                               bc_length)
  }

  if (extra_elements) {
    # extra elements are the TG & GGC Namrata added to CD36 MPRA for some reason
    tot_construct_length = tot_construct_length + 5
  }

  if (flip_RV & grepl('RV', snp$INFO)) {
    notes = c(notes, 'The alleles for this SNP were flipped from the input VCF because of the presence of the RV tag in the INFO field. This feature is new so please double-check the result. ')

    if (isSNV){

      snp$REF = snp$REF %>%
        Biostrings::DNAString() %>%
        Biostrings::complement() %>%
        Biostrings::toString()

      snp$ALT = snp$ALT %>%
        Biostrings::DNAString() %>%
        Biostrings::complement() %>%
        Biostrings::toString()
    } else if (isINS) {

      # If it's an insertion, REF = '-', so you only need to flip ALT
      snp$ALT = snp$ALT %>%
        Biostrings::DNAString() %>%
        Biostrings::complement() %>%
        Biostrings::toString()
    } else {
      # else it's a deletion
      # if it's a deletion, ALT = '-', so you only need to flip REF
      snp$REF = snp$REF %>%
        Biostrings::DNAString() %>%
        Biostrings::complement() %>%
        Biostrings::toString()
    }
  }

  if (!is.null(max_construct_size) && tot_construct_length > max_construct_size) {
    excess = tot_construct_length - max_construct_size
    amount_to_remove = ceiling(excess / 2)

    upstreamContextRange = upstreamContextRange - amount_to_remove
    downstreamContextRange = downstreamContextRange - amount_to_remove
    notes = c(notes, paste0(' Shortened context by ', amount_to_remove, ' bp on each side to account for input maximum construct size. '))
  }

  #### Start generating the constructs depending on SNP type
  if (isSNV) {

    rangestart = snp$POS - upstreamContextRange
    rangeend = snp$POS + downstreamContextRange

    if (snp$reverseGene) {
      rangestart = snp$POS - downstreamContextRange
      rangeend = snp$POS + upstreamContextRange
    }

    snpseq = subseq(genome[[paste0('chr', as.character(snp$CHROM))]], # the chrom field needs to be only digits
                    start = rangestart,
                    end = rangeend)
    if (snp$reverseGene) {
      snpseq %<>% reverseComplement()
    }

    ndigsite = countDigSites(snpseq, enzyme1, enzyme2, enzyme3)
    ndigsite_in_context = ndigsite

    if (ndigsite > 0) {
      if (alter_aberrant) {

        # the digestion sites get randomly fixed later, so this is all that's done here
        dig_site_locations = purrr::map(dig_patterns, Biostrings::matchPattern, subject = snpseq, fixed = FALSE)
      } else {
        failureRes = data_frame(ID = snp$ID,
                                CHROM = snp$CHROM,
                                POS = snp$POS,
                                REF = snp$REF,
                                ALT = snp$ALT,
                                result = 'Failed - Context contained a digestion site')
        return(failureRes)
      }

    }

    refseq = toString(snpseq)
    altseq = toString(replaceLetterAt(snpseq, upstreamContextRange + 1, snp$ALT))

    res = data_frame(ID = snp$ID,
                     CHROM = snp$CHROM,
                     snpIndex = 1:(nper*2),
                     type = rep(c('ref', 'alt'), each = nper),
                     mid = ifelse(type == 'ref', snp$REF, snp$ALT),
                     barcodes = sample(snp$bcPools %>% unlist,
                                       2*nper,
                                       replace = FALSE),
                     constrseq = type %>% purrr::map_chr(~ifelse(.x == 'ref',
                                                          refseq,
                                                          altseq)))

    if (extra_elements) {
      res %<>% mutate(sequence = paste0(fwprimer,
                                        'TG',
                                        constrseq,
                                        enzyme1,
                                        enzyme2,
                                        barcodes,
                                        'GGC',
                                        revprimer),
                      ndigSites = sequence %>% purrr::map_int(~countDigSites(DNAString(.x), enzyme1, enzyme2, enzyme3)))
    } else {
      res %<>% mutate(sequence = paste0(fwprimer,
                                        constrseq,
                                        enzyme1,
                                        enzyme2,
                                        barcodes,
                                        revprimer),
                      ndigSites = sequence %>% purrr::map_int(~countDigSites(DNAString(.x), enzyme1, enzyme2, enzyme3)))
    }

    #If all of the sequences contained > 3 digestion sites, there's probably
    #some location at the context/other parts boundary that generates a site.
    #This is too complicated to fix automatically, so just fail the SNP
    if (all(res$ndigSites >= 3)) {

      if (alter_aberrant & ndigsite_in_context > 0) {

        if (!exists('dig_site_locations')) {

          dig_site_locations = purrr::map(dig_patterns, Biostrings::matchPattern, subject = snpseq, fixed = FALSE)

          # dig_site_locations = tidyr::crossing(dig_pattern = dig_patterns, constr_seq = res$constrseq) %>%
          #   mutate(pattern_loc = purrr::map2(dig_pattern, constr_seq,
          #                             ~Biostrings::matchPattern(.x, subject = DNAString(.y), fixed = FALSE))) %>%
          #   pull(pattern_loc) %>%
          #   unique
        }

        multiple_aberrant_dig_sites = sum(purrr::map_int(dig_site_locations,
                                                  ~length(BiocGenerics::width(.x)))) > 1

        if (multiple_aberrant_dig_sites) {
          failureRes = data_frame(ID = snp$ID,
                                  CHROM = snp$CHROM,
                                  POS = snp$POS,
                                  REF = snp$REF,
                                  ALT = snp$ALT,
                                  result = 'Failed - SNP sequence could not be generated without multiple aberrant digestion sites')
          return(failureRes)
        }

        if (extra_elements) {
          res = randomly_fix(snp,
                             res,
                             dig_patterns,
                             dig_site_locations) %>%
            mutate(sequence = paste0(fwprimer,
                                     'TG',
                                     constrseq_fixed,
                                     enzyme1,
                                     enzyme2,
                                     barcodes,
                                     'GGC',
                                     revprimer),
                   ndigSites = sequence %>% purrr::map_int(~countDigSites(DNAString(.x), enzyme1, enzyme2, enzyme3))) %>%
            tidyr::nest_legacy(aberrant_pattern:constrseq_fixed,
                        .key = 'site_fix_info')
        } else {
          res = randomly_fix(snp,
                             res,
                             dig_patterns,
                             dig_site_locations) %>%
            mutate(sequence = paste0(fwprimer,
                                     constrseq_fixed,
                                     enzyme1,
                                     enzyme2,
                                     barcodes,
                                     revprimer),
                   ndigSites = sequence %>% purrr::map_int(~countDigSites(DNAString(.x), enzyme1, enzyme2, enzyme3))) %>%
            tidyr::nest_legacy(aberrant_pattern:constrseq_fixed,
                        .key = 'site_fix_info')
        }
      } else {
        failureRes = data_frame(ID = snp$ID,
                                CHROM = snp$CHROM,
                                POS = snp$POS,
                                REF = snp$REF,
                                ALT = snp$ALT,
                                result = 'Failed - SNP sequence could not be generated without an aberrant digestion site in all constructs')
        return(failureRes)
      }
    }

    if (any(res$ndigSites > 3)) {
      #divide the results into broken and working sequences
      working = res %>% filter(ndigSites < 3)
      broken = res %>% filter(ndigSites >= 3)

      #remove the barcodes for the sequences that work
      brokenPool = unlist(snp$bcPools)[!(unlist(snp$bcPools) %in% working$barcodes)]

      resample_attempts = 0
      while (any(res$ndigSites > 3)) {
        resample_attempts = resample_attempts + 1

        if (resample_attempts > 40) {

          failureRes = data_frame(ID = snp$ID,
                                  CHROM = snp$CHROM,
                                  POS = snp$POS,
                                  REF = snp$REF,
                                  ALT = snp$ALT,
                                  result = 'Failed - SNP sequence could not be generated without an aberrant digestion site even after barcode resampling')
          return(failureRes)
        }

        #For the subset of sequences that don't work, resample the barcodes and try again.
        if (extra_elements) {
          fixed = broken %>% mutate(barcodes = sample(brokenPool,
                                                      nrow(broken)),
                                    sequence = paste0(fwprimer,
                                                      'TG',
                                                      constrseq,
                                                      enzyme1,
                                                      enzyme2,
                                                      barcodes,
                                                      'GGC',
                                                      revprimer),
                                    ndigSites = sequence %>% purrr::map_int(~countDigSites(DNAString(.x), enzyme1, enzyme2, enzyme3)))
        } else {
          fixed = broken %>% mutate(barcodes = sample(brokenPool,
                                                      nrow(broken)),
                                    sequence = paste0(fwprimer,
                                                      constrseq,
                                                      enzyme1,
                                                      enzyme2,
                                                      barcodes,
                                                      revprimer),
                                    ndigSites = sequence %>% purrr::map_int(~countDigSites(DNAString(.x), enzyme1, enzyme2, enzyme3)))
        }
        res = rbind(working, fixed)
      }
    }
  } else if (isINS) {

    #If the snp is an insertion, the range of context to get is the same, but the middle allele (variable 'mid') is different
    rangestart = snp$POS - upstreamContextRange
    rangeend = snp$POS + downstreamContextRange

    if (snp$reverseGene) {
      rangestart = snp$POS - downstreamContextRange
      rangeend = snp$POS + upstreamContextRange
    }

    snpseq = subseq(genome[[paste0('chr', as.character(snp$CHROM))]], # the chrom field needs to be only digits
                    start = rangestart,
                    end = rangeend)

    if (snp$reverseGene) {
      snpseq %<>% reverseComplement()
    }

    ndigsite = countDigSites(snpseq, enzyme1, enzyme2, enzyme3)
    ndigsite_in_context = ndigsite

    if (ndigsite > 0) {
      if (alter_aberrant) {
        dig_patterns = c(enzyme1,
                         enzyme2,
                         enzyme3,
                         enzyme1 %>% reverse,
                         enzyme2 %>% reverse,
                         enzyme3 %>% reverse)

        dig_site_locations = purrr::map(dig_patterns, Biostrings::matchPattern, subject = snpseq, fixed = FALSE)
      } else {
        failureRes = data_frame(ID = snp$ID,
                                CHROM = snp$CHROM,
                                POS = snp$POS,
                                REF = snp$REF,
                                ALT = snp$ALT,
                                result = 'Failed - Context contained a digestion site')
        return(failureRes)
      }
    }

    altseq = generateInsConstruct(snpseq, DNAString(snp$ALT), snp$reverseGene, upstreamContextRange, downstreamContextRange)

    if (extra_elements) {
      res = data_frame(ID = snp$ID,
                       CHROM = snp$CHROM,
                       snpIndex = 1:(nper*2),
                       type = rep(c('ref', 'alt'), each = nper),
                       mid = ifelse(type == 'ref', '', snp$ALT), #This line is line is unique to the isINS block
                       barcodes = sample(snp$bcPools %>% unlist,
                                         2*nper),
                       constrseq = purrr::map2_chr(mid, type, ~ifelse(.y == 'ref',
                                                                      toString(snpseq),
                                                                      altseq)),
                       sequence = paste0(fwprimer,
                                         'TG',
                                         constrseq,
                                         enzyme1,
                                         enzyme2,
                                         barcodes,
                                         'GGC',
                                         revprimer),
                       ndigSites = sequence %>% purrr::map_int(~countDigSites(DNAString(.x), enzyme1, enzyme2, enzyme3)))
    } else {
      res = data_frame(ID = snp$ID,
                       CHROM = snp$CHROM,
                       snpIndex = 1:(nper*2),
                       type = rep(c('ref', 'alt'), each = nper),
                       mid = ifelse(type == 'ref', '', snp$ALT), #This line is line is unique to the isINS block
                       barcodes = sample(snp$bcPools %>% unlist,
                                         2*nper),
                       constrseq = purrr::map2_chr(mid, type, ~ifelse(.y == 'ref',
                                                                      toString(snpseq),
                                                                      altseq)),
                       sequence = paste0(fwprimer,
                                         constrseq,
                                         enzyme1,
                                         enzyme2,
                                         barcodes,
                                         revprimer),
                       ndigSites = sequence %>% purrr::map_int(~countDigSites(DNAString(.x), enzyme1, enzyme2, enzyme3)))
    }

    #If all of the sequences contained > 3 digestion sites, there's probably
    #some location at the context/other parts boundary that generates a site.
    #This is too complicated to fix automatically, so just fail the SNP
    if (all(res$ndigSites >= 3)) {
      if (alter_aberrant & ndigsite_in_context > 0) {

        if (!exists('dig_site_locations')) {

          dig_site_locations = purrr::map(dig_patterns, Biostrings::matchPattern, subject = snpseq, fixed = FALSE)

          # dig_site_locations = tidyr::crossing(dig_pattern = dig_patterns, constr_seq = res$constrseq) %>%
          #   mutate(pattern_loc = purrr::map2(dig_pattern, constr_seq,
          #                             ~Biostrings::matchPattern(.x, subject = DNAString(.y), fixed = FALSE))) %>%
          #   pull(pattern_loc) %>%
          #   unique
        }

        multiple_aberrant_dig_sites = sum(purrr::map_lgl(dig_site_locations,
                                                  ~length(BiocGenerics::width(.x)) != 0)) > 1

        if (multiple_aberrant_dig_sites) {
          failureRes = data_frame(ID = snp$ID,
                                  CHROM = snp$CHROM,
                                  POS = snp$POS,
                                  REF = snp$REF,
                                  ALT = snp$ALT,
                                  result = 'Failed - SNP sequence could not be generated without multiple aberrant digestion sites')
          return(failureRes)
        }

        if (extra_elements) {
          res = randomly_fix(snp,
                             res,
                             dig_patterns,
                             dig_site_locations) %>%
            mutate(sequence = paste0(fwprimer,
                                     'TG',
                                     constrseq_fixed,
                                     enzyme1,
                                     enzyme2,
                                     barcodes,
                                     'GGC',
                                     revprimer),
                   ndigSites = sequence %>% purrr::map_int(~countDigSites(DNAString(.x), enzyme1, enzyme2, enzyme3))) %>%
            tidyr::nest_legacy(aberrant_pattern:constrseq_fixed,
                        .key = 'site_fix_info')
        } else {
          res = randomly_fix(snp,
                             res,
                             dig_patterns,
                             dig_site_locations) %>%
            mutate(sequence = paste0(fwprimer,
                                     constrseq_fixed,
                                     enzyme1,
                                     enzyme2,
                                     barcodes,
                                     revprimer),
                   ndigSites = sequence %>% purrr::map_int(~countDigSites(DNAString(.x), enzyme1, enzyme2, enzyme3))) %>%
            tidyr::nest_legacy(aberrant_pattern:constrseq_fixed,
                        .key = 'site_fix_info')
        }
      } else {
        failureRes = data_frame(ID = snp$ID,
                                CHROM = snp$CHROM,
                                POS = snp$POS,
                                REF = snp$REF,
                                ALT = snp$ALT,
                                result = 'Failed - SNP sequence could not be generated without an aberrant digestion site in all constructs')
        return(failureRes)
      }
    }

    if (any(res$ndigSites > 3)) {
      #divide the results into broken and working sequences
      working = res %>% filter(ndigSites < 3)
      broken = res %>% filter(ndigSites >= 3)

      #remove the barcodes for the sequences that work
      brokenPool = unlist(snp$bcPools)[!(unlist(snp$bcPools) %in% working$barcodes)]

      resample_attempts = 0
      while (any(res$ndigSites > 3)) {
        resample_attempts = resample_attempts + 1

        if (resample_attempts > 40) {
          failureRes = data_frame(ID = snp$ID,
                                  CHROM = snp$CHROM,
                                  POS = snp$POS,
                                  REF = snp$REF,
                                  ALT = snp$ALT,
                                  result = 'Failed - SNP sequence could not be generated without an aberrant digestion site')
          return(failureRes)
        }

        #For the subset of sequences that don't work, resample the barcodes and try again.
        if (extra_elements) {
          fixed = broken %>% mutate(barcodes = sample(brokenPool,
                                                      nrow(broken)),
                                    sequence = paste0(fwprimer,
                                                      'TG',
                                                      constrseq,
                                                      enzyme1,
                                                      enzyme2,
                                                      barcodes,
                                                      'GGC',
                                                      revprimer),
                                    ndigSites = sequence %>% purrr::map_int(~countDigSites(DNAString(.x), enzyme1, enzyme2, enzyme3)))
        } else {
          fixed = broken %>% mutate(barcodes = sample(brokenPool,
                                                      nrow(broken)),
                                    sequence = paste0(fwprimer,
                                                      constrseq,
                                                      enzyme1,
                                                      enzyme2,
                                                      barcodes,
                                                      revprimer),
                                    ndigSites = sequence %>% purrr::map_int(~countDigSites(DNAString(.x), enzyme1, enzyme2, enzyme3)))
        }
        res = rbind(working, fixed)
      }
    }
  } else if (isDEL) {

    genome = BSgenome.Hsapiens.UCSC.hg38
    refwidth = nchar(snp$REF)

    rangestart = snp$POS - upstreamContextRange
    rangeend = snp$POS + downstreamContextRange

    if (snp$reverseGene) {
      rangestart = snp$POS - downstreamContextRange
      rangeend = snp$POS + upstreamContextRange
    }

    snpseq = subseq(genome[[paste0('chr', as.character(snp$CHROM))]], # the chrom field needs to be only digits
                    start = rangestart,
                    end = rangeend)
    if (snp$reverseGene) {
      snpseq %<>% reverseComplement()
    }

    ndigsite = countDigSites(snpseq, enzyme1, enzyme2, enzyme3)
    ndigsite_in_context = ndigsite

    if (ndigsite > 0) {
      if (alter_aberrant) {
        dig_patterns = c(enzyme1,
                         enzyme2,
                         enzyme3,
                         enzyme1 %>% reverse,
                         enzyme2 %>% reverse,
                         enzyme3 %>% reverse)

        dig_site_locations = purrr::map(dig_patterns, Biostrings::matchPattern, subject = snpseq, fixed = FALSE)
      } else {
        failureRes = data_frame(ID = snp$ID,
                                CHROM = snp$CHROM,
                                POS = snp$POS,
                                REF = snp$REF,
                                ALT = snp$ALT,
                                result = 'Failed - Context contained a digestion site')
        return(failureRes)
      }
    }

    delUpstreamRange = ifelse(snp$reverseGene,
                              downstreamContextRange,
                              upstreamContextRange)

    altseq = generateDelConstruct(snpseq, refwidth, delUpstreamRange)

    res = data_frame(ID = snp$ID,
                     CHROM = snp$CHROM,
                     snpIndex = 1:(nper*2),
                     type = rep(c('ref', 'alt'), each = nper),
                     mid = ifelse(type == 'ref', snp$REF, snp$ALT),
                     barcodes = sample(snp$bcPools %>% unlist,
                                       2*nper),
                     constrseq = purrr::map(type, ~if (.x == 'ref') {snpseq} else {altseq}))

    res %<>% mutate(constrseq = constrseq %>% purrr::map_chr(~toString(DNAString(.x))))

    if (extra_elements) {
      res %<>% mutate(sequence = paste0(fwprimer,
                                        'TG',
                                        constrseq,
                                        enzyme1,
                                        enzyme2,
                                        barcodes,
                                        'GGC',
                                        revprimer),
                      ndigSites = sequence %>% purrr::map_int(~countDigSites(DNAString(.x), enzyme1, enzyme2, enzyme3)))
    } else {
      res %<>% mutate(sequence = paste0(fwprimer,
                                        constrseq,
                                        enzyme1,
                                        enzyme2,
                                        barcodes,
                                        revprimer),
                      ndigSites = sequence %>% purrr::map_int(~countDigSites(DNAString(.x), enzyme1, enzyme2, enzyme3)))
    }

    #If all of the sequences contained > 3 digestion sites, there's probably
    #some location at the context/other parts boundary that generates a site.
    #This is too complicated to fix automatically, so just fail the SNP
    if (all(res$ndigSites >= 3)) {
      if (alter_aberrant & ndigsite_in_context > 0) {
        dig_site_locations = purrr::map(dig_patterns, Biostrings::matchPattern, subject = snpseq, fixed = FALSE)

        multiple_aberrant_dig_sites = sum(purrr::map_lgl(dig_site_locations,
                                                         ~length(BiocGenerics::width(.x)) != 0)) > 1

        if (multiple_aberrant_dig_sites) {
          failureRes = data_frame(ID = snp$ID,
                                  CHROM = snp$CHROM,
                                  POS = snp$POS,
                                  REF = snp$REF,
                                  ALT = snp$ALT,
                                  result = 'Failed - SNP sequence could not be generated without multiple aberrant digestion sites')
          return(failureRes)
        }

        if (extra_elements) {
          res = randomly_fix(snp,
                             res,
                             dig_patterns,
                             dig_site_locations) %>%
            mutate(sequence = paste0(fwprimer,
                                     'TG',
                                     constrseq_fixed,
                                     enzyme1,
                                     enzyme2,
                                     barcodes,
                                     'GGC',
                                     revprimer),
                   ndigSites = sequence %>% purrr::map_int(~countDigSites(DNAString(.x), enzyme1, enzyme2, enzyme3))) %>%
            tidyr::nest_legacy(aberrant_pattern:constrseq_fixed,
                               .key = 'site_fix_info')
        } else {
          res = randomly_fix(snp,
                             res,
                             dig_patterns,
                             dig_site_locations) %>%
            mutate(sequence = paste0(fwprimer,
                                     constrseq_fixed,
                                     enzyme1,
                                     enzyme2,
                                     barcodes,
                                     revprimer),
                   ndigSites = sequence %>% purrr::map_int(~countDigSites(DNAString(.x), enzyme1, enzyme2, enzyme3))) %>%
            tidyr::nest_legacy(aberrant_pattern:constrseq_fixed,
                               .key = 'site_fix_info')
        }
      } else {
        failureRes = data_frame(ID = snp$ID,
                                CHROM = snp$CHROM,
                                POS = snp$POS,
                                REF = snp$REF,
                                ALT = snp$ALT,
                                result = 'Failed - SNP sequence could not be generated without an aberrant digestion site in all constructs')
        return(failureRes)
      }
    }

    if (any(res$ndigSites > 3)) {
      #divide the results into broken and working sequences
      working = res %>% filter(ndigSites < 3)
      broken = res %>% filter(ndigSites >= 3)

      #remove the barcodes for the sequences that work
      brokenPool = unlist(snp$bcPools)[!(unlist(snp$bcPools) %in% working$barcodes)]

      resample_attempts = 0
      while (any(res$ndigSites > 3)) {
        resample_attempts = resample_attempts + 1

        if(resample_attempts > 40) {
          failureRes = data_frame(ID = snp$ID,
                                  CHROM = snp$CHROM,
                                  POS = snp$POS,
                                  REF = snp$REF,
                                  ALT = snp$ALT,
                                  result = 'Failed - SNP sequence could not be generated without an aberrant digestion site')
          return(failureRes)
        }

        if (extra_elements) {
          fixed = broken %>% mutate(barcodes = sample(brokenPool,
                                                      nrow(broken)),
                                    sequence = paste0(fwprimer,
                                                      'TG',
                                                      constrseq,
                                                      enzyme1,
                                                      enzyme2,
                                                      barcodes,
                                                      'GGC',
                                                      revprimer),
                                    ndigSites = sequence %>% purrr::map_int(~countDigSites(DNAString(.x), enzyme1, enzyme2, enzyme3)))
        } else {
          fixed = broken %>% mutate(barcodes = sample(brokenPool,
                                                      nrow(broken)),
                                    sequence = paste0(fwprimer,
                                                      constrseq,
                                                      enzyme1,
                                                      enzyme2,
                                                      barcodes,
                                                      revprimer),
                                    ndigSites = sequence %>% purrr::map_int(~countDigSites(DNAString(.x), enzyme1, enzyme2, enzyme3)))
        }
        res = rbind(working, fixed)
      }
    }

  } else {
    failureRes = data_frame(ID = snp$ID,
                            CHROM = snp$CHROM,
                            POS = snp$POS,
                            REF = snp$REF,
                            ALT = snp$ALT,
                            result = 'Failed - Not identifiable as SNV, insertion, or deletion.')
    return(failureRes)
  }

  # Do some checks that all the barcodes are unique and otherwise return the result.
  if (length(unique(res$barcodes)) != nrow(res)) { # This should never happen.
    stop('Sequence generation finished but barcodes are nonunique')
  }

  if (!('site_fix_info' %in% names(res))) {
    res$site_fix_info = NA
  }

  if (!is.null(notes)) {
    res$notes = notes
  } else {
    res$notes = NA
  }

  return(res)
}



#' Process VCF into MPRA sequences
#'
#' \code{processVCF} takes a VCF of SNPs (preferably from dbSNP) and turns them
#' into a set of labeled MPRA sequences barcoded with inert n-mers
#' @param vcf the path to the input VCF
#' @param nper the number of barcoded sequences to be generated per allele per
#'   SNP
#' @param upstreamContextRange the amount of sequence context to acquire
#'   upstream of the SNP
#' @param downstreamContextRange the amount of sequence context to acquire
#'   downstream of the SNP
#' @param fwprimer a string containing the forward PCR primer to be used
#' @param revprimer a string containing the reverse PCR primer to be used
#' @param enzyme1 a string containing the pattern for the first restriction
#'   enzyme. Defaults to KpnI.
#' @param enzyme2 a string containing the pattern for the second restriction
#'   enzyme. Defaults to XbaI.
#' @param enzyme3 a string containing the pattern for the third restriction
#'   enzyme. Defaults to SfiI.
#' @param filterPatterns a character vector of patterns to filter out of the
#'   barcode pool (along with their reverse complements)
#' @param outPath an optional path stating where to write a .tsv of the results
#' @param alter_aberrant under development - logical indicating whether to
#'   randomly alter aberrant digestion sites across barcodes
#' @param extra_elements under development - logical indicating whether to
#'   include the extra TG / GGC as shown on the sequence diagram on the shiny
#'   app
#' @param max_construct_size under development - integer indicating the maximum
#'   construct size to generate. If provided, constructs that end up longer than
#'   this have sequence context evenly removed from both sides until
#'   sufficiently short.
#' @param barcode_set string - indicating the barcode set to use. Alternatively
#'   a vector containing custom barcodes. See below for details.
#' @param ensure_all_4_nuc logical -- if true, barcodes are filtered to only
#'   those containaing all four nucleotides.
#' @param flip_RV logical - if true, take the reverse complement of any alleles
#'   with "RV" in the INFO field. This is to account for \href{https://www.ncbi.nlm.nih.gov/variation/docs/oldglossary_dbSNP1_vcf/}{SNPs that are encoded in terms of the reverse strand alleles in dbSNP}.
#' @details The \code{"filterPatterns"} argument is used to remove barcodes
#'   containing patterns that may perform badly in a MPRA setting. For example,
#'   the default, 'AATAAA', corresponds to a sequence required for cleavage and
#'   polyadenylation of pre-mRNAs in eukaryotic cells.
#'
#'   The \code{upstreamContextRange} and \code{downstreamContextRange} arguments
#'   are handled intuitively for minus strand SNPs (i.e. those that have the
#'   MPRAREV tag). So for a minus strand SNP you get the complement of
#'   \code{downstreamContextRange} - SNP - \code{upstreamContextRange} as the
#'   genomic context.
#'
#'   The three \code{enzyme} arguments may contain ambiguous nucleotides by
#'   including an N character at the appropriate base (for example the 5 N's in
#'   the SfiI default).
#'
#'   The sequence for \code{enzyme3} does not show up in the output sequences,
#'   however it is necessary to check for it's presence in the output sequences
#'   as it is used when preparing the plasmid library. Aberrant \code{enzyme3}
#'   sites could cause the library preparation to fail.
#'
#'   Alternative barcode sets may be used by specifying the \code{barcode_set}
#'   argument to \code{processVCF} one of the following values. The first number
#'   indicates the length of the barcodes in basepairs, the second indicates the
#'   number of errors correctable while still being able to identify the
#'   original barcode. These are provided by the freebarcodes package, detailed
#'   at the publication below and available from the subsequently listed github
#'   repository. The original barcode set provided with mpradesigntools is
#'   available as the \code{twelvemers} barcode set. See the README on github
#'   for a listing of the number of barcodes available per set. The freebarcodes
#'   sets only meet the traditional MPRA barcode requirements to varying degree.
#'   \itemize{ \item contains all four nucleotides \item doesn't contain runs of
#'   4 or more of the same nucleotide \item doesn't contain miR seed sequences }
#'   Alternatively, \code{barcode_set} can be a character vector containing a
#'   custom set of all barcodes you'd like to use.
#' @return A list of two data_frames. The first, named 'result', is a data_frame
#'   containing the labeled MPRA sequences. The second, named 'failed', is a
#'   data_frame listing the SNPs that are not able to have MPRA sequences
#'   generated and the reason why.
#' @export
#' @importFrom dplyr select
#' @importFrom dplyr %>%
#' @importFrom dplyr filter
#' @importFrom dplyr rename
#' @importFrom dplyr mutate
#' @importFrom dplyr rowwise
#' @importFrom dplyr do
#' @importFrom dplyr as.tbl
#' @importFrom dplyr ungroup
#' @importFrom magrittr %<>%
#' @importFrom purrr map
#' @importFrom purrr map_int
#' @importFrom purrr map_chr
#' @importFrom purrr map2
#' @importFrom purrr map2_chr
#' @importFrom purrrlyr by_row
#' @importFrom readr read_tsv
#' @importFrom readr write_tsv
#' @importFrom stringr str_split
#' @importFrom stringr str_locate
#' @importFrom tidyr unnest_legacy
#' @importFrom tibble rownames_to_column
#' @importFrom Biostrings DNAStringSet
#' @importFrom Biostrings reverseComplement
#' @importFrom Biostrings toString
processVCF = function(vcf,
                      nper,
                      upstreamContextRange,
                      downstreamContextRange,
                      fwprimer,
                      revprimer,
                      enzyme1 = 'GGTACC',
                      enzyme2 = 'TCTAGA',
                      enzyme3 = 'GGCCNNNNNGGCC',
                      filterPatterns = 'AATAAA',
                      alter_aberrant = FALSE,
                      extra_elements = FALSE,
                      max_construct_size = NULL,
                      barcode_set = 'twelvemers',
                      ensure_all_4_nuc = TRUE,
                      flip_RV = TRUE,
                      outPath = NULL){

  # kpn = 'GGTACC' #KpnI
  # xba = 'TCTAGA' #XbaI
  # sfi = 'GGCCNNNNNGGCC' #SfiI

  #skip metadata lines
  vcf_lines = readLines(con = vcf)
  skipNum = length(grep(x = vcf_lines, pattern = '^##'))

  # Check that the header doesn't have spaces in place of tabs. If it does (why
  # dbSNP, why?), replace the spaces with tabs and create a new col_names
  # variable
  vcfColumns = vcf_lines[skipNum + 1] %>%
    gsub('#', '', .) %>%
    gsub('[ ]+', '\t', .) %>% #replace spaces with tabs if applicable
    str_split('\t') %>%
    unlist

  vcf = readr::read_tsv(vcf,
                 skip = skipNum + 1,
                 col_names = vcfColumns,
                 col_types = readr::cols(.default = readr::col_character(),
                                         POS = readr::col_integer()))

  #select = dplyr::select

  vcf %<>%
    purrrlyr::by_row(spreadAllelesAcrossRows) %>%
    .$.out %>%
    Reduce('rbind', .)


  #mers = twelvemers
  # list.files('data') %>% gsub('\\.RData|\\.rda', '', .) %>% dput()
  available_sets = c("barcodes10-1", "barcodes10-2", "barcodes11-1", "barcodes11-2",
                     "barcodes12-1", "barcodes12-2", "barcodes13-1", "barcodes13-2",
                     "barcodes14-1", "barcodes14-2", "barcodes15-1", "barcodes15-2",
                     "barcodes16-1", "barcodes16-2", "barcodes17-2", "barcodes3-1",
                     "barcodes4-1", "barcodes5-1", "barcodes5-2", "barcodes6-1", "barcodes6-2",
                     "barcodes7-1", "barcodes7-2", "barcodes8-1", "barcodes8-2", "barcodes9-1",
                     "barcodes9-2", "twelvemers")

  if (length(barcode_set) == 1 & barcode_set %in% available_sets){
    mers = get(barcode_set)
  } else if (length(barcode_set) > 1 & class(barcode_set) == 'character') {
    print('Note: Using custom barcode set. Users recommended to ensure barcode set meets appropriate barcode parameters for length, inert properties, etc....')
    mers = barcode_set
  } else {
    stop('Input barcode set improperly formatted. Use the names of one of the available sets or a character vector of same-length DNA n-mers.')
  }

  if (nrow(vcf)*2*nper > length(mers)) {
    if (barcode_set == 'barcodes16-1') {
      stop('Your design requires more barcodes than is possible with the largest available barcode_set.\n\nPoke the developer about integrating the freebarcodes barcode-concatenation trick.')
    } else {
      stop('Your design requires more barcodes than is possible with the selected barcode_set. Try a bigger set.')
    }
  }

  if (ensure_all_4_nuc) {
    print('Filtering barcode set to ensure that all barcodes contain all four barcodes...')
    start_amount = length(mers)
    mers = mers[purrr::map_lgl(mers, matches_all_nucleotides)]
    end_amount = length(mers)

    print(paste0('Removed ', start_amount - end_amount, ' barcodes out of ', start_amount, ' (', round((start_amount - end_amount)/start_amount * 100, digits = 2), '%)'))
  }


  filterRegex = paste(c(filterPatterns, # the patterns
                        filterPatterns %>% DNAStringSet %>% reverseComplement() %>% toString %>% str_split(', ') %>% unlist), # and their reverse complements
                      collapse = '|')

  print('Filtering undesired barcode patterns...')
  undesired_start_n = length(mers)
  barcodeFilter = mers %>%
    str_locate(filterRegex) %>%
    as.data.frame() %>%
    as.tbl() %>%
    tibble::rownames_to_column('removeIndex') %>%
    mutate(removeIndex = as.integer(removeIndex)) %>%
    na.omit

  if (nrow(barcodeFilter) > 1) {
    mers %<>% .[-barcodeFilter$removeIndex]
  }
  undesired_end_n = length(mers)

  print(paste0('Removed ',
               undesired_start_n - undesired_end_n,
               ' barcodes from the usable pool out of the original ',
               undesired_start_n,
               ' (', round(100*(undesired_start_n - undesired_end_n)/undesired_start_n, digits = 3), '%)'))


  #Create a pool of barcodes for each snp
  shuffled_mers = mers[sort(runif(length(mers)), index.return = TRUE)$ix]

  vcf %<>% mutate(bcPools = split(shuffled_mers,
                                  ceiling(1:length(shuffled_mers) %% nrow(vcf))),
                  reverseGene = grepl('MPRAREV', INFO),
                  snpNums = 1:nrow(vcf),
                  snpTot = nrow(vcf))

  print('Processing SNPs...')

  processed = vcf %>%
    rowwise %>%
    do(seqs = processSnp(.,
                         nper = nper,
                         upstreamContextRange = upstreamContextRange,
                         downstreamContextRange = downstreamContextRange,
                         fwprimer,
                         revprimer,
                         enzyme1,
                         enzyme2,
                         enzyme3,
                         alter_aberrant = alter_aberrant,
                         extra_elements = extra_elements,
                         max_construct_size = max_construct_size,
                         flip_RV = flip_RV)) %>%
    mutate(dataNames = names(seqs) %>% list,
           failed = any(grepl('result', dataNames)))

  if (!any(processed$failed)) { # if none failed

    successes = processed %>%
      filter(!failed) %>%
      .$seqs %>%
      Reduce('rbind', .) %>%
      mutate(.,
             constrseq = constrseq %>% unlist, # not sure how this got turned into a list
             totIndex = 1:nrow(.)) %>%
      rename(allele = mid,
             barcode = barcodes) %>%
      select(ID, type, allele, snpIndex, totIndex, barcode, sequence, site_fix_info, notes)

    res = list(result = successes, failed = NA)

    if (!is.null(outPath)) {
      outPath %<>% gsub('.tsv', '', .) %>% paste0(., '.tsv')
      write_tsv(successes, path = outPath)
    }

  } else if (all(processed$failed)) { # if they ALL failed :(

    failures = processed %>%
      filter(failed) %>%
      select(seqs) %>%
      ungroup %>%
      tidyr::unnest_legacy() %>%
      dplyr::rename(reason = result)

    res = list(result = NA, failed = failures)

  } else {

    failures = processed %>%
      filter(failed) %>%
      select(seqs) %>%
      ungroup %>%
      tidyr::unnest_legacy() %>%
      dplyr::rename(reason = result)

    successes = processed %>%
      filter(!failed) %>%
      .$seqs %>%
      Reduce('rbind', .) %>%
      mutate(.,
             constrseq = constrseq %>% unlist, # not sure how this got turned into a list
             totIndex = 1:nrow(.)) %>%
      dplyr::rename(allele = mid,
             barcode = barcodes) %>%
      select(ID, type, allele, snpIndex, totIndex, barcode, sequence, site_fix_info, notes)

    res = list(result = successes, failed = failures)

    if (!is.null(outPath)) {
      outPath %<>% gsub('\\.tsv', '', .) %>% paste0(., '.tsv')
      save(res,
           file = outPath %>% gsub('\\.tsv', '\\.RData', .))

      if (alter_aberrant) {

        output = res$result

        fixed_sites = res$result %>%
          dplyr::filter(!purrr::map_lgl(site_fix_info,
                                        ~all(class(.x) == 'logical')))

        if (nrow(fixed_sites) > 0) {
          tmp = fixed_sites %>%
            dplyr::pull(site_fix_info) %>%
            .[[1]]

          tmp[2,] = NA
          empty_fix = tmp[2,]

          output$site_fix_info = purrr::map(output$site_fix_info,
                                            fix_site_fix_info,
                                            empty = empty_fix)

          output %>%
            tidyr::unnest_legacy() %>%
            write_tsv(path = outPath)
        } else {
          res$result %>%
            write_tsv(path = outPath)
        }


      } else {
        write_tsv(successes, path = outPath)
      }
    }

  }
  print('Output construct size distribution:')
  res$result %>%
    dplyr::mutate(n_bp = nchar(sequence)) %>%
    dplyr::count(n_bp) %>%
    {print(.,
           n = nrow(.))}

  return(res)
}

fix_site_fix_info = function(.x, empty){
  if (any(class(.x) == 'data.frame')){
    return(.x)
  } else {
    return(empty)
  }
}

matches_all_nucleotides = function(bc){
  # https://stackoverflow.com/questions/469913/regular-expressions-is-there-an-and-operator
  grepl('(?=.*A)(?=.*C)(?=.*G)(?=.*T)', bc,
        perl = TRUE)
}

fix_indels = function(REF, ALT){
  # fixes indels to put "-" in the appropriate allele
  #returns REF_ALT
  # to be called like mutate(fixed_allele = map2(REF, ALT, fix_indels)) %>% separate(col = fixed_allele, into = c('REF', 'ALT'))
  if (nchar(REF) == 1 & nchar(ALT) == 1) {
    return(paste0(REF, '_', ALT))
  }

  if (nchar(REF) == 1 & nchar(ALT) > 1) {
    # it's an insertion

    if (REF == substr(ALT, 1,1)) {

      fixed_ins = paste0('-_', substr(ALT, 2, nchar(ALT)))
      return(fixed_ins)
    } else {
      stop('detected insertion but initial bases do not match')
    }
  }

  if (nchar(REF) > 1 & nchar(ALT) == 1) {
    # it's a deletion

    if (substr(REF, 1,1) == ALT) {

      fixed_del = paste0(substr(REF, 2, nchar(REF)), '_-')
      return(fixed_del)
    } else {
      return('detected deletion but initial bases do not match')
    }
  }

  else {
    stop('cannot repair two multibase alleles')
  }
}

#' Prepare input VCF
#'
#' Prepare an input VCF for sequence generation
#'
#' @details This function takes an input vcf, reads it in, spreads any SNPs with
#'   alternate alleles across multiple rows, and fixes any indels that are
#'   improperly formatted. This means if REF and ALT are listed as "A" and "ATC"
#'   they will be replaced with "-" and "TC". If they're "A" and "T,C" this will
#'   be spread into two entries of "A" and "T" and "A" and "C".
#'
#'   The output is written to the same directory as the input named
#'   "*_fixed.vcf". This is ready to be fed into processVCF()
#' @param vcf_path path to a vcf to be fixed
#' @return a data_frame of the fixed VCF
#' @export
#' @importFrom stringr str_split
#' @importFrom readr read_tsv
#' @importFrom purrrlyr by_row
#' @importFrom dplyr pull
#' @importFrom dplyr bind_rows
#' @importFrom tidyr separate
#' @importFrom magrittr %>%
#' @importFrom magrittr %<>%
spread_and_fix_indels = function(vcf_path){

  vcf_lines = readLines(con = vcf_path)
  skipNum = length(grep(x = vcf_lines, pattern = '^##'))

  # Check that the header doesn't have spaces in place of tabs. If it does (why
  # dbSNP, why?), replace the spaces with tabs and create a new col_names
  # variable
  vcfColumns = vcf_lines[skipNum + 1] %>%
    gsub('#', '', .) %>%
    gsub('[ ]+', '\t', .) %>% #replace spaces with tabs if applicable
    stringr::str_split('\t') %>%
    unlist

  vcf = readr::read_tsv(vcf_path,
                 skip = skipNum + 1,
                 col_names = vcfColumns,
                 col_types = readr::cols(.default = readr::col_character(),
                                         POS = readr::col_integer()))

  vcf %<>%
    purrrlyr::by_row(spreadAllelesAcrossRows) %>%
    dplyr::pull(.out) %>%
    dplyr::bind_rows() %>%
    dplyr::mutate(fixed_allele = purrr::map2_chr(REF, ALT, fix_indels)) %>%
    tidyr::separate(col = fixed_allele, into = c('REF', 'ALT'), sep = '_')

  names(vcf)[names(vcf) == 'CHROM'] = '#CHROM'
  vcf %>%
    write_tsv(gsub('.vcf', '_fixed.vcf', vcf_path))
  vcf
}

