#Function to take an input vcf and generate an output vcf with the specified up/downstream context ranges and # of barcodes per snp
# And now let's do it quickly using dplyr and purrr. processVCF

#' @importFrom stringr str_split
spreadAllelesAcrossRows = function(snp){
  #snp is a row from a vcf data_frame
  # if the 'ALT' column has commas in it, spread those out across otherwise identical rows

  if (!grepl(',', snp$ALT)) {
    return(snp)
  } else {
    altAlleles = snp$ALT %>% str_split(',', simplify = TRUE)
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

  if (reverseGene) {
    toString(c(subseq(snpseq, 1, downstreamContextRange), #The insertion goes to the right of the position given
               complement(mid),
               subseq(snpseq, downstreamContextRange + 1, length(snpseq))))
  } else {
    toString(c(subseq(snpseq, 1, upstreamContextRange + 1), #The insertion goes to the right of the position given
               mid,
               subseq(snpseq, upstreamContextRange + 2, length(snpseq))))
  }
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
#' @return a data_frame of labeled sequences or a data_frame listing the SNP and
#'   why it failed
#' @import BSgenome.Hsapiens.UCSC.hg38
#' @importFrom Biostrings reverseComplement
#' @importFrom Biostrings replaceLetterAt
#' @importFrom tibble data_frame
processSnp = function(snp, nper, upstreamContextRange, downstreamContextRange, fwprimer, revprimer, enzyme1, enzyme2, enzyme3){
  # snp is one row from the expanded vcf, including the reverseGene column which
  # indicates whether or not to use the reverse complement genomic context. It
  # also has a dedicated pool of barcodes to select from

  genome = BSgenome.Hsapiens.UCSC.hg38
  # kpn = 'GGTACC' #KpnI
  # xba = 'TCTAGA' #XbaI
  # sfi = 'GGCCNNNNNGGCC' #SfiI

  isSNV = (snp$REF %in% c('A', 'C', 'G', 'T') && snp$ALT %in% c('A', 'C', 'G', 'T')) #look at me interchanging between SNV and SNP willy-nilly
  isINS = snp$REF == '-'
  isDEL = snp$ALT == '-'

  #There are three code blocks below for each of these cases. Use RStudio's code folding to open up only the one of interest.
  # This could probably be made to be 1/3 the length by writing a function that
  # adapts to the type of SNP when generating the sequences, but that might be
  # hard. There are subtle differences in each of the types of variants.

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

    ndigsite = countDigSites(snpseq, enzyme1, enzyme2, enzyme3)

    if (ndigsite > 0) {
      failureRes = data_frame(ID = snp$ID,
                              CHROM = snp$CHROM,
                              POS = snp$POS,
                              REF = snp$REF,
                              ALT = snp$ALT,
                              result = 'Failed - Context contained a digestion site')
      return(failureRes)
    }

    refseq = toString(snpseq)
    altseq = toString(replaceLetterAt(snpseq, upstreamContextRange + 1, snp$ALT))

    res = data_frame(ID = snp$ID,
                     CHROM = snp$CHROM,
                     snpIndex = 1:(nper*2),
                     type = rep(c('ref', 'alt'), each = nper),
                     mid = ifelse(type == 'ref', snp$REF, snp$ALT),
                     barcodes = sample(snp$bcPools %>% unlist,
                                       2*nper),
                     constrseq = type %>% map_chr(~ifelse(.x == 'ref',
                                                          refseq,
                                                          altseq)))

    if (snp$reverseGene) {
      res %<>% mutate(constrseq = constrseq %>% map_chr(~toString(reverseComplement(.x))))
    }

    res %<>% mutate(sequence = paste0(fwprimer,
                                      'TG',
                                      constrseq,
                                      enzyme1,
                                      enzyme2,
                                      barcodes,
                                      'GGC',
                                      revprimer),
                    ndigSites = sequence %>% map_int(~countDigSites(DNAString(.x), enzyme1, enzyme2, enzyme3)))

    #If all of the sequences contained > 3 digestion sites, there's probably some location at the context/other parts boundary that generates a site. This is too complicated to fix automatically, so just fail the SNP
    if (all(res$ndigSites > 3)) {
      failureRes = data_frame(ID = snp$ID,
                              CHROM = snp$CHROM,
                              POS = snp$POS,
                              REF = snp$REF,
                              ALT = snp$ALT,
                              result = 'Failed - SNP sequence could not be generated without an aberrant digestion site')
      return(failureRes)
    }

    if (any(res$ndigSites > 3)) {
      #divide the results into broken and working sequences
      working = res %>% filter(ndigSites <= 3)
      broken = res %>% filter(ndigSites > 3)

      #remove the barcodes for the sequences that work
      brokenPool = snp$bcPools[!(unlist(snp$bcPools) %in% working$barcodes)] %>% unlist

      while (any(res$ndigSites > 3)) {
        #For the subset of sequences that don't work, resample the barcodes and try again.
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
                                  ndigSites = sequence %>% map_int(~countDigSites(DNAString(.x), enzyme1, enzyme2, enzyme3)))
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

    if (ndigsite > 0) {
      failureRes = data_frame(ID = snp$ID,
                              CHROM = snp$CHROM,
                              POS = snp$POS,
                              REF = snp$REF,
                              ALT = snp$ALT,
                              result = 'Failed - Context contained a digestion site')
      return(failureRes)
    }

    altseq = generateInsConstruct(snpseq, DNAString(snp$ALT), snp$reverseGene, upstreamContextRange, downstreamContextRange)

    res = data_frame(ID = snp$ID,
                     CHROM = snp$CHROM,
                     snpIndex = 1:(nper*2),
                     type = rep(c('ref', 'alt'), each = nper),
                     mid = ifelse(type == 'ref', '', snp$ALT), #This line is line is unique to the isINS block
                     barcodes = sample(snp$bcPools %>% unlist,
                                       2*nper),
                     constrseq = map2_chr(mid, type, ~ifelse(.y == 'ref',
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
                     ndigSites = sequence %>% map_int(~countDigSites(DNAString(.x), enzyme1, enzyme2, enzyme3)))

    #If all of the sequences contained > 3 digestion sites, there's probably some location at the context/other parts boundary that generates a site. This is too complicated to fix automatically, so just fail the SNP

    if (all(res$ndigSites > 3)) {
      failureRes = data_frame(ID = snp$ID,
                              CHROM = snp$CHROM,
                              POS = snp$POS,
                              REF = snp$REF,
                              ALT = snp$ALT,
                              result = 'Failed - SNP sequence could not be generated without an aberrant digestion site')
      return(failureRes)
    }

    if (any(res$ndigSites > 3)) {
      #divide the results into broken and working sequences
      working = res %>% filter(ndigSites <= 3)
      broken = res %>% filter(ndigSites > 3)

      #remove the barcodes for the sequences that work
      brokenPool = snp$bcPools[!(unlist(snp$bcPools) %in% working$barcodes)] %>% unlist

      while (any(res$ndigSites > 3)) {
        #For the subset of sequences that don't work, resample the barcodes and try again.
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
                                  ndigSites = sequence %>% map_int(~countDigSites(DNAString(.x), enzyme1, enzyme2, enzyme3)))
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

    ndigsite = countDigSites(snpseq, enzyme1, enzyme2, enzyme3)

    if (ndigsite > 0) {
      failureRes = data_frame(ID = snp$ID,
                              CHROM = snp$CHROM,
                              POS = snp$POS,
                              REF = snp$REF,
                              ALT = snp$ALT,
                              result = 'Failed - Context contained a digestion site')
      return(failureRes)
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
                     constrseq = map(type, ~if (.x == 'ref') {snpseq} else {altseq}))

    if (snp$reverseGene) {
      res %<>% mutate(constrseq = constrseq %>% map_chr(~toString(reverseComplement(.x))))
    } else {
      res %<>% mutate(constrseq = constrseq %>% map_chr(~toString(.x)))
    }

    res %<>% mutate(sequence = paste0(fwprimer,
                                      'TG',
                                      constrseq,
                                      enzyme1,
                                      enzyme2,
                                      barcodes,
                                      'GGC',
                                      revprimer),
                    ndigSites = sequence %>% map_int(~countDigSites(DNAString(.x), enzyme1, enzyme2, enzyme3)))

    #If all of the sequences contained > 3 digestion sites, there's probably some location at the context/other parts boundary that generates a site. This is too complicated to fix automatically, so just fail the SNP

    if (all(res$ndigSites > 3)) {
      failureRes = data_frame(ID = snp$ID,
                              CHROM = snp$CHROM,
                              POS = snp$POS,
                              REF = snp$REF,
                              ALT = snp$ALT,
                              result = 'Failed - SNP sequence could not be generated without an aberrant digestion site')
      return(failureRes)
    }

    if (any(res$ndigSites > 3)) {
      #divide the results into broken and working sequences
      working = res %>% filter(ndigSites <= 3)
      broken = res %>% filter(ndigSites > 3)

      #remove the barcodes for the sequences that work
      brokenPool = snp$bcPools[!(unlist(snp$bcPools) %in% working$barcodes)] %>% unlist

      while (any(res$ndigSites > 3)) {
        #For the subset of sequences that don't work, resample the barcodes and try again.
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
                                  ndigSites = sequence %>% map_int(~countDigSites(DNAString(.x), enzyme1, enzyme2, enzyme3)))
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

  return(res)
}



#' Process VCF into MPRA sequences
#'
#' \code{processVCF} takes a VCF of SNPs (preferably from dbSNP) and turns them
#' into a set of labeled MPRA sequences barcoded with inert twelvemers
#' @param vcf the path to the input VCF
#' @param nper the number of barcoded sequences to be generated per allele per
#'   SNP
#' @param upstreamContextRange the amount of sequence context to acquire
#'   upstream of the SNP
#' @param downstreamContextRange the amount of sequence context to acquire
#'   downstream of the SNP
#' @param fwprimer a string containing the forward PCR primer to be used
#' @param revprimer a string containing the reverse PCR primer to be used
#' @param enzyme1 a string containing the pattern for the first restriction enzyme. Defaults to KpnI.
#' @param enzyme2 a string containing the pattern for the second restriction enzyme. Defaults to XbaI.
#' @param enzyme3 a string containing the pattern for the third restriction enzyme. Defaults to SfiI.
#' @param filterPatterns a character vector of patterns to filter out of the
#'   barcode pool (along with their reverse complements)
#' @param outPath an optional path stating where to write a .tsv of the results
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
#' @importFrom tidyr unnest
#' @importFrom tibble rownames_to_column
#' @importFrom Biostrings DNAStringSet
#' @importFrom Biostrings reverseComplement
#' @importFrom Biostrings toString
processVCF = function(vcf, nper, upstreamContextRange, downstreamContextRange, fwprimer, revprimer, enzyme1 = 'GGTACC', enzyme2 = 'TCTAGA', enzyme3 = 'GGCCNNNNNGGCC', filterPatterns = 'AATAAA', outPath = NULL){

  # kpn = 'GGTACC' #KpnI
  # xba = 'TCTAGA' #XbaI
  # sfi = 'GGCCNNNNNGGCC' #SfiI

  #skip metadata lines
  skipNum = system(paste0('grep ^## ', vcf, ' | wc -l'), intern = TRUE) %>% as.numeric

  # Check that the header doesn't have spaces in place of tabs. If it does (why
  # dbSNP, why?), replace the spaces with tabs and create a new col_names
  # variable
  vcfColumns = system(paste0('head -', skipNum + 1, ' ', vcf, ' | tail -1'),
                      intern = TRUE) %>%
    gsub('#', '', .) %>%
    gsub('[ ]+', '\t', .) %>% #replace spaces with tabs if applicable
    str_split('\t') %>%
    unlist

  vcf = read_tsv(vcf,
                 skip = skipNum + 1,
                 col_names = vcfColumns)

  #select = dplyr::select

  vcf %<>%
    by_row(spreadAllelesAcrossRows) %>%
    .$.out %>%
    Reduce('rbind', .)

  if (nrow(vcf)*2*nper > 1140292) {
    stop('Your design requests requires more barcodes than is possible')
  }

  #load('outputs/inertTwelveMersChar.RData')
  mers = twelvemers

  filterRegex = paste(c(filterPatterns, # the patterns
                        filterPatterns %>% DNAStringSet %>% reverseComplement() %>% toString %>% str_split(', ') %>% unlist), # and their reverse complements
                      collapse = '|')


  print('Filtering undesired barcode patterns...')
  barcodeFilter = mers %>%
    str_locate(filterRegex) %>%
    as.data.frame() %>%
    as.tbl() %>%
    rownames_to_column('removeIndex') %>%
    mutate(removeIndex = as.integer(removeIndex)) %>%
    na.omit

  mers %<>% .[-barcodeFilter$removeIndex]
  print(paste0('Removed ',
               nrow(barcodeFilter),
               ' barcodes from the usable pool out of the original ',
               length(twelvemers),
               ' (', round(100*nrow(barcodeFilter)/length(twelvemers), digits = 3), '%)'))


  #Create a pool of barcodes for each snp
  vcf %<>% mutate(bcPools = split(mers, ceiling(1:length(mers)/(length(mers) / nrow(vcf)))),
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
                         enzyme3)) %>%
    mutate(dataNames = names(seqs) %>% list,
           failed = any(grepl('result', dataNames)))

  failures = processed %>%
    filter(failed) %>%
    select(seqs) %>%
    unnest %>%
    rename(reason = result)

  successes = processed %>%
    filter(!failed) %>%
    .$seqs %>%
    Reduce('rbind', .) %>%
    mutate(.,
           constrseq = constrseq %>% unlist, # not sure how this got turned into a list
           totIndex = 1:nrow(.)) %>%
    rename(allele = mid,
           barcode = barcodes) %>%
    select(ID, type, allele, snpIndex, totIndex, barcode, sequence)

  if (!is.null(outPath)) {
    outPath %<>% gsub('.tsv', '', .) %>% paste0(., '.tsv')
    write_tsv(successes, path = outPath)
  }

  res = list(result = successes, failed = failures)
  return(res)
}
