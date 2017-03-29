#Function to take an input vcf and generate an output vcf with the specified context width and # of barcodes per snp
# And now let's do it quickly using dplyr and purrr. processVCF
library(stringr)

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

countDigSites = function(biostring) {
  #Count the number of times KpnI, XbaI, and SfiI sites occur in a given biostring
  
  kpn = DNAString('GGTACC') #KpnI
  xba = DNAString('TCTAGA') #XbaI
  sfi = DNAString('GGCCNNNNNGGCC') #SfiI
                  
  sum(countPattern(kpn, biostring),
      countPattern(xba, biostring),
      countPattern(sfi, biostring, fixed = FALSE),
      countPattern(kpn %>% reverse, biostring),
      countPattern(xba %>% reverse, biostring),
      countPattern(sfi %>% reverse, biostring, fixed = FALSE))
}

generateInsConstruct = function(snpseq, mid, reverseGene, seqwidth){
  #This function generates a mutant construct sesquence based on 
  # snpseq - the genomic context 
  # mid - the insertion allele
  # reverseGene - a logical indicating whether or not the SNP is for a gene that's reversed
  # seqwidth - the width of the context
  
  # If the insertion is for a gene that's transcribed from the reverse strand, it needs to be the COMPLEMENT of the allele to the LEFT of the position
  # Otherwise it's the mutant allele to the right of the position
  
  if (reverseGene) {
    toString(c(subseq(snpseq, 1, seqwidth), #The insertion goes to the right of the position given
               complement(mid),
               subseq(snpseq, seqwidth + 1, length(snpseq))))
  } else {
    toString(c(subseq(snpseq, 1, seqwidth + 1), #The insertion goes to the right of the position given
               mid,
               subseq(snpseq, seqwidth + 2, length(snpseq))))
  }
}

generateDelConstruct = function(snpseq, refwidth, seqwidth) {
  c(subseq(snpseq,
           1, 
           seqwidth),
    subseq(snpseq,
           seqwidth + refwidth + 1,
           length(snpseq)))
}

processSnp = function(snp, nper, seqwidth, fwprimer, revprimer, updateProgress = NULL){
  # snp is one row from the expanded vcf, including the reverseGene column which
  # indicates whether or not to use the reverse complement genomic context. It
  # also has a dedicated pool of barcodes to select from
  
  genome = BSgenome.Hsapiens.UCSC.hg38
  kpn = 'GGTACC' #KpnI
  xba = 'TCTAGA' #XbaI
  sfi = 'GGCCNNNNNGGCC' #SfiI
  
  isSNV = (snp$REF %in% c('A', 'C', 'G', 'T') && snp$ALT %in% c('A', 'C', 'G', 'T'))
  isINS = snp$REF == '-'
  isDEL = snp$ALT == '-'
  
  #There are three code blocks below for each of these cases. Use RStudio's code folding to open up only the one of interest.
  # This could probably be made to be 1/3 the length by writing a function that
  # adapts to the type of SNP when generating the sequences, but that might be
  # hard. There are subtle differences in each of the types of variants.
  
  if (isSNV) {
    
    rangestart = snp$POS - seqwidth
    rangeend = snp$POS + seqwidth
    
    snpseq = subseq(genome[[paste0('chr', as.character(snp$CHROM))]], # the chrom field needs to be only digits
                    start = rangestart, 
                    end = rangeend)
    
    ndigsite = countDigSites(snpseq)
    
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
    altseq = toString(replaceLetterAt(snpseq, seqwidth + 1, snp$ALT))
    
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
                                      kpn,
                                      xba,
                                      barcodes,
                                      'GGC',
                                      revprimer),
                    ndigSites = sequence %>% map_int(~countDigSites(DNAString(.x))))
    
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
                                                    kpn,
                                                    xba,
                                                    barcodes,
                                                    'GGC',
                                                    revprimer),
                                  ndigSites = sequence %>% map_int(~countDigSites(DNAString(.x))))
        res = rbind(working, fixed)
      }
    }
  } else if (isINS) {
    
    #If the snp is an insertion, the range of context to get is the same, but the middle allele (variable 'mid') is different
    rangestart = snp$POS - seqwidth
    rangeend = snp$POS + seqwidth
    
    snpseq = subseq(genome[[paste0('chr', as.character(snp$CHROM))]], # the chrom field needs to be only digits
                    start = rangestart, 
                    end = rangeend)
    
    if (snp$reverseGene) {
      snpseq %<>% reverseComplement()
    }
    
    ndigsite = countDigSites(snpseq)
    
    if (ndigsite > 0) {
      failureRes = data_frame(ID = snp$ID,
                              CHROM = snp$CHROM,
                              POS = snp$POS,
                              REF = snp$REF,
                              ALT = snp$ALT,
                              result = 'Failed - Context contained a digestion site')
      return(failureRes)
    }
    
    altseq = generateInsConstruct(snpseq, DNAString(snp$ALT), snp$reverseGene, seqwidth)
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
                                       kpn,
                                       xba,
                                       barcodes,
                                       'GGC',
                                       revprimer),
                     ndigSites = sequence %>% map_int(~countDigSites(DNAString(.x))))
    
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
                                                    kpn,
                                                    xba,
                                                    barcodes,
                                                    'GGC',
                                                    revprimer),
                                  ndigSites = sequence %>% map_int(~countDigSites(DNAString(.x))))
        res = rbind(working, fixed)
      }
    }
  } else if (isDEL) {
    
    genome = BSgenome.Hsapiens.UCSC.hg38
    refwidth = nchar(snp$REF)
    rangestart = snp$POS - seqwidth
    rangeend = snp$POS + seqwidth
    
    snpseq = subseq(genome[[paste0('chr', as.character(snp$CHROM))]], # the chrom field needs to be only digits
                    start = rangestart, 
                    end = rangeend)
    
    ndigsite = countDigSites(snpseq)
    
    if (ndigsite > 0) {
      failureRes = data_frame(ID = snp$ID,
                              CHROM = snp$CHROM,
                              POS = snp$POS,
                              REF = snp$REF,
                              ALT = snp$ALT,
                              result = 'Failed - Context contained a digestion site')
      return(failureRes)
    }
    
    altseq = generateDelConstruct(snpseq, refwidth, seqwidth)
    
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
                                      kpn,
                                      xba,
                                      barcodes,
                                      'GGC',
                                      revprimer),
                    ndigSites = sequence %>% map_int(~countDigSites(DNAString(.x))))
    
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
                                                    kpn,
                                                    xba,
                                                    barcodes,
                                                    'GGC',
                                                    revprimer),
                                  ndigSites = sequence %>% map_int(~countDigSites(DNAString(.x))))
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
  
  # http://shiny.rstudio.com/gallery/progress-bar-example.html
  if (is.function(updateProgress)) {
    text <- paste0(snp$ID, ' completed -- ', snp$snpNums, ' / ', snp$snpTot)
    updateProgress(value = snp$snpNums / snp$snpTot, detail = text)
  }
  
  return(res)
}

processVCF = function(vcf, nper, seqwidth, fwprimer, revprimer, updateProgress = NULL){
  library(BSgenome.Hsapiens.UCSC.hg38)
  library(magrittr)
  library(tidyverse)
  
  #expand = S4Vectors::expand
  select = dplyr::select
  
  vcf %<>% 
    by_row(spreadAllelesAcrossRows) %>% 
    .$.out %>% 
    Reduce('rbind', .)
  
  if (nrow(vcf)*2*nper > 1140292) {
    stop('Your design requests requires more barcodes than is possible')
  }
  
  load('outputs/inertTwelveMersChar.RData')
  mers = twelvemers
  
  # kpn = 'GGTACC' #KpnI
  # xba = 'TCTAGA' #XbaI
  # sfi = 'GGCCNNNNNGGCC' #SfiI
  # 
  # maxbc = nrow(vcf)*2*nper
  
  #Create a pool of barcodes for each snp
  vcf %<>% mutate(bcPools = split(mers, ceiling(1:length(mers)/(length(mers) / nrow(vcf)))),
                  reverseGene = grepl('MPRAREV', INFO),
                  snpNums = 1:nrow(vcf),
                  snpTot = nrow(vcf))
  
  print(sessionInfo())
  #writeLines(capture.output(sessionInfo()), "/mnt/labhome/andrew/designMPRA/sessionInfo.txt")
  
  processed = vcf %>% 
    rowwise %>% 
    do(seqs = processSnp(., nper = nper, seqwidth = seqwidth, fwprimer, revprimer, updateProgress)) %>% 
    mutate(dataNames = names(seqs) %>% list,
           failed = any(grepl('result', dataNames)))
  
  #writeLines(capture.output(sessionInfo()), "/mnt/labhome/andrew/designMPRA/sessionInfoEnd.txt")
  
  failures = processed %>% 
    filter(failed) %>% 
    select(seqs) %>% 
    unnest %>% 
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
    select(ID, type, allele, snpIndex, totIndex, barcode, sequence)

  res = list(result = successes, failed = failures)
  #write_tsv(successes, paste0('outputs/seqFileOutputs/', Sys.Date() %>% gsub('-', '_', .), '.tsv'))
  return(res)
}
