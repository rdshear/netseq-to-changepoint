# TODO - convert from package to script
#' HRseqDataFromBAM
#' 
#' Create a list of HRseqData objects from a list of BAM files
#' 
#' TODO - more description
#' 
#' @import dplyr
#' @import purrr
#' @import Rsamtools
#' @import GenomicAlignments
#' @import GenomicRanges
#' @importMethodsFrom rtracklayer import
#' @importClassesFrom GenomeInfoDb Seqinfo
#' @exportClass HRseqData
#' @exportMethod HRseqDataFromBAM
#' 

#' @param sampleId String
#' @param bam_file Location of BAM file
#' @param gene_list Filename of gff3 file containing regions of interest
#' @param seqinfo a SeqInfo object for the appropriate genome
#' @return A HRseqData object
#' 
setMethod("HRseqDataFromBAM", signature = c("character"), 
  function(sampleId, bam_file, gene_list, seqinfo) {
    if (!isa(seqinfo, "Seqinfo") || length(seqinfo) < 1) {
      seqinfo = seqinfo(gene_list)
    }
    if (!isa(seqinfo, "Seqinfo") || length(seqinfo) < 1) {
      stop("valid Seqinfo object required")
    }

    # This is a 'rough' filter. Gets rid of anti-sense reads and reads far away from target gene
    # Just for efficiency. There is no provision for strand-aware filtering in readGAlignments
    bam_read_mask <- gene_list
    strand(bam_read_mask) <- "*"
    bam_read_mask <- GenomicRanges::reduce(bam_read_mask)
    
    # The reads are from cDNA, therefore the 5'-end is the 3'-end of the nascent RNA
    # which is the last base exposed from the elongation complex.
    # Therefore we will declare the occupancy to be the 5'-end. And we will
    # swap the strand information when we actually du=o the counts.
    dplyr::tibble(sample_id = sampleId, bam_file = bam_file)  %>%
      dplyr::mutate(bamreads = map(bam_file, function(u) {
        u <- GenomicRanges::GRanges(GenomicAlignments::readGAlignments(u, 
               param = Rsamtools::ScanBamParam(tag = c("NH", "HI"),
                                    what = c("qname", "cigar", "qwidth"), 
                                    which = bam_read_mask)), 
             seqinfo = seqinfo)
        glst <- gene_list
        strand(glst) <- .complementStrand(glst)
        u <- .OverlappedRanges(glst, u)
        u$cigstart <- explodeCigarOps(u$cigar) %>% 
          map(paste0, collapse="") %>% 
          unlist %>%
          stringi::stri_sub(., 
            if_else(as.character(strand(u)) == "+", 1, -1), length = 1)
        u <- u[u$cigstart == "M"]
        strand(u) <- .complementStrand(u)
        u <- split(u, u$HI)
        tibble(n_multi = as.integer(names(u)), gmask = as.list(u))
      })) %>%
      mutate(gsignal = map(bamreads, function(u) {
        GenomicRanges::resize(u$gmask[[1]], width = 1, fix = "end")  %>%
          sort
      })) %>%
      mutate(gsignal = map(gsignal, function(u)
        SamToScore(.OverlappedRanges(gene_list, u)))) %>% 
      mutate(gmask = map(bamreads, function(masktab) {
        cum <- GRanges(seqinfo = seqinfo(masktab$gmask[[1]]))
        result <- list()
        masktab <- masktab[-1, ]
        for (i in masktab$gmask) {
          cum <- GenomicRanges::reduce(.OverlappedRanges(gene_list, c(cum, i)))
          seqinfo(cum) <- seqinfo
          result <- append(result, list(cum))
        }
        masktab$gmask <- result
        masktab
      })) %>%
      mutate(topmask = map(gmask, function(u) 
        {
          w <- which.max(u$n_multi)
          if (!is_empty(w)) {
            u$gmask[which.max(u$n_multi)][[1]]
          } else
          {
            GRanges()
          }
        })) %>% #TODO shrink this down!
      mutate(g_scores = map2(gsignal, topmask, function(s, m){
        masked_scores <- .OverlappedRanges(gaps(m) %>%
                                            .[as.character(strand(.)) != "*"], s)
        ov <- findOverlaps(gene_list, masked_scores)
        result <- split(masked_scores[subjectHits(ov)], names(gene_list)[queryHits(ov)])
      })) %>% 
      mutate(r = pmap(list(gsignal, sample_id, topmask),
                         function(gsignal, sample_id, topmask) {
        HRseqData(scores = gsignal[gsignal$score != 0],
                   sampleId = sample_id,
                   seqinfo = seqinfo,
                   mask = topmask)
      })) -> nsdl
    nsdl$r
})
