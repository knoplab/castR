# FUNCTIONS FOR GENOME PROCESSING ==============================================

#' Select Genomic Features from a  GFF File for CASTING
#'
#' Creates a CASTING lookup table from a genome assembly and its corresponding 
#' genomic feature annotation file. This serves as a starting point for further 
#' processing.
#' @details The CASTING lookup table is an expanded form of a GFF table with the
#' last \code{range_width} nucleotides upstream and the first \code{range_width}
#' nucleotides downstream the inteded site of integration/feature manipulation,
#' referred to as landing sequences \code{ls_5_full} and \code{ls_3_full}. The
#' \code{feature_site} and \code{offset} are kept as (mandatory) attributes of
#' this \code{data.table}. Note that the GFF “attributes” column is renamed to 
#' “\code{attribute}” in order to avoid potential confustion with base 
#' \code{attributes()} function.
#' In addition, the GFF feature attributes “locus_tag” and “gene” are kept as
#' human-readable identifiers in column if present.
#' @param genome_assembly A \code{DNAStringSet} of the genome assembly to be 
#' processed. It must be structured by the same landmarks used in the coordinate 
#' system of the genome annotation. 
#' @param genome_annotation A \code{data.frame} or \code{data.table} with nine 
#' columns that match the General Feature Format (GFF3) specification.
#' @param landmark_pattern A regular expression that can be used to extract the
#' landmark identifier from the \code{genome_assembly}'s range name. Defaults
#' to the first word in \code{@ranges@NAMES}. This must correspond to the 
#' \code{genome_annotation}'s \code{seqid} values (first column).
#' @param feature_type A character vector with the names of all feature types 
#' to be processed. These must correspond to the \code{genome_annotation}'s 
#' \code{type} values (third column).
#' @param feature_site \code{TRUE} if the feature should be manipulated at its 
#' annotated downstream site (e. g. N-terminal tagging); \code{FALSE} if the
#' upstream site should be modified (e. g. C-terminal tagging). 
#' Defaults \code{FALSE}.
#' @param offset Integer that indicates how many nucleotides should be skipped
#' after the feature start (\code{feature_site = TRUE}) or before the feature
#' end (\code{feature_site = FALSE}) with respect to the intended site of tag
#' integration/feature maniplation. This is useful for example to tag in-frame. 
#' On default, the first (last) codon is skipped.
#' @param range_width Number of nucleotides that should be kept up- and 
#' downstream of the site of modification for further processing such as NGS
#' reference genome assembly. 
#' Defaults \code{200L}.
#' @return A \code{data.table} compatible with CASTING. All sequences are
#' returned in sense target_orientation.
#' @import data.table
#' @importFrom stringr regex str_match
#' @importFrom Biostrings DNAString reverseComplement subseq
#' @export
#' @examples
#' my_lookup_table <- find_assembly(my_assembly, my_annotation, 
#'                 feature_type = "gene")

find_features <- function(genome_assembly, genome_annotation, 
                          landmark_pattern = "([^\\s]+)",
                          feature_type = NULL, 
                          feature_site = FALSE, 
                          offset = 3L, 
                          range_width = 200L
) {
  
  # PRE-PROCESSING OF GENOME ANNOTATION --------------------------------------
  
  GFF.locus_tag <- function(description) stringr::str_match(description, 
      pattern = stringr::regex("locus_tag=([^;]+)", ignore_case = TRUE))[, 2]
  GFF.locus_name <- function(description) stringr::str_match(description, 
      pattern = stringr::regex("gene=([^;]+)", ignore_case = TRUE))[, 2]
  
  # ensure column names comply with GFF3 standard (http://gmod.org/wiki/GFF3)
  colnames(genome_annotation) <- c("seqid", "source", "type", 
     "start", "end", "score", "strand", "phase", "attribute")
  
  genome_subset <- as.data.table(genome_annotation)[type %in% feature_type]
  
  # make identifiers explicit
  genome_subset[, locus_tag  := GFF.locus_tag(attribute)]
  genome_subset[, locus_name := GFF.locus_name(attribute)]
  
  # PRE-PROCESSING OF GENOME ASSEMBLY ------------------------------------------
  # 
  # time spent on the following operations can be greatly enhanced by
  # transfroming a DNAStringSet into a data.table; although we depend
  # on Biostrings' reverseComplement(...) method for DNA conversion.
  #
  # The procedure is about three times faster than other approaches.
  
  genome_subset <- merge(genome_subset, data.table(
    seqid   = stringr::str_match(genome_assembly@ranges@NAMES,
                                 pattern = landmark_pattern)[, 2],
    genomic = as.character(genome_assembly),
    max_end = genome_assembly@ranges@width
  ))
  
  # EXTRACTION OF LANDING SEQUENCES --------------------------------------------
  
  genome_subset[, ls_5_full := character(0)] # avoid coercion warnings
  genome_subset[, ls_3_full := character(0)]
  
  if (feature_site) {
	  
	# around feature start
    
    genome_subset[strand == "+",
                  ls_5_full := mapply(Biostrings::subseq, x = genomic,
                                end   = start + offset - 1, 
                                width = pmin(range_width, start + offset - 1))]
    genome_subset[strand == "+", 
                  ls_3_full := mapply(Biostrings::subseq, x = genomic,
                                start = start + offset, 
                                width = pmin(range_width, max_end - start 
                                             - offset))]
    
    genome_subset[strand == "-",
                  ls_5_full := mapply(Biostrings::subseq, x = genomic,
                                start = end - offset + 1, 
                                width = pmin(range_width, max_end - end 
                                             + offset))]
    genome_subset[strand == "-", 
                  ls_3_full := mapply(Biostrings::subseq, x = genomic,
                                end   = end - offset, 
                                width = pmin(range_width, end - offset))]
    
  } else {
	  
  	# around feature end
    
	genome_subset[strand == "+",
	              ls_5_full := mapply(Biostrings::subseq, x = genomic,
	                            end   = end - offset, 
	                            width = pmin(range_width, end - offset))]
	genome_subset[strand == "+", 
	              ls_3_full := mapply(Biostrings::subseq, x = genomic,
	                            start = end - offset + 1, 
	                            width = pmin(range_width, max_end - end 
	                                         + offset))]

	genome_subset[strand == "-",
	              ls_5_full := mapply(Biostrings::subseq, x = genomic,
	                            start = start + offset,
	                            width = pmin(range_width, max_end - start 
	                                         - offset + 1))]
	genome_subset[strand == "-", 
	              ls_3_full := mapply(Biostrings::subseq, x = genomic,
	                            end   = start + offset - 1, 
	                            width = pmin(range_width, start + offset - 1))]
    
  }
  
  # remove helper columns
  genome_subset[, genomic := NULL]
  genome_subset[, max_end := NULL]
  
  # make reverse complement from freatures on negative strand
  genome_subset[strand == "-",
                ls_3_full := as.character(Biostrings::reverseComplement(
                  Biostrings::DNAStringSet(genome_subset[strand == "-", 
                                                         ls_3_full])))
                ]
  
  genome_subset[strand == "-",
                ls_5_full := as.character(Biostrings::reverseComplement(
                  Biostrings::DNAStringSet(genome_subset[strand == "-", 
                                                         ls_5_full])))
                ]
  
  # remember attributes
  data.table::setattr(genome_subset, "feature_site", feature_site)
  data.table::setattr(genome_subset, "offset", offset)
  
  return(genome_subset)
  
}

#' Identify Candidate CRISPR Targets for CASTING
#'
#' This function will expand a CASTING lookup table with CRISPR targets that
#' can be used with CASTING.
#' @details The CASTING lookup table will be expanded by the following columns:
#' \code{target_orientation} and \code{target_position}, which describe the position of the 
#' first nucleotide of the target before/after the PAM with respect to the feature 
#' start (\code{start}). The meaning of \code{target_orientation} is dependent on 
#' \code{target_upstream} and refers to the (+)-strand if \code{target_upstream
#' == FALSE} (\code{TRUE}) and \code{target_orientation < 0} (\code{> 0}) or the 
#' (—)-strand if \code{target_upstream == FALSE} (\code{TRUE}) and 
#' \code{target_orientation > 0} (\code{< 0}). The columns \code{target_sequence}, \code{PAM_verbose},
#' and \code{PAM_generic} contain the respective sequence information about the
#' crRNA traget. \code{target_penalty} equals \code{TRUE} if a pre-mature RNA-polymerase
#' terminator (TTTTT) is present in the crRNA sequence. The leavage distance
#' (\code{cleavage_distance}) and whether the target spans the landing
#' sequence junction or not (\code{target_removed}) are indicated.
#' On default, unsuitable targets will not be removed from the lookup table, but 
#' are annotated in such a way that the user can exclude them manually by simple
#' subset operations.
#' @param subject A CASTING lookup table (\code{data.table} object).
#' @param pattern The protospacer-adjacent motif (PAM) of the selected CRISPR
#' endonuclease. This can contain IUPAC ambiguity letters.
#' @param target_length The length of the crRNA or its genomic target respectively
#' (without counting the PAM nucleotides). Defaults 20 nt.
#' @param target_core The seed length on the crRNA in which mismatches between
#' the crRNA and the genomic target are not well tolerated. Defaults 8 nt.
#' @param target_upstream Whether the PAM precedes the genomic target or not.
#' If \code{TRUE} this is 5'-NNNNNNPAM-3', and 5'-PAMNNNNNN-3' otherwise.
#' Default is \code{FALSE}.
#' @param cleave_at The distance at which the CRISPR endonucleases cleaves the
#' genomic target beyond the PAM. Defaults 15 nt.
#' @param look_around The maximum distance between the site of modification and
#' the cleavage of the CRISPR endonuclease. Defaults 20 nt.
#' @param reference_genome A \code{DNAStringSet} to test CRISPR off-targets
#' against. This should be a complete genome assembly.
#' @param mismatches An integer vector with all number of mismatches that should
#' be allowed for an off-target outside it's seed sequence. 
#' @return A \code{data.table} compatible with CASTING. All sequences are
#' returned in sense target_orientation.
#' @import data.table
#' @importFrom stringr str_match str_count str_length
#' @importFrom Biostrings DNAString reverseComplement subseq PDict vcountPDict vmatchPattern 
#' @export
#' @examples
#' my_processed_lookup_table <- find_targets(my_lookup_table, "TTV", 
#'            reference_genome = my_assembly, mismatches = c(1, 3))

find_targets <- function(subject, pattern, 
                       look_around   = 20L,
                       target_length = 20L,   target_core =  8L, 
                       target_upstream = FALSE, cleave_at = 15L,
                       reference_genome = NULL, mismatches = integer(0)
) {
  
  # PARAMETER EXTRACTION -------------------------------------------------------
  
  if (is.null(attr(subject, "feature_site"))) stop(
    "could not identify intended site of genomic feature manipulation 
     that 'subject' was created for, attribute 'feature_site' missing"
  )
  
  # sometimes attributes get cached as character ...
  feature_site <- as.logical(attr(subject, "feature_site"))
  
  if (is.null(attr(subject, "offset"))) {
    
    warning ("could not identify offset of the annotated feature
             start/end and the 'subject' sequence information, 
             attribute 'offset' missing; treated as zero")
    
    offset <- 0
    
  } else {
    
    # sometimes attributes get cached as character ...
    offset <- as.integer(attr(subject, "offset"))
    
  }
  
  # PROTOSPACER CREATION
  
  protospacer.ele <- c(
    # contains PAM site
    pattern,
    # contains guide sequence of defined length; this will prevent
    # matches beyond the borders of the query sequence
    paste(rep("N", target_length), collapse = "")
  )
  
  protospacer.pos <- DNAString(paste(c(
    # either PAM[NNNN] (target downstream) or [NNNN]PAM (target upstream)
    protospacer.ele[( target_upstream) + 1],
    protospacer.ele[(!target_upstream) + 1]
  ), collapse = ""
  ))
  protospacer.neg <- Biostrings::reverseComplement(protospacer.pos)
  
  protospacer_length <- protospacer.pos@length
  
  # EXTRACTION OF QUERY SEQUENCE -----------------------------------------------
  
  # The sequence material provided is narrowed further to limit computation time 
  # spent on pattern matching. The frame is chosen (broadended) so that cleavage
  # must within this range when valid. As the validity of any particular target 
  # is strand-dependent, this will be checked at a later point.
  #
  # The case of cleavage outside the protospacer was considered obsolete.
  #
  # Query sequence framing is performed on data.table level to preserve infor-
  # mation on possible irregular feature ends at distal scaffold sites.
  
  # combine both landing sequences
  subject[, subject.seq := paste0(ls_5_full, ls_3_full)]
  
  # subset
  subject[, subject.sub := mapply(Biostrings::subseq, x = subject.seq,
                           start = nchar(ls_5_full) - pmin(nchar(ls_5_full), 
                                   look_around + protospacer_length) + 1, 
                           end   = nchar(ls_5_full) + pmin(nchar(ls_3_full), 
                                   look_around + protospacer_length)
                                 )]
  
  subject.query <- Biostrings::DNAStringSet(subject[, subject.sub])
  subject.query@ranges@NAMES <- subject[, locus_tag] # in case needed
  
  # PROTOSPACER MATCHING -------------------------------------------------------
  
  # Initial operations are performed based on the position of the right-most 
  # nucleotide beloning to the target sequence irrespective of this being 3' 
  # or 5' of the PAM.
  
  target.end.pos <- c(0, nchar(pattern))[target_upstream + 1]
  target.end.neg <- c(nchar(pattern), 0)[target_upstream + 1]
  
  all.PAM.pos <- lapply(Biostrings::vmatchPattern(protospacer.pos, subject.query,
                        fixed = FALSE) # false allows N to match any base
                        @ ends, function(sites) sites - target.end.pos)
  all.PAM.neg <- lapply(Biostrings::vmatchPattern(protospacer.neg, subject.query,
                        fixed = FALSE)
                        @ ends, function(sites) sites - target.end.neg)
  
  # FILTER 1: remove guides in where a PAM precedes another PAM ----------------

  f_1.arg.pos <- c(0, 1)[target_upstream + 1]
  f_1.arg.neg <- c(1, 0)[target_upstream + 1]
  
  # Since moving the diff()-map by +1 will always lead to loss of the first 
  # entry, create a pseudo-first element that can safely be lost.
  
  if (f_1.arg.pos == 1) lapply(all.PAM.pos, function(sites)
    if (length(sites) > 0) c(sites[[1]] - nchar(pattern), sites) else 0)
  if (f_1.arg.neg == 1) lapply(all.PAM.neg, function(sites)
    if (length(sites) > 0) c(sites[[1]] - nchar(pattern), sites) else 0)
  
  f_1.PAM.pos <- lapply(all.PAM.pos, function(sites) if (length(sites) > 0) c(
    sites[[1]] * f_1.arg.pos, 
    sites[which(diff(sites) > nchar(pattern)) + f_1.arg.pos], 
    sites[[length(sites)]] * f_1.arg.neg))
  f_1.PAM.neg <- lapply(all.PAM.neg, function(sites) if (length(sites) > 0) c(
    sites[[1]] * f_1.arg.neg,
    sites[which(diff(sites) > nchar(pattern)) + f_1.arg.neg],
    sites[[length(sites)]] * f_1.arg.pos))
  
  # FILTER 2: remove guides cleaved outside the region of interest -------------
  
  f_2.arg.pos <- c(-1, +1)[target_upstream + 1]
  f_2.arg.neg <- c(+1, -1)[target_upstream + 1]
  
  f_2.len.pos <- target_length * f_1.arg.neg + cleave_at * f_2.arg.pos
  f_2.len.neg <- target_length * f_1.arg.pos + cleave_at * f_2.arg.neg
  
  f_2.PAM.pos <- lapply(seq_along(subject.query), function(i)
    f_1.PAM.pos[[i]][ # subset by gene:
      f_1.PAM.pos[[i]] - f_2.len.pos > protospacer_length - 1 &
      f_1.PAM.pos[[i]] - f_2.len.pos < subject.query@ranges@width[[i]] 
      - protospacer_length + 1]
  )

  f_2.PAM.neg <- lapply(seq_along(subject.query), function(i)
    f_1.PAM.neg[[i]][ # subset by gene:
      f_1.PAM.neg[[i]] - f_2.len.neg > protospacer_length - 1 &
      f_1.PAM.neg[[i]] - f_2.len.neg < subject.query@ranges@width[[i]] 
      - protospacer_length + 1]
  )
  
  # MERGING TARGET CANDIDATES INTO SUBJECT -------------------------------------
  
  # This was formerly keyed by "locus_tag", however we can not rely on this
  # value being present in all annotations, so "attribute" is used instead.
  
  subject.result <- merge(subject, data.table(
    target_positions   =     c(f_2.PAM.pos, f_2.PAM.neg),
    target_orientation = rep(c(f_2.arg.neg, f_2.arg.pos), 
                        each = length(subject.query)), # sic!
    attribute   = subject[, attribute]
  )[lengths(target_positions) > 0][, .(target_position = unlist(target_positions)),
                            keyby = list(attribute, target_orientation)],
  by = "attribute")
  
  # EXTRACTION OF TARGET SEQUENCE ----------------------------------------------
  
  subject.result[, target_sequence:= mapply(Biostrings::subseq, x = subject.sub, 
                                  width = target_length, end = target_position)]
  
  subject.result[target_orientation == c(-1, +1)[target_upstream + 1], 
                 target_sequence:= as.character(Biostrings::reverseComplement(
                   Biostrings::DNAStringSet(
                 subject.result[target_orientation == c(-1, +1)[target_upstream + 1], 
                                target_sequence])))]

  # EXTRACTION OF PROTOSPACER ADJACENT MOTIF -----------------------------------
  
  # For small subject.results, subsets may be empty. If so, mapply gets mixed 
  # zero/non-zero length input and throws and error.
  
  if (subject.result[target_orientation < 0, .N] > 0) 
      subject.result[target_orientation < 0, PAM_verbose := mapply(Biostrings::subseq, 
                 x = subject.sub, start = target_position + 1, 
                 width = nchar(pattern))]
  if (subject.result[target_orientation > 0, .N] > 0)
      subject.result[target_orientation > 0, PAM_verbose := mapply(Biostrings::subseq, 
                 x = subject.sub, end = target_position - target_length, 
                 width = nchar(pattern))]
  
  subject.result[target_orientation == c(-1, +1)[target_upstream + 1], 
                 PAM_verbose := as.character(
                   Biostrings::reverseComplement(Biostrings::DNAStringSet(
                 subject.result[target_orientation == c(-1, +1)[target_upstream + 1], 
                                PAM_verbose])))]
   
  subject.result[, PAM_generic := pattern]
  
  subject.result[, subject.seq := NULL]
  subject.result[, subject.sub := NULL]
  
  # ANNOTATION 1: EARLY TERMINATION SIGNAL -------------------------------------
  
  # faster than stringr
  subject.result[, target_penalty := mapply(grepl, x = target_sequence, pattern = "TTTTT")]
  
  # ANNOTATION 2: GENOME-WIDE MISMATCHES ---------------------------------------
  
  if (length(mismatches) > 0 & !is.null(reference_genome)) {
  
    # PAM proximal to 'target_position'
    targets.prox <- Biostrings::DNAStringSet(subject.result[target_orientation < 0, target_sequence])
    # PAM distal from 'target_position'
    targets.dist <- Biostrings::DNAStringSet(subject.result[target_orientation > 0, target_sequence])
    
    targets.prox.fwd <- Biostrings::PDict(targets.prox, 
                              tb.start = target_length - target_core,
                              tb.width = target_core)
    targets.prox.rev <- Biostrings::PDict(Biostrings::reverseComplement(targets.prox),
                              tb.start = 1,
                              tb.width = target_core)
    
    targets.dist.fwd <- Biostrings::PDict(targets.dist, 
                              tb.start = 1,
                              tb.width = target_core)
    targets.dist.rev <- Biostrings::PDict(Biostrings::reverseComplement(targets.dist),
                              tb.start = target_length - target_core,
                              tb.width = target_core)
  
  for (i in mismatches) {
    
    matches.prox <- Biostrings::vcountPDict(targets.prox.fwd, reference_genome, 
                                collapse = 1, max.mismatch = i) +
                    Biostrings::vcountPDict(targets.prox.rev, reference_genome, 
                                collapse = 1, max.mismatch = i) - 1
    matches.dist <- Biostrings::vcountPDict(targets.dist.fwd, reference_genome, 
                                collapse = 1, max.mismatch = i) +
                    Biostrings::vcountPDict(targets.dist.rev, reference_genome, 
                                collapse = 1, max.mismatch = i) - 1
    
    subject.result[target_orientation < 0, paste0("off_targets_", i) := matches.prox]
    subject.result[target_orientation > 0, paste0("off_targets_", i) := matches.dist]
    
  }}
  
  # ANNOTATION 3: TARGET REMOVAL BY INTEGRATION --------------------------------
  
  subject.result[, target_removed :=
                   target_position - look_around - protospacer_length < target_length &
                   target_position - look_around - protospacer_length > 0
                   ]
  
  # ANNOTATION 4: CLEAVAGE DISTANCE --------------------------------------------
  
  subject.result[target_orientation < 0, cleavage_distance :=
                   target_position - cleave_at -
                   look_around - protospacer_length
                 ]
  
  subject.result[target_orientation > 0, cleavage_distance :=
                   target_position - target_length + cleave_at -
                   look_around - protospacer_length
                 ]
  
  # ANNOTATION 5: POSITION IN RELATION TO FEATURE START ------------------------
  
  subject.result[target_orientation < 0, target_position := 
                   target_position +
                   (end - start) * c(1, 0)[feature_site + 1] # F: add feature width
                 - look_around - protospacer_length 
                 + offset * c(-1, +1)[feature_site + 1]      # F: substract offset
                 ]
  
  subject.result[target_orientation > 0, target_position :=
                   target_position - target_length +
                   (end - start) * c(1, 0)[feature_site + 1] # F: add feature width
                 - look_around - protospacer_length 
                 + offset * c(-1, +1)[feature_site + 1]      # F: substract offset
                 ]
  
  subject.result[target_orientation == c(+1, -1)[feature_site + 1],
                 target_position := target_position + c(+1, -1)[feature_site + 1]
                 ]
  
  subject.result[target_position >= 0, target_position := target_position + 1]
  
  # ANNOTATION 6: TARGET AND TARGET SEED GC CONTENT ----------------------------
  
  gc_content <- function(x) stringr::str_count(x, pattern = "[GCgc]") / stringr::str_length(x)
  
  # TODO target_gc; target_gc_seed
  
  data.table::setattr(subject.result, "feature_site", feature_site)
  data.table::setattr(subject.result, "offset", offset)

  return(subject.result)
  
}
