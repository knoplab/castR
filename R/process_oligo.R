# FUNCTIONS FOR OLIGO PROCESSING ===============================================

#' Landing Sequence Adjustment Table
#' 
#' This \code{data.frame} is used to infer the length of landing sequence homology 
#' provided in dependence of the absulte cleavage distance from the intended site
#' of integration/feature manipulation for CASTING.
#' @param hom.in Homology towards the feature site, e. g. 5\' landing sequence
#' for C-terminal tagging.
#' @param hom.ex Homology away from the feature site, e. g. 3\' landing sequence
#' for C-terminal tagging.
#' @param ID Positional argument (not used).
#' @export

LANDING_SEQUENCE_ADJUSTMENT <- data.frame(
  # homology lengths towards the ORF
  hom.in = c(40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 
             39, 38, 37, 36, 35, 34, 33, 32, 31, 
             rep(30, 27)),
  # homology lengths outwards the ORF
  hom.ex = c(40, 39, 38, 37, 36, 35, 34, 33, 32, 31, 30, 29, 28, 27, 
             27, 27, 27, 27, 27, 27, 27, 27, 27, 
             seq(27, 1)),
  # absolute amount of cleavage distance
  ID = seq(0, 49)
)

#' Generate pre-CP Element Oligo Pool Sequences for CASTING
#'
#' This function will expand a (processed) CASTING lookup table which contains 
#' CRISPR target sequences with the respective oligo pool sequences that must be
#' ordered as pre-CP elements according to the pre-CP pool template schemed.
#' @details The CASTING lookup table must contain at least the \code{ls_5_full}, 
#' the \code{ls_3_full}, and the \code{target_sequence} sequences, as well as the 
#' inferred (signed/directional) \code{cleavage_distance} from the landing 
#' sequence junction.
#' @param subject A processed CASTING lookup table (\code{data.table} object);
#' see details.
#' @param template A string template of the pre-CP element, in which “\code{str_C}” 
#' is replaced by the crRNA sequence, “\code{str_3}” by the adjusted 3\' landing 
#' sequence, and “\code{str_5}” by the adjusted 5\' landing sequence. Landing 
#' sequence adjustment is defined in \code{LANDING_SEQUENCE_ADJUSTMENT}. 
#' The \code{str_X} parameters can be used to force the reverse complement
#' of the respective element.
#' @param avoid_pattern A character vector of any pattern to avoid during oligo
#' pool design. This can contain IUPAC ambiguity letters from base 3 onwards. 
#' @return A \code{data.table} compatible with CASTING. All sequences are
#' returned in sense orientation.
#' @import data.table
#' @importFrom stringr str_replace
#' @importFrom Biostrings DNAString reverseComplement PDict vcountPDict 
#' @export
#' @examples
#' generate_oligo_sequence(my_processed_lookup_table, 
#'                        "aaastr_Ctttstr_5gagcgggcstr_3cggc",
#'                         avoid_pattern = c("GATATC", "AAGCTT"))

generate_oligo_sequence <- function(subject, template, str_C.RC = FALSE,
                                    str_5.RC = FALSE,  str_3.RC = FALSE,
                                    avoid_pattern = character(0)) {
  
  if (is.null(attr(subject, "feature_site"))) stop(
    "could not identify intended site of genomic feature manipulation 
     that 'subject' was created for, attribute 'feature_site' missing"
  )
  
  # LANDING SEQUENCE ADJUSTMENT/TRIMMING ---------------------------------------
  
  if (as.logical(attr(subject, "feature_site"))) {
    
    subject[, ls_5_homology := LANDING_SEQUENCE_ADJUSTMENT$hom.ex[
      abs(cleavage_distance) + 1]]
    subject[, ls_3_homology := LANDING_SEQUENCE_ADJUSTMENT$hom.in[
      abs(cleavage_distance) + 1]]
    
  } else {
    
    subject[, ls_5_homology := LANDING_SEQUENCE_ADJUSTMENT$hom.in[
      abs(cleavage_distance) + 1]]
    subject[, ls_3_homology := LANDING_SEQUENCE_ADJUSTMENT$hom.ex[
      abs(cleavage_distance) + 1]]
    
  }
  
  subject[cleavage_distance < 1, ls_5_seq_homology := substr(
    ls_5_full,
    start = nchar(ls_5_full) - ls_5_homology + cleavage_distance + 1,
    stop  = nchar(ls_5_full)
  )]
  subject[cleavage_distance < 1, ls_3_seq_homology := substr(
    ls_3_full,
    start = 1,
    stop  = ls_3_homology
  )]

  subject[cleavage_distance > 0, ls_5_seq_homology := substr(
    ls_5_full,
    start = nchar(ls_5_full) - ls_5_homology + 1,
    stop  = nchar(ls_5_full)
  )]
  subject[cleavage_distance > 0, ls_3_seq_homology := substr(
    ls_3_full,
    start = 1,
    stop  = ls_3_homology + cleavage_distance
  )]
  
  # OLIGO SEQUENCE ASSEMBLY ----------------------------------------------------
  
  if (str_C.RC) subject[, target_sequence := as.character(
    Biostrings::reverseComplement(Biostrings::DNAStringSet(
                                  subject[, target_sequence])))]
  if (str_3.RC) subject[, ls_3_seq_homology := as.character(
    Biostrings::reverseComplement(Biostrings::DNAStringSet(
                                  subject[, ls_3_seq_homology])))]
  if (str_5.RC) subject[, ls_5_seq_homology := as.character(
    Biostrings::reverseComplement(Biostrings::DNAStringSet(
                                  subject[, ls_5_seq_homology])))]

  subject[, oligo_sequence := template]
  
  subject[, oligo_sequence := mapply(stringr::str_replace, oligo_sequence, 
                                     "str_3", ls_3_seq_homology)]
  subject[, oligo_sequence := mapply(stringr::str_replace, oligo_sequence, 
                                     "str_5", ls_5_seq_homology)]
  subject[, oligo_sequence := mapply(stringr::str_replace, oligo_sequence, 
                                     "str_C", target_sequence)]
  
  # remove helper columns
  subject[, ls_5_seq_homology := NULL]
  subject[, ls_3_seq_homology := NULL]
  
  if (str_C.RC) subject[, target_sequence := as.character(Biostrings::reverseComplement(
                      Biostrings::DNAStringSet(subject[, target_sequence])))]
  
  # ANNOTATION OF INCOMPATIBILITIES --------------------------------------------
  
  if (length(avoid_pattern) > 0) {
    
    # All sequences are supposed to be 5' to 3', no need to reverse complement.
    tmp.avoid_pattern <- Biostrings::vcountPDict(Biostrings::PDict(avoid_pattern, 
                                     tb.start = 1, tb.width = 3),
                         Biostrings::DNAStringSet(subject[, oligo_sequence]),
                                     fixed = FALSE, collapse = 2)
    
    subject[, oligo_penalty := tmp.avoid_pattern]
    
  } else {
    
    subject[, oligo_penalty := 0]
    
  }
  
}

#' Generate a Reference Genome for post-CASTING NGS Analysis
#'
#' This function will return a (processed) CASTING lookup table which contains 
#' CRISPR target sequences with the respective reference genome sequences that
#' should be used for Vectorette NGS.
#' @details The CASTING lookup table must contain at least the \code{ls_5_full}, 
#' the \code{ls_3_full}, and the \code{target_sequence} sequences, as well as the 
#' inferred (signed/directional) \code{cleavage_distance} from the landing 
#' sequence junction.
#' @param subject A processed CASTING lookup table (\code{data.table} object);
#' see details.
#' @param template A string template of the pre-CP element, in which “\code{str_C}” 
#' is replaced by the crRNA sequence, “\code{str_3}” by the adjusted 3\' landing 
#' sequence, and “\code{str_5}” by the adjusted 5\' landing sequence. Overlapping
#' sequence material with the \code{tagging_cassette} will be removed. 
#' @param tagging_cassette The sequence of the tagging cassette intended to be
#' integrated at the genomic site.
#' @param primer The Vectorette initiating primer for NGS. Use the annealing part
#' only.
#' @return A \code{data.table} compatible with CASTING. All sequences are
#' returned in sense orientation.
#' @import data.table
#' @importFrom utils head tail
#' @importFrom stringr str_locate str_replace
#' @importFrom Biostrings DNAString reverseComplement matchPattern
#' @export
#' @examples
#' my_ref_table <- generate_reference_genome(my_processed_lookup_table, 
#'                          "aaastr_Ctttstr_5gagcgggcstr_3cggc",
#'                           tagging_cassette = my_cassette,
#'                           primer = "gtcgacctgcagcgtacg")
#' 
#' data.table::fwrite(my_ref_table[, .(oligo_id, junction)], sep = "\n",
#'        row.names = FALSE, col.names = FALSE, file = my_file)
                    
generate_reference_genome <- function(subject, template, tagging_cassette, 
                                      primer) {
  
  if (is.null(attr(subject, "feature_site"))) stop(
    "could not identify intended site of genomic feature manipulation 
     that 'subject' was created for, attribute 'feature_site' missing"
  )
  
  feature_site <- as.logical(attr(subject, "feature_site"))
  
  # IDENTIFICATION OF PRIMER BINDING -------------------------------------------
  
  # cast character vector to DNAString object
  dna.prm <- Biostrings::DNAString(primer)
  dna.cas <- Biostrings::DNAString(tagging_cassette)
  
  # the 3' LS will be sequenced when:
  primer_bind.fwd <- Biostrings::matchPattern(dna.prm, 
                                              dna.cas)@ranges
  # the 5' LS will be sequenced when:
  primer_bind.rev <- Biostrings::matchPattern(reverseComplement(dna.prm), 
                                              dna.cas)@ranges
  
  if (!any(primer_bind.rev@start, primer_bind.fwd@start)) 
    stop("primer does not bind in tagging cassette")
  
  # IDENTIFICATION OF TEMPLATE/CASSETTE OVERLAP --------------------------------
  
  find_prefix <- function(str_A, str_B) for (i in 1:nchar(str_A)) 
    if (startsWith(tolower(str_B), substr(tolower(str_A), i, nchar(str_A)))) return(i)
  find_suffix <- function(str_A, str_B) for (i in nchar(str_A):1) 
    if (endsWith(tolower(str_B), substr(tolower(str_A), 1, i))) return(i)
  
  # REPLACEMENT OF TARGET SEQUENCE ---------------------------------------------
  
  replace_str_C <- function(x, str_C.RC = FALSE) {
    
    if (str_C.RC) x[, target_sequence := as.character(
      Biostrings::reverseComplement(Biostrings::DNAStringSet(x[, target_sequence])))]
      
    x[, junction := mapply(stringr::str_replace, junction, 
                                         "str_C", target_sequence)]      
    
    if (str_C.RC) x[, target_sequence := as.character(
      Biostrings::reverseComplement(Biostrings::DNAStringSet(x[, target_sequence])))]
    
    x[, oligo_id := paste(oligo_id, target_sequence, sep = "|")]
    
  }
  
  # ASSEMBLY OF REFERENCE GENOME SEQUENCE --------------------------------------
  
  # add information on potential sequence redundancy
  sub.subject <- merge(subject, subject[, .N, by =  ORF_name.long], 
                                              by = "ORF_name.long")
  
  sub.subject[, oligo_id := paste(">lcl", ORF_name.long, sep = "|")]
  
  # The primer binds on forward strand of tagging cassette, the 5' site of the
  # template becomes relevant up to the 3' LS.
  
  if (any(primer_bind.fwd@start)) {
    
    seq_from.cas <- substr(tagging_cassette, start = primer_bind.fwd@start, 
                           stop = nchar(tagging_cassette))
    
    seq_from.tmp <- utils::head(unlist(stringr::str_split(template, "str_3")), 1)
    
    # check/correct overlap
    
    tmp.pos <- find_suffix(seq_from.tmp, seq_from.cas)
    if (!is.null(tmp.pos)) seq_from.tmp <- substr(seq_from.tmp, 
                                 start = tmp.pos + 1, stop = nchar(seq_from.tmp))
    
    sub.subject[, junction := paste0(seq_from.cas, seq_from.tmp, ls_3_full)]
    
    if (grepl("str_C", seq_from.tmp)) {
      
      replace_str_C(sub.subject, feature_site)
      
    } else {

      sub.subject <- sub.subject[sub.subject[, .I[which.max(ls_3_homology
      + ifelse(cleavage_distance > 0, cleavage_distance, 0))], by = ORF_name.long]$V1]        
      
    }
    
    sub.subject[, oligo_id := paste(oligo_id, nchar(junction) - nchar(ls_3_full) 
                  + ls_3_homology + ifelse(cleavage_distance > 0, 
                                                 cleavage_distance, 0) + 1,
                  nchar(junction), N, sep = "|")]
    
  }
  
  # The primer binds on reverse strand of tagging cassette, the 3' site of the
  # template becomes relevant up to the 5' LS.
  
  if (any(primer_bind.rev@start)) {

    seq_from.cas <- substr(tagging_cassette, start = 1,
                           stop = primer_bind.rev@start + primer_bind.rev@width - 1)
    
    seq_from.tmp <- utils::tail(unlist(stringr::str_split(template, "str_5")), 1)
    
    # check/correct overlap
    
    tmp.pos <- find_prefix(seq_from.tmp, seq_from.cas)
    if (!is.null(tmp.pos)) seq_from.tmp <- substr(seq_from.tmp, 
                                 start = 1, stop = tmp.pos - 1)
    
    sub.subject[, junction := paste0(ls_5_full, seq_from.tmp, seq_from.cas)]
    
    if (grepl("str_C", seq_from.tmp)) {
      
      replace_str_C(sub.subject, feature_site)
    
    } else {

      sub.subject <- sub.subject[sub.subject[, .I[which.max(ls_5_homology
      - ifelse(cleavage_distance < 0, cleavage_distance, 0))], by = ORF_name.long]$V1]        
      
    }
    
    sub.subject[, oligo_id := paste(oligo_id, 1, nchar(ls_5_full) 
                  - ls_5_homology + ifelse(cleavage_distance < 0, 
                                                 cleavage_distance, 0),
                  N, sep = "|")]    
    
  }
  
  # re-sample
  sub.subject <- sub.subject[sample(1:nrow(sub.subject)), ]
  
  return(sub.subject)

}
