# INTERNAL DATA RESOURCES ======================================================

# ENDONUCLEASE PROPERTIES ------------------------------------------------------

endonucleases <- data.table(
    name    = c("FnCpf1", "FnCpf1", "LbCpf1", "AsCpf1", "MbCpf1"),
    pattern = c("TTV", "TYN", "HYBV", "TTTV", "TYYV"),
    handle  = c(
        "taatttctactgttgtagat",  # FnCpf1
        "taatttctactgttgtagat",  # - " -
        "taatttctactaagtgtagat", # LbCpf1
        "taatttctactcttgtagat",  # AsCpf1
        "aaatttctactgtttgtagat"  # MbCpf1
    ),
    key = "pattern"
)

# SEQUENCE JUNKS FOR OLIGO POOL DESIGN -----------------------------------------

oligo.homologies_FC <- list(
    "SNR52" = "gataaatgatc"
)

oligo.homologies_RC <- list(
    "SNR52" = "gatcatttatc"
)

oligo.linkers_FC <- list(
    "S1/S3 linker"  = "cgtacgctgcaggtcgac",
    "S2 linker"     = "atcgatgaattcgagctcg",
    "S4 linker"     = "catcgatgaattctctgtcg",
    "NgoMIV linker" = "gccggcaatacttctact"
)

oligo.linkers_RC <- list(
    "S1/S3 linker"  = "gtcgacctgcagcgtacg",
    "S2 linker"     = "cgagctcgaattcatcgat",
    "S4 linker"     = "cgacagagaattcatcgatg",
    "NgoMIV linker" = "agtagaagtattgccggc"
)

oligo.terminators_FC <- list(
    "SUP4" = "tttttttt"          # minimal Pol(III) terminator
)

oligo.terminators_RC <- list(
    "SUP4" = "aaaaaaaa"          # minimal Pol(III) terminator
)

# It is usually unnecassary to reverse complement the restriction sites as well
# as the landing sequences in any scenario.

oligo.restriction_sites <- list(
    "BstXI A" = "ccatcaacttgg",  # BstXI site cutting only
    "BstXI B" = "accacggcgctggt" # BstXI site cutting and amplification
)

primers.vectorette <- list(
  "mCherry (5'—CTCCTCGCCCTTGCTCACC—3')" = "ctcctcgcccttgctcacc",
  "S1/S3 (5'—GTCGACCTGCAGCGTACG—3')" = "gtcgacctgcagcgtacg",
  "SNR52 (5'—GCAGTGAAAGATAAATGATCTAATTTCTACTG—3')" = "gcagtgaaagataaatgatctaatttctactg"
)

# DEFAULTS ---------------------------------------------------------------------

default.endonuclease  <- "FnCpf1"
default.look_around   <- 200L
default.target_length <-  20L
default.target_core   <-   8L
default.cleave_at     <-  15L

# SHINY HELPER FUNCTIONS =======================================================

textAreaInput2 <- function (inputId, label, value = "", width = NULL, height = NULL, 
    cols = NULL, rows = NULL, placeholder = NULL, resize = NULL) 
{
    value <- restoreInput(id = inputId, default = value)
    if (!is.null(resize)) {
        resize <- match.arg(resize, c("both", "none", "vertical", 
            "horizontal"))
    }
    style <- paste("max-width: 100%; font-family: monospace;", if (!is.null(width)) 
        paste0("width: ", validateCssUnit(width), ";"), if (!is.null(height)) 
        paste0("height: ", validateCssUnit(height), ";"), if (!is.null(resize)) 
        paste0("resize: ", resize, ";"))
    if (length(style) == 0) 
        style <- NULL
    div(class = "form-group", 
        tags$label(label, `for` = inputId), tags$textarea(id = inputId, 
        class = "form-control", placeholder = placeholder, style = style, 
        rows = rows, cols = cols, value))
}

`%AND%` <- function(x, y) {
  if (!is.null(x) && !is.na(x))
    if (!is.null(y) && !is.na(y))
      return(y)
  return(NULL)
}

textInput2 <- function (inputId, label, value = "", width = NULL, placeholder = NULL) 
{
    value <- restoreInput(id = inputId, default = value)
    div(class = "form-group shiny-input-container", style = paste("font-family: monospace;", 
        if (!is.null(width)) paste0("width: ", validateCssUnit(width), ";")), 
        label %AND% tags$label(label, `for` = inputId), tags$input(id = inputId, 
        type = "text", class = "form-control", value = value, 
        placeholder = placeholder))
}

disableDownloadButton <- function(id, session) {
  
  session$sendCustomMessage(type = "jsCode", list(code = 
                    paste0("$('#", id, "').attr('disabled', true)")))
  
}

enableDownloadButton <- function(id, session) {
  
  session$sendCustomMessage(type = "jsCode", list(code = 
                    paste0("$('#", id, "').attr('disabled', false)")))
  
}

disableActionButton <- function(id, session) {
  
  session$sendCustomMessage(type = "jsCode", list(code = 
                    paste0("$('#", id, "').prop('disabled', true)")))
  
}

enableActionButton <- function(id, session) {
  
  session$sendCustomMessage(type = "jsCode", list(code = 
                    paste0("$('#", id, "').prop('disabled', false)")))

}
