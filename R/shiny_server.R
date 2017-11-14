# SHINY SERVER =================================================================

CASTLING_server <- function(session, input, output) {
  
  # RESPONSIVE UI ELEMENTS =====================================================
  
  # responsive on taging cassette ----------------------------------------------
  
  output$inp_vec_cass_text <- renderUI({
    
    if (nchar(input$inp_vec_cass) == 0) em("Enter feature cassette") else
      em(paste0("Feature cassette (", nchar(input$inp_vec_cass), " nt)"))
    
  })
  
  # genome settings summary ----------------------------------------------------
  
  output$inp_genome_summary.1 <- renderUI({
    
    return(HTML(paste0(
      "<p><table style='width:100%'>
  <tr>
    <td>GFF feature type</td>
    <td>", 
      ifelse(is.null(input$inp_feature_type), 
             "(upload a genome annotation)",
             paste(input$inp_feature_type, collapse = ", ")
      ), "</td> 
  </tr>
  <tr>
    <td>site to engineer</td>
    <td>",
      ifelse(is.null(input$inp_feature_dist),
             def_feature_dist,
             input$inp_feature_dist
      ),
      ifelse(is.null(input$inp_feature_site),
             ifelse(def_feature_site, " nt after START", " nt before STOP"),
             ifelse(input$inp_feature_site, " nt after START", " nt before STOP")
      ), "</td> 
  </tr>
  <tr>
    <td>area to target</td>
    <td>± ", 
      ifelse(is.null(input$inp_feature_look), 
             def_feature_look, input$inp_feature_look
      ),
      ifelse(is.null(input$inp_feature_site),
             ifelse(def_feature_site, " nt around START", " nt around STOP"),
             ifelse(input$inp_feature_site, " nt around START", " nt around STOP")
      ), "</td> 
  </tr>
</table></p>")))
    
  })
  
  output$inp_genome_summary.2 <- renderUI({
    
    return(HTML(paste0(
      "<p><table style='width:100%'>
  <tr>
    <td>PAM selectivity</td>
    <td>", 
      ifelse(is.null(input$inp_endonuclease_pam), 
             endonucleases[name %in% input$inp_endonuclease, pattern],
             paste(input$inp_endonuclease_pam, collapse = ", ")
      ), 
      ifelse(ifelse(is.null(input$inp_target_is_upstream),
             def_target_is_upstream, as.logical(input$inp_target_is_upstream)) == TRUE,
             " (downstream of protospacer)",
             " (upstream of protospacer)"
             ), "</td> 
  </tr>
  <tr>
    <td>CRISPR spacer</td>
    <td>",
      ifelse(is.null(input$inp_target_totlen),
             def_target_totlen, input$inp_target_totlen
      ), " nt (", 
      ifelse(is.null(input$inp_target_corlen),
             def_target_corlen, input$inp_target_corlen
      ), " nt seed)</td> 
  </tr>
  <tr>
    <td>target cleavage</td>
    <td>",
      ifelse(is.null(input$inp_target_cleave),
             def_target_cleave, input$inp_target_cleave
      ),
      ifelse(ifelse(is.null(input$inp_target_is_upstream),
             def_target_is_upstream, as.logical(input$inp_target_is_upstream)) == TRUE,
             " nt upstream of PAM",
             " nt downstream of PAM"
      ), "</td> 
  </tr>
</table></p>")))
    
  })

  
  # responsive on inp_target_is_upstream ------------------------------------------
  
  output$inp_cleave_at_text <- renderUI({
    
    tmp.text <- paste("The target is cleaved after …", 
                      ifelse(as.logical(input$inp_target_is_upstream), "upstream", "downstream"),
                      "of the PAM.")
    
    return(helpText(tmp.text))
    
  })
  
  # responsive on inp_target_totlen --------------------------------------------
  
  output$inp_target_corlen <- renderUI({
    
    sliderInput("inp_target_corlen", label = NULL, 
                round = TRUE, post = " nt",
                min = 5, max = input$inp_target_totlen, step = 1, 
                value = def_target_corlen)
  
  })
  
  output$inp_cleave_at <- renderUI({
    
    sliderInput("inp_cleave_at", label = NULL,
                round = TRUE, post = " nt",
                min = 0, max = input$inp_target_totlen, step = 1, 
                value = def_target_cleave)
  
  })
  
  # responsive on inp_feature_site ----------------------------------------------
  
  output$inp_feature_dist_text <- renderUI({
    
    tmp.text <- "To engineer the loci of the selected GFF features, the ‘self-integrating cassette’ (SIC) will integrate itself"
    
    if (input$inp_feature_site) {
      
      if (input$inp_feature_dist == 0) {
        
        tmp.text <- paste(tmp.text, "right at the feature start.")
        
      } else {
        
        tmp.text <- paste(tmp.text, "after the first ", input$inp_feature_dist, 
                          "leading bases.")
        
      }
      
    } else {
      
      if (input$inp_feature_dist == 0) {
        
        tmp.text <- paste(tmp.text, "straight at the feature end.")
        
      } else {
        
        tmp.text <- paste(tmp.text, "before the last", input$inp_feature_dist, 
                          "trailing bases.")
        
      }
      
    }
    
    return(helpText(tmp.text))
    
  })
  
  # others ---------------------------------------------------------------------
  
  output$out_oligo_tmp <- renderUI({
    
    HTML(paste0("<p>5'—<tt>", 
                stringr::str_replace_all(
                stringr::str_replace_all(parse_oligo_template(), 
                               "str_(.{1})", "</tt> (\\1'-HA) <wbr><tt>"), "C'-HA", "CRISPR spacer"), 
                "</tt>—3'</p>"))
    
  })
  
  # OBSERVERS ==================================================================
    
  # input FEATURE TYPES --------------------------------------------------------
  
  observe({
    
    if (is.null(input$inp_genome_ann)) {
      
      feature_type_choices <- character(0)

    } else {
      
      feature_type_choices <- unique(get_genome_annotation()$type)
      
    }
    
    updateSelectizeInput(session, "inp_feature_type", label = NULL,
                      choices  = feature_type_choices,
                      selected = def_feature_type
    )
    
  })

  # input PAM ------------------------------------------------------------------

  observe({

    updateSelectInput(session, "inp_endonuclease_pam",
       selected = endonucleases[name %in% input$inp_endonuclease, pattern])

  })

  # input OLIGO STRUCTURE TEMPLATE ---------------------------------------------

  observe({

    if (input$inp_feature_site) {

      updateSelectizeInput(session, "inp_oligo_tmp",
        choices = list(
          "Linker" = oligo.linkers_RC,
          "Sequence wildcards"  = list("3'-HA" = "str_3"),
          "Restriction Sites"   = oligo.restriction_sites,
          "Sequence wildcards"  = list("5'-HA" = "str_5"),
          "Minimal Terminators" = oligo.terminators_RC,
          "Sequence wildcards"  = list("CRISPR spacer" = "str_C",
                                      "handle" = "str_H"),
          "Promoter Homologies" = oligo.homologies_RC
        ),
        selected = c(
          oligo.linkers_RC["S4 linker"],
          "str_3",
          oligo.restriction_sites[[1]],
          "str_5",
          oligo.terminators_RC[[1]],
          "str_C", "str_H",
          oligo.homologies_RC["SNR52"]
        )
      )

    } else {

      updateSelectizeInput(session, "inp_oligo_tmp",
        choices  = list(
          "Promoter Homologies" = oligo.homologies_FC,
          "Sequence wildcards"  = list("handle" = "str_H",
                                       "CRISPR spacer"  = "str_C"),
          "Minimal Terminators" = oligo.terminators_FC,
          "Sequence wildcards"  = list("3'-HA" = "str_3"),
          "Restriction Sites"   = oligo.restriction_sites,
          "Sequence wildcards"  = list("5'-HA" = "str_5"),
          "Linker" = oligo.linkers_FC
        ),
        selected = c(
          oligo.homologies_FC["SNR52"],
          "str_H", "str_C",
          oligo.terminators_FC[[1]],
          "str_3",
          oligo.restriction_sites[[1]],
          "str_5",
          oligo.linkers_FC["S1/S3 linker"]
        )

      )

    }

  })

  # input ACTION BUTTONS -------------------------------------------------------

  observe({

    if (is.null(input$inp_genome_ass) |
        is.null(input$inp_genome_ann) |
        length( input$inp_feature_type) < 1 |
        length( input$inp_endonuclease_pam) < 1
        ) {

      disableActionButton("run_genome", session)
      disableDownloadButton("run_genome_download", session)

    } else {

      enableActionButton("run_genome", session)
      if (input$run_genome > 0) enableDownloadButton("run_genome_download", 
                                                     session)

    }

  })

  observe({

    if (exists(get_lookup_table_name()) | input$run_genome > 0) {

      enableActionButton("run_oligos", session)
      enableActionButton("run_vectorette", session)
      
      if (input$run_oligos > 0) {
        
        enableDownloadButton("run_oligos_download", session)
        enableDownloadButton("run_vectorette_download", session)
        
      }

    } else {

      disableActionButton("run_oligos", session)
      disableActionButton("run_vectorette", session)
      disableDownloadButton("run_oligos_download", session)
      disableDownloadButton("run_vectorette_download", session)

    }

  })
  
  # REACTIVES ==================================================================
  
  # PARSE GENOME ASSEMBLY/ANNOTATION -------------------------------------------
  
  get_genome_assembly <- reactive({
    
    if (is.null(input$inp_genome_ass)) return(NULL)
    
    tmp.genome_assembly <- input$inp_genome_ass
    
    file.rename(tmp.genome_assembly$datapath,
      file.path(dirname(tmp.genome_assembly$datapath), 
      tmp.genome_assembly$name))
    
    # returns a DNAStringSet object:
    Biostrings::readDNAStringSet(filepath = file.path(
      dirname(tmp.genome_assembly$datapath), tmp.genome_assembly$name))
    
  })
  
  get_genome_annotation <- reactive({
    
    if (is.null(input$inp_genome_ann)) return(NULL)

    tmp.genome_annotation <- input$inp_genome_ann
    
    file.rename(tmp.genome_annotation$datapath,
      file.path(dirname(tmp.genome_annotation$datapath), 
      tmp.genome_annotation$name))
    
    # returns a data.frame object:
    tmp.table <- readr::read_delim(file = file.path(
      dirname(tmp.genome_annotation$datapath), tmp.genome_annotation$name), 
      delim = "\t", col_names = FALSE, comment = "#", col_types = "ccciicccc")
    
    colnames(tmp.table) <- c("seqid", "source", "type", 
     "start", "end", "score", "strand", "phase", "attribute")
    
    return(tmp.table)
    
  })
  
  get_lookup_table_name <- reactive({
    
    return(paste0("cast-", paste(input$inp_endonuclease_pam, collapse = "_")))
    
  })
  
  # parse OLIGO STRUCTURE/TEMPLATE ---------------------------------------------
  
  parse_oligo_template <- reactive({

    tmp.template <- input$inp_oligo_tmp # will be modified
    
    pos_3 <- which(tmp.template == "str_3")
    pos_5 <- which(tmp.template == "str_5")
    pos_H <- which(tmp.template == "str_H")
    pos_C <- which(tmp.template == "str_C")
    
    # replace placeholder with correct handle
    if (length(pos_H) > 0) tmp.template[pos_H] <- unique(endonucleases[name %in%
                                               input$inp_endonuclease, handle])
    
    # expected element order (not tested for at this stage)
    #
    # upstream modification:   (5')-...-str_3-...-str_5-...-str_C-str_H-(3')
    # downstream modification: (5')-str_H-str_C-...-str_3-...-str_5-...-(3')
    
    # For upstream feature modification, it is necessary to reverse complement
    # each element of the template, but the landing sequences and restriction
    # sites. Note that the crRNA needs to be reverse complemented afterwards.

    if (input$inp_feature_site) {
      
      wrong.order <- is.unsorted(c(pos_3, pos_5, pos_C, pos_H))
      
      # The following code can be used to reverse complement the ends (and other
      # portions thought to be required) of the oligo template automatically.
      #
      # tmp.order <- c(0, which(startsWith(tmp.template, "str_")), 
      #                length(tmp.template) + 1)
      # tmp.valid <- which(diff(c(1, which(startsWith(tmp.template, "str_")), 
      #                           length(tmp.template))) > 0)
      # 
      # if (length(tmp.valid) > 1) {
      #   
      #   # at least the first element is not "str_X" and can be replaced
      #   tmp.range <- (tmp.order[tmp.valid[1]] + 1):(
      #     tmp.order[tmp.valid[1 + 1]] - 1)
      #   tmp.template[tmp.range] <- lapply(tmp.template[tmp.range], function(x)
      #     tolower(reverseComplement(DNAString(x))))
      #   
      # }
      # 
      # if (length(tmp.valid) > 1) {
      #   
      #   # at least the last element is not "str_X" and can be replaced
      #   tmp.range <- (tmp.order[tmp.valid[length(tmp.valid)]] + 1):(
      #     tmp.order[tmp.valid[length(tmp.valid)] + 1] - 1)
      #   tmp.template[tmp.range] <- lapply(tmp.template[tmp.range], function(x)
      #     tolower(reverseComplement(DNAString(x))))
      #   
      # }
      # 
      # if (length(tmp.valid) > 2) {
      # 
      #   # at least the before-last element is not "str_X" and can be replaced
      #   tmp.range <- (tmp.order[tmp.valid[length(tmp.valid) - 1]] + 1):(
      #     tmp.order[tmp.valid[length(tmp.valid)]] - 1)
      #   tmp.template[tmp.range] <- lapply(tmp.template[tmp.range], function(x)
      #     tolower(reverseComplement(DNAString(x))))
      # 
      # }

    } else {
      
      wrong.order <- is.unsorted(c(pos_H, pos_C, pos_3, pos_5))
      
    }
    
    if (wrong.order) warning("oligo template may have wrong order of elements")
    
    return(paste0(tmp.template, collapse = ""))

  })
  
  # create CASTLING LOOK-UP TABLE -----------------------------------------------
  
  create_look_up_table <- eventReactive(input$run_genome, {
    
    # In order to have a nice progress bar, the set is split up into smaller
    # portions which are processed sequentially.
    
    tmp.subsets_1 <<- list() # will contain look-up after find_features()
    tmp.subsets_2 <<- list() # will contain look-up after find_targets()
    
    steps.max <- 4 * length(get_genome_assembly()) # slice by landmarks
    steps.mrk <- rep(1:length(get_genome_assembly()), each = 4)
    steps.val <- round(seq(from = 1, to = nrow(get_genome_annotation()), 
                           length.out = steps.max)) # equivalent portions
    
    # PROGRESS FOR FEATURE SELECITON
    
    withProgress(message = NULL, value = 0, {
      
      for (i in 1:(steps.max - 1)) {
        
        tmp.subsets_1[[i]] <- find_features(
          genome_assembly   = get_genome_assembly(),
          genome_annotation = get_genome_annotation()[(steps.val[i]):(
                                                       steps.val[i + 1] - 1), ],
          feature_type = input$inp_feature_type,
          feature_site = input$inp_feature_site,
          offset = input$inp_feature_dist,
          range_width = max_feature_look
        )
        
        incProgress(amount  = 0.5 * 1 / steps.max, 
                    message = "Creating genomic lookup table:",
                    detail  = paste("Scaffold", steps.mrk[i]))
        
      }
      
      # PROGRESS FOR crRNA DETERMINATION
      
      for (p in input$inp_endonuclease_pam) {
        
      for (i in 1:(steps.max - 1)) {
        
        tmp.subsets_2[[p]][[i]] <- find_targets(
          subject = tmp.subsets_1[[i]],
          pattern = p,
          look_around = input$inp_feature_look,
          target_upstream = as.logical(input$inp_target_is_upstream),
          # If the following panels were not selected, those values might not be
          # set to default, so treat them as default:
          target_length = ifelse(is.null(input$inp_target_totlen),
                                 def_target_totlen, input$inp_target_totlen),
          target_core   = ifelse(is.null(input$inp_target_corlen),
                                 def_target_corlen, input$inp_target_corlen),
          cleave_at     = ifelse(is.null(input$inp_cleave_at),
                                 def_target_cleave, input$inp_cleave_at),
          reference_genome = get_genome_assembly(),
          mismatches = unique(na.omit(as.numeric(input$inp_mismatches)))
        )
        
        incProgress(amount  = 0.5 * 1 / (steps.max * length(input$inp_endonuclease_pam)),
                    message = paste0("Identifying targets for ", p, ":"),
                    detail  = HTML(paste("Scaffold", steps.mrk[i])))
        
      }}
      
      tmp.look_up_table <- data.table::rbindlist(
        unlist(tmp.subsets_2, recursive = FALSE)
      )
      
      data.table::setattr(tmp.look_up_table, "feature_site", 
                          as.logical(input$inp_feature_site))
      data.table::setattr(tmp.look_up_table, "offset",
                          as.integer(input$inp_feature_dist))
      
      # make this result globally available
      assign(get_lookup_table_name(), tmp.look_up_table, envir = .GlobalEnv)
      
      return(tmp.look_up_table)

    })
    
  })
  
  # process CASTLING LOOK-UP TABLE ----------------------------------------------
  
  process_look_up_table <- eventReactive(input$run_oligos, {
    
    if(exists(get_lookup_table_name())) {
      
      tmp.look_up_table <- eval(parse(text = paste0("`", 
                                                get_lookup_table_name(), "`")))

    } else {
      
      tmp.look_up_table <- create_look_up_table()
      
    }

    tmp.processed_table <- generate_oligo_sequence(
      subject  = tmp.look_up_table,
      template = parse_oligo_template(),
      str_C.RC = attr(tmp.look_up_table, "feature_site"),
      avoid_pattern = input$inp_oligo_rem
        )
    
    # modify this table globally
    assign(get_lookup_table_name(), tmp.processed_table, envir = .GlobalEnv)
    
    return(tmp.processed_table)
    
  })
  
  # create REFERENCE GENOME ----------------------------------------------------
  
  create_reference_genome <- eventReactive(input$run_vectorette, {
    
    if(exists(get_lookup_table_name())) {
      
      tmp.look_up_table <- eval(parse(text = paste0("`", 
                                                get_lookup_table_name(), "`")))

    } else {
      
      tmp.look_up_table <- create_look_up_table()
      
    }    
    
    tmp.reference_genome <- generate_reference_genome(
      subject  = tmp.look_up_table, 
      template = parse_oligo_template(), 
      tagging_cassette = input$inp_vec_cass, 
      primer = input$inp_vec_prim)
    
    return(tmp.reference_genome)

  })
  
  # OUTPUTS --------------------------------------------------------------------
  
  output$run_genome_download <- downloadHandler(
    
    filename = function() { paste0(get_lookup_table_name(), "-raw.tar.gz") },
    
    content  = function(file) { 
      
      if (length(input$set_input) < 1) return(NULL)
      
      tmp.dir <- file.path(tempfile(pattern = "cast", tmpdir = "."))
      
      dir.create(tmp.dir)
      
      if ("RDS" %in% input$set_input) saveRDS(create_look_up_table(), 
                file = file.path(tmp.dir, paste0(get_lookup_table_name(), ".rds")))
      
      if ("CSV" %in% input$set_input) data.table::fwrite(create_look_up_table(), 
                file = file.path(tmp.dir, paste0(get_lookup_table_name(), ".csv")))
      
      tmp.status <- tar(tarfile = file, files = tmp.dir)
      
      unlink(tmp.dir, recursive = TRUE)
      
      return(file)
      
    }
    
  )
  
  output$run_oligos_download <- downloadHandler(
    
    filename = function() { paste0(get_lookup_table_name(), "-annotated.tar.gz") },
    
    content  = function(file) { 
      
      if (length(input$set_input) < 1) return(NULL)
      
      tmp.dir <- file.path(tempfile(pattern = "cast", tmpdir = "."))
      
      dir.create(tmp.dir)
      
      if ("RDS" %in% input$set_input) saveRDS(process_look_up_table(), 
                file = file.path(tmp.dir, paste0(get_lookup_table_name(), ".rds")))
      
      if ("CSV" %in% input$set_input) data.table::fwrite(process_look_up_table(), 
                file = file.path(tmp.dir, paste0(get_lookup_table_name(), ".csv")))
      
      tmp.status <- tar(tarfile = file, files = tmp.dir)
      
      unlink(tmp.dir, recursive = TRUE)
      
      return(file)
      
    }
    
  )
  
  output$run_vectorette_download <- downloadHandler(
    
    filename = function() { paste0(get_lookup_table_name(), "-reference.tar.gz") },
    
    content  = function(file) { 
      
      if (length(input$set_input) < 1) return(NULL)
      
      tmp.dir <- file.path(tempfile(pattern = "cast", tmpdir = "."))
      
      dir.create(tmp.dir)
      
      if ("RDS" %in% input$set_input) saveRDS(create_reference_genome(), 
                file = file.path(tmp.dir, paste0(get_lookup_table_name(), "-reference.rds")))
      
      if ("CSV" %in% input$set_input) data.table::fwrite(
        create_reference_genome()[, .(seq.ID, junction)], sep = "\n",
        row.names = FALSE, col.names = FALSE,
        file = file.path(tmp.dir, paste0(get_lookup_table_name(), "-reference.fasta"))
        )
      
      tmp.status <- tar(tarfile = file, files = tmp.dir)
      
      unlink(tmp.dir, recursive = TRUE)
      
      return(file)
      
    }
    
  )
  
}
