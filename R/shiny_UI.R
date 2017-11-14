# SHINY UI =====================================================================

# SETTINGS ---------------------------------------------------------------------

panel.settings <- wellPanel(
  
  checkboxGroupInput("set_input", label = em("Save/Export data as …"),
                      choices = c("CSV and/or FASTA" = "CSV", 
                                  "Single R data objects (RDS)" = "RDS"),
                     selected = c("CSV", "RDS")
  )
  
)

# NGS UTILITY PANEL ------------------------------------------------------------

panel.vectorette <- fluidPage(

  
  uiOutput("inp_vec_cass_text"),
  helpText("The sequence information should be given in sense with respect to 
            the genomic features to be manipulated. Overlaps with the oligo 
            template will be removed from the reference."),
  
  textAreaInput2("inp_vec_cass", label = NULL,
                 width = "100%", rows = 7, resize = "vertical",
                 placeholder = "No cassette sequence provided."),
  
  em("Select or add a vectorette initiating primer"),
  helpText("Include only the annealed primer part in this field."),
  
  # textInput2("inp_vec_prim", label = NULL,
  #           width = "100%",
  #           placeholder = "DNA sequence expected …")
  selectizeInput("inp_vec_prim", label = NULL,
                 width = "100%", 
                 options = list(create = TRUE),
                 choices = primers.vectorette,
                 selected = primers.vectorette[[1]]
  )
  
)

# OLIGO POOL CREATION PANELS ---------------------------------------------------

panel.oligo_str <- fluidRow(

  column(12,  

  tags$div(title = "Wildcards for homology arms and CRISPR spacer will be replaced and reverse complemented if applicable.",
           
           em("Design of locus-specific oligo template"),
           helpText("The structure of this template must be given in sense with respect
           to the genomic features to be manipulated."),
           
           selectizeInput("inp_oligo_tmp", label = NULL, 
                          width = "100%", multiple = TRUE, 
                          options = list(create = TRUE, maxItems = 8),
                          choices = list()
                          
           )
  ),
  
  htmlOutput("out_oligo_tmp")
           
))

panel.oligo_rem <- fluidRow(
  
  column(12,
  
  em("Select or add sequence motifs to avoid"),
  helpText("Exclude restriction sites required for wet lab pool processing."),
  
  selectizeInput("inp_oligo_rem", label = NULL, 
                 width = "100%", multiple = TRUE,
                 options = list(create = TRUE, placeholder = "No motif will be excluded"),
                 choices  = c("CCANNNNNNTGG", "GCCGGC", "ACNNNNGTAYC"),
                 selected = "CCANNNNNNTGG"
  )
                                 
))

# GENOME PROCESSING PANELS -----------------------------------------------------

panel.genome_ref <- fluidRow(
        
    column(
      
      width = 6,
           
      em("Host genome sequence"),
      
      tags$div(title = "Get the full genomic sequence of your host organism in (compressed) FASTA format and upload it here. This file is used to extract target sequences, homology arms, and to identify potential off-target hits.",
               
               helpText("Upload in FASTA format."),
               
               fileInput("inp_genome_ass", label = NULL, width = "100%", accept = )
      )
      
    ),
    
    column(
      
      width = 6,
      
      em("Host genome annotation"),
      
      tags$div(title = "Get the corresponding genome annotation that specifies where individual features such as genes are located. You can also upload a custom-tailored annotation. Please upload in (compressed) ‘generic feature format’ (GFF3).",
               
               helpText("Upload in generic feature fromat (GFF3)."),
               
               fileInput("inp_genome_ann", label = NULL, width = "100%")
               
      )
           
    )
    
)

panel.genome_sum <- tabPanel(
  
  title = "Summary",
  
  em("Properties of loci to engineer"),
  
  uiOutput("inp_genome_summary.1"),
  
  em("Properties of CRISPR system"),
  
  uiOutput("inp_genome_summary.2"),
  
  fluidRow(
        
        column(6, align = "left",
               actionButton("run_genome", 
                            "Apply changes …")
               ),

        column(6, align = "right",
               downloadButton("run_genome_download", width = "100%",
                              "Look for targets …")
        )
  )
                             
)

panel.genome_sel <- tabPanel(
  
  title = "… Loci to Engineer",

  tags$div(title = "The loci of these features (as appearing in the host genome annotation) will be screened for candidate CRISPR targets.",                         
                             
           em("Which GFF features to engineer?"),
           helpText(" Beware of discontinuous feature annotations 
           in the GFF. Not all combinations might be meaningful."),
  
           selectizeInput("inp_feature_type", label = NULL, 
                          multiple = TRUE,
                          choices = character(0),
                          options = list(placeholder = "No feature type selected")
                          )
           
  ),
  
  em("Where to engineer the GFF features?"),
  
  radioButtons("inp_feature_site", label = NULL,
               choices  = list("At feature START (e.g. N-terminal tagging)" = TRUE, 
                               "At feature STOP (e.g. C-terminal tagging)"  = FALSE), 
               selected = def_feature_site
  ),
  
  uiOutput("inp_feature_dist_text"),
  
  sliderInput("inp_feature_dist", label = NULL,
              round = TRUE, post = " nt",
              min = 0, max = 21, step = 1, value = def_feature_dist), 

  em("Maximal cleavage distance from site to engineer"),
  helpText("Allow a potential DNA double-strand break to be introduced no 
            farther away from the integration site than …"),
  
  sliderInput("inp_feature_look", label = NULL,
              round = TRUE, post = " nt", pre = "± ",
              min = 0, max = max_feature_look, step =  1, value = def_feature_look)
         
)

panel.genome_endo <- tabPanel(
  
  title = "… CRISPR Endonuclease",
  
  fluidRow(
  column(width = 6,
   
  selectizeInput("inp_endonuclease_pam", 
                 label = em("Select or add a (generic) PAM"), 
                 multiple = TRUE,
                 options = list(create = TRUE),
                 choices = c("", sort(unique(endonucleases$pattern)))
  )),   
  
  column(width = 6,
                              
  selectizeInput("inp_endonuclease", 
                 label = em("Or: RNA-guided endonuclease"),
                 options = list(create = TRUE),
                 choices = c("", sort(unique(endonucleases$name)), "custom"), 
                 select  = def_endonuclease
  ))
  
  ),
  
  em("Where is the protospacer located?"),
  
  radioButtons("inp_target_is_upstream", label = NULL,
               choices  = list("5'-[PAM]target-3'" = FALSE, 
                               "5'-target[PAM]-3'"= TRUE), 
               selected = def_target_is_upstream
  ),
  
  em("Target sequence structure"),
  helpText("The genomic target totals … without the PAM."),
  
  sliderInput("inp_target_totlen", label = NULL, 
              round = TRUE, post = " nt",
              min = 5, max = 25, step = 1, value = def_target_totlen),
  
  helpText("The first … (seed region) do not tolerate any mismatch between the
           crRNA and the genomic target."),
  uiOutput("inp_target_corlen"),
  
  uiOutput("inp_cleave_at_text"),
  uiOutput("inp_cleave_at")
  
)

panel.genome_misc <- tabPanel(
  
  title = "Miscellaneous",
  
  tags$div(title = "Perform a genome-wide screening for potential off-targets also out of a PAM context.",
                                 
           em("Count potential off-target sites"),
           helpText("Allow … substitutions outside the crRNA\'s seed region for an off-target. 
                     Note, this can be computational intense."),
  
           selectizeInput("inp_mismatches", label = NULL,
                          multiple = TRUE,
                          options = list(create = TRUE, 
                                         placeholder = "No off-targets to be counted"),
                          choices = c("", 1, 3)
           )
  )
  
)

# COMPLETE UI ==================================================================

CASTLING_ui <- fixedPage(title = NULL,
  
  # this allows to handle additional JavaScript functions
  tags$head(

    tags$script(
      HTML("
      Shiny.addCustomMessageHandler('jsCode',
        function(message) {
          console.log(message)
          eval(message.code);
        }
      );
    ")),
    
    tags$style(
      HTML("
      table, th, td {
        border: 1px solid lightgray;
        border-collapse: collapse;
        border-spacing: 50px;
        padding: 6px;
      };
      ul {
        list-style-type: square;
        padding-left: 14pt;
      };
      "))
  ),
  
  tags$style(type = "text/css",
    "label { font-weight: normal; }"
  ),
                   
  navbarPage(
    
    title = "CASTLING castR",
    
    tabPanel(
      
      title = "Oligo Pool Design",
      
      # ---- HOST GENOME PROCESSING --------------------------------------------
      
      fluidPage(           
        
        column(12,
        
        strong("STEP 1 — Genome of Host Organism"),
        
        helpText(
          "Annotated genome assemblies are publicly available from", 
          a("NCBI", href="https://www.ncbi.nlm.nih.gov/genome/?term=txid4932[orgn]",
            target="_blank"),
          "and other repositories. Please note, the size for file upload is limited to",
          ifelse(is.null(getOption("shiny.maxRequestSize")), 5, getOption("shiny.maxRequestSize") / (1024 ^ 2)),
          "MB on this server."
  
        ),   
        
        panel.genome_ref
        
      ),
      
      column(12, hr()),
      
      column(12,
        
        strong("STEP 2 — Select Genomic Loci to Engineer and Properties of Candidate CRISPR Targets"),
        
        helpText(" "),
        

        navlistPanel(
          
          panel.genome_sum,
          
          panel.genome_sel,
          
          panel.genome_endo,
          
          panel.genome_misc
          
        )
        
      ),
      
      column(12, hr()),
      
      column(12,
        
        strong("STEP 3 — Select Genomic Loci to Engineer and Properties of Candidate CRISPR Targets"),
        
        helpText(" "),
        
        panel.oligo_str,
      
        panel.oligo_rem),
      
      column(7, align = "left",
               actionButton("run_oligos", 
                            "Apply changes …")
      ),

      column(5, align = "right",
             downloadButton("run_oligos_download", width = "100%", align = "right",
                            "Create oligo sequences …")
             
      ),

      hr()
                 
    )),
    
    tabPanel(title = "NGS Utilities",
             
             panel.vectorette,
             
             fluidPage(fluidRow(
                 
               column(7,
                      actionButton("run_vectorette", 
                                   "Apply changes …"),
               # ),
               # 
               # column(4, 
                      downloadButton("run_vectorette_download", 
                                     "Download reference …")
               ),

               column(5
                      # helpText("This may take several minutes.")
               )
                 
             ))

             ),
    
    tabPanel(title = "Settings",
             
             panel.settings
             
             ),
    
    # verbatimTextOutput("placeholder", placeholder = TRUE),
    
    hr()
    
  )
              
)
