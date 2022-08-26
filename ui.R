#USER INTERFACE
#NB:File -> Reopen with encoding UTF-8
#NB2: peut ne pas marcher en UTF8 en local mais marcher online
#
#

# Obtention des librairies et mise en place de la fonction pour creer le site ----

library(shiny)
library(shinyjs)
library(bslib)
library(httr)   

shinyUI(
  fluidPage(
    theme = bs_theme(base_font = "Helvetica"),
    useShinyjs(),
    # Mise en place de l'interface du site ----
    titlePanel(
      title = div(img(src="methylation.jpg",width="120",height="60"),"MemiR: A tool for miRNA methylation analysis"),
      windowTitle  = "MemiR"
    ),
    singleton(tags$head(
      tags$link(rel = "stylesheet", type = "text/css", href = "css_perso.css")
    )),
    singleton(tags$head(
      tags$script(src = "javascript_perso.js")
    )),
    # theme = "flatly.min.css",
    
    navbarPage(
      title = NULL,
      id = "tabs",
      header = NULL,
      footer = NULL,
      fluid = T,
      # Creation de la page d'accueil ----
      tabPanel("Home",
               tags$hr(),
               tags$div(HTML('<div id="container-home" style="height:100%;
    margin:0;
    padding:0;
    text-align:center;">
    
    <h4 id="text-home" style="width:49%;
    display:inline-block;
    text-align:center;
    min-width:450px;
    padding-left:10px;
    vertical-align:middle;
    -moz-box-sizing:border-box;
    box-sizing:border-box;">
    Welcome to MemiR<br><br>
    <span style="font-weight: bold;">Me</span><span style="font-weight: bold; color:orange;">miR</span>
    is a platform of bioinformatics tools aimed at predicting the presence of different types of methylations affecting miRNAs (<span style="font-weight: bold;">Me</span><span style="font-weight: bold; color:red;">miR</span><span style="font-weight: bold;">Pred</span>),
the impact that these methylations could have on the structure and exosomal export of miRNAs (<span style="font-weight: bold;">Me</span><span style="font-weight: bold; color:green;">miR</span><span style="font-weight: bold;">Struct</span> and
<span style="font-weight: bold;">Me</span><span style="font-weight: bold; color:maroon;">miR</span><span style="font-weight: bold;">Exo</span>)
but also on the presence of therapeutic options (<span style="font-weight: bold; color:blue;">miR</span><span style="font-weight: bold;">4TO</span>) and/or resistance to anti-cancer therapies 
(<span style="font-weight: bold; color:violet;">miR</span><span style="font-weight: bold;">4Resist</span>).
</h4>
</div>')),
      ),
      # Creation onglet MemiRPred ----
      navbarMenu(title = "MemiRPred", 
#                  # Home ----
#                  tabPanel("Home",
#                           tags$hr(),
#                           tags$div(HTML('<div id="container-home" style="height:100%;
#     margin:0;
#     padding:0;
#     text-align:center;">
#     
# 
#     <img src="image_explication.jpg" id="img-home" style="padding:0;
#     margin:0;
#     width:49%;
#     min-width:450px;
#     max-width:650px;   
#     vertical-align:top;"/>
#     
#     <h4 id="text-home" style="width:49%;
#     display:inline-block;
#     text-align:left;
#     min-width:450px;
#     padding-left:10px;
#     vertical-align:middle;
#     -moz-box-sizing:border-box;
#     box-sizing:border-box;">
#     <span style="font-weight: bold;">miR4TO</span>  is an updated web tool.<br>
# The evidence for methylated miRNA is based on the selection of publications available on Pubmed.<br>
# The link between miR and targeted protein is based on the use of the miRTarBase website.<br>
# The target treatments list included in our tool is constantly read in the daytime according to molecules used in clinical trials and discussions with expert oncologists.
# </h4>
# </div>')),
#                  ),
                 # Prediction of methylation(1) ----
                 tabPanel(
                   "Prediction of methylation for a single miR",
                   tags$hr(),
                   
                   fluidPage(
                     h5("MemiRPred predicts the presence of adenosine (m6), cytosine (m5) and guanosine (m7G) methylation 
                        of miRNAs based on the calculation of a methylation score. This calculation is based on the addition 
                        of two coefficients reflecting two parameters: the experimental evidence of the considered methylation  
                        (experimental presence coefficient, EPC) and the search for homology of presence of methylation motifs or 
                        methylmotifs (methylmotif coefficient, MMC). The methylmotifs considered by MemiRPred are CG, DRACH, RAGGU and UUBD.
                        Both parameters can be weighted in a range of 0 to 10. By default, a weighting of 1:1 is applied."),
                     sidebarLayout(
                       # Sidebar with a slider input
                       sidebarPanel(
                         h3("Select the type of miRNA"),
                         radioButtons("long", label = "", choices = c("mature", "stem loop")),
                         h3("Enter the miRNA name"),
                         selectizeInput(inputId = 'rna', selected = '', multiple = FALSE, 
                                        choices = NULL, label = "", 
                                        options = list(maxOptions = 5, placeholder = "miR name")),
                         h3("Sequence"),
                         uiOutput("rna.correspondant"),
                         h3("Weigthing of methylation score calculation"),
                         h5("The prediction score can be weighted in function to importance given to the experimental evidence
                 and the presence of methylated consensus sequences"),
                         sliderInput("poids_motif", "MMC", 0, 10, 1),
                         sliderInput("poids_liste", "EPC", 0, 10, 1)
                       ),
                       
                       mainPanel(
                         plotOutput("plot",width = "80%",height = "350px"),
                         br(),
                         uiOutput("download_plot")
                       )
                     ),
                     tags$hr(),
                     h3("Experimental evidence"),
                     
                     dataTableOutput("rna_listes"),
                     tags$a(
                       href = "https://memir.shinyapps.io/MemiR_sheets/",
                       "Modify (admin only)",
                       target = "_blank"
                     ),
                     
                     tags$hr(),
                     h3("Predictive evidence based on methyl motifs presence"),
                     dataTableOutput("rna_res")
                   ),
                 ),
                 # Prediction of methylation(2) ----
                 tabPanel(
                   "Prediction of methylation for a list of miR",
                   tags$hr(),
                   fluidPage(
                     h5("MemiRPred predicts the presence of adenosine (m6), cytosine (m5) and guanosine (m7G) methylation 
                        of miRNAs based on the calculation of a methylation score. This calculation is based on the addition 
                        of two coefficients reflecting two parameters: the experimental evidence of the considered methylation 
                        (experimental presence coefficient, EPC) and the search for homology of presence of methylation motifs or 
                        methylmotifs (methylmotif coefficient, MMC). The methylmotifs considered by MemiRPred are CG, DRACH, RAGGU and UUBD.
                        Both parameters can be weighted in a range of 0 to 10. By default, a weighting of 1:1 is applied."),
                     sidebarLayout(
                       # Sidebar with a slider input
                       sidebarPanel(
                         h3("List of miRNAs "),
                         h5("One miRNA per line"),
                         textAreaInput(inputId="liste_miR",value ='',
                                       label = "",
                                       resize = 'both',
                                       cols = 80,
                                       rows = 10),
                         h3("Weigthing of methylation score calculation"),
                         h5("The prediction score can be weighted in function to importance given to the experimental evidence
                 and the presence of methylated consensus sequences"),
                         sliderInput("poids_motif_2", "MMC", 0, 10, 1),
                         sliderInput("poids_liste_2", "EPC", 0, 10, 1)
                       ),
                       
                       mainPanel(
                         h3("Predictive evidence based on methyl motifs presence"),
                         h5("Please, don't use any parenthesis and reset all the input"),
                         dataTableOutput("liste_miR_methy"),
                         uiOutput("download_liste_methy")
                       )
                     )
                   )
                 ),
                 # # Recup pattern-----
                 # tabPanel(
                 #   "Methylation list recovery (pattern only)",
                 #   fluidPage(
                 #     sidebarLayout(
                 #       sidebarPanel(
                 #         radioButtons("long2", label = "Type of miRNA", choices = c("mature", "stem loop")),
                 #         h3("Input your methylation"),
                 #         selectInput("methyl", "Methylation", choices = NULL),
                 #         selectInput(inputId = 'seuil_score_liste', 
                 #                     label = "Target score threshold for methylated miR", 
                 #                     choices = NULL),
                 #         actionButton("go", "Submit", class = "btn-success"),
                 #       ),
                 #       mainPanel(
                 #         uiOutput("download_miR_liste"),
                 #         dataTableOutput("listes"),
                 #       )
                 #     )
                 #   )
                 # ),
                 # Recup Score ----
                 tabPanel(
                   "Methylation list recovery",
                   fluidPage(
                     h5("MemiRPred predicts the presence of adenosine (m6), cytosine (m5) and guanosine (m7G) methylation 
                        of miRNAs based on the calculation of a methylation score. This calculation is based on the addition 
                        of two coefficients reflecting two parameters: the experimental evidence of the considered methylation 
                        (experimental presence coefficient, EPC) and the search for homology of presence of methylation motifs or 
                        methylmotifs (methylmotif coefficient, MMC). The methylmotifs considered by MemiRPred are CG, DRACH, RAGGU and UUBD.
                        Both parameters can be weighted in a range of 0 to 10. By default, a weighting of 1:1 is applied."),
                     sidebarLayout(
                       sidebarPanel(
                         h3("Select the type of miRNA"),
                         radioButtons("long3", label = "", choices = c("mature", "stem loop")),
                         selectizeInput(inputId = 'methyl2', selected = '', 
                                        multiple = FALSE, choices = NULL, 
                                        label = "Select the methylation type"),
                         h3("Weigthing of methylation score calculation"),
                         h5("Set a threshold for miRNA methylation (according to MemiRPred)"),
                         numericInput("seuil_score2", "Target score threshold for methylated miR", value = 2, min = 1, max = 20),
                         h5("The prediction score can be weighted in function to importance given to the experimental evidence
                 and the presence of methylated consensus sequences"),
                         sliderInput("poids_motif2", "MMC", 0, 10, 1),
                         sliderInput("poids_liste2", "EPC", 0, 10, 1)
                       ),
                       mainPanel(
                         uiOutput("download_miR_liste2"),
                         dataTableOutput("listes2"),
                       )
                     )
                   )
                 )
                 ),
      # Creation onglet MemiRStruct ---- 
      navbarMenu("MemiRStruct",
                 # Boucle ----
                 tabPanel(
                   "Impact of methylation on the secondary structure",
                   fluidPage(
                     h5("MemiRStruct allows the visualization of loops and homodulexes that can form within miRNAs while considering 
                        the presence of different methylations that can take place in miRNAs"),
                     sidebarLayout(
                       sidebarPanel(
                         h3("Select the type of miRNA"),
                         radioButtons("long4", label = "", choices = c("mature", "stem loop")),
                         h3("Enter miRNA name"),
                         selectizeInput("rna_ss", selected = '', 
                                        multiple = FALSE, choices = NULL, 
                                        label = "",options = list(maxOptions = 10)),
                         h3("Select the type of methylation"),
                         selectizeInput("choix_methy", selected = "", 
                                        multiple = FALSE, choices = NULL, 
                                        label = "The different types of methylation proposed are based on the use of MemiRPred",
                                        options = list(maxOptions = 10)),
                         actionButton("go2", "Submit", class = "btn-success")
                       ),
                       mainPanel(
                         dataTableOutput("met"),
                         plotOutput("plot_boucle",width = "700px",height = "450px"),
                         # dataTableOutput("boucle"),
                         br(),
                         br(),
                         textOutput("Warning"),
                         actionButton("go3", "De-twisted", class = "btn-success"),
                         uiOutput("download_impact")
                         # br(),
                         # h3("Impact of the methylation"),
                         # selectInput(inputId = 'seuil_score_loop', 
                         #             label = "Target score threshold for methylated miR", 
                         #             choices = c('25%', '50%', "75%", "100%")),
                         # dataTableOutput("Impact")
                       )
                     )
                   )
                 )),
      # Creation onglet MemiRExo ----
      navbarMenu("MemiRExo",
                 # Research of exomotif and overlapping----
                 tabPanel(
                   "Research of exomotif",
                   tags$hr(),
                   fluidPage(
                     h5("MemiRExo predicts the existence of potential overlap 
                        between the presence of exomotif and methylmotif in miRNAs."),
                     h5("The exomotif term here refers to a nucleotide motif existing in miRNA that controls their loading into exosomes/extracellular vesicles. 
                        The list of exomotifs considered here is taken from the literature."),
                     sidebarLayout(position="left",
                                   sidebarPanel(
                                     h3("Select the type of miRNA"),
                                     radioButtons("long1", label = "", choices = c("mature", "stem loop")),
                                     h3("Enter miRNA name"),
                                     selectizeInput(inputId = 'rna_exo', selected = '', multiple = FALSE, 
                                                    choices = NULL,label = '' ,
                                                    options = list(maxOptions = 5, placeholder = "miR name")),
                                   ),
                                   mainPanel(
                                     
                                   )),
                     h3("Exomotif"),
                     dataTableOutput ("rna_exomes"),
                     uiOutput("download_exomes"),
                     tags$a(
                       href = "https://memir.shinyapps.io/MemiR_sheets/",
                       "Modify (admin only)",
                       target = "_blank"
                     )
                   )
                 )),
      # Creation de l'onglet miR4TO ----
      navbarMenu("miR4TO",
                 # MT ----
                 tabPanel(
                   "Single miRNA analyses",
                   tags$hr(),
                   fluidPage(
                     tags$div(HTML('<img src="image_explication.jpg" id="img-home" style="padding:0;
                      display:block;
                      width:35%;
                      min-width:250px;
                      max-width:650px;
                      vertical-align:top;
                      margin-left: auto;
                      margin-right: auto;"/>')),
                     h5('Via its integrated miRNA analysis mode, miR4TO predicts the presence or loss of a therapeutic option 
                        using the results of integrated expression and methylation analysis of miRNAs. For this purpose, 
                        miR4TO integrates the triple hypothesis that a) the under-expression or hypermethylation of a miRNA 
                        can induce the increase of expression of a protein that can be the target of a therapy (gain of therapeutic option), 
                        b) under-expression or hypomethylation can induce the loss of expression of a protein that could be the target of a therapy (loss of therapeutic option),
                        and c) these 2 phenomena can "cancel each other out" when an under-expressed or hypermethylated miRNA and an over-expressed or hypomethylated miRNA 
                        have the same target. By integrating these 3 hypotheses, using the miRTarBase (for the identification of miRNA/target interactions) and our targeted therapy database 
                        (continuously updated), the integrated miRNA analysis mode of miR4TO delivers triads associating a miRNA, a target and a targeted therapy molecule associated to 
                        a context of gain or loss of therapeutic option.'),
                     sidebarLayout(
                       sidebarPanel(
                         h3("Enter miRNA name"),
                         selectizeInput(inputId = 'rna_mtt', selected = '', multiple = FALSE, 
                                        choices = NULL, label = "", 
                                        options = list(maxOptions = 5, placeholder = "miR name")),
                         h3("Select the amount of evidence validating the miRNA/target interaction to be used to predict the miRNA/target/treatment option triads."),
                         h6("Selecting the amount of evidence validating the miRNA/target interaction to use to predict miRNA/target/treatment option triads increases 
                            the robustness of miRNA/target/treatment option triad predictions."),
                         numericInput('nombre_evidence_mini', label = "Target evidence threshold for triads", value = 4, min = 1, max = 20)
                       ),
                       mainPanel(
                         tags$hr(),
                         h2("miR4TO: miR-Target-Treatment triads"),
                         h3("Weight the relative importance of experiments validating the miRNA/target interaction."),
                         h6("Modifying the weighting coefficients allows to modulate the relative importance of the experiments validating the 
                            miRNA/target interaction in accordance with the miRTarBase site (strong evidence: Reporter assay, Western blot and qPCR; 
                            weaker evidence: Microassay, NGS, pSILAC, CLIP-seq and other methods)"),
                         sliderInput("poids_fort", "Weighting coefficient for strong evidence", 0, 10, 1),
                         sliderInput("poids_faible", "Weighting coefficient for less strong evidence", 0, 10, 1),
                         uiOutput("download_triade"),
                         dataTableOutput("triade"),
                         tags$a(
                           href = "https://memir.shinyapps.io/MemiR_sheets/",
                           "Modify (admin only)",
                           target = "_blank"
                         )
                       )
                     )
                   )
                 ), 
                 # Activation/Inhibition ----
                 tabPanel(
                   title = 'Integrated miRNA analyses',
                   tags$hr(),
                   fluidPage(
                     tags$div(HTML('<img src="image_explication.jpg" id="img-home" style="padding:0;
                      display:block;
                      width:35%;
                      min-width:250px;
                      max-width:650px;
                      vertical-align:top;
                      margin-left: auto;
                      margin-right: auto;"/>')),
                     h5('Via its integrated miRNA analysis mode, miR4TO predicts the presence or loss of a therapeutic option 
                        using the results of integrated expression and methylation analysis of miRNAs. For this purpose, 
                        miR4TO integrates the triple hypothesis that a) the under-expression or hypermethylation of a miRNA 
                        can induce the increase of expression of a protein that can be the target of a therapy (gain of therapeutic option), 
                        b) under-expression or hypomethylation can induce the loss of expression of a protein that could be the target of a therapy (loss of therapeutic option),
                        and c) these 2 phenomena can "cancel each other out" when an under-expressed or hypermethylated miRNA and an over-expressed or hypomethylated miRNA 
                        have the same target. By integrating these 3 hypotheses, using the miRTarBase (for the identification of miRNA/target interactions) and our targeted therapy database 
                        (continuously updated), the integrated miRNA analysis mode of miR4TO delivers triads associating a miRNA, a target and a targeted therapy molecule associated to 
                        a context of gain or loss of therapeutic option.'),
                     sidebarLayout(
                       sidebarPanel(
                         h3("Enter the list of hypomethylated or overexpressed miRNA"),
                         h6("One miRNA per line"),
                         textAreaInput(inputId = "up_regul",value ='',
                                       label = "",
                                       resize = 'both',
                                       cols = 80,
                                       rows = 10),
                         h3("Enter the list of hypermethylated or downexpressed miRNA"),
                         h6("One miRNA per line"),
                         textAreaInput(inputId = "down_regul",value ='',
                                       label = "",
                                       resize = 'both',
                                       cols = 80,
                                       rows = 10),
                         h3("Weight the relative importance of experiments validating the miRNA/target interaction."),
                         h6("Modifying the weighting coefficients allows to modulate the relative importance of the experiments 
                                        validating the miRNA/target interaction in accordance with the miRTarBase site (strong evidence: 
                                        Reporter assay, Western blot and qPCR; weaker evidence: Microassay, NGS, pSILAC, CLIP-seq and other methods)"),
                         sliderInput("poids_fort_ai", "Weighting coefficient for strong evidence", 0, 10, 1),
                         sliderInput("poids_faible_ai", "Weighting coefficient for less strong evidence", 0, 10, 1),
                         h3("Select the amount of evidence validating the miRNA/target interaction to be used to predict the miRNA/target/treatment option triads."),
                         h6("Selecting the amount of evidence validating the miRNA/target interaction to use to predict miRNA/target/treatment option triads increases the 
                                        robustness of miRNA/target/treatment option triad predictions."),
                         numericInput('seuil_preuves_acti_repr', "Target evidence threshold for activated proteines", value = 4, min = 1, max = 20)
                       ),
                       
                       mainPanel(
                         h3("Gain of therapeutic option"),
                         dataTableOutput("acti"),
                         uiOutput("download_acti"),
                         tags$hr(),
                         h3("Loss of therapeutic option"),
                         dataTableOutput("repr"),
                         uiOutput("download_repr")
                       )
                     ),
                   ))),
      # Creation onglet miR4Resist ----
      navbarMenu("miR4Resist",
                 # Analysis : Treatment to miR----
                 tabPanel(
                   "Analysis : Treatment to miR",
                   tags$hr(),
                   
                   fluidPage(
                     h5("miR4Resist predicts miRNAs and methylations of these miRNAs that may 
                        cause resistance to a targeted therapy of interest. For this, it is possible 
                        to weight the relative importance of experiments validating the miRNA/target 
                        interaction, the amount of evidence validating the miRNA/target interaction 
                        to be used to predict miRNA/target/treatment option triads, and to set a 
                        threshold for miRNA methylation (according to MemiRPred)."),
                     sidebarLayout(position="left",
                                   # Sidebar with a slider input
                                   sidebarPanel(
                                     h3("Select a therapeutic molecule"),
                                     selectizeInput(inputId = 'therapie', selected = '', 
                                                    multiple = FALSE, choices = NULL, 
                                                    label = "", 
                                                    options = list(maxOptions = 5, placeholder = "Treatment")),
                                     h3("Weight the relative importance of experiments validating the miRNA/target interaction."),
                                     h6("Modifying the weighting coefficients allows to modulate the relative importance of the experiments 
                                        validating the miRNA/target interaction in accordance with the miRTarBase site (strong evidence: 
                                        Reporter assay, Western blot and qPCR; weaker evidence: Microassay, NGS, pSILAC, CLIP-seq and other methods)"),
                                     sliderInput("poids_fort_ttm", "Weighting coefficient for strong evidence", 0, 10, 1),
                                     sliderInput("poids_faible_ttm", "Weighting coefficient for less strong evidence", 0, 10, 1),
                                     h3("Select the amount of evidence validating the miRNA/target interaction to be used to predict the miRNA/target/treatment option triads."),
                                     h6("Selecting the amount of evidence validating the miRNA/target interaction to use to predict miRNA/target/treatment option triads increases the 
                                        robustness of miRNA/target/treatment option triad predictions."),
                                     numericInput('seuil_preuves', "Target evidence threshold for triads", value = 4, min = 1, max = 20),
                                     h3("Set a threshold for miRNA methylation (according to MemiRPred)"),
                                     selectizeInput(inputId = 'seuil_score', selected = '25%', 
                                                    label = "Target score threshold for methylated miR", 
                                                    options = list(maxOptions = 4), 
                                                    multiple = FALSE, choices = c('25%', '50%', "75%", "100%"))),
                                   mainPanel(
                                     dataTableOutput("mir"),
                                     uiOutput("download_treat"),
                                   )
                     ),
                   ))),
      # # Creation onglet miRFunction ----
      # navbarMenu('miRFunction',
      #            # Tam_miR ----
      #            tabPanel(
      #              "miR to function",
      #              fluidPage(
      #                sidebarLayout(
      #                  sidebarPanel(
      #                    h3("Input your miR(Works only for the Homo Sapiens"),
      #                    selectizeInput(inputId = 'rna_tam', selected = '', multiple = FALSE, 
      #                                   choices = NULL, label = "Input the miR:", 
      #                                   options = list(maxOptions = 5, placeholder = "miR name")),
      #                  ),
      #                  mainPanel(
      #                    h3("Function"),
      #                    dataTableOutput("Tam_fon")
      #                  )
      #                )
      #              )
      #            ),
      #            # Tam_autre ----
      #            tabPanel(
      #              "Function to miR",
      #              fluidPage(
      #                sidebarLayout(
      #                  sidebarPanel(
      #                    h3("Input your Title"),
      #                    selectizeInput("Tit", selected = '', 
      #                                   multiple = FALSE, choices = NULL, 
      #                                   label = "Title",options = list(maxOptions = 10)),
      #                    h3("Input your category"),
      #                    selectizeInput(
      #                      inputId = "Cate", 
      #                      label = "Category", 
      #                      choices = NULL, 
      #                      selected = "",
      #                      multiple = FALSE
      #                    )
      #                  ),
      #                  mainPanel(
      #                    h3("miR"),
      #                    dataTableOutput("mir_tam")
      #                  )
      #                )
      #              )
      #            )),
      # 
      # Creation onglet About ----
      tabPanel("About",
               tags$hr(),
               tags$div(HTML('<div id="container-home" style="height:100%;
    margin:0;
    padding:0;
    text-align:center;">
    
    <h4 id="text-home" style="width:49%;
    display:inline-block;
    text-align:center;
    min-width:450px;
    padding-left:10px;
    vertical-align:middle;
    -moz-box-sizing:border-box;
    box-sizing:border-box;">
    <span style="font-weight: bold;">CONTACT :</span> <br><br>
    <span style="font-weight: bold;">MemiR</span>  platform team is composed by:<br><br>

Dr Pierre-François CARTRON<br>
Quentin ARBERET<br>
Axel PREVOT<br><br>
Any questions about the MemiR platform should be directed at MemiRTools@gmail.com<br>
MemiR platform use is free only for the academic public. For commercial usage, please contact MemiR platform at MemiRTools@gmail.com<br>
The different tools proposed in the MemiR platform echo the research work of our team.<br>
</h4>
</div>')),
               h3("The visit number :"),
               uiOutput("logs_table"),
               h3("References : "),
               dataTableOutput("Ref")
               
               
      ),
    )
  ))
