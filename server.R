# Obtention des listes et des fichiers

library(shiny)
library(fmsb)
library(data.table)

#https://shiny.rstudio.com/tutorial/written-tutorial/lesson5/
source("util.R")

shinyServer(function(input, output, session) {
  
  ############## obtention des listes et creation des choix ----
  # Obtention listes des pattern de methylation ----
  
  data_pattern <- list()
  nom_pattern <<- sheet_names(ss = sheet_id_p)
  for (name in nom_pattern) {
    data <- data.frame(read_sheet(ss = sheet_id_p,
                                  sheet = name))
    # names(data) <- name
    data_pattern <- append(data_pattern, data)
  }
  data_pattern <<- data_pattern
  
  # Obtention listes de presence de methylation ----
  
  data_liste <- c()
  nom_listes <<- sheet_names(ss = sheet_id)
  for (name in nom_listes) {
    data <- data.frame(read_sheet(ss = sheet_id,
                                  sheet = name))
    names(data)[1] <- name
    data_liste <- append(data_liste, data[1])
  }
  data_liste <<- data_liste
  
  # Obtention des therapies et de leur cibles ----
  target_therapie <<- data.frame(read_sheet(ss = sheet_id_2,
                                            sheet = "target_therapie"))
  # Obtention des exomotif ----
  Liste_Exomotif <<- data.frame(read_sheet(ss = sheet_id_3,
                                           sheet = "Exomotif"))
  # Obtention des Reférences ----
  Liste_ref <<- data.frame(read_sheet(ss = sheet_id_4,
                                           sheet = "References"))
  # Mis a jour des choix d'input ----
  # Input Mir----
  
  observeEvent(input$long, {
    updateSelectizeInput(
      session,
      'rna',
      selected = '' ,
      choices = Base_exploit_pred()$Nom,
      server = TRUE
    )
  })
  
  updateSelectizeInput(
    session,
    'rna_mtt',
    selected = '' ,
    choices = ExtracNom,
    server = TRUE
  )
  
  observeEvent(input$long1, {
    updateSelectizeInput(
      session,
      'rna_exo',
      selected = '' ,
      choices = Base_exploit_exo()$Nom,
      server = TRUE
    )
  })
  
  updateSelectizeInput(
    session,
    'rna_tam',
    selected = '' ,
    choices = Base_pour_Tam$Nom,
    server = TRUE
  )
  
  observeEvent(input$long4, {
    updateSelectizeInput(
      session,
      'rna_ss',
      selected = '' ,
      choices = Base_exploit_ss()$Nom,
      server = TRUE
    )
  })
  
  # Methyl -----
  # observe({
  #   if (input$methyl == "m5C"){
  #     updateSelectInput(session, "seuil_score_liste", choices = c("100%"), selected = c("100%"))
  #   }
  #   if(input$methyl == "m7G"){
  #     updateSelectInput(session, "seuil_score_liste", choices = c("50%", "100%"), selected = c('50%'))
  #   }
  #   if(input$methyl == "m6A"){
  #     updateSelectInput(session, "seuil_score_liste", choices = c('25%', '50%', "75%", "100%"), selected = c('25%'))
  #   }
  # })
  updateSelectInput(session, "methyl", choices = nom_pattern, selected = "")
  
  observeEvent(input$methyl,{
    choix <- c()
    ind <- grep(input$methyl, names(data_pattern))
    for (i in ind){
      choix <- append(choix, paste0(fct_lire(names(data_pattern)[i], "_", T)[2], "%"))
    }
    updateSelectInput(session, "seuil_score_liste", choices = choix, selected = "")
  }
               )
  
  updateSelectInput(session, "methyl2", choices = nom_pattern, selected = "")
  
  # Here, update allows to select nothing during the first use
  # updateSelectizeInput(
  #   session,
  #   'methyl2',
  #   selected = '' ,
  #   choices = c("m5C", "m6A", "m7G"),
  #   server = TRUE
  # )
  # Autre ----
  updateSelectizeInput(
    session,
    'therapie',
    selected = '' ,
    choices = target_therapie$Therapie,
    server = TRUE
  )
  updateSelectizeInput(
    session,
    'Tit',
    selected = '' ,
    choices = ExtracTam,
    server = TRUE
  )
  observeEvent(input$rna_ss, {
    updateSelectizeInput(
      session,
      'choix_methy',
      selected = '' ,
      choices = Liste_up(),
      server = TRUE
    )
  })
  ############## outputs_button_dowloading ----
  # Graphique ----
  output$download_plot <- renderUI({
    req(input$rna)
    downloadButton(outputId = "down_plot",
                   label = "",
                   style = "float:right; margin-right:20%;")
  })
  # Liste_pred ----
  output$download_liste_methy<- renderUI({
    req(input$liste_miR)
    downloadButton(outputId = "down_liste_methy",
                   label = "",
                   style = "float:right; margin-right:20%;")
  })
  # Exomes -----
  output$download_exomes<- renderUI({
    req(input$rna_exo)
    downloadButton(outputId = "down_exomes",
                   label = "",
                   style = "float:right; margin-right:20%;")
  })
  # Triades ----
  output$download_triade <- renderUI({
    req(input$rna_mtt)
    downloadButton(outputId = "down_triade",
                   label = "",
                   style = "float:right;")
  })
  # Acti/Repr ----
  output$download_acti <- renderUI({
    req(input$down_regul)
    req(input$up_regul)
    downloadButton(outputId = "down_acti",
                   label = "",
                   style = "float:right;")
  })
  output$download_repr <- renderUI({
    req(input$down_regul)
    req(input$up_regul)
    downloadButton(outputId = "down_repr",
                   label = "",
                   style = "float:right;")
  })
  # Treatment ----
  output$download_treat <- renderUI({
    req(input$therapie)
    downloadButton(outputId = "down_treat",
                   label = "",
                   style = "float:right;")
  })
  # Liste pattern ----
  output$download_mir_liste <- renderUI({
    req(input$go)
    downloadButton(outputId = "down_mir_liste",
                   label = "",
                   style = "float:right;")
  })
  # Liste Score ----
  output$download_mir_liste2 <- renderUI({
    req(input$methyl2)
    downloadButton(outputId = "down_mir_liste2",
                   label = "",
                   style = "float:right;")
  })
  # Secondary structure -----
  output$download_impact <- renderUI({
    req(input$rna_ss)
    req(input$choix_methy)
    req(input$go2)
    downloadButton(outputId = "down_impact",
                   label = "",
                   style = "float:left;")
  })
  
  ############## Outputs_dowloading ---- 
  # Graphique ----
  # downloadHandler contains 2 arguments as functions, namely filename, content 
  output$down_plot <- downloadHandler(
    filename =  function() {
      paste(paste("MemiRPred - Predictive score - ", input$rna, sep = ""),
            "png",
            sep = ".")
    },
    # content is a function with argument file. content writes the plot to the device
    content = function(file) {
      png(file) # open the png device
      if (input$rna == ""){return("")}
      fct_plot_scores(input$rna,input$poids_motif, input$poids_liste, Base_util)
      dev.off()  # turn the device off
      
    }
  )
  # Liste_pred----
  output$down_liste_methy <- downloadHandler(
    filename =  function() {
      paste("MemiRPred - Prediction of methylation for your list of miRNA",
            "csv",
            sep = ".")
    },
    content = function(file) {
      write.csv(liste_pred(), file)
    }
  )
  # Exomes----
  output$down_exomes <- downloadHandler(
    filename =  function() {
      paste(paste("MemiRExo - Exomotif and overlapping with methylated pattern for the miRNA ", input$rna_exo, sep = " "),
            "csv",
            sep = ".")
    },
    content = function(file) {
      write.csv(exome(), file)
    }
  )

  # Triades ----
  output$down_triade <- downloadHandler(
    filename =  function() {
      paste(paste("miR4TO - Triades mir-Target-Treatment for the mir ", input$rna_mtt, sep = " "),
            "csv",
            sep = ".")
    },
    content = function(file) {
      write.csv(triade(), file)
    }
  )
  # Acti/repr----
  output$down_acti <- downloadHandler(
    filename =  function() {
      paste("miR4TO - Integrated miRNA analyses - Gain of therapeutic option ",
            "csv",
            sep = ".")
    },
    content = function(file) {
      write.csv(acti(), file)
    }
  )
  output$down_repr <- downloadHandler(
    filename =  function() {
      paste("miR4TO - Integrated miRNA analyses - Loss of therapeutic option ",
            "csv",
            sep = ".")
    },
    content = function(file) {
      write.csv(repr(), file)
    }
  )
  # Treatment ----
  output$down_treat <- downloadHandler(
    filename =  function() {
      paste(paste("miR4Resist-List of mirs with susceptibility to be methylated correlated with the treatment",input$therapie, sep = " "),
            "csv",
            sep = ".")
    },
    content = function(file) {
      write.csv(mir(), file)
    }
  )
  # Liste pattern ----
  output$down_mir_liste <- downloadHandler(
    filename =  function() {
      paste(paste("miR4TO-List of mirs with a pattern of", input$methyl, "with a score of",input$seuil_score_liste, sep = " "),
            "csv",
            sep = ".")
    },
    content = function(file) {
      write.csv(liste(), file)
      
    }
  )
  # Liste Score ----
  output$down_mir_liste2 <- downloadHandler(
    filename =  function() {
      paste(paste("miR4TO-List of mirs with a pattern of", input$methyl2,
                  "with a score of", input$seuil_score2,
                  "with a cofficient motif ",input$poids_motif2, 
                  'and a coefficient liste',input$poids_liste2, sep = " "),
            "csv",
            sep = ".")
    },
    content = function(file) {
      write.csv(listes2(), file)
      
    }
  )
  # Secondary structure -----
  output$down_impact <- downloadHandler(
    filename =  function() {
      paste(paste("miR4TO-Structure of",input$rna_ss, "for the methylation", input$choix_methy, sep = " "),
            "png",
            sep = ".")
    },
    content = function(file) {
      png(file, width = 800, height = 480, units = "px")
      if(input$rna_ss == ""){return("")}
      i <- grep(req(input$choix_methy), Liste_up())
      Seq <- fct_lire(req(methy()[ i, "Sequence"]), "[AUCGaucg]", FALSE)
      liste <- fct_calcul(tolower(Seq))
      calcul_boucle <- list()
      if (min(liste[[1]]) < 0){
        Sep <- fct_struc(liste, "Yes")
        calcul_boucle <-list(Seq, Sep)
      }
      else{
        Sep1 <- fct_struc(fct_calcul(tolower(Seq)), "No")
        Sep2 <- fct_struc(fct_calcul(tolower(Seq)), "Yes")
        calcul_boucle <-list(Seq, Sep1, Seq, Sep2)
      }
      if(length(calcul_boucle())==2){
        Seq <- calcul_boucle[[1]]
        Sep <- calcul_boucle[[2]]
        fct_graph(Seq, Sep, req(methy())[i ,c("Position")])
      }
      else{
        par(mar = c(4, 5, 0.5, 0.5), mfrow = c(1, 2))
        Seq <- calcul_boucle[[1]]
        Sep <- calcul_boucle[[2]]
        fct_graph(Seq, Sep, req(methy())[i ,c("Position")])
        Seq <- calcul_boucle[[3]]
        Sep <- calcul_boucle[[4]]
        fct_graph(Seq, Sep, req(methy())[i ,c("Position")])
      }
      dev.off()
    }
  )
  ############## Methylation ----
  
  Base_exploit_pred <- reactive({
    if(input$long=="stem loop"){test <- T}
    else{test <- F}
    fct_Choix_Base("hsa", test)
  })
  
  # Rechercher sequence ----
  
  output$rna.correspondant <- renderUI({
    inFile <- input$rna
    return(toupper(fct_nom_to_seq(req(inFile), Base_exploit_pred())))
  })
  # Creation tableau prediction ----
  output$rna_res <- renderDataTable({
    inFile <- input$rna
    #permet de ne rien retourner si rien n'est entrÃ©
    rna.res <- fct_Matchtotal(req(inFile), Base_exploit_pred())
    if((rna.res[1, 1]== "No methylation")){
      rna.res.final <- rna.res
    }
    else{
      rna.res.final <- subset(rna.res, select = c('Methylation', 'miRNA name', 'Sequence', 'Sequence_mutee', 'Position', 'Score'))
      colnames(rna.res.final)[4] <- "Methylated sequence"
      return(rna.res.final)
    }
  }, options = list(
    paging = FALSE,
    searching = FALSE, 
    info = FALSE)
  )
  # Recherche presence liste ----
  output$rna_listes <- renderDataTable({
    inFile <- input$rna
    #permet de ne rien retourner si rien n'est entrÃ©
    rna.listes <- fct_presence_listes(req(inFile))
    return(rna.listes)
  },
  options = list(
    paging = FALSE,
    searching = FALSE,
    info = FALSE,
    drawCallback = I(
      'function() {var x = document.getElementsByTagName("td");
            for (var i= x.length; i-->0;){
              if (x[i].innerText=="Present"){
                x[i].style.background= "#3BE6C4";
}}}'
    )))
  # Creation graphique ----
  plot_score <- reactive({
    return(fct_plot_scores(req(input$rna), input$poids_motif, input$poids_liste, Base_exploit_pred()))
  })
  
  output$plot <- renderPlot ({
    return(plot_score())
  })
  
  # Liste methy -----
  
  liste_pred <- reactive({
    mirs <- fct_lire(req(input$liste_miR), "\n", TRUE)
    return(fct_liste_methyl(tolower(mirs), input$poids_motif_2, input$poids_liste_2))
  })

  output$liste_miR_methy <- renderDataTable({
    mirs <- fct_lire(req(input$liste_miR), "\n", TRUE)
    return(fct_stats_methy(tolower(mirs), input$poids_motif_2, input$poids_liste_2, 4))
  }, options = list(
    paging = FALSE,
    searching = FALSE,
    info = FALSE,
    drawCallback= I(
      'function() {
            var x = liste_miR_methy.getElementsByTagName("td");
            var y = liste_miR_methy.getElementsByTagName("tr");
            var z = liste_miR_methy.getElementsByTagName("th");
            for (var i= x.length; i-->0;){
              let l  = (x.length/(y.length - 1) -1)/2
              var k = i - l
              if (x[i].innerText=="0%-25%"){
                x[k].style.background= "#ff0000";
                x[i].style.display="none";
              }
              if (x[i].innerText=="25%-50%"){
                x[k].style.background= "#fe8200";
                x[i].style.display="none";
              }
              if (x[i].innerText=="50%-75%"){
                x[k].style.background= "#ffff00";
                x[i].style.display="none";
              }
              if (x[i].innerText=="75%-100%"){
                x[k].style.background= "#00ff00";
                x[i].style.display="none";
              }
            }
            for (var j= z.length; j-->0;){
              if (z[j].innerText=="Stat"){
                z[j].style.display="none";
              }
            }
      }'
    )
  )
  )
  
  ############## Exomotif ----
  
  Base_exploit_exo <- reactive({
    if(input$long1=="stem loop"){test <- T}
    else{test <- F}
    fct_Choix_Base("hsa", test)
  })
  exome <- reactive({
    inFile <- input$rna_exo
    exo <- fct_exomes(req(inFile), Base_exploit_exo())
    return(exo)
  })
  
  output$rna_exomes <- renderDataTable({
    exome()
  },
  options = list(
    paging = TRUE,
    searching = FALSE,
    info = FALSE,
    pageLength = 10)
  )
  
  ############## Target and Triades ----
  # Base -----
  
  Base_targ <- reactive({
    req(input$rna_mtt)
    req(input$nombre_evidence_mini)
    nom <- "Base_targ/hsa_MTI.csv"
    Base_ex <- read.csv(nom, header = TRUE)
    incides_fun <- grep("Non-Functional MTI", Base_ex$Support.Type, invert = TRUE)
    miRNA <- Base_ex$miRNA[incides_fun]
    Target.gene <- Base_ex$Target.Gene[incides_fun]
    experi <- Base_ex$Experiments[incides_fun]
    Clean_target <- data.frame(Mir = miRNA, Target = Target.gene, Experiment = experi)
    return(Clean_target)
  })
  
  # Triades----
  triade <- reactive({
    inFile <- req(input$rna_mtt)
    triade <- fct_triade_m_t_t_g(inFile, fct_Target_General(inFile, Base_targ()),
                          input$poids_fort, input$poids_faible, 
                          req(input$nombre_evidence_mini))
    return(triade)
  })
  
  output$triade <- renderDataTable({
    return(triade())
  },
  options = list(
    paging = TRUE,
    searching = FALSE,
    info = FALSE,
    pageLength = 10
  ))
  ############## Treatment to mir ----
  
  Base_targ_2 <- reactive({
    req(input$therapie)
    nom <- "Base_targ/hsa_MTI.csv"
    Base_ex <- read.csv(nom, header = TRUE)
    incides_fun <- grep("Non-Functional MTI", Base_ex$Support.Type, invert = TRUE)
    miRNA <- Base_ex$miRNA[incides_fun]
    Target.gene <- Base_ex$Target.Gene[incides_fun]
    experi <- Base_ex$Experiments[incides_fun]
    Clean_target <- data.frame(Mir = miRNA, Target = Target.gene, Experiment = experi)
    return(Clean_target)
  })
  
  mir <- reactive({
    inFile <- input$therapie
    mirs <- fct_therapie_to_methyl(req(inFile),input$poids_fort_ttm, input$poids_faible_ttm,
                                   input$seuil_preuves, input$seuil_score, 
                                   Base_targ_2(), Base_util)
    return(mirs)
  })
  
  output$mir <- renderDataTable({
    inFile<- input$therapie
    seuil_preuves <- input$seuil_preuves
    seuil_score <- input$seuil_score
    return(fct_therapie_to_methyl(req(inFile),input$poids_fort_ttm, input$poids_faible_ttm,
                                  seuil_preuves, seuil_score, 
                                  Base_targ_2(), Base_util))
  },
  options = list(
    paging = TRUE,
    searching = FALSE,
    info = FALSE,
    pageLength = 10)
  )
  
  ############## Activation and Inhibition ----
  # Activation ----
  # output$acti <- renderDataTable({
  #   inFile <- input$up_regul
  #   inFile2 <- input$down_regul
  #   up <- fct_lire(req(inFile), "\n", TRUE)
  #   down <- fct_lire(req(inFile2), "\n", TRUE)
  #   return(fct_activation(up, down, input$seuil_preuves_acti, 
  #                         input$poids_fort_ai, input$poids_faible_ai))
  #   
  # },
  # options = list(
  #   paging = TRUE,
  #   searching = FALSE,
  #   info = FALSE,
  #   pageLength = 10)
  # )
  
  acti_repr <- reactive({
  inFile <- input$up_regul
  inFile2 <- input$down_regul
  up <- fct_lire(req(inFile), "\n", TRUE)
  down <- fct_lire(req(inFile2), "\n", TRUE)
  return(fct_acti_repr(up, down))
})
  
  acti <- reactive({
    activa <- acti_repr()[[1]]
    repres <- acti_repr()[[2]]
    # selection ----
    if (is.null(activa[1, 1])){activa_dt <- data.frame()}
    else{
      activa_dt <- activa[as.numeric(activa$`Nb of strong evidence`)*input$poids_fort_ai+
                            as.numeric(activa$`Nb of less strong evidence`)*input$poids_faible_ai 
                          >= input$seuil_preuves_acti_repr,c("Mir","Target")]
    }
    if(is.null(repres[1, 1])){repres_vect <- c()}
    else{
      repres_vect <- repres[as.numeric(repres$`Nb of strong evidence`)*input$poids_fort_ai+
                              as.numeric(repres$`Nb of less strong evidence`)*input$poids_faible_ai 
                            >= input$seuil_preuves_acti_repr,c("Target")]
    }

    activa_unique <- fct_mir_unique(activa_dt)
    
    go <- T
    if (is.null(activa_unique[1, 1])){go <- F}
    else{
      final_dt <- activa_unique[!sapply(activa_unique$Target, is.element, repres_vect), c("Mir", "Target")]
    }
    
    # Therapie ----
    target_vect <- c()
    therapie_vect <- c()
    mir_vect <- c()
    
    if (go){
      for (target in final_dt$Target) {
        target_exact <- paste0('^', target, '$')
        mir <- final_dt$Mir[grep(target_exact, final_dt$Target)]
        indices <- grep(target_exact, target_therapie$Target)
        if (length(indices)==0){
          therapie_vect <- append(therapie_vect, "No therapy available")
          target_vect <-append(target_vect, target)
          mir_vect <-append(mir_vect, mir)
        }
        else{
          therapie <- target_therapie$Therapie[indices]
          therapie_vect <- append(therapie_vect, therapie)
          target_vect <-append(target_vect, c(rep(target, length(indices))))
          mir_vect <-append(mir_vect, c(rep(mir, length(indices))))
        }
      }
    }
    # Creation du dataframe -----
    final <- data.frame(miRNA = mir_vect, Prot =target_vect, Treatment = therapie_vect)
    if (is.null(final[1, 1])){
      return(data.frame(Issue = "No results"))
    }
    else{
      colnames(final)[2] <- "Target"
      return(final)
    }
  })
  output$acti <- renderDataTable({
    acti()
    },options = list(
      paging = TRUE,
      searching = FALSE,
      info = FALSE,
      pageLength = 10))
  
  # Inhibition ----
  repr <- reactive({
    activa <- acti_repr()[[1]]
    repres <- acti_repr()[[2]]
    # selection ----
    if (is.null(activa[1, 1])){activa_vect <- c()}
    else{
    activa_vect <- activa[as.numeric(activa$`Nb of strong evidence`)*input$poids_fort_ai+
                            as.numeric(activa$`Nb of less strong evidence`)*input$poids_faible_ai 
                          >= input$seuil_preuves_acti_repr,c("Target")]
    }
    if (is.null(repres[1, 1])){repres_dt <- data.frame()}
    else{
    repres_dt <- repres[as.numeric(repres$`Nb of strong evidence`)*input$poids_fort_ai+
                            as.numeric(repres$`Nb of less strong evidence`)*input$poids_faible_ai 
                          >= input$seuil_preuves_acti_repr,c("Mir","Target")]
    }
    
    repres_unique <- fct_mir_unique(repres_dt)
    go <- T
    
    if (is.null(repres_unique[1, 1])){go<- F}
    else{
      final_dtr<- repres_unique[!sapply(repres_unique$Target, is.element, activa_vect), c("Mir", "Target")]
    }
    
    # Therapie ----
    target_vect <- c()
    therapie_vect <- c()
    mirs_vect <- c()
    
    if (go){
      for (target in final_dtr$Target) {
        target_exact <- paste0('^', target, '$')
        mirs <- final_dtr$Mir[grep(target_exact, final_dtr$Target)]
        indices <- grep(target_exact, target_therapie$Target)
        if (length(indices)==0){
          therapie_vect <- append(therapie_vect, "No therapy available")
          target_vect <-append(target_vect, target)
          mirs_vect <-append(mirs_vect, mirs)
        }
        else{
          therapie <- target_therapie$Therapie[indices]
          therapie_vect <- append(therapie_vect, therapie)
          target_vect <-append(target_vect, c(rep(target, length(indices))))
          mirs_vect <-append(mirs_vect, c(rep(mirs, length(indices))))
        }
      }
    }
    # Creation du dataframe -----
    final <- data.frame(miRNA = mirs_vect, Prot =target_vect, Treatment = therapie_vect)
    if (is.null(final[1, 1])){
      return(data.frame(Issue = "No results"))
    }
    else{
      colnames(final)[2] <- "Target"
      return(final)
    }
  })
  output$repr <- renderDataTable({
    repr()
  },options = list(
    paging = TRUE,
    searching = FALSE,
    info = FALSE,
    pageLength = 10))
  
  # repr <- eventReactive(input$up_regul,
  #                       input$down_regul,{
  #                         notify <- function(msg, id = NULL) {
  #                           showNotification(msg, id = id, duration = NULL, closeButton = FALSE,type="message")
  #                         }
  #                         req(input$up_regul)
  #                         req(input$down_regul)
  #                         id <- notify("Processing...")
  #                         inFile <- input$up_regul
  #                         inFile2 <- input$down_regul
  #                         up <- fct_lire(req(inFile), "\n", T)
  #                         down <- fct_lire(req(inFile2), "\n", T)
  #                         on.exit(removeNotification(id), add = TRUE)
  #                         return(fct_inhibition(up, down, input$seuil_preuves_repr,
  #                                               input$poids_fort_ai, input$poids_faible_ai))
  #   
  # })
  # 
  # output$repr <- renderDataTable(
  #   repr(),
  #   options = list(
  #     paging = TRUE,
  #     searching = FALSE,
  #     info = FALSE,
  #     pageLength = 10)
  # )
  # 
  ############## Liste des mirs ----
  # Partie motif ----  
  Base_exploit_mot <- reactive({
    if(input$long2=="stem loop"){test <- T}
    else{test <- F}
    fct_Choix_Base("hsa", test)
  })
  
  liste <- eventReactive(input$go, {
    notify <- function(msg, id = NULL) {
      showNotification(msg, id = id, duration = NULL, closeButton = FALSE,type="message")
    }
    req(input$go)
    id <- notify("Processing...")
    tempo <- fct_listes_mirs(input$methyl, input$seuil_score_liste, Base_exploit_mot())
    on.exit(removeNotification(id), add = TRUE)
    return(tempo)
  })
  
  output$listes <- renderDataTable(
    liste(),
    options = list(
      paging = TRUE,
      pageLength = 10,
      searching = FALSE,
      info = FALSE,
      info = FALSE)
  )
  # Partie Score -----
  
  Base_exploit_sco <- reactive({
    if(input$long3=="stem loop"){test <- T}
    else{test <- F}
    fct_Choix_Base("hsa", test)
  })
  
  listes2 <- eventReactive(input$methyl2,{
    notify <- function(msg, id = NULL) {
      showNotification(msg, id = id, duration = NULL, closeButton = FALSE,type="message")
    }
    req(input$methyl2)
    id <- notify("Processing...")
    tempo <- fct_listes_mirs_2(req(input$methyl2), 
                               Base_exploit_sco())
    on.exit(removeNotification(id), add = TRUE)
    return(tempo)
    
  })
  
  output$listes2 <- renderDataTable(
    return(fct_select_mirs_score(listes2(), req(input$seuil_score2),
                                 input$poids_motif2, input$poids_liste2)),
    
    options = list(
      paging = TRUE,
      pageLength = 10,
      searching = FALSE,
      info = FALSE)
  )
  ############## Tam_Mir----
  output$Tam_fon <- renderDataTable({
    return(fct_Tam_mir(req(input$rna_tam)))
  },
  options = list(
    paging = TRUE,
    searching = FALSE,
    info = FALSE,
    pageLength = 10)
  ) 
  ############## Tam_fun----
  observe({
    dt <- Base_Tam$Category[Base_Tam$Title %in% input$Tit]
    updateSelectizeInput(session, "Cate", choices = dt, selected = '', server = TRUE)
  })
  output$mir_tam <- renderDataTable({
    return(fct_Tam_other(req(input$Tit), req(input$Cate)))
  },
  options = list(
    paging = TRUE,
    searching = FALSE,
    info = FALSE,
    pageLength = 10)
  )
  ############## Secondary_structure ----
  
  Base_exploit_ss <- reactive({
    if(input$long4=="stem loop"){test <- T}
    else{test <- F}
    fct_Choix_Base("hsa", test)
  })
  
  methy <- reactive({
    return(fct_methyl_boucle(req(input$rna_ss), Base_exploit_ss()))
  })
  
  Liste_up <- reactive({
    return(fct_update(req(methy())))
  })
  
  calcul_boucle <-eventReactive(input$go2, {
    i <- grep(req(input$choix_methy), Liste_up())
    Seq <- fct_lire(req(methy()[ i, "Sequence"]), "[AUCGaucg]", FALSE)
    liste <- fct_calcul(tolower(Seq))
    if (min(liste[[1]]) < 0){
      Sep <- fct_struc(liste, "Yes")
      return(list(Seq, Sep))
    }
    else{
      Sep1 <- fct_struc(fct_calcul(tolower(Seq)), "No")
      Sep2 <- fct_struc(fct_calcul(tolower(Seq)), "Yes")
      return(list(Seq, Sep1, Seq, Sep2))
    }

  })
  
  plot_boucle_t <- eventReactive(c(input$go2, input$go3), {
    i <- grep(req(input$choix_methy), Liste_up())
    if (length(calcul_boucle())==2){
      Seq <- calcul_boucle()[[1]]
      Sep <- calcul_boucle()[[2]]
      fct_graph(Seq, Sep, req(methy())[i ,c("Position")])
    }
    else{
      par(mar = c(4, 5, 0.5, 0.5), mfrow = c(1, 2))
      Seq <- calcul_boucle()[[1]]
      Sep <- calcul_boucle()[[2]]
      fct_graph(Seq, Sep, req(methy())[i ,c("Position")])
      Seq <- calcul_boucle()[[3]]
      Sep <- calcul_boucle()[[4]]
      fct_graph(Seq, Sep, req(methy())[i ,c("Position")])
    }

  })
  
  affi_methy <- eventReactive(input$go2, {
    req(input$rna_ss)
    i <- grep(req(input$choix_methy), Liste_up())
    return(methy()[i,])
  })
  output$met <- renderDataTable({
    methy <-affi_methy()
    colnames(methy)[1] <- c("miRNA name")
    methy
  },options = list(
    paging = FALSE,
    searching = FALSE,
    info = FALSE,
    pageLength = 10)
  )
  
  output$plot_boucle <- renderPlot (
    plot_boucle_t(), 
    width = "auto",
    height = 450,
  )

  
  output$Warning <- renderText({
    req(plot_boucle_t())
    return("Click on De-twisted to obtain a structure devoid of twist")
  })
  # Ancien programme ----
  # boucle <- reactive({
  #   return(fct_methyl_boucle(req(input$rna_ss), Base_util))
  # })
  # boucle_downlo <- reactive({
  #   boucle_final <- boucle()
  #   #boucle_final <- subset(boucle_final, select = c("Mir", "Taille_pont", "Dessin"))
  #   return(boucle_final)
  # })
  # 
  # output$boucle <- renderDataTable({
  #   boucle_final <- boucle()
  #   #boucle_final <- subset(boucle_final, select = c("Mir", "Taille_pont", "Taille_boucle","Dessin"))
  #   return(boucle_final)
  #   
  # }, options = list(
  #   paging = TRUE,
  #   searching = FALSE,
  #   info = FALSE,
  #   pageLength = 5,
  #   autoWidth = TRUE,
  #   order = list(list(3, 'desc'))
  # ))
  
  # output$Impact <- renderDataTable({
  #   impact <- fct_match_boucle(req(input$rna_ss), Base_util, input$seuil_score_loop)
  #   return(impact)
  # },options = list(
  #   paging = TRUE,
  #   searching = FALSE,
  #   info = FALSE,
  #   pageLength = 5)
  # )
  ############## About -----
  # Creation compteur visite ----
  values <- reactiveValues(
    logs = NULL # Initialisation
  )    
  observe({
    if (file.exists("logs.csv")) {
      values$logs <- read.table(
        "logs.csv", 
        sep = ",", 
        colClasses = "character",
        col.names = c("clic", "time")
      )
    } else {
      NULL
    }
  })    
  observeEvent(input$methyl, 
               {
                 log <- data.frame(
                   clic = "Viste",
                   time = as.character(Sys.time())
                 )
                 write.table(log, 
                             file = "logs.csv", 
                             sep = ",",
                             append = TRUE, 
                             row.names = FALSE,
                             col.names = FALSE)
                 if (is.null(values$logs)) {
                   values$logs <- log
                 } else {
                   values$logs <- rbind(
                     values$logs,
                     log
                   )
                 }
               })
  output$logs_table <- renderUI({
    if (!is.null(values$logs)) {
      test<-values$logs
      return(nrow(test))
    }
  })
  # Ref -----
  output$Ref <- renderDataTable({
    
    return(data.frame(References = Liste_ref$References))
  },options = list(
    paging = FALSE,
    searching = FALSE,
    info = FALSE,
    pageLength = 10) 
  )
  
  
  
  
})
