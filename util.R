# Obtention des librairies ----

library(shiny)
library(googledrive)
library(googlesheets4)
library(pracma)
library(stringr)
library(igraph)


# Obtention du drive pour l'obtention et la modification des listes ----
options(# whenever there is one account token found, use the cached token
  gargle_oauth_email = TRUE,
  # specify auth tokens should be stored in a hidden directory ".secrets"
  gargle_oauth_cache = ".secrets")

  # Obtention des listes du drive ----

    # Liste des pattern ---- 
sheet_name_p <- "Pattern"
sheet_id_p <- googledrive::drive_get(sheet_name_p)$id

    # Liste de presence ----
sheet_name <- "custom_lists"
sheet_id <- googledrive::drive_get(sheet_name)$id

    # Liste des therapies et de leur cibles ----
sheet_name_2 <- "target_therapie"
sheet_id_2 <- googledrive::drive_get(sheet_name_2)$id

  # Liste des exomotifs ----
sheet_name_3 <- "Liste_Exomotif"
sheet_id_3 <- googledrive::drive_get(sheet_name_3)$id
  # Liste des references ----
sheet_name_4 <- "References"
sheet_id_4<- googledrive::drive_get(sheet_name_4)$id


# Fonction Passe-Partout(de Fort-boyard)----

fct_lire <- function(cara, pattern, inversion){
  
  # Permet de lire  
  
  new_cara <- c()
  m <- gregexpr(pattern, cara)
  new_cara_2 <- regmatches(cara, m, invert = inversion)
  new_cara_3 <- new_cara_2[[1]]
  for (i in new_cara_3){new_cara <- append(new_cara, i)}
  #NC <- data.frame(nc = new_cara) #test pour vérifier si la fonction marchait
  return(new_cara)
}

fct_unique <- function(vect, cara){
  
  # N'ajoute un element que s'il n'est pas present dans la liste
  
  if (!is.element(cara, vect)){
    vect <- append(vect, cara)
  }
  return(vect)
}

# Obtention de la base final, des couples Mir-Target et des listes de methylation ----

  # Creation de la base finale ------ 

  # Recuperation de la base complete -----
Base <-
  read.csv(file = "Base_Totale.csv",
           col.names = c("Nom", "Seq"),
           header = FALSE)
Base$Seq <- tolower(Base$Seq)

Base_long <-
  read.csv(file = "miRNA_tige_boucle.csv",
           col.names = c("Nom", "Seq"),
           header = FALSE)
Base_long$Seq <- tolower(Base_long$Seq)

  # Extraction des mir humain ----

fct_Choix_Base <- function(type, long){
  if (long){
    indices <- grep(type, Base_long$Nom, value = FALSE)
    ExtracNom <- Base_long$Nom[indices]
    ExtracSeq <- Base_long$Seq[indices]
  }
  else{
    indices <- grep(type, Base$Nom, value = FALSE)
    ExtracNom <- Base$Nom[indices]
    ExtracSeq <- Base$Seq[indices]
  }
  BaseFinale <- data.frame(Nom = ExtracNom, Seq = ExtracSeq)
  return(BaseFinale)
}

Base_pour_Tam<-fct_Choix_Base('hsa', F)

Base_util <- fct_Choix_Base("hsa", F)

Base_util_long <- fct_Choix_Base("hsa", T)

ExtracNom <- unique(Base_util$Nom)

ExtracNom_long <- unique(Base_util_long$Nom)

  # Obtention des couples Mir-Target ----
Mirtargets <- read.csv(
  file = "Mirtargets.csv",
  col.names = c("mir", "target", "evidences"),
  header = FALSE,
  na.strings = "x"
)

nom <- "Base_targ/hsa_MTI.csv"
Base_ex <- read.csv(nom, header = TRUE)
indices_fun <- grep("Non-Functional MTI", Base_ex$Support.Type, invert = TRUE)
miRNA <- Base_ex$miRNA[indices_fun]
Target.gene <- Base_ex$Target.Gene[indices_fun]
experi <- Base_ex$Experiments[indices_fun]
Clean_target <- data.frame(Mir = miRNA, Target = Target.gene, Experiment = experi)


  # Obtention des listes de methylation ----

  # Liste m5C -----
# liste_m5C <- read.csv(
#   "Pattern/cg.csv",
#   header = FALSE,
#   na.strings = "x"
# )
# 
# colnames(liste_m5C) = c("100%")

  # Liste m6A ----
# liste_m6A <- read.csv(
#   "Pattern/drach.csv",
#   header = FALSE,
#   na.strings = "x"
# )
# 
# colnames(liste_m6A) = c("100%", "75%", "50%", "25%")

  # Liste m7G ----
# liste_m7G <- read.csv(
#   "Pattern/raggu.csv",
#   header = FALSE,
#   na.strings = "x"
# )
# 
# colnames(liste_m7G) = c("100%", "50%")

  # Obtention de la base Tam ----

Base_Tam <- read.csv("Mir_Tam_Finale.csv", header = FALSE)
colnames(Base_Tam) <- c("Title", "Category", "Mir")
extratam <- Base_Tam$Title
ExtracTam <- unique(extratam)


# Function for the prediction of methylation ----

fct_nom_to_seq <- function(arn_nom, BaseFinale) { 
  
  # Fais correspondre le nom du mir a sa sequence
  
  # indice <- which(BaseFinale$Nom == arn_nom)
  # if(length(indice)==0){return("")}
  
  #Autre manière pour avoir qu'un seul match (vlookup)
  arn_nom_exact <- paste0('^', arn_nom, '$')
  indice <- grep(arn_nom_exact, BaseFinale$Nom, ignore.case = T)
  if(length(indice)==0){return("")}
  
  arn_seq <- BaseFinale$Seq[indice]
  return(arn_seq)
  
}

fct_seq_to_nom <- function(arn_seq, BaseFinale) { 
  
  # Fais correspondre la sequence du mir a son nom
  
  arn_seq <- tolower(arn_seq)
  indice <- which(BaseFinale$Seq == arn_seq)
  arn_nom <- BaseFinale$Nom[indice[1]]
  return(arn_nom)
  
}

fct_methylation <- function(rnaseq, nom_muta, Score, col_de_liste, BaseFinale) {
  # Initialisation ----
  
  methylations_df <- data.frame()
  Mutation_vect <- c()
  Nom_vect <- c()
  Sequence_vect <- c()
  Sequence_mutee_vect <- c()
  Position_reele <- c()
  Position_debut <- c()
  Position_fin <- c()
  Score_vect <- c()
  
  
  # Recherche des sequences de methylation ----
  # g3<- c("uugc", "uuga", "uugu")
  
  test <- sapply(col_de_liste, grepl, rnaseq, ignore.case = T)
  mutas_found <- names(test)[which(test == TRUE)]
  if (!is.element(TRUE, test)) {return(methylations_df)}
  seq_mutee <- c()
  for (tofind in mutas_found) {
    x <- str_locate_all(tolower(rnaseq), pattern = tolower(tofind))
    all <- do.call(rbind, x)
    start <- all[, 1]
    stop <- all[, 2]
    # if(nom_muta == "m7G"){
    #   Position_debut <- append(Position_debut, start)
    #   if (tofind %in% g3){Position_reele <- append(Position_reele, start+2)}
    #   else{Position_reele <- append(Position_reele, start+3)}
    #   Position_fin <- append(Position_fin, stop)
    # }
    # if(nom_muta == "m5C"){
    #   Position_debut <- append(Position_debut, start)
    #   Position_reele <- append(Position_reele, start)
    #   Position_fin <- append(Position_fin, stop)
    # }
    # if(nom_muta == "m6A"){
    #   Position_debut <- append(Position_debut, start)
    #   Position_reele <- append(Position_reele, start+2)
    #   Position_fin <- append(Position_fin, stop)
    # }
    pos <- str_locate(tofind, "[AUCG]")[1,1] - 1
    Position_debut <- append(Position_debut, start)
    Position_reele <- append(Position_reele, start+pos)
    Position_fin <- append(Position_fin, stop)
    for (i in 1:length(start)) {
      Mutation_vect <- append(Mutation_vect, nom_muta)
      Nom_vect <- append(Nom_vect, fct_seq_to_nom(rnaseq, BaseFinale))
      Sequence_vect <- append(Sequence_vect, toupper(rnaseq))
      Score_vect <- append(Score_vect, Score)
      seq_mutee <- append(seq_mutee, tofind)
    }
  }
  Sequence_mutee_vect <-
    append(Sequence_mutee_vect, toupper(seq_mutee))
  
  # Creation du dataFrame final ----
  methylations_df <-
    data.frame(
      Mutation = Mutation_vect,
      Mirname = Nom_vect,
      Sequence = Sequence_vect,
      Sequence_mutee = Sequence_mutee_vect,
      Position_Reele = Position_reele,
      Position_Debut = Position_debut,
      Position_Fin = Position_fin,
      Score = Score_vect, 
      row.names = NULL
      
    )
  #on renomme les colonnes
  colnames(methylations_df)[1] <- "Methylation"
  colnames(methylations_df)[2] <- "miRNA name"
  colnames(methylations_df)[3] <- "Sequence"
  colnames(methylations_df)[5] <- "Position"
  rownames(methylations_df) <- NULL
  return(methylations_df)
  
  
}

fct_Matchtotal <- function(rna, BaseFinale) {
  
  # Cree un dataframe avec toutes les informations sur les motifs de methylation
  
  # Initialisation ----
  rnaseq <- fct_nom_to_seq(rna, BaseFinale)
  rna.res <- data.frame()
  # Ancienne versionn ----
  # # Motif de m5C 
  # col_de_liste = data_pattern$m5C_100
  # score = "100%"
  # nom_muta = "m5C"
  # rna.res <-
  #   rbind(rna.res,
  #         fct_methylation(rnaseq, nom_muta, score, col_de_liste, BaseFinale))
  # 
  # # Motif de m7G 
  # col_de_liste = data_pattern$m7G_50
  # score = "50%"
  # nom_muta = "m7G"
  # rna.res <-
  #   rbind(rna.res,
  #         fct_methylation(rnaseq, nom_muta, score, col_de_liste, BaseFinale))
  # 
  # 
  # col_de_liste = data_pattern$m7G_100
  # score = "100%"
  # nom_muta = "m7G"
  # rna.res <-
  #   rbind(rna.res,
  #         fct_methylation(rnaseq, nom_muta, score, col_de_liste, BaseFinale))
  # 
  # # Motif de m6A
  # col_de_liste = data_pattern$m6A_25
  # score = "25%"
  # nom_muta = "m6A"
  # rna.res <-
  #   rbind(rna.res,
  #         fct_methylation(rnaseq, nom_muta, score, col_de_liste, BaseFinale))
  # 
  # col_de_liste = data_pattern$m6A_50
  # score = "50%"
  # nom_muta = "m6A"
  # rna.res <-
  #   rbind(rna.res,
  #         fct_methylation(rnaseq, nom_muta, score, col_de_liste, BaseFinale))
  # 
  # col_de_liste = data_pattern$m6A_75
  # score = "75%"
  # nom_muta = "m6A"
  # rna.res <-
  #   rbind(rna.res,
  #         fct_methylation(rnaseq, nom_muta, score, col_de_liste, BaseFinale))
  # 
  # 
  # col_de_liste = data_pattern$m6A_100
  # score = "100%"
  # nom_muta = "m6A"
  # rna.res <-
  #   rbind(rna.res,
  #         fct_methylation(rnaseq, nom_muta, score, col_de_liste, BaseFinale))
  # Nouvelle version ----
  for (nom in names(data_pattern)){
    ind <- grep(nom, names(data_pattern))
    col_de_liste <- data_pattern[[ind]]
    lect <- fct_lire(nom, "_", T)
    nom_muta <- lect[1]
    score <- paste0(lect[2], "%")
    rna.res <-
        rbind(rna.res,
              fct_methylation(rnaseq, nom_muta, score, col_de_liste, BaseFinale))
  }
  
  # Retourne le dataframe finale ----
  if (is.null(rna.res[1, 1])){
    return(data.frame(Methylation = "No methylation", Mir = rna, Sequence = rnaseq))
  } 
  else{
    return(rna.res)
  }
}

fct_presence_listes <- function(rna) {
  
  # Cree la matrice de presence d'un mir dans une liste 
  
  # Initialisation ----
  Taille <- length(nom_listes)
  rna.liste <- data.frame(matrix(ncol = Taille, nrow = 1))
  x <-
    c(nom_listes)
  colnames(rna.liste) <- x
  
  # Verification de presence ----
  for (nom in nom_listes) {
    test <- grepl(rna, data_liste[nom], ignore.case = TRUE)
    if (is.element(TRUE, test)) {
      rna.liste[nom] <- 'Present'
    } else{
      rna.liste[nom] <- 'Absent'
    }
  }
  
  return(rna.liste)
}

fct_scores <- function(rna.res){
  
  # Calcule le score de methylation totale avec le nombre de pattern present et leur score
  
  # Initialisation ----
  
  # scores_m5C <- 0
  # scores_m7G <- 0
  # scores_m6A <- 0
  score_total <- c()
  nom_mutas <- c()
  
  # Calcul du score ----
  
  if (!(rna.res[1, 1]== "No methylation")){ # On verifie si le mir possede des patterns
    Mut <- rna.res$Methylation
    Sco <- rna.res$Score
    l<- length(Mut)
    
    for (i in 1:l){
      # if (Mut[i] == 'm5C'){
      #   scores_m5C <- scores_m5C + fct_transfo_pourcen_nb(Sco[i])
      # }
      # if (Mut[i] == 'm7G'){
      #   scores_m7G <- scores_m7G + fct_transfo_pourcen_nb(Sco[i])
      # }
      # if (Mut[i] == 'm6A'){
      #   scores_m6A <- scores_m6A + fct_transfo_pourcen_nb(Sco[i])
      # }
      ind  <- grep(Mut[i], nom_mutas)
      if (length(ind)==0){
        nom_mutas <- append(nom_mutas, Mut[i])
        score_total <- append(score_total, fct_transfo_pourcen_nb(Sco[i]))
      }
      else{
        score_total[ind] <- score_total[ind] + fct_transfo_pourcen_nb(Sco[i])
      }
    }
    
    
  }
  if (length(grep("m5C", nom_mutas))==0){
    nom_mutas <- append(nom_mutas, "m5C")
    score_total <- append(score_total, 0)
  }
  if (length(grep("m6A", nom_mutas))==0){
    nom_mutas <- append(nom_mutas, "m6A")
    score_total <- append(score_total, 0)
  }
  if (length(grep("m7G", nom_mutas))==0){
    nom_mutas <- append(nom_mutas, "m7G")
    score_total <- append(score_total, 0)
  }
  # Creation du dataframe resumant le score ----
  
  score_motif <- data.frame()
  score_motif <- rbind(score_motif, score_total)
  colnames(score_motif) <- nom_mutas
  
  return(score_motif)
}

fct_scores_listes <- function(rna.liste){
  
  # Calcule le score de presence du mir dans une liste
  
  # Initialisation ----
  
  Taille <- length(nom_listes)
  scores_listes <- data.frame(matrix(0,ncol = Taille, nrow = 1))
  x <-
    c(nom_listes)
  colnames(scores_listes) <- x
  name_lis <- c()
  
  # Calcul du score ----
  
  if (!is.null(rna.liste[1,1])){
    for (names in nom_listes){
      name_lis <- append(name_lis, fct_lire(names, "-", T)[2])
      if (rna.liste[names] == 'Present'){
        scores_listes[names] <- scores_listes[names] + 1
      }
    }
  }
  # Creation du dataframe resumant le score ----
  
  # colnames(scores_listes)[1]<- 'm5C'
  # colnames(scores_listes)[2]<- 'm6A'
  # colnames(scores_listes)[3]<- 'm7G'
  colnames(scores_listes) <- name_lis
  
  
  
  return(scores_listes)
}

fct_plot_scores <- function(rna, poids_motif, poids_liste, BaseFinale){
  
  # Cree le graphique resumant la prediction
  
  # Initialisation ----
  rna.res <- fct_Matchtotal(rna, BaseFinale)
  rna.liste <- fct_presence_listes(rna)
  scores_motif <- fct_scores(rna.res)
  scores_listes <- fct_scores_listes(rna.liste)
  
  data_inter <- data.frame(m5C = 0, m6A = 0, m7G = 0)
  
  for (i in 1:length(scores_motif)){
    if (names(scores_motif)[i] %in% names(data_inter)){
      if (names(scores_motif)[i] == "m5C"){
        data_inter$m5C[1] <- scores_motif[1, i]*poids_motif
      }
      if (names(scores_motif)[i] == "m6A"){
        data_inter$m6A[1] <- scores_motif[1, i]*poids_motif
      }
      if (names(scores_motif)[i] == "m7G"){
        data_inter$m7G[1] <- scores_motif[1, i]*poids_motif
      }
    }
    else{
      data_inter <- cbind(data_inter, scores_motif[1, i]*poids_motif)
      names(data_inter)[length(data_inter)] <- names(scores_motif)[i]
    }
  }
  for (nom in nom_pattern){
    if (!(nom %in% names(data_inter))){
      data_inter <- cbind(data_inter, 0)
      names(data_inter)[length(data_inter)] <- nom
    }
  }
  for (i in 1:length(scores_listes)){
    if (names(scores_listes)[i] %in% names(data_inter)){
      ind <- grep(names(scores_listes)[i], names(data_inter))
      data_inter[1, ind] <- data_inter[1, ind] + poids_liste*scores_listes[1, i]
    }
    else{
      data_inter <- cbind(data_inter, scores_listes[i]*poids_liste)
      names(data_inter)[length(data_inter)] <- names(scores_listes)[i]
    }
  }
  
  # data_inter <- (poids_liste*scores_listes)+(poids_motif*scores_motif)
  
  # Creation du graph ----
  if(!is.null(data_inter[1, 1])){ # On verifie que l'on a des données
    min <- 0
    max <-rowSums(data_inter)
    data <- rbind(rep(max, length(data_inter)) , rep(min, length(data_inter)) , data_inter)
    titre <-
      paste("MemiRPred: Predictive methylation score for", rna)
    par(mar = c(0, 2, 4, 2))
    plot <- radarchart(
      data,
      centerzero = TRUE,
      #custom polygon
      pcol = rgb(0.2, 0.5, 0.5, 0.9) ,
      pfcol = rgb(0.2, 0.5, 0.5, 0.5) ,
      plwd = 3 ,
      
      #custom the grid
      cglcol = "grey",
      cglty = 1,
      axislabcol = "grey",
      # caxislabels = seq(0, max, max / 3),
      cglwd = 0.8,
      
      #custom labels
      vlcex = 1.3,
      title = titre
    )
    
    box(col = "black",
        which = "figure",
        lwd = 2)
    return(plot)
  }
}

fct_liste_methyl <- function(mirs, poids_motif, poids_liste){
  dt <- data.frame()
  dt_final <- data.frame()
  nom_final <- c("m5C", "m6A", "m7G")
  
  for (mir in mirs){
    
    test <- c(F, F)
    if (fct_nom_to_seq(mir, Base_util)!= ""){test[1] <- T}
    if (fct_nom_to_seq(mir, Base_util_long)!= ""){test[2] <- T}
    
    
    if (test[1]){
      scores_motif <- fct_scores(fct_Matchtotal(mir, Base_util))
      scores_listes <- fct_scores_listes(fct_presence_listes(mir))
      # m5c <- as.numeric(score_motif[1,1]*coef_motif+score_liste[1,1]*coef_liste)
      # m6a <- as.numeric(score_motif[1,2]*coef_motif+score_liste[1,2]*coef_liste)
      # m7g <- as.numeric(score_motif[1,3]*coef_motif+score_liste[1,3]*coef_liste)
      data_inter <- data.frame(m5C = 0, m6A = 0, m7G = 0)
      
      for (i in 1:length(scores_motif)){
        if (names(scores_motif)[i] %in% names(data_inter)){
          if (names(scores_motif)[i] == "m5C"){
            data_inter$m5C[1] <- scores_motif[1, i]*poids_motif
          }
          if (names(scores_motif)[i] == "m6A"){
            data_inter$m6A[1] <- scores_motif[1, i]*poids_motif
          }
          if (names(scores_motif)[i] == "m7G"){
            data_inter$m7G[1] <- scores_motif[1, i]*poids_motif
          }
        }
        else{
          data_inter <- cbind(data_inter, scores_motif[1, i]*poids_motif)
          names(data_inter)[length(data_inter)] <- names(scores_motif)[i]
          nom_final <- fct_unique(nom_final, names(scores_motif)[i])
        }
      }
      for (nom in nom_pattern){
        if (!(nom %in% names(data_inter))){
          data_inter <- cbind(data_inter, 0)
          names(data_inter)[length(data_inter)] <- nom
          nom_final <- fct_unique(nom_final, nom)
        }
      }
      for (i in 1:length(scores_listes)){
        if (names(scores_listes)[i] %in% names(data_inter)){
          ind <- grep(names(scores_listes)[i], names(data_inter))
          data_inter[1, ind] <- data_inter[1, ind] + poids_liste*scores_listes[1, i]
        }
        else{
          data_inter <- cbind(data_inter, scores_listes[i]*poids_liste)
          names(data_inter)[length(data_inter)] <- names(scores_listes)[i]
          nom_final <- fct_unique(nom_final, names(scores_listes)[i])
        }
      }
      
      data_final <- cbind(mir, data_inter)
      
      dt_final <- rbind(dt_final, data_final)
    }
    if (test[2]){
      scores_motif <- fct_scores(fct_Matchtotal(mir, Base_util_long))
      scores_listes <- fct_scores_listes(fct_presence_listes(mir))
      data_inter <- data.frame(m5C = 0, m6A = 0, m7G = 0)
      
      for (i in 1:length(scores_motif)){
        if (names(scores_motif)[i] %in% names(data_inter)){
          if (names(scores_motif)[i] == "m5C"){
            data_inter$m5C[1] <- scores_motif[1, i]*poids_motif
          }
          if (names(scores_motif)[i] == "m6A"){
            data_inter$m6A[1] <- scores_motif[1, i]*poids_motif
          }
          if (names(scores_motif)[i] == "m7G"){
            data_inter$m7G[1] <- scores_motif[1, i]*poids_motif
          }
        }
        else{
          data_inter <- cbind(data_inter, scores_motif[1, i]*poids_motif)
          names(data_inter)[length(data_inter)] <- names(scores_motif)[i]
          nom_final <- fct_unique(nom_final, names(scores_motif)[i])
        }
      }
      for (nom in nom_pattern){
        if (!(nom %in% names(data_inter))){
          data_inter <- cbind(data_inter, 0)
          names(data_inter)[length(data_inter)] <- nom
          nom_final <- fct_unique(nom_final, nom)
        }
      }
      for (i in 1:length(scores_listes)){
        if (names(scores_listes)[i] %in% names(data_inter)){
          ind <- grep(names(scores_listes)[i], names(data_inter))
          data_inter[1, ind] <- data_inter[1, ind] + poids_liste*scores_listes[1, i]
        }
        else{
          data_inter <- cbind(data_inter, scores_listes[i]*poids_liste)
          names(data_inter)[length(data_inter)] <- names(scores_listes)[i]
          nom_final <- fct_unique(nom_final, names(scores_listes)[i])
        }
      }
      
      data_final <- cbind(paste("stem loop", mir, " "), data_inter)
      colnames(data_final) <- names(dt_final)
      
      dt_final <- rbind(dt_final, data_final)
    }
    if ((!test[1])&&(!test[2])){
      if (!is.null(dt_final[1, 1])){
        Taille <- length(dt_final) 
        data_inter <- matrix(data= 0, ncol = Taille, nrow = 1)
        data_inter[1, 1] <- mir
        colnames(data_inter) <- names(dt_final)
        
        dt_final <- rbind(dt_final, data_inter)
      }
      else{
        return(data.frame(Issue = "No mir matching your input"))
      }
    }
  }
  nom <- c("miRNA")
  for(nm in nom_final){
    nom <- append(nom, paste("Score", nm, " ")) 
  }
  colnames(dt_final) <- nom
  return(dt_final)
}

fct_stats_methy <- function(mirs, poids_motif, poids_liste, nb){
  dt <- data.frame()
  dt_final <- data.frame()
  nom_final <- c("m5C", "m6A", "m7G")
  
  for (mir in mirs){
    
    test <- c(F, F)
    if (fct_nom_to_seq(mir, Base_util)!= ""){test[1] <- T}
    if (fct_nom_to_seq(mir, Base_util_long)!= ""){test[2] <- T}
    
    
    if (test[1]){
      scores_motif <- fct_scores(fct_Matchtotal(mir, Base_util))
      scores_listes <- fct_scores_listes(fct_presence_listes(mir))
      # m5c <- as.numeric(score_motif[1,1]*coef_motif+score_liste[1,1]*coef_liste)
      # m6a <- as.numeric(score_motif[1,2]*coef_motif+score_liste[1,2]*coef_liste)
      # m7g <- as.numeric(score_motif[1,3]*coef_motif+score_liste[1,3]*coef_liste)
      data_inter <- data.frame(m5C = 0, m6A = 0, m7G = 0)
      
      for (i in 1:length(scores_motif)){
        if (names(scores_motif)[i] %in% names(data_inter)){
          if (names(scores_motif)[i] == "m5C"){
            data_inter$m5C[1] <- scores_motif[1, i]*poids_motif
          }
          if (names(scores_motif)[i] == "m6A"){
            data_inter$m6A[1] <- scores_motif[1, i]*poids_motif
          }
          if (names(scores_motif)[i] == "m7G"){
            data_inter$m7G[1] <- scores_motif[1, i]*poids_motif
          }
        }
        else{
          data_inter <- cbind(data_inter, scores_motif[1, i]*poids_motif)
          names(data_inter)[length(data_inter)] <- names(scores_motif)[i]
          nom_final <- fct_unique(nom_final, names(scores_motif)[i])
        }
      }
      for (nom in nom_pattern){
        if (!(nom %in% names(data_inter))){
          data_inter <- cbind(data_inter, 0)
          names(data_inter)[length(data_inter)] <- nom
          nom_final <- fct_unique(nom_final, nom)
        }
      }
      for (i in 1:length(scores_listes)){
        if (names(scores_listes)[i] %in% names(data_inter)){
          ind <- grep(names(scores_listes)[i], names(data_inter))
          data_inter[1, ind] <- data_inter[1, ind] + poids_liste*scores_listes[1, i]
        }
        else{
          data_inter <- cbind(data_inter, scores_listes[i]*poids_liste)
          names(data_inter)[length(data_inter)] <- names(scores_listes)[i]
          nom_final <- fct_unique(nom_final, names(scores_listes)[i])
        }
      }
      
      data_final <- cbind(mir, data_inter)
      
      dt_final <- rbind(dt_final, data_final)
    }
    if (test[2]){
      scores_motif <- fct_scores(fct_Matchtotal(mir, Base_util_long))
      scores_listes <- fct_scores_listes(fct_presence_listes(mir))
      data_inter <- data.frame(m5C = 0, m6A = 0, m7G = 0)
      
      for (i in 1:length(scores_motif)){
        if (names(scores_motif)[i] %in% names(data_inter)){
          if (names(scores_motif)[i] == "m5C"){
            data_inter$m5C[1] <- scores_motif[1, i]*poids_motif
          }
          if (names(scores_motif)[i] == "m6A"){
            data_inter$m6A[1] <- scores_motif[1, i]*poids_motif
          }
          if (names(scores_motif)[i] == "m7G"){
            data_inter$m7G[1] <- scores_motif[1, i]*poids_motif
          }
        }
        else{
          data_inter <- cbind(data_inter, scores_motif[1, i]*poids_motif)
          names(data_inter)[length(data_inter)] <- names(scores_motif)[i]
          nom_final <- fct_unique(nom_final, names(scores_motif)[i])
        }
      }
      for (nom in nom_pattern){
        if (!(nom %in% names(data_inter))){
          data_inter <- cbind(data_inter, 0)
          names(data_inter)[length(data_inter)] <- nom
          nom_final <- fct_unique(nom_final, nom)
        }
      }
      for (i in 1:length(scores_listes)){
        if (names(scores_listes)[i] %in% names(data_inter)){
          ind <- grep(names(scores_listes)[i], names(data_inter))
          data_inter[1, ind] <- data_inter[1, ind] + poids_liste*scores_listes[1, i]
        }
        else{
          data_inter <- cbind(data_inter, scores_listes[i]*poids_liste)
          names(data_inter)[length(data_inter)] <- names(scores_listes)[i]
          nom_final <- fct_unique(nom_final, names(scores_listes)[i])
        }
      }
      
      data_final <- cbind(paste("stem loop", mir, " "), data_inter)
      colnames(data_final) <- names(dt_final)
      
      dt_final <- rbind(dt_final, data_final)
    }
    if ((!test[1])&&(!test[2])){
      if (!is.null(dt_final[1, 1])){
        Taille <- length(dt_final) 
        data_inter <- matrix(data= 0, ncol = Taille, nrow = 1)
        data_inter[1, 1] <- mir
        colnames(data_inter) <- names(dt_final)
        
        dt_final <- rbind(dt_final, data_inter)
      }
      else{
        return(data.frame(Issue = "No mir matching your input"))
      }
    }
  }
  liste_quant <- list()
  
  for (i in 2:length(dt_final)){
    liste_quant[[i-1]] <- quantile(as.numeric(dt_final[,i]), linspace(0, 1, nb+1))
  }

  for (i in 1:length(dt_final[,1])){
    # dt <- rbind(dt, c(fct_quantile(C, as.numeric(dt_final[i, 2])),
    #                   fct_quantile(A, as.numeric(dt_final[i, 3])),
    #                   fct_quantile(G, as.numeric(dt_final[i, 4]))))
    # dt <- rbind(dt, c(fct_quantile_general(C, as.numeric(dt_final[i, 2])),
    #                   fct_quantile_general(A, as.numeric(dt_final[i, 3])),
    #                   fct_quantile_general(G, as.numeric(dt_final[i, 4]))))
    quan <- c()
    for (j in 1:length(liste_quant)){
      quan <- append(quan, fct_quantile_general(liste_quant[[j]], as.numeric(dt_final[i, j+1])))
    }
    dt <- rbind(dt, quan)
  }
  dt_final <- cbind(dt_final, dt)
  nom <- c("miRNA")
  for(nm in nom_final){
    nom <- append(nom, paste("Score", nm, " ")) 
  }
  for(nm in nom_final){
    nom <- append(nom, "Stat") 
  }
  colnames(dt_final) <- nom
  return(dt_final)
}

fct_quantile <- function(Quan, n){
  
  q1 <- as.numeric(Quan[2])
  m <- as.numeric(Quan[3])
  q3 <- as.numeric(Quan[4])
  max <- as.numeric(Quan[5])
  
  q1_m <- FALSE
  m_q3 <- FALSE
  q3_max <- FALSE
  
  if (q1 == m){q1_m <- TRUE}
  if (q3 == m){m_q3 <- TRUE}
  if (q3 == max){q3_max <- TRUE}
  
  if((n <= as.numeric(Quan[2]))&&!q1_m){return("0-25%")}
  if((n >=as.numeric(Quan[2]))&&(n<=as.numeric(Quan[3]))&&!m_q3){return("25-50%")}
  if((n >=as.numeric(Quan[3]))&&(n<=as.numeric(Quan[4]))&&!q3_max){return("50-75%")}
  if(n >= as.numeric(Quan[4])){return("75-100%")}
  
}

fct_quantile_general <- function(Quan, n){
  
  if (n == 0){
    return("0%-25%")
  }
  
  names <- row.names(as.data.frame(Quan))
  Inter <- as.data.frame(Quan)$Quan
  
  indice <- findInterval(n, Inter)
  
  if (indice != length(names)){ final <- paste(names[indice], names[indice + 1], sep="-")}
  else{final <- paste(names[indice-1], names[indice], sep="-")}
  
  return(final)
  
}

# Function for the exomotif ----

fct_exomes <- function(arn_nom, BaseFinale){
  
  # Cree un dataframe resumant les informations sur les exomotifs
  
  # Recuperation des donnees ----
  
  rna.res <- fct_Matchtotal(arn_nom, BaseFinale)
  if (!(rna.res[1, 1]== "No methylation")){
    arn_seq <- fct_nom_to_seq(arn_nom, BaseFinale)
    methylation <- rna.res$Position
    methy <- rna.res$Sequence_mutee
    
    # Initialisation ----
    
    Nom_vect <- c()
    Sequence_vect <- c()
    Position_debut_ex <- c()
    Position_fin_ex <- c()
    Position_ex <- c()
    Position_methy <- c()
    Motif_methylation <- c()
    l<- length(methylation)
    
    # Recherche d'exomotif et d'overlapping ----
    
    test <- sapply(Liste_Exomotif$Exomotif, grepl, arn_seq)
    exom_found <- names(test)[which(test == TRUE)]
    if (is.element(TRUE, test)) {
      presence <- c()
      overlapping <- c()
      for (tofind in exom_found) {
        x <- str_locate_all(arn_seq, pattern = tofind)
        all <- do.call(rbind, x)
        start <- all[, 1]
        stop <- all[, 2]
        Position_debut_ex <- append(Position_debut_ex, start)
        Position_fin_ex <- append(Position_fin_ex, stop)
        m <- length(start)
        
        for(i in 1:l){
          for (p in 1:m){
            if ((methylation[i]<=stop[p])&&(methylation[i]>=start[p])){
              overlapping <- append(overlapping, 'Overlapping')
            }
            else{
              overlapping <- append(overlapping, 'No overlapping')
            }
            Position_methy <- append(Position_methy, methylation[i])
            Motif_methylation <- append(Motif_methylation, methy[i])
            Position_ex <- append(Position_ex, start[p])
            presence <- append(presence, toupper(tofind))
            Nom_vect <- append(Nom_vect, arn_nom)
            Sequence_vect <- append(Sequence_vect, toupper(arn_seq))
          } 
        }
      }
      
      # Creation du dataframe resumant les informations ----
      Exomes_df <-
        data.frame(
          Mirname = Nom_vect,
          Sequence = Sequence_vect,
          Motif_Methy = Motif_methylation,
          Position_Methylation = Position_methy,
          Motif_Exomes = presence,
          Position_Debut_ex = Position_ex,
          Overlapping = overlapping
        )
      
      colnames(Exomes_df)[1] <- "miRNA name"
      colnames(Exomes_df)[3] <- "Methylation Pattern"
      colnames(Exomes_df)[4] <- "Position of the methylation nucleoside"
      colnames(Exomes_df)[5] <- "Exomotif"
      colnames(Exomes_df)[6] <- "Position of the exomotif"
      return(Exomes_df)
    }
    else{
      Exomes_df <-
        data.frame(
          Mirname = arn_nom,
          Sequence = toupper(arn_seq),
          Motif_Exomes = c("No exomotif")
        )
      
      colnames(Exomes_df)[1] <- "miRNA name"
      colnames(Exomes_df)[3] <- "Exomotif"
      return(Exomes_df)
    }
  }
  else{
    arn_seq <- fct_nom_to_seq(arn_nom, BaseFinale)
    # Initialisation ----
    
    Nom_vect <- c()
    Sequence_vect <- c()
    Position_debut_ex <- c()
    Position_fin_ex <- c()
    Position_ex <- c()
    presence <- c()
    
    # Recherche d'exomotif et d'overlapping ----
    
    test <- sapply(Liste_Exomotif$Exomotif, grepl, arn_seq)
    exom_found <- names(test)[which(test == TRUE)]
    if (is.element(TRUE, test)) {
      for (tofind in exom_found) {
        x <- str_locate_all(arn_seq, pattern = tofind)
        all <- do.call(rbind, x)
        start <- all[, 1]
        stop <- all[, 2]
        Position_debut_ex <- append(Position_debut_ex, start)
        Position_fin_ex <- append(Position_fin_ex, stop)
        m <- length(start)
        
        for(i in 1:m){
            presence <- c(presence, toupper(tofind))
            Nom_vect <- append(Nom_vect, arn_nom)
            Sequence_vect <- append(Sequence_vect, toupper(arn_seq))
          } 
        }
      }
      
      # Creation du dataframe resumant les informations ----
      if(length(presence)==0){
        Exomes_df <-
          data.frame(
            Mirname = arn_nom,
            Sequence = toupper(arn_seq),
            Motif_Exomes = c("No exomotif")
          )
        
        colnames(Exomes_df)[1] <- "miRNA name"
        colnames(Exomes_df)[3] <- "Exomotif"
        return(Exomes_df)
      }
      Exomes_df <-
        data.frame(
          Mirname = Nom_vect,
          Sequence = Sequence_vect,
          Motif_Exomes = presence,
          Position_Debut_ex = Position_debut_ex
        )
      
      colnames(Exomes_df)[1] <- "miRNA name"
      colnames(Exomes_df)[3] <- "Exomotif"
      colnames(Exomes_df)[4] <- "Position of the exomotif"
      return(Exomes_df)
  }
 
}

# Function for the targets and triades ----

fct_score_tar_fort <- function(ex){
  score <- 0
  match_exp <- c("Reporter assay", "Western", "PCR")
  exp_pres <- sapply(match_exp, grepl, ex, ignore.case = TRUE)
  if (length(ex)<2){
    if (is.element(TRUE, exp_pres["Reporter assay"])){score <- score +1}
    if (is.element(TRUE, exp_pres["Western"])){score <- score +1}
    if (is.element(TRUE, exp_pres["PCR"])){score <- score +1}
  }
  else{
    if (is.element(TRUE, exp_pres[,"Reporter assay"])){score <- score +1}
    if (is.element(TRUE, exp_pres[,"Western"])){score <- score +1}
    if (is.element(TRUE, exp_pres[,"PCR"])){score <- score +1}
  }
  return (score)
}

fct_score_tar_faible <- function(ex){
  score <- 0
  match_exp <- c("Microarray", "NGS", "pSILAC", "CLIP")
  exp_pres <- sapply(match_exp, grepl, ex, ignore.case = TRUE)
  if (length(ex)<2){
    if (is.element(TRUE, exp_pres["Microarray"])){score <- score +1}
    if (is.element(TRUE, exp_pres["NGS"])){score <- score +1}
    if (is.element(TRUE, exp_pres["pSILAC"])){score <- score +1}
    if (is.element(TRUE, exp_pres["CLIP"])){score <- score +1}
  }
  else{
    if (is.element(TRUE, exp_pres[,"Microarray"])){score <- score +1}
    if (is.element(TRUE, exp_pres[,"NGS"])){score <- score +1}
    if (is.element(TRUE, exp_pres[,"pSILAC"])){score <- score +1}
    if (is.element(TRUE, exp_pres[,"CLIP"])){score <- score +1}
  }
  return (score)
}

fct_Targets <- function(arn_nom) {
  
  # Recherhe les targets d'un mir donnees avec le nombre de preuves associés
  
  arn_nom_exact <- paste0('^', arn_nom, '$')
  indices <- grep(arn_nom_exact, Mirtargets$mir)
  targets_vect <- Mirtargets$target[indices]
  evidences_vect <- Mirtargets$evidences[indices]
  
  if (length(targets_vect) == 0) {
    return('')
  }
  targets <-
    data.frame(Targets = targets_vect, Evidences = evidences_vect)
  return(targets)
}

fct_triade_m_t_t <- function(rna, nombre_evidence_mini) {
  
  # Cree les triades Mir-Target-Treatment
  
  # Recuperation de donnees --------
  
  mirtargets <- fct_Targets(rna)
  
  # test pour verifier s'il existe des targets --------
  
  targets_found <-mirtargets[mirtargets[2] >= nombre_evidence_mini, ]
  if (is.na(targets_found[1, 1])) {
    return(data.frame(Mir = rna, Issue = "No target with enough evidence"))
  }
  
  # Initialisation ------
  
  target_vect <- c()
  therapie_vect <- c()
  mir_vect <- c()
  
  # Recherche de treatment -------
  
  for (target in targets_found$Targets) {
    target_exact <- paste0('^', target, '$')
    indices <- grep(target_exact, target_therapie$Target)
    if (length(indices)==0){
      therapie_vect <- append(therapie_vect, "No therapie available")
      mir_vect <- append(mir_vect, rna)
      target_vect <-
        append(target_vect, target)
    }
    else{
      therapie <- target_therapie$Therapie[indices]
      therapie_vect <- append(therapie_vect, therapie)
      mir_vect <- append(mir_vect, c(rep(rna, length(indices))))
      target_vect <-
        append(target_vect, c(rep(target, length(indices))))
    }
  }
  
  # Creation des triades -----
  
  m_t_t <-
    data.frame(miR = mir_vect,
               Target = target_vect,
               Therapie = therapie_vect)
  colnames(m_t_t)[3] <- "Treatment"
  return(m_t_t)
}

fct_Target_General <- function(mir, Base){
  
  
  target <- data.frame()
  if (mir == ""){return(target)}
  mir_f <- paste0('^', mir, '$') 
  indices <- grep(mir_f, Base$Mir)
  if( is.na(indices[1]) ){return(target)}
  mir_tar <- Base[indices, c('Target', 'Experiment')]
  prot_unique <- c()
  prot_unique <- unique(mir_tar$Target)
  for (prot in prot_unique){
    score_fort <- 0
    score_faible <- 0
    indices2<- grep(prot, mir_tar$Target)
    ex<-mir_tar$Experiment[indices2]
    clean_ex <- sapply(ex, fct_lire, "//|;", TRUE)
    exp_unique <- c()
    for (g_exp in clean_ex){
      exp_unique <- append(exp_unique, unique(g_exp))
    }
    score_fort <- fct_score_tar_fort(exp_unique)
    score_faible <- fct_score_tar_faible(exp_unique)
    if (score_faible+score_fort < length(exp_unique)){
      score_faible <- score_faible+1
    }
    target<- rbind(target, c(prot, score_fort, score_faible))
  }
  if (!is.null(target[1, 1])){
    colnames(target)<- c('Target', "Nb of strong evidence", "Nb of less strong evidence")
  }
  return(target)
}

fct_triade_m_t_t_g <- function(mir, targ, coeff_fort, coeff_faible, seuil_evid){
  
  prot_seuil <- targ$Target[as.numeric(targ$`Nb of strong evidence`)*coeff_fort + as.numeric(targ$`Nb of less strong evidence`)*coeff_faible >= seuil_evid]
  
  if (is.na(prot_seuil[1])) {
    return(data.frame(Mir = mir, Issue = "No target with enough evidence"))
  }
  
  target_vect <- c()
  therapie_vect <- c()
  mir_vect <- c()
  
  # Recherche de treatment -------
  
  for (target in prot_seuil) {
    target_exact <- paste0('^', target, '$')
    indices <- grep(target_exact, target_therapie$Target)
    if (length(indices)>0){
      therapie <- target_therapie$Therapie[indices]
      therapie_vect <- append(therapie_vect, therapie)
      mir_vect <- append(mir_vect, c(rep(mir, length(indices))))
      target_vect <-
        append(target_vect, c(rep(target, length(indices))))
    }
    else{
      therapie_vect <- append(therapie_vect, "No therapy available")
      mir_vect <- append(mir_vect, mir)
      target_vect <- append(target_vect, target)
    }
  }
  
  # Creation des triades -----
  m_t_t <-
    data.frame(miRNA = mir_vect,
               Target = target_vect,
               Therapie = therapie_vect)
  colnames(m_t_t)[3] <- "Treatment"
  return(m_t_t)
  
}

# Function for the Treatment to mir ----

fct_therapie_to_target <- function(therapie){
  
  # Sort les targets d'une therapie donnees
  
  target_vect <- c()
  indices <- grep(therapie, target_therapie$Therapie)
  target <- target_therapie$Target[indices]
  target_vect<- append(target_vect, target)
  
  return(target_vect)
}

fct_target_to_mir <- function(target) {
  
  # Sort les mirs impactant une target donnee
  
  target_exact <- paste0('^', target, '$')
  indices <- grep(target_exact, Mirtargets$target)
  mir_vect <- Mirtargets$mir[indices]
  evidences <- Mirtargets$evidences[indices]
  if (length(mir_vect) == 0) {
    return(data.frame())
  }
  mirs <-
    data.frame(Proteine = target, Mirs = mir_vect, Evidences = evidences)
  return(mirs)
}

fct_target_to_mir_general <- function(target, Base){
    
    mirs <- data.frame()
    indices <- grep(target, Base$Target)
    if( is.na(indices[1]) ){return(mirs)}
    base_tar <- Base[indices, c('Mir', 'Experiment')]
    mir_unique <- c()
    mir_unique <- unique(base_tar$Mir)
    for (mir in mir_unique){
      score_fort <- 0
      score_faible <- 0
      indices2<- grep(mir, base_tar$Mir)
      ex<-base_tar$Experiment[indices2]
      clean_ex <- sapply(ex, fct_lire, "//|;", TRUE)
      exp_unique <- c()
      for (g_exp in clean_ex){
        exp_unique <- append(exp_unique, unique(g_exp))
      }
      score_fort <- fct_score_tar_fort(exp_unique)
      score_faible <- fct_score_tar_faible(exp_unique)
      if (score_faible+score_fort < length(exp_unique)){
        score_faible <- score_faible+1
      }
      mirs<- rbind(mirs, c(mir, target, score_fort, score_faible))
    }
    if (!is.null(mirs[1, 1])){
      colnames(mirs)<- c("Mir","Target", "Nb of strong evidence", "Nb of less strong evidence")
    }
    return(mirs)
  
  
}

fct_therapie_to_mir<- function(therapie, Base){
  
  # Assosie une therapie a des mirs ayant une target soignee par la therapie
  
  mirs_final<- data.frame()
  mirs_dt<-data.frame()
  target_vect<-fct_therapie_to_target(therapie)
  
  for (target in target_vect){
    #mirs_dt<- fct_target_to_mir(target)3
    mirs_dt<- fct_target_to_mir_general(target, Base)
    if (!is.null(mirs_dt[1,1])){
      mirs_final<- rbind(mirs_final, mirs_dt)
    }
  }
  
  if (length(mirs_final)==0){
    return(data.frame())
  }
  return(mirs_final)
}

fct_therapie_to_methyl <- function(therapie, coeff_fort, coeff_faible, 
                                   seuil_preuves, seuil_score, Base_targ, BaseFinale){
  
  # Associe une therapie a des mirs ayant des pattern de methylation 
  
  # Recuperation des donnees -----
  
  target_mir<- fct_therapie_to_mir(therapie, Base_targ)
  
  # Verification des donnees -----
  
  if (is.null(target_mir[1,1])){
    return(data.frame(Issue = "No miRNA found"))
  }
  else{
    final <- data.frame()
    inter <- data.frame()
    tempo <- data.frame()
  }
  l <- length(target_mir$Mir)
  
  # Recherche des mirs ayant un pattern de methylation ----
  
  
  for (i in 1:l){
    inter <- fct_Matchtotal(target_mir$Mir[i], BaseFinale)
    if (!(inter[1, 1] == "No methylation")){
      l_sc <- length(inter$Score)
      for (k in 1:l_sc){
        if (fct_transfo_pourcen_nb(inter$Score[k]) >= fct_transfo_pourcen_nb(seuil_score)){
          tempo <- data.frame(Proteine = target_mir$Target[i],
                              Evidences = as.numeric(target_mir$`Nb of strong evidence`[i])*coeff_fort
                              +as.numeric(target_mir$`Nb of less strong evidence`[i])*coeff_faible, 
                              Mir_name = target_mir$Mir[i], 
                              Methylation = inter$Methylation[k], 
                              Position = inter$Position[k], 
                              Score = inter$Score[k])
          tempo <- tempo[tempo$Evidences >= seuil_preuves ,]
          final<- rbind(final, tempo)
        }
      }
    }
  }
  if (is.na(final[1,1])){
    return(data.frame(Issue = "No miRNA found"))
  }
  else{
    colnames(final) <- c("Proteine", "Evidence", "miRNA name", "Methylation", "Position of the methylated nucleoside", "Score")
    return(final)
  }
}

fct_transfo_pourcen_nb<- function(cara){
  
  # Permet de passer d'un score en % à un nombre
  
  score<-0
  if (cara =='100%'){score<-1}
  if (cara =='75%'){score<-0.75}
  if (cara =='50%'){score<-0.5}
  if (cara =='25%'){score<-0.25}
  return(score)
}

fct_transfo_nb_pourcent <- function(nb){
  score<-""
  if (nb ==1){score<-"100%"}
  if (nb ==0.75){score<-"75%"}
  if (nb ==0.5){score<-"50%"}
  if (nb ==0.25){score<-"25%"}
  if (nb == 0){score<- "0%"}
  return(score)
}

# Function for the Activation and Inhibition ----

fct_acti_repr <- function(up_reg, down_reg){
  
  # Cree la liste des proteines activees
  
  # Initialisation ----
  
  activa <- data.frame()
  activa_vect <- c()
  activa_unique <- c()
  repres <- data.frame()
  repres_vect <- c()
  final_vect <- c()
  final <- data.frame()
  
  # Recuperation des donnees ----
  
  for (mir_u in up_reg){
    if (mir_u != ""){
      tar <- fct_Target_General(mir_u, Clean_target)
      mir <- c(rep(mir_u, length(tar[,1])))
      inter <- cbind(mir, tar)
      if (!is.null(inter[1, 1])){colnames(inter)[1] <- "Mir"}
      activa <- rbind(activa, inter)
    }
  }
  for (mir_d in down_reg){
    if (mir_d != ""){
      tar <- fct_Target_General(mir_d, Clean_target)
      mir <- c()
      if(!is.null(tar[1, 1])){mir <- c(rep(mir_d, length(tar[,1])))}
      inter <- cbind(mir, tar)
      if (!is.null(inter[1, 1])){colnames(inter)[1] <- "Mir"}
      repres <- rbind(repres, inter)
    }
  }
  
  return(list(activa, repres))
  

}

fct_mir_unique <- function(dt){
  
  final <- data.frame()
  tar_uni <- c()
  for (tar in dt$Target){
    if (!(tar %in% tar_uni)){
      tar_uni <- append(tar_uni, tar)
      tot <- dt$Mir[grep(paste0('^', tar, "$"), dt$Target)]
      mirs <- tot[1]
      for(mir in tot[-1]){
        mirs <- paste(mirs, mir, sep = ", ")
      }
      inter <- data.frame(Mir = mirs, Target = tar)
      final <- rbind(final, inter) 
    }
    
  }
  return(final)
}

# Function for the Liste pattern ----

fct_meilleurs_scores <- function(rna.res){
  
  # Cree un dataframe avec le meilleur score d'un pattern de methylation
  
  # Initialisation ----
  
  Taille <- length(nom_pattern)
  score_meilleur <- data.frame(matrix(0,ncol = Taille, nrow = 1))
  x <-
    c(nom_pattern)
  colnames(score_meilleur) <- x
  
  # Recherche du score ----
  
  if (!(rna.res[1, 1]== "No methylation")){
    Mut <- rna.res$Methylation
    Sco <- rna.res$Score
    l<- length(Mut)
    for (i in 1:l){
      # if (Mut[i]=='m5C' && fct_transfo_pourcen_nb(Sco[i]) > mm5C){
      #   mm5C <- fct_transfo_pourcen_nb(Sco[i])
      # }
      # if (Mut[i]=='m6A' && fct_transfo_pourcen_nb(Sco[i]) > mm6A){
      #   mm6A <- fct_transfo_pourcen_nb(Sco[i])
      # }
      # if (Mut[i]=='m7G' && fct_transfo_pourcen_nb(Sco[i]) > mm7G){
      #   mm7G <- fct_transfo_pourcen_nb(Sco[i])
      # }
      
      ind <- grep(Mut[i], nom_pattern)
      score_meilleur[1, ind] <- max(score_meilleur[1, ind], fct_transfo_pourcen_nb(Sco[i]))
    }
  }

  return(score_meilleur)
}

fct_listes_mirs <- function(methylation, seuil_score, BaseFinale){
  
  # Cree la liste des mirs ayant un score au-dessus du seuil
  
  # Initialisation ----
  
  list <- data.frame()
  list_vect <- c()
  
  # Extraction des mirs----
  
  for (mir in BaseFinale$Nom){
    scores_max <- fct_meilleurs_scores(fct_Matchtotal(mir, BaseFinale))
    if (!is.null(scores_max[1, 1])){
      if (scores_max[methylation]>= fct_transfo_pourcen_nb(seuil_score)){
        list_vect <- append(list_vect, mir)
      }
    }
  }
  
  # Creation du dataframe ----
  list<- data.frame(Mirs = list_vect)
  return(list)
}

fct_decodage <- function(methylation){
  
  # Fonction palliative pour lire la methylation(N'est plus utilise)
  
  if (methylation == 'm5C'){
    return(1)
  }
  if (methylation == 'm6A'){
    return(2)
  }
  if (methylation == 'm7G'){
    return(3)
  }
}

# Function for the Liste Score ----

fct_listes_mirs_2 <- function(methylation, BaseFinale){
  
  # Cree un dataframe regroupant les scores de pattern et de liste de chaque mir
  
  # Initialisation ----
  
  list <- data.frame()
  list_vect <- c()
  score_vect <- c()
  score_vect_liste <- c()
  
  # Calcul des scores ----
  
  for (mir in BaseFinale$Nom){
    list_vect <- append(list_vect, mir)
    score_vect <- append(score_vect, max(fct_scores(fct_Matchtotal(mir, BaseFinale))[1, methylation], 0))
    score_vect_liste <- append(score_vect_liste, max(fct_scores_listes(fct_presence_listes(mir))[1, methylation], 0))
  }
  
  # Creation du dataframe ----
  
  list<- data.frame(Mirs = list_vect, Score_motif = score_vect, Score_liste = score_vect_liste)
  return(list) 
}

fct_select_mirs_score <- function(liste, seuil_score, coeff_mot, coeff_liste){
  
  # Ressort la liste des mirs ayant un score total superieur au seuil
  if (is.null(liste[1, 1])){return(NULL)}
  list_final <- subset(liste[liste$Score_motif * coeff_mot + liste$Score_liste * coeff_liste >= seuil_score ,], 
                       select = c("Mirs"))
  colnames(list_final) <- "miRNA"
  return(list_final)
}

# Function for the Tam ----

fct_Tam_mir <- function(mir){
  
  pattern <- "-3p"
  pattern2 <- "-5p"
  m <- gregexpr(pattern, mir)
  new_cara_2 <- regmatches(mir, m, invert = TRUE)
  mir <- new_cara_2[[1]][1]
  m <- gregexpr(pattern2, mir)
  new_cara_2 <- regmatches(mir, m, invert = TRUE)
  mir <- new_cara_2[[1]][1]
  
  
  indices <- grep(mir, Base_Tam$Mir)
  Tit <- Base_Tam$Title[indices]
  Cat <- Base_Tam$Category[indices]
  Fin <- data.frame(Title = Tit, Category = Cat)
  return(Fin)
}

fct_Tam_other <- function(Title, Ctg){
  indices1 <- grep(Title, Base_Tam$Title)
  Mirinter <- Base_Tam$Mir[indices1]
  cata <- Base_Tam$Category[indices1]
  Inter <- data.frame(Cata = cata, Mir = Mirinter)
  indices2 <- grep(Ctg, Inter$Cata)
  Mir_fi <- data.frame(Mir = Inter$Mir[indices2])
  return(Mir_fi)
}

# Function for the structure ----

fct_complementary <- function(arn){
  A<-chartr("u", "A", arn)
  U<-chartr("a", "U", A)
  C<-chartr("g", "C", U)
  G<-chartr("c", "G", C)
  FI<-tolower(G)
  return (FI)
}

fct_equal <- function(arn1, arn2){
  l1 <-length(arn1)
  l2 <- length(arn2)
  tarn1 <- tolower(arn1)
  tarn2 <- tolower(arn2) 
  compte <- 0 
  true_l <- min(l1, l2)
  if(true_l==0){return(compte)}
  for (i in 1:true_l){
    if (!(((tarn1[i]=="a")&&(tarn2[i]=="g"))|| 
          ((tarn1[i]=="c")&&(tarn2[i]=="u"))||
          (tarn1[i]==tarn2[i]))){return(compte)}
    compte <- compte + 1
  }
  return(compte)
}

Match_Boucle <- c("uucg", "gaaa", "gcaa", "gaga", "guga", "ggaa", "cuug", "uuug")

fct_boucle <-function(mir, Base_util){
  taille_boucle <- c()
  taille_pont <- c()
  position_debut <- c()
  position_fin <- c()
  pattern_boucle_final <- c()
  pattern_final <- c()
  brin_1 <- c()
  brin_2 <- c()
  dessin_vect <- c()
  lect_mir <- fct_lire(mir, "[aucgAUCG]", FALSE)
  lmir <- length(lect_mir)
  for (lp in 2:4){
    true_l <- lmir-lp
    for (p in 2:(true_l)){
      true_lp <- p+lp
      cop_mir<-lect_mir
      pattern <- ""
      for (i in p:true_lp){
        pattern <- paste0(pattern, lect_mir[i])
        cop_mir[i] <- "X"
      }
      mir_tempo <- ""
      for (nuc in cop_mir){mir_tempo <- paste(mir_tempo, nuc, sep = "")}
      arn1 <- fct_lire(mir_tempo, "X{1,}", TRUE)[1]
      arn2 <- fct_lire(mir_tempo, "X{1,}", TRUE)[2]
      lect_arn1 <- fct_complementary(rev(fct_lire(arn1, "[aucgAUCG]", FALSE)))
      lect_arn2 <- fct_lire(arn2, "[aucgAUCG]", FALSE)
      pont<-fct_equal(lect_arn1, lect_arn2)
      if (pont > 1){
        amont <- ""
        for(k in pont:1){amont <- paste(amont, lect_mir[p-k], sep = "")}
        aval <- ""
        for(k in 1:pont){aval <- paste(aval, lect_mir[true_lp+k], sep = "")}
        pattern_final <- append(pattern_final, 
                                paste(amont, toupper(pattern), aval, sep = ""))
        brin_amont <- ""
        t_l <- p-pont-1 
        if (t_l == 0){brin_amont <- 'Pas de brin amont'}
        else{
          for (m in 1:t_l){
            brin_amont <- paste(brin_amont, lect_mir[m], sep = "")
          }
        }
        brin_aval <- ""
        t_m <- true_lp+pont+1
        if (t_m == lmir+1){brin_aval <- 'Pas de brin aval'}
        else{
          for (q in t_m:lmir){
            brin_aval <- paste(brin_aval, lect_mir[q], sep = "")
          } 
        }
        dessin <- fct_dessin_boucle(lect_mir, pont, p-pont, true_lp+pont, pattern)
        dessin_vect <- append(dessin_vect, dessin)
        taille_pont <- append(taille_pont, pont)
        taille_boucle <- append(taille_boucle, lp+1)
        position_debut <- append(position_debut, p-pont)
        position_fin <- append(position_fin, true_lp+pont)
        pattern_boucle_final <- append(pattern_boucle_final, toupper(pattern))
        brin_1 <- append(brin_1, brin_amont)
        brin_2 <- append(brin_2, brin_aval)
      }
    }
  }
  if(is.null(taille_boucle[1])){return(data.frame())}
  dt_final <- data.frame(Mir = fct_seq_to_nom(mir, Base_util), Nb_liasion = taille_pont,
                         Taille_boucle = taille_boucle, 
                         Motif_Boucle = pattern_final,
                         Motif_Turn = pattern_boucle_final,
                         Position_debut = position_debut,
                         Position_fin = position_fin,
                         Brin_amont = brin_1, Brin_aval = brin_2, 
                         Dessin = dessin_vect)
  return(dt_final)
}

fct_dessin_boucle <-function(lect_mir, pont, position_deb, position_fin, pattern){
  
  dessin <- ""
  dessin2 <-""
  dessin1 <- ""
  longueur_pattern <- length(fct_lire(pattern, "[AUGCaugc]", inversion = FALSE))
  lmi <- length(lect_mir)
  true_length <- max((position_deb - 1), (lmi - position_fin))
  lp2 <- floor(longueur_pattern/2)
  pos_pat <- position_deb - 1 + pont
  avm <- pos_pat+lp2
  if ((position_deb - 1) != 0){
    for (i in 1:(position_deb - 1)){
      dessin <- paste(dessin, lect_mir[i], sep = "")
    }
  }
  for (i in position_deb:pos_pat){
    dessin <- paste(dessin, lect_mir[i], sep = "")
  }
  if(longueur_pattern%%2 != 0){
    milieu <- lect_mir[pos_pat+1+lp2]
  }else{milieu <- ")"}
  
  for(i in (pos_pat+1):avm){
    dessin <- paste(dessin, lect_mir[i], sep = "")
  }
  ######################
  for (i in position_deb:pos_pat){
    dessin1 <- paste(dessin1, "| " , sep ="" )
  }
  if ((lp2 - 1)!=0){for (i in 1:(lp2 - 1)){dessin1 <- paste(dessin1, "--", sep = "")}}
  dessin1 <- paste(dessin1, milieu, sep = "")
  #############################################
  inv_ind <- position_fin - pont - lp2 + 1
  
  if(lmi > (position_fin+1)){
    for (i in lmi : (position_fin+1)){
      dessin2 <- paste(dessin2, lect_mir[i], sep = "")
    }
  }
  for (i in position_fin:inv_ind){
    dessin2 <- paste(dessin2, lect_mir[i], sep = "")
  }
  
  Dessin <- c("", "\n", dessin, "\n", dessin1, "\n", dessin2, "\n", "")
  Dessin_t <- toString(Dessin)
  Dessin_finale <- str_remove_all(toString(Dessin_t), ",")
  return(Dessin_finale)
}

fct_match_boucle <- function(Mir, Base_util, seuil_sco){
  
  methyl <- fct_Matchtotal(Mir, Base_util)
  if (is.null(methyl[1, 1])){return(data.frame())}
  
  boucle <- fct_boucle(fct_nom_to_seq(Mir, Base_util), Base_util)
  if (is.null(boucle[1, 1])){return(data.frame())}
  
  lb<- length(boucle$Mir)
  lm <- length(methyl$Position)
  
  methylation <- c()
  Pos_methyl <- c()
  Pos_pont <- c()
  sco <- c()
  dt_final <- data.frame()
  
  for (i in 1:lb){
    for (j in 1:lm){
      if(((methyl$Position[j] >= boucle$Position_debut[i])&&(methyl$Position[j] <= boucle$Position_fin[i]))||
         ((methyl$Position_Fin[j] >= boucle$Position_debut[i])&&(methyl$Position_Fin[j] <= boucle$Position_fin[i]))){
        methylation <- append(methylation, methyl$Methylation[j])
        Pos_methyl <- append(Pos_methyl, methyl$Position_Reele[j])
        Pos_pont <- append(Pos_pont, boucle$Position_debut[i])
        sco <- append(sco, methyl$Score[j])
      }
    }
  }
  
  ls <- length(sco)
  if (ls == 0){return(dt_final)}
  for (t in 1:ls){
    if (fct_transfo_pourcen_nb(sco[t]) >= fct_transfo_pourcen_nb(seuil_sco)){
      dt_final <- rbind(dt_final, c(Mir, Pos_pont[t], methylation[t], Pos_methyl[t], sco[t]))
    }
  }
  
  colnames(dt_final) <- c("Mir", "Position of the beginning of the loop",
                          "Methylation", "Position of the methylation", "Score")
  
  return(dt_final)
}

fct_methyl_boucle <- function(Mir, Base_util){
  
  methyl <- fct_Matchtotal(Mir, Base_util)
  #if (is.null(methyl[1, 1])){return(fct_boucle(fct_nom_to_seq(Mir, Base_util), Base_util))}
  if ((methyl[1, 1]== "No methylation")){
    return(data.frame(Sequence = toupper(fct_nom_to_seq(Mir, Base_util)), 
                      Position = 0, 
                      Methylation = c("No methylation")))
  }
  
  boucle_methyl <- data.frame()
  
  l <- length(methyl$Sequence)
  
  for (i in 1:l){
    mir_seq <- fct_lire(tolower(methyl$Sequence), "[aucg]", FALSE)
    mir_seq_methy <- ""
    if (methyl$Position_Debut[i] != 1){
      for (j in 1:(methyl$Position_Debut[i]-1)){
        mir_seq_methy <- paste(mir_seq_methy, mir_seq[j], sep = "")
      }
    }
    for (j in methyl$Position_Debut[i]:methyl$Position_Fin[i]){
      mir_seq_methy <- paste(mir_seq_methy, toupper(mir_seq[j]), sep="")
    }
    if (methyl$Position_Fin[i] != length(mir_seq)){
      for(j in (methyl$Position_Fin[i] + 1):length(mir_seq)){
        mir_seq_methy <- paste(mir_seq_methy, mir_seq[j], sep="")
      }
    }
    
    #boucle <- fct_boucle(mir_seq_methy, Base_util)
    
    tempo <- data.frame(Mir_name = Mir, Sequence = mir_seq_methy,
                        Position = methyl$Position[i],
                        Position_debut = methyl$Position_Debut[i],
                        Methylation = methyl$Methylation[i], 
                        Score = methyl$Score[i]) 
    #Nb_liasion = boucle$Nb_liasion)
    #Dessin = boucle$Dessin)
    
    boucle_methyl <- rbind(boucle_methyl, tempo)
    
  }
  colnames(boucle_methyl)[4] <- "Position Start"
  return(boucle_methyl)
}


# Structure Zuker ----
  # Donnée ----
    # Donnée sur l'empilement des liaisons ----
empile_98 <- matrix(data = -1*c(0.9, 1.7, 2.1, 0.9, 0.5, 1.0, 
                             1.8, 2.9, 2.0, 1.7, 1.2, 1.9,
                             2.3, 3.4, 2.9, 2.1, 1.4, 2.1, 
                             1.1, 2.3, 1.8, 0.9, 0.8, 1.1, 
                             1.1, 2.1, 1.9, 1.0, 0.4, 1.5, 
                             0.8, 1.4, 1.2, 0.5, 0.2, 0.4 ),
                 nrow = 6, ncol = 6)
colnames(empile_98)<-c("A/U", "C/G","G/C","U/A","G/U","U/G")
rownames(empile_98)<-c("A/U", "C/G","G/C","U/A","G/U","U/G")

empile_04 <- matrix(data = -1*c(0.9, 1.3, 2.4, 2.1, 1.3, 1.0,
                                1.1, 0.9, 2.2, 2.1, 1.4, 0.6, 
                                2.1, 2.1, 3.3, 2.4, 2.1, 1.4, 
                                2.2, 2.4, 3.4, 3.3, 2.5, 1.5, 
                                0.6, 1.0, 1.5, 1.4, 0.5, -0.3, 
                                1.4, 1.3, 2.5, 2.1, -1.3, 0.5),
                    nrow = 6, ncol = 6)

colnames(empile_04)<-c("A/U", "U/A", "C/G", "G/C", "G/U", "U/G")
rownames(empile_04)<-c("A/U", "U/A", "C/G", "G/C", "G/U", "U/G")
    # Donnée et interpolation sur le desequilibre des boucles ----
x_h <- c(5, 10, 20, 30)
y_h <- c(4.4, 5.3, 6.1, 6.5)
model_h <- lm(y_h~log(x_h))
a_h <- model_h$coefficients[2]
b_h <- model_h$coefficients[1]

x_b <- c(1, 5, 10, 20, 30)
y_b <- c(3.9, 4.8, 5.5, 6.3, 6.7)
model_b <- lm(y_b~log(x_b))
a_b <- model_b$coefficients[2]
b_b <- model_b$coefficients[1]

x_i <- c(5, 10, 20, 30)
y_i <- c(5.3, 6.6, 7.0, 7.4)
model_i <- lm(y_i~log(x_i))
a_i <- model_i$coefficients[2]
b_i <- model_i$coefficients[1]

bou_perte <- matrix(data = NA, nrow = 3, ncol = 120)
rownames(bou_perte)<-c("Hairpin", "Bulge", "Interne_loop")
for (i in 1:120){
  if(i <= 30){
    bou_perte["Hairpin", i] <- b_h+a_h*log(i)
    bou_perte["Bulge", i] <- b_b+a_b*log(i)
    bou_perte["Interne_loop", i] <- b_i+a_i*log(i) 
  }
  else{
    bou_perte["Hairpin", i] <- 1000
    bou_perte["Bulge", i] <- 1000
    bou_perte["Interne_loop", i] <- 1000
  }

}
bou_perte["Hairpin", 1] <- NA # Pas d'Hairpin d'un longueur d'une base
bou_perte["Hairpin", 2] <- NA # Pas d'Hairpin d'un longueur de deux base
bou_perte["Interne_loop", 1] <- NA # Pas de boucle interne de longueur de une base

bou_perte_04 <- read.table(file = "loop.txt")
colnames(bou_perte_04)<- c("Internal loop", "Bulge", "Hairpin")

  # Fonction ----

fct_calcul <- function(Seq){
  
  # Calcul 3 matrices. La premiere represente l'energie minimale entre deux paires de bases. 
  # La deuxieme correspond a la "direction" de repliement du miARN
  # La troisieme correspond a l'energie minimale si il y a appariement entre i et j
  
  # Creation des matrices
  len <- length(Seq)
  
  E <- matrix(data = 1000, nrow = len, ncol = len) # On le met à 1000 pour etre sur de ne pas fausser les calculs
  V <- matrix(data = 1000, nrow = len, ncol = len) # On le met à 1000 pour etre sur de ne pas fausser les calculs
  
  
  # Petit point sur Dir[i, j]
  # 1 = Meilleur structure entre i et j-1
  # 2 = Appariement entre i et j
  # 3 = Meilleur structure entre i+1 et j
  # 5,XX = Meilleur structure est la somme de (i, XX) et (XX + 1, j)
  # 8 = Hairpin entre i et j
  # 7,XX = Renfoncement de XX base en i
  # 9,XX = Renfoncement de XX base en i
  Dir <- matrix(data = 0, nrow = len, ncol = len)
  
  # L'algorithme
  for (k in 4:(len - 1)){ # On a besoin d'au moins 5 base pour faire une structure coherente
    for (i in 1:(len - k)){
      j <- i+k
      
      brin_1 <- fct_complementary(c(Seq[i], Seq[i+1]))
      brin_2 <- c(Seq[j], Seq[j-1])
      com <- fct_equal(brin_1, brin_2) # les resultats possibles sont 0 si il n'y a aucun appariement, 1 si i et j s'apparient et 2 si tout est apparie
      
      # Energie de l'assemblage de deux structures entre i et j
      El_final <- 10000
      l_final<- 0
      if (j-i > 2){
        for (l in ((i+1):(j-1))){
          El <- E[i, l]+E[l+1, j]
          if (El < El_final){
            El_final <- El
            l_final <- l
          }
        }
      }
      
      # Cas où il n'y a aucun appariement 
      if (com == 0){
        E[i, j] <- min(E[i+1, j], E[i, j-1], El_final, na.rm = TRUE)
        if(E[i, j]== El_final){Dir[i, j] <- 1 + 0.001*l_final}
        if(E[i, j] == E[i+1, j]){Dir[i, j] <- 3}
        if(E[i, j] == E[i, j-1]){Dir[i, j] <- 1}
      }
      
      # Calcul des renfoncement possible
      bul_final_i<- 10000
      mfinal_i <- 0
      bul_final_j<- 10000
      mfinal_j <- 0
      if(j-i > 6){
        for(m in 1:(j-i-6)){
          if(fct_equal(fct_complementary(c(Seq[i], Seq[i+1+m])), c(Seq[j], Seq[j-1]))==2){
            # li <- paste(toupper(Seq[j]), toupper(Seq[i]), sep = "/")
            # co <- paste(toupper(Seq[j-1]), toupper(Seq[i+m+1]), sep = "/")
            # alpha <- empile_98[li, co]
            li <- paste(toupper(Seq[i]), toupper(Seq[j]), sep = "/")
            co <- paste(toupper(Seq[i+m+1]), toupper(Seq[j-1]),  sep = "/")
            alpha <- empile_04[li, co]
            if(m == 1){bul <- alpha+V[i+m+1, j-1]+bou_perte["Bulge", m]}
            else{bul <- V[i+m+1, j-1]+bou_perte["Bulge", m]}
            # bul <- alpha+V[i+m+1, j-1]+bou_perte["Bulge", m]
            #bul <- alpha+V[i+m+1, j-1]+as.numeric(bou_perte_04$Bulge[m])

            if (bul < bul_final_i){
              mfinal_i <- m
              bul_final_i <- bul
            }
          }
          if(fct_equal(fct_complementary(c(Seq[i], Seq[i+1])), c(Seq[j], Seq[j-1-m]))==2){
            # li <- paste(toupper(Seq[j]), toupper(Seq[i]), sep = "/")
            # co <- paste(toupper(Seq[j-m-1]), toupper(Seq[i+1]), sep = "/")
            # alpha <- empile_98[li, co]
            li <- paste(toupper(Seq[i]), toupper(Seq[j]), sep = "/")
            co <- paste(toupper(Seq[i+1]), toupper(Seq[j-m-1]), sep = "/")
            alpha <- empile_04[li, co]
            if(m == 1){bul <- alpha+V[i+1, j-m-1]+bou_perte["Bulge", m]}
            else{bul <- V[i+1, j-m-1]+bou_perte["Bulge", m]}
            # bul <- alpha+V[i+1, j-m-1]+bou_perte["Bulge", m]
            #bul <- alpha+V[i+1, j-m-1]+as.numeric(bou_perte_04$Bulge[m])
            if (bul < bul_final_j){
              mfinal_j <- m
              bul_final_j <- bul
            }
          }
        }
      }
      # Calcul des boucles internes possible
      int_loop_final <- 10000
      k_final <- 0
      if (j-i >= 8){
        for (k1 in 1:(j-i-7)){
          for (k2 in 2:(j-i-7)){
            if (j-i-6-k1-k2 >= 0){
              if(fct_equal(fct_complementary(c(Seq[i], Seq[i+k1+1])), c(Seq[j], Seq[j-1-k2]))==2){
                # li <- paste(toupper(Seq[j]), toupper(Seq[i]), sep = "/")
                # co <- paste(toupper(Seq[j-k2-1]), toupper(Seq[i+1+k1]), sep = "/")
                # alpha <- empile_98[li, co]
                li <- paste(toupper(Seq[i]), toupper(Seq[j]), sep = "/")
                co <- paste(toupper(Seq[i+1+k1]), toupper(Seq[j-k2-1]), sep = "/")
                alpha <- empile_04[li, co]
                kp <- (k1+k2)
                if (kp == 2){int_loop <- alpha+V[i+1+k1, j-k2-1]+bou_perte["Interne_loop", kp]}
                else{int_loop <- V[i+1+k1, j-k2-1]+bou_perte["Interne_loop", kp]}
                # int_loop <- alpha+V[i+1+k1, j-k2-1]+bou_perte["Interne_loop", kp]
                #int_loop <- alpha+V[i+1+k1, j-k2-1]+as.numeric(bou_perte_04$`Internal loop`[kp])
                if (int_loop < int_loop_final){
                  k_final<- k1+k2*0.001
                  int_loop_final <- int_loop
                }
              }
            }
          }
        }
      }
      
      # Cas où i et j s'apparient et comparaison des structures
      
      if (com == 1){
        #V[i, j] <- min(as.numeric(bou_perte_04$Hairpin[(j-i-1)]), bul_final_i, bul_final_j, int_loop_final)
        V[i, j] <- min(bou_perte["Hairpin",(j-i-1)], bul_final_i, bul_final_j, int_loop_final)
        E[i, j] <- min(E[i+1, j], E[i, j-1], V[i, j], El_final, na.rm = TRUE)
        if(E[i, j] == El_final){Dir[i, j] <- 5 + 0.001*l_final}
        if(E[i, j] == E[i+1, j]){Dir[i, j] <- 3}
        if(E[i, j] == E[i, j-1]){Dir[i, j] <- 1}
        #if(V[i, j] == bou_perte_04$Hairpin[(j-i-1)]){Dir[i, j] <- 8}
        if(V[i, j] == bou_perte["Hairpin",(j-i-1)]){Dir[i, j] <- 8}
        if(V[i, j] == bul_final_i){Dir[i, j] <- 7+mfinal_i*0.001}
        if(V[i, j] == bul_final_j){Dir[i, j] <- 9+mfinal_j*0.001}
        if(V[i, j] == int_loop_final){Dir[i, j] <- 4+k_final*0.001}
      }
      # Cas où tout est apparié et comparaison des structures
      if (com == 2){
        # li <- paste(toupper(Seq[j]), toupper(Seq[i]), sep = "/")
        # co <- paste(toupper(Seq[j-1]), toupper(Seq[i+1]), sep = "/")
        # alpha <- empile_98[li, co]
        li <- paste(toupper(Seq[i]), toupper(Seq[j]), sep = "/")
        co <- paste(toupper(Seq[i+1]), toupper(Seq[j-1]), sep = "/")
        alpha <- empile_04[li, co]
        #hairpin <- as.numeric(bou_perte_04$Hairpin[(j-i-1)])
        hairpin <- bou_perte["Hairpin", (j-i-1)]
        
        V[i, j] <- min(alpha+V[i+1, j-1], hairpin,bul_final_j, bul_final_i, int_loop_final, na.rm = TRUE)
        E[i, j] <- min(E[i+1, j], E[i, j-1], V[i, j], El_final, na.rm = TRUE)
        if(E[i, j] == El_final){Dir[i, j] <- 5 + 0.001*l_final}
        if(E[i, j] == E[i+1, j]){Dir[i, j] <- 3}
        if(E[i, j] == E[i, j-1]){Dir[i, j] <- 1}
        if(V[i, j] == (alpha+V[i+1, j-1])){Dir[i, j] <- 2}
        if(V[i, j] == hairpin){Dir[i, j] <- 8}
        if(V[i, j] == bul_final_i){Dir[i, j] <- 7+mfinal_i*0.001}
        if(V[i, j] == bul_final_j){Dir[i, j] <- 9+mfinal_j*0.001}
        if(V[i, j] == int_loop_final){Dir[i, j] <- 4+k_final*0.001}
      }
      # Changement pour permettre la bonne création de structure dans les algorithmes prochains
      if (j-i <=5 ){
        Dir[i, j] <- 8
      }
    }
  }
  # # Remise à 0 des Energies
  # for(i in 1:length(E)){
  #   if (round(E[i]) == 1000){
  #     E[i] <- 0
  #   }
  #   if (round(V[i]) == 1000){
  #     V[i] <- 0
  #   }
  # }
  return(list(E, Dir, V))
}

fct_struc <- function(list, fold){
  
  # Amorce de la création de la structure secondaire sous forme de parenthèse
  
  i <- 1
  j <- sqrt(length(list[[1]]))
  seq <- c()
  for(l in 1:sqrt(length(list[[1]]))){
    seq <- append(seq, c("."))
  }
  
  # On cherche le minimum d'Energie
  min_E<- min(list[[1]])
  
  if ((min_E > 0)&&(fold == "No")){return(list(seq, 0))}

  
  # On parcourt la matrice pour trouver l'endroit 
  # avec l'energie minimale avec le minimum de bases 
  while(list[[1]][i+1, j] == min_E){
    i<- i+1
  }
  while(list[[1]][i, j-1] == min_E){
    j<- j-1
  }
  
  # Creation de la structure parenthésé
  
  seq <- fct_re(i, j, seq, list[[2]])
  return(list(seq, min_E))
}

fct_re <- function(i, j, seq, Dir){
  
  # Creation de la structure secondaire la structure secondaire
  # Avec la fonction d'avant, on s'assure de commencer sur un 2, 5, 7, 9
  
  
  if(Dir[i, j]==2){
    seq[i]<- "("
    seq[j]<- ")"
    seq <- fct_re(i+1, j-1, seq, Dir)
  }

  if(Dir[i, j]==8){
    seq[i]<- "("
    seq[j]<- ")"
    return(seq)
  }
  if(floor(Dir[i, j])==5){
    jonc <-(Dir[i, j]-5 )*1000
    i_1 <- i
    j_1 <- jonc
    i_2 <- jonc + 1
    j_2 <- j
    seq <- fct_re(i_1, j_1, seq, Dir)
    seq <- fct_re(i_2, j_2, seq, Dir)
  }
  if(floor(Dir[i, j])==7){
    seq[i]<- "("
    seq[j]<- ")"
    decale <- round(((Dir[i, j])-7)*1000)
    seq <- fct_re((i+1+decale), j-1, seq, Dir)
  }
  if(floor(Dir[i, j])==9){
    seq[i]<- "("
    seq[j]<- ")"
    decale <- round(((Dir[i, j])-9)*1000)
    seq <- fct_re(i+1, (j-1-decale), seq, Dir)
  }
  if(floor(Dir[i, j])==4){
    seq[i]<- "("
    seq[j]<- ")"
    decale_i <- round(((Dir[i, j])-4)*1000)
    decale_j <- round(((Dir[i, j]-4)*1000-decale_i)*1000)
    seq <- fct_re((i+1+decale_i), (j-1-decale_j), seq, Dir)
  }
  return(seq)
}

fct_graph<- function(Seq, l_Seqp, Pos_met){
  
  # Creation de la structure secondaire sous format graphique.
  # Etant donné que l'on passe par des graphes dirigé, il est possible que la 
  # chaine se croise de temps en temps.
  Seqp <- l_Seqp[[1]]
  
  l <- length(Seqp)
  from <- c(1:(l-1))
  to <- c(2:l)
  tmp <- matrix(100, nrow = 1, ncol = (l-1))
  poids <- tmp[1,]
  pile <- c()
  for(i in 1:l){
    if (Seqp[i]=="("){
      pile <- append(pile, i)
    }
    if(Seqp[i]==")"){
      from <- append(from, pile[length(pile)])
      to <- append(to, i)
      pile <- pile[-length(pile)]
      poids <- append(poids, 5)
    }
  }
  Edges <- data.frame(From = from, To = to, Weight = poids)
  
  
  id <- c(1:l)
  if (Pos_met == 0){
    methylated <- c()
    for (i in id){
      methylated <- append(methylated, 1)
    }
    Nodes <- data.frame(id = id, Base = toupper(Seq), Methylated = methylated)
    coul <- c("gray50")
  }
  else {
    methylated <- c()
    for (i in id){
      if(grepl("[AUCG]", Seq[i])){
        methylated <- append(methylated, 2)
      }
      else{
        methylated <- append(methylated, 1)
      }
    }
    methylated[Pos_met] <- 3
    Nodes <- data.frame(id = id, Base = toupper(Seq), Methylated = methylated)
    coul <- c("gray50", "tomato", "gold")
  }
  
  net <- graph_from_data_frame(d=Edges, vertices=Nodes, directed=F)
  V(net)$color <- coul[V(net)$Methylated]
  #V(net)$frame.color <- "red"
  l <- layout_with_fr(net)
  l <- norm_coords(l, ymin=-2, ymax=2, xmin=-5, xmax=5)
  plot(net, vertex.label = Nodes$Base, layout = l, rescale=T)
  if(Pos_met !=0){
    legend(x=-1.5, y=-1.1, c("Nucleoside","Methylation motif nucleoside", "Methylated nucleoside"), pch=21,
           
           col="#777777", pt.bg=coul, pt.cex=2, cex=1, bty="n", ncol=1, title = paste("Prediction of the folding of the miRNA \n with an energy of", l_Seqp[[2]], "kcal/mole", sep = " "))
  }
  else{
    legend(x=-1.5, y=-1.1, c("No methylation"), pch=21, 
           
           col="#777777", pt.bg=coul, pt.cex=2, cex=1, bty="n", ncol=1, title = paste("Prediction of the folding of the miRNA with an energy of", l_Seqp[[2]], "kcal/mole", sep = " "))
  }
  
}

fct_update <- function(methy){
  Pos <- methy$Position
  Methyla <- methy$Methylation
  Liste <- c()
  for(i in 1:(length(Pos))){
    if (Methyla[i] == 'm7G'){
      Liste <- append(Liste, paste(Methyla[i], Pos[i], methy$`Position Start`[i], 
                                   sep = "__"))
    }
    else{
      Liste <- append(Liste, paste(Methyla[i], Pos[i], sep = "__"))
    }
  }
  return(Liste)
}