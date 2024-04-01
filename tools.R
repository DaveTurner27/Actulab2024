###
### Fonctions R pour le projet d'Actulab
###

# Paquetages
liste.paquetage <- c("ggplot2", "dplyr", "tidyverse", "factoextra", "plotly", "reshape2", "cluster", "mice", "foreign", "nnet", "gbm", "xgboost")
data.init <- data
inst <- liste.paquetage %in% installed.packages()
if(length(liste.paquetage[!inst]) > 0) install.packages(liste.paquetage[!inst])
paquetages <- lapply(liste.paquetage, require, character.only = TRUE)

# Tester la validité des primes
test.primes <- function(soumissions) {
    primes_upper_limit <- 3770
    primes_lower_limit <- 860
    primes_moy_upper <- 1666
    primes_moy_lower <- 1633

    if(!"PRIME" %in% colnames(soumissions)) stop('La colonne "PRIME" doit se trouver dans votre csv')
    if(!"POLICE_ID" %in% colnames(soumissions)) stop('La colonne "POLICE_ID" doit se trouver dans votre csv')
    if(any(duplicated(soumissions$POLICE_ID))) {
        warning(paste0('Il y a ', sum(duplicated(soumissions$POLICE_ID)), ' doublons dans "POLICE_ID". Pour chaque doublon, le premier sera conservé'))
        soumissions <- soumissions[!duplicated(soumissions$POLICE_ID), ]
    }
    if(!is.numeric(soumissions$PRIME)) stop('La colonne "PRIME" doit est numérique en R en utilisant read.csv')

    tmp <- !as.character(format(10000000 + 1:1000000, scientific = FALSE)) %in% soumissions$POLICE_ID
    if(any(tmp)) warning(paste0('Il vous manque ', sum(tmp), ' polices à votre fichier.'))
    if(any(is.na(soumissions$PRIME))) {
        keep <- !is.na(soumissions$PRIME)
        message(paste0(sum(!keep), " primes manquantes (NA), ces assurés n'auront pas de soumissions"))
        soumissions <- soumissions[keep, ]
    }

    if(any(soumissions$PRIME < primes_lower_limit)) {
        warning(paste0(sum(soumissions$PRIME < primes_lower_limit),
                       " primes sont sous le seuil minimal. Ces soumissions seront ignorées"))
        soumissions <- soumissions[soumissions$PRIME >= primes_lower_limit, ]
    }
    if(any(soumissions$PRIME > primes_upper_limit)) {
        warning(paste0(sum(soumissions$PRIME > primes_upper_limit),
                       " primes sont au-dessus du seuil maximal, ces assurés n'accepteront pas ces offres"))
        soumissions <- soumissions[soumissions$PRIME <= primes_upper_limit, ]
    }
    if(mean(soumissions$PRIME, na.rm=TRUE) < primes_moy_lower || mean(soumissions$PRIME, na.rm=TRUE) > primes_moy_upper) stop(paste0(
        "ERREUR! la moyenne de vos primes valides est de ", round(mean(soumissions$PRIME, na.rm=TRUE), 2), ". Cette moyenne doit être entre ", primes_moy_lower, " et ", primes_moy_upper, "."))

    if(mean(abs(soumissions$PRIME - round(soumissions$PRIME, 2))) > 0.001) message("Vérifiez si vous avez arrondi vos primes à 0.01$ près.")

    message("Données valides!")
}

# impact analysis graph
impact_analysis <- function(actual, ref, pred, mod.names=c("REF", "PRED"), k=4){

    # Préparation
    idx.order <- order(ref)
    n <- length(ref)
    # normalisé
    ref <- ref * mean(actual) / mean(ref)
    pred <- pred * mean(actual) / mean(pred)

    # Sort
    ref <- ref[idx.order]
    pred <- pred[idx.order]
    actual <- actual[idx.order]
    data <- data.frame(
        FIGURE=c(ref, pred, actual),
        Model=as.factor(c(rep(mod.names[1], n), rep(mod.names[2], n), rep("ACTUAL", n)))
    )

    # Cut
    data$GROUPE <- rep(cut(1:n, k, labels = F), 3)

    # Moyenne
    data <- data %>%
        group_by(Model, GROUPE) %>%
        summarise(MOYENNE = mean(FIGURE))

    # Graphique
    p <- ggplot(data, aes(x=GROUPE, y=MOYENNE, col=Model)) +
        geom_point(size=2.25) + geom_line(linewidth=1) + labs(title = "Impact analysis", y="LC", x="Group") + theme_bw()
    print(p)
}

# double lift curve graph
double_lift <- function(actual, ref, pred, mod.names=c("REF", "PRED"), k=4){

    # Préparation
    idx.order <- order(pred/ref)
    n <- length(ref)
    # normalisé
    ref <- ref * mean(actual) / mean(ref)
    pred <- pred * mean(actual) / mean(pred)

    # Sort
    ref <- ref[idx.order]
    pred <- pred[idx.order]
    actual <- actual[idx.order]
    data <- data.frame(
        FIGURE=c(ref, pred),
        Model=as.factor(c(rep(mod.names[1], n), rep(mod.names[2], n)))
    )

    # Cut
    data$GROUPE <- rep(cut(1:n, k, labels = F), 2)

    # LR
    sum_actual <- sapply(split(actual, cut(1:n, k, labels = F)), sum)
    data <- data %>%
        group_by(Model, GROUPE) %>%
        summarise(LR = sum(FIGURE))
    data$LR <- sum_actual / data$LR

    # Graphique
    p <- ggplot(data, aes(x=GROUPE, y=LR, col=Model)) +
        geom_point(size=2.25) + geom_line(linewidth=1) +
        geom_hline(yintercept = 1, linewidth=0.3, col="black", alpha=0.4)+
        labs(title = "Double lift chart", y="LR", x="Group") + theme_bw()
    print(p)
}

# Calculateur de primes
predict_df <- function(df.prep, freq_data, test_premium=FALSE){

    # Prédictions par modèle
    frequence <- predict(freq, freq_data, reshape = TRUE)

    sev <- cbind(
        exp(predict(eau, df.prep, type="response") + eau$var/2),
        exp(predict(feu, df.prep, type="response") + feu$var/2),
        exp(predict(met, df.prep, type="response") + met$var/2)
    )

    # Primes pures
    PP_sans_infl <- unlist(rowSums(frequence[, -1] * sev))
    if(test_premium)
        return(PP_sans_infl)

    # Ajout d'inflation
    ajustement <- (exp(0.063613*2) - 1) / (2*0.063613)
    PP <- PP_sans_infl * ajustement
    print(paste("Inflation", ajustement))

    # Ajout de profit
    idx <- which(PP > 3600)
    trop_risque <- which(df.prep$MONTANT_COUVERTURE > 2E6)
    idx <- sort(unique(c(idx, trop_risque)))
    print(length(idx))
    if(length(idx) != 0)
        PP[idx] <- NA
    ajust_profit <- 1664/mean(PP, na.rm=TRUE)
    GROSS_P <- ajust_profit * PP
    GROSS_P <- pmax(GROSS_P, 860)
    GROSS_P <- pmin(GROSS_P, 3770)
    print(paste("Moyenne des primes :", round(mean(GROSS_P, na.rm=TRUE), 2)))
    print(paste("Profit approx. :", round(100 * (ajust_profit - 1), 2), "%"))

    # petit graphique des primes
    hist(GROSS_P, main="Distributions des primes", xlab="Prime ($)", ylab="", breaks=100, col="lightblue", lwd=1.5)
    abline(v=mean(PP, na.rm=TRUE), lwd=3, col="blue")
    abline(v=mean(GROSS_P, na.rm=TRUE), lwd=3, col="red")

    # Retour
    data.frame(
        POLICE_ID=df.prep$POLICE_ID,
        PRIME=round(GROSS_P, 2)
    )
}

# Préparer le data raw
prepare_df <- function(df.init){

    # Initialisation
    df.init$OCCUPATION_ASSURE <- as.factor(df.init$OCCUPATION_ASSURE)
    df.init$NOMMUNI <- as.factor(df.init$NOMMUNI)
    df.init$TYPE_CHAUFFAGE_PRINCIPAL <- as.factor(df.init$TYPE_CHAUFFAGE_PRINCIPAL)
    df.init$ALIM_CHAUFFAGE_PRINCIPAL <- as.factor(df.init$ALIM_CHAUFFAGE_PRINCIPAL)
    df.init$IND_CHAUF_COMBUSTION_SECONDAIRE <- as.factor(df.init$IND_CHAUF_COMBUSTION_SECONDAIRE)
    df.init$INDICATEUR_MULTIPLE_LOGIS <- as.factor(df.init$INDICATEUR_MULTIPLE_LOGIS)

    # Occupation
    cat1 <- c("Employé de bureau - Personnel administratif", "Propriétaire d'un commerce ou d'une entreprise", "Administrateur - Cadre", "Enseignant / professeur", "Journalier", "Comptable", "Ingénieur", "Assurance (Agent, Courtier, Analyste de risque, etc)", "Informaticien", "Technicien de la santé", "Finance (agent, conseiller, courtier financier)", "Entrepreneur général", "Représentant", "Policier", "Agent-Courtier immobilier", "Médecin (incluant médecin spécialiste)", "Avocat", "Pompier", "Préposé aux bénéficiaires", "Science / Professionnel", "Technicien de recherche - laboratoire", "Ressources humaines", "Pharmacien", "Marketing - Communication (technicien, administrateur, publiciste, relationniste)", "Architecte", "Notaire", "Dentiste - Orthodontiste", "Travailleur social", "Expert en sinistre", "Arpenteur-géomètre", "Physiothérapeute", "Inspecteur", "Vétérinaire", "Estimateur - Évaluateur", "Actuaire", "Chiropraticien", "Personnel aérien - Pilote", "Optométriste", "Traducteur", "Denturologiste", "Diététicien - Nutritionniste", "Ostéopathe", "Fiscaliste", "Opticien", "Orthothérapeute", "Kinésithérapeute", "Urbaniste", "Audioprothésiste", "Économiste")
    cat2 <- c("Infirmier / infirmier auxiliaire", "Cultivateur-Fermier-Agriculteur", "Conducteur - (camions, public, taxi, livreur, etc.)", "Ouvrier de la construction", "Mécanicien", "Éducatrice en garderie", "Restauration (Serveur, Cuisinier, Plongeur, Traiteur)", "Électricien", "Menuisier - ébéniste", "Vendeur", "Coiffeur", "Entretien", "Plombier", "Dessinateur, Désigner, graphiste, imprimeur, typographe", "Commerce au détail - Employé", "Étudiant à temps plein, âgé de 24 ans ou moins", "Conducteur - Transport public", "Agronome", "Militaire", "Caissier", "Agent de sécurité", "Massothérapeute", "Ambulancier", "Horticulteur (botaniste, fleuriste, jardinier, paysagiste)", "Agent correctionnel", "Étudiant", "Boucher", "Couturier", "Artiste", "Photographe", "Agent de voyage", "Réceptionniste", "Pâtissier - Boulanger", "Bijoutier - Horloger", "Acupuncteur", "Personnalité publique (radio, télévision , journaux, etc)", "Préventionniste", "Bibliothécaire", "Libraire", "Pêcheur", "Cordonnier", "Athlète professionnel")
    autre <- c("Autre occupation", "Aucune occupation", "Occupation inconnue", "Personne au foyer")
    levels(df.init$OCCUPATION_ASSURE)[levels(df.init$OCCUPATION_ASSURE) %in% autre] <- "autre"
    levels(df.init$OCCUPATION_ASSURE)[levels(df.init$OCCUPATION_ASSURE) %in% cat1] <- "cat1"
    levels(df.init$OCCUPATION_ASSURE)[levels(df.init$OCCUPATION_ASSURE) %in% cat2] <- "cat2"
    df.init$OCCUPATION_ASSURE <- droplevels(df.init$OCCUPATION_ASSURE)
    df.init$OCCUPATION_ASSURE <- relevel(df.init$OCCUPATION_ASSURE, ref = "autre")

    # COTE_CREDIT
    df.init <- df.init %>% mutate(
        COTE_CREDIT = case_when(
            COTE_CREDIT < 700 ~ "<700",
            COTE_CREDIT < 800 ~ "700-800",
            COTE_CREDIT >= 800 ~ ">=800",
            TRUE ~ "NonDisp"
        )
    )
    df.init$COTE_CREDIT <- as.factor(df.init$COTE_CREDIT)

    # Expérience
    df.init <- df.init %>% mutate(
        NB_SINISTRE_0_5ANS = case_when(
            NB_SINISTRE_0_5ANS < 2 ~ as.character(NB_SINISTRE_0_5ANS),
            TRUE ~ ">=2"
        )
    )
    df.init$NB_SINISTRE_0_5ANS <- as.factor(df.init$NB_SINISTRE_0_5ANS)

    # Municipalité
    muni <- c("AUTRE", "Gatineau", "Laval", "Lévis", "Longueuil", "Montréal", "Québec", "Saguenay", "Sherbrooke", "Trois-Rivières")
    levels(df.init$NOMMUNI)[! (levels(df.init$NOMMUNI) %in% muni)] <- "AUTRE"
    df.init$NOMMUNI <- droplevels(df.init$NOMMUNI)

    # Nombre chauff.
    df.init <- df.init %>% mutate(
        NOMBRE_CHAUFFAGES_EMPLACEMENT = case_when(
            NOMBRE_CHAUFFAGES_EMPLACEMENT < 3 ~ as.character(NOMBRE_CHAUFFAGES_EMPLACEMENT),
            TRUE ~ ">=3"
        )
    )
    df.init$NOMBRE_CHAUFFAGES_EMPLACEMENT <- as.factor(df.init$NOMBRE_CHAUFFAGES_EMPLACEMENT)

    # TYPE_CHAUFFAGE_PRINCIPAL
    type <- c("Central", "Plinthes_fixes", "Thermopompe")
    levels(df.init$TYPE_CHAUFFAGE_PRINCIPAL)[! (levels(df.init$TYPE_CHAUFFAGE_PRINCIPAL) %in% type)] <- "AUTRE"
    df.init$TYPE_CHAUFFAGE_PRINCIPAL <- droplevels(df.init$TYPE_CHAUFFAGE_PRINCIPAL)

    # ALIM_CHAUFFAGE_PRINCIPAL
    alim <- c("Bois", "Electricite", "Mazout")
    levels(df.init$ALIM_CHAUFFAGE_PRINCIPAL)[! (levels(df.init$ALIM_CHAUFFAGE_PRINCIPAL) %in% alim)] <- "NON_ELECTRIQUE"
    df.init$ALIM_CHAUFFAGE_PRINCIPAL <- droplevels(df.init$ALIM_CHAUFFAGE_PRINCIPAL)
    df.init$ALIM_CHAUFFAGE_PRINCIPAL <- relevel(df.init$ALIM_CHAUFFAGE_PRINCIPAL, ref = "Bois")

    # Retourner le df pret
    df.init
}
