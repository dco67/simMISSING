# ============================================================================
#  simManquant — Application pédagogique sur les DONNÉES MANQUANTES
#  Mécanismes (MCAR, MAR, NMAR), impacts, méthodes de traitement
#  Auteur : Daniel Coulombe  |  dcstats.org
#  Paquets : shiny, bslib, mice, ggplot2, MASS, DT
# ============================================================================

library(shiny)
library(bslib)
library(mice)
library(ggplot2)
library(MASS)
library(DT)

# Éviter le masquage de validate()/need() par jsonlite (attaché par certains environnements)
validate <- shiny::validate
need     <- shiny::need

# ----------------------------------------------------------------------------
#  PALETTE ET THÈME
# ----------------------------------------------------------------------------
COUL_OBS <- "#2c7fb8"   # observé
COUL_NA  <- "#d95f0e"   # manquant
COUL_ACC <- "#31a354"   # accent vérité

theme_app <- bs_theme(
  version = 5,
  bootswatch = "flatly",
  primary = "#2c7fb8",
  base_font = font_google("Source Sans Pro"),
  heading_font = font_google("Source Sans Pro")
)

theme_gg <- function() {
  theme_minimal(base_size = 14) +
    theme(panel.grid.minor = element_blank(),
          plot.title = element_text(face = "bold"),
          legend.position = "bottom")
}

# Tableaux DT en français (sans dépendance à un CDN)
options(DT.options = list(language = list(
  search = "Rechercher :", lengthMenu = "Afficher _MENU_ lignes",
  info = "Lignes _START_ à _END_ sur _TOTAL_", infoEmpty = "Aucune ligne",
  infoFiltered = "(filtré sur _MAX_ lignes)", zeroRecords = "Aucun résultat",
  emptyTable = "Aucune donnée", decimal = ",",
  paginate = list(first = "Premier", previous = "Précédent", `next` = "Suivant", last = "Dernier"))))

# Nombres à la française : virgule décimale (textes, tableaux, axes)
fmt_fr <- function(x, d = 2, signe = FALSE) {
  out <- formatC(x, format = "f", digits = d, decimal.mark = ",",
                 flag = if (signe) "+" else "")
  out[is.na(x)] <- "—"
  out
}
pct_fr <- function(x, d = 1) paste0(fmt_fr(x, d), " %")   # x déjà en %
p_fr <- function(p) ifelse(is.na(p), "—",
                           ifelse(p < 0.001, "< 0,001", fmt_fr(p, 3)))
lab_fr <- function(x) format(x, decimal.mark = ",", big.mark = " ",
                             scientific = FALSE, trim = TRUE, drop0trailing = TRUE)
sx_fr <- function(...) scale_x_continuous(labels = lab_fr, ...)
sy_fr <- function(...) scale_y_continuous(labels = lab_fr, ...)

# ============================================================================
#  1. FONCTIONS DE CALCUL (moteur)  — testables hors de Shiny
# ============================================================================

# ---- 1.1 Simulation de données complètes multinormales corrélées ----------
simuler_complet <- function(N = 300L, p = 4L, rho = 0.5,
                            moyennes = NULL, ecarts = NULL, graine = 2024L) {
  N <- as.integer(N); p <- as.integer(p)
  if (is.null(moyennes)) {
    moyennes <- rep(0, p)
  }
  if (is.null(ecarts)) {
    ecarts <- rep(1, p)
  }
  moyennes <- rep_len(moyennes, p)
  ecarts   <- rep_len(ecarts, p)
  set.seed(graine)
  R <- matrix(rho, p, p)
  diag(R) <- 1
  D <- diag(ecarts, p, p)
  Sigma <- D %*% R %*% D
  X <- MASS::mvrnorm(N, mu = moyennes, Sigma = Sigma)
  X <- as.data.frame(X)
  names(X) <- paste0("X", seq_len(p))
  X
}

# ---- 1.2 Amputation selon un mécanisme (MCAR / MAR / MNAR) -----------------
# NB. Dans mice::ampute, prop = proportion de CAS incomplets (bycases = TRUE),
#     et non le pourcentage de cellules manquantes. On renvoie les deux.
amputer <- function(df, prop = 0.30, mech = "MCAR", graine = 2024L,
                    var_cible = NULL, var_pilote = NULL, direction = "RIGHT",
                    patron = "ponctuel") {
  mech <- match.arg(mech, c("MCAR", "MAR", "MNAR"))
  vn <- names(df); p <- length(vn)
  # Garde-fous : la cible peut être une ou plusieurs variables.
  cibles <- unique(var_cible[var_cible %in% vn])
  if (!length(cibles)) cibles <- NULL
  if (!is.null(var_pilote) && !(var_pilote %in% vn)) var_pilote <- NULL
  patron <- match.arg(patron, c("ponctuel", "bloc"))
  set.seed(graine)
  if (is.null(cibles)) {
    # mode diffus (rétrocompatible) : schéma ampute standard
    res <- mice::ampute(df, prop = prop, mech = mech, bycases = TRUE)
  } else {
    # Un patron par variable = omissions réparties; un seul patron = bloc omis.
    if (length(cibles) == 1L || identical(patron, "bloc")) {
      patt <- matrix(as.integer(!(vn %in% cibles)), nrow = 1,
                     dimnames = list(NULL, vn))
    } else {
      patt <- do.call(rbind, lapply(cibles, function(v)
        setNames(as.integer(vn != v), vn)))
      colnames(patt) <- vn
    }
    freq <- rep(1 / nrow(patt), nrow(patt))
    if (mech == "MCAR") {
      res <- mice::ampute(df, prop = prop, mech = "MCAR", bycases = TRUE,
                          patterns = patt, freq = freq)
    } else {
      w <- matrix(0, nrow = nrow(patt), ncol = p,
                  dimnames = list(NULL, vn))
      if (mech == "MAR") {
        # Un prédicteur MAR doit rester observé dans le jeu incomplet.
        if (is.null(var_pilote) || var_pilote %in% cibles) {
          var_pilote <- setdiff(vn, cibles)[1]
        }
        if (is.null(var_pilote) || is.na(var_pilote)) {
          stop("Sous MAR, laissez au moins une variable observée comme prédicteur du manque.")
        }
        w[, match(var_pilote, vn)] <- 1
      } else {                       # MNAR : piloté par la cible elle-même
        if (nrow(patt) == 1L) {
          w[1, match(cibles, vn)] <- 1 / length(cibles)
        } else {
          for (j in seq_len(nrow(patt))) w[j, match(cibles[j], vn)] <- 1
        }
        var_pilote <- cibles
      }
      res <- mice::ampute(df, prop = prop, mech = mech, bycases = TRUE,
                          patterns = patt, freq = freq, weights = w,
                          type = rep(direction, nrow(patt)))
    }
  }
  amp <- res$amp
  list(amp = amp,
       pct_cellules = mean(is.na(amp)),
       pct_cas = mean(!complete.cases(amp)),
       var_cible = cibles,
       var_pilote = if (!is.null(cibles) && mech != "MCAR") var_pilote else NULL)
}

# ---- 1.3 Comparaison des méthodes : estimation d'une MOYENNE ---------------
# ---- FIML (maximum de vraisemblance à information complète), moyenne --------
fiml_moyenne <- function(ampute, var) {
  if (!requireNamespace("lavaan", quietly = TRUE)) return(NA_real_)
  vars <- names(ampute)
  moyennes <- paste0(vars, " ~ 1")                       # modèle saturé
  covs <- if (length(vars) >= 2)
    apply(utils::combn(vars, 2), 2,
          function(p) paste0(p[1], " ~~ ", p[2])) else character(0)
  modele <- paste(c(moyennes, covs), collapse = "\n")
  fit <- tryCatch(
    suppressWarnings(lavaan::lavaan(modele, data = ampute, missing = "fiml",
                                    auto.var = TRUE)),
    error = function(e) NULL)
  if (is.null(fit) || !lavaan::lavInspect(fit, "converged")) return(NA_real_)
  co <- lavaan::coef(fit); key <- paste0(var, "~1")
  if (key %in% names(co)) unname(co[key]) else NA_real_
}

# ---- FIML : coefficients de régression -------------------------------------
fiml_regression <- function(ampute, cible, predicteurs, termes) {
  if (!requireNamespace("lavaan", quietly = TRUE)) return(rep(NA_real_, length(termes)))
  modele <- paste0(cible, " ~ ", paste(predicteurs, collapse = " + "))
  fit <- tryCatch(
    suppressWarnings(lavaan::sem(modele, data = ampute, missing = "fiml",
                                 fixed.x = FALSE, meanstructure = TRUE)),
    error = function(e) NULL)
  if (is.null(fit) || !lavaan::lavInspect(fit, "converged"))
    return(rep(NA_real_, length(termes)))
  co <- lavaan::coef(fit)
  vapply(termes, function(t) {          # (Intercept) -> cible~1 ; Xk -> cible~Xk
    key <- if (t == "(Intercept)") paste0(cible, "~1") else paste0(cible, "~", t)
    if (key %in% names(co)) unname(co[key]) else NA_real_
  }, numeric(1))
}

# ---- FIML : corrélation (covariance ajustée, standardisée) ------------------
fiml_correlation <- function(ampute, v1, v2) {
  if (!requireNamespace("lavaan", quietly = TRUE)) return(NA_real_)
  vars <- names(ampute)
  moyennes <- paste0(vars, " ~ 1")
  covs <- if (length(vars) >= 2)
    apply(utils::combn(vars, 2), 2,
          function(p) paste0(p[1], " ~~ ", p[2])) else character(0)
  modele <- paste(c(moyennes, covs), collapse = "\n")
  fit <- tryCatch(
    suppressWarnings(lavaan::lavaan(modele, data = ampute, missing = "fiml",
                                    auto.var = TRUE)),
    error = function(e) NULL)
  if (is.null(fit) || !lavaan::lavInspect(fit, "converged")) return(NA_real_)
  S <- lavaan::lavInspect(fit, "cov.ov")
  if (!all(c(v1, v2) %in% rownames(S))) return(NA_real_)
  unname(S[v1, v2] / sqrt(S[v1, v1] * S[v2, v2]))
}

# ---- Imputation unique par régression (m = 1) : déterministe ou stochastique --
#  déterministe  : valeur prédite par la régression (mice « norm.predict »)
#  stochastique  : valeur prédite + résidu aléatoire (mice « norm.nob »)
imputer_regression <- function(ampute, type = c("deterministe", "stochastique"),
                               graine = 2024L) {
  type <- match.arg(type)
  set.seed(graine)
  meth <- if (type == "deterministe") "norm.predict" else "norm.nob"
  mice::complete(mice::mice(ampute, m = 1L, method = meth, printFlag = FALSE), 1)
}

imputer_moyenne <- function(ampute) {
  for (v in names(ampute)) ampute[[v]][is.na(ampute[[v]])] <- mean(ampute[[v]], na.rm = TRUE)
  ampute
}

METHODES_COMP <- c("Vérité (données complètes)", "Cas complets (listwise)",
                   "Valeurs disponibles (pairwise)", "Imputation par la moyenne",
                   "Imputation par régression (déterministe)",
                   "Imputation par régression (stochastique)",
                   "Imputation multiple (mise en commun)",
                   "Maximum de vraisemblance (FIML)")

comparer_moyenne <- function(complet, ampute, var, m = 10L, graine = 2024L) {
  y_vrai <- mean(complet[[var]])
  cc <- ampute[complete.cases(ampute), ]
  y_cc <- mean(cc[[var]])                                   # cas complets (listwise)
  y_disp <- mean(ampute[[var]], na.rm = TRUE)               # valeurs disponibles
  y_moy <- mean(imputer_moyenne(ampute)[[var]])             # imputation par la moyenne
  y_det <- mean(imputer_regression(ampute, "deterministe", graine)[[var]])
  y_sto <- mean(imputer_regression(ampute, "stochastique", graine)[[var]])
  set.seed(graine)                                          # imputation multiple
  mi <- mice::mice(ampute, m = as.integer(m), printFlag = FALSE)
  ests <- vapply(seq_len(as.integer(m)),
                 function(k) mean(mice::complete(mi, k)[[var]]), numeric(1))
  y_mi <- mean(ests)
  y_fiml <- fiml_moyenne(ampute, var)
  out <- data.frame(Methode = METHODES_COMP,
                    Estimation = c(y_vrai, y_cc, y_disp, y_moy, y_det, y_sto, y_mi, y_fiml),
                    stringsAsFactors = FALSE)
  out$Biais <- out$Estimation - y_vrai
  attr(out, "mids") <- mi
  out
}

# ---- 1.4 Comparaison des méthodes : CORRÉLATION entre deux variables -------
comparer_correlation <- function(complet, ampute, v1, v2, m = 10L, graine = 2024L) {
  r_of <- function(d) cor(d[[v1]], d[[v2]])
  r_vrai <- r_of(complet)
  r_cc <- r_of(ampute[complete.cases(ampute), ])
  r_disp <- cor(ampute[[v1]], ampute[[v2]], use = "pairwise.complete.obs")
  r_moy <- r_of(imputer_moyenne(ampute))
  r_det <- r_of(imputer_regression(ampute, "deterministe", graine))
  r_sto <- r_of(imputer_regression(ampute, "stochastique", graine))
  set.seed(graine)
  mi <- mice::mice(ampute, m = as.integer(m), printFlag = FALSE)
  rs <- vapply(seq_len(as.integer(m)), function(k) r_of(mice::complete(mi, k)), numeric(1))
  r_mi <- tanh(mean(atanh(rs)))          # mise en commun sur l'échelle z de Fisher
  r_fiml <- fiml_correlation(ampute, v1, v2)
  out <- data.frame(Methode = METHODES_COMP,
                    Estimation = c(r_vrai, r_cc, r_disp, r_moy, r_det, r_sto, r_mi, r_fiml),
                    stringsAsFactors = FALSE)
  out$Biais <- out$Estimation - r_vrai
  attr(out, "mids") <- mi
  out
}

# ---- 1.5 Comparaison des méthodes : RÉGRESSION -----------------------------
COLS_REG <- c(Verite = "Vérité", Cas_complets = "Cas complets",
              Imp_moyenne = "Imp. moyenne", Imp_reg_det = "Imp. régr. déterministe",
              Imp_reg_sto = "Imp. régr. stochastique", Imp_multiple = "Imp. multiple",
              FIML = "FIML")

comparer_regression <- function(complet, ampute, cible, predicteurs,
                                m = 10L, graine = 2024L) {
  form <- as.formula(paste(cible, "~", paste(predicteurs, collapse = " + ")))
  b_vrai <- coef(lm(form, data = complet))
  trm <- names(b_vrai)
  b_cc   <- coef(lm(form, data = ampute))          # lm applique listwise
  b_moy  <- coef(lm(form, data = imputer_moyenne(ampute)))
  b_det  <- coef(lm(form, data = imputer_regression(ampute, "deterministe", graine)))
  b_sto  <- coef(lm(form, data = imputer_regression(ampute, "stochastique", graine)))
  set.seed(graine)                                 # imputation multiple + Rubin
  mi <- mice::mice(ampute, m = as.integer(m), printFlag = FALSE)
  fits <- lapply(seq_len(as.integer(m)),
                 function(k) lm(form, data = mice::complete(mi, k)))
  pooled <- summary(mice::pool(mice::as.mira(fits)))
  b_mi <- setNames(pooled$estimate, as.character(pooled$term))[trm]
  out <- data.frame(
    Terme = trm,
    Verite = as.numeric(b_vrai),
    Cas_complets = as.numeric(b_cc[trm]),
    Imp_moyenne = as.numeric(b_moy[trm]),
    Imp_reg_det = as.numeric(b_det[trm]),
    Imp_reg_sto = as.numeric(b_sto[trm]),
    Imp_multiple = as.numeric(b_mi),
    FIML = fiml_regression(ampute, cible, predicteurs, trm),
    stringsAsFactors = FALSE, check.names = FALSE
  )
  attr(out, "mids") <- mi
  out
}

# ---- 1.6 Effet sur la dispersion (écart-type) -----------------------------
#  mids : imputation multiple déjà calculée (on réutilise les m jeux)
effet_dispersion <- function(complet, ampute, var, mids = NULL, graine = 2024L) {
  sd_mi <- if (!is.null(mids)) mean(vapply(seq_len(mids$m),
             function(k) sd(mice::complete(mids, k)[[var]]), numeric(1))) else NA_real_
  data.frame(
    Situation = c("Données complètes (vérité)", "Valeurs disponibles",
                  "Après imputation par la moyenne",
                  "Après régression déterministe",
                  "Après régression stochastique",
                  "Imputation multiple (moyenne des m jeux)"),
    Ecart_type = c(sd(complet[[var]]), sd(ampute[[var]], na.rm = TRUE),
                   sd(imputer_moyenne(ampute)[[var]]),
                   sd(imputer_regression(ampute, "deterministe", graine)[[var]]),
                   sd(imputer_regression(ampute, "stochastique", graine)[[var]]),
                   sd_mi),
    stringsAsFactors = FALSE
  )
}

# ============================================================================
#  2. FONCTIONS GRAPHIQUES
# ============================================================================

# ---- 2.1 Carte des valeurs manquantes (heatmap cellule par cellule) -------
carte_manquants <- function(df, max_lignes = 200L) {
  n <- nrow(df)
  idx <- if (n > max_lignes) {
    round(seq(1, n, length.out = as.integer(max_lignes)))
  } else {
    seq_len(n)
  }
  d <- df[idx, , drop = FALSE]
  long <- data.frame(
    ligne = rep(seq_along(idx), times = ncol(d)),
    variable = rep(names(d), each = nrow(d)),
    manquant = as.vector(is.na(as.matrix(d)))
  )
  long$variable <- factor(long$variable, levels = names(d))
  ggplot(long, aes(x = variable, y = ligne, fill = manquant)) +
    geom_tile() +
    scale_fill_manual(values = c(`FALSE` = COUL_OBS, `TRUE` = COUL_NA),
                      labels = c("Observé", "Manquant"), name = NULL) +
    scale_y_reverse() +
    labs(title = "Carte des valeurs manquantes",
         subtitle = sprintf("%d observations affichées (sur %d)", length(idx), n),
         x = "Variable", y = "Observation") +
    theme_gg()
}

# ---- 2.2 Fréquence des patrons de manquants -------------------------------
graphe_patrons <- function(df) {
  mp <- mice::md.pattern(df, plot = FALSE)
  # dernière ligne = nb de NA par variable ; dernière colonne = nb de NA du patron
  freq <- as.integer(rownames(mp)[-nrow(mp)])
  n_manq <- as.integer(mp[-nrow(mp), ncol(mp)])
  d <- data.frame(patron = factor(seq_along(freq)),
                  freq = freq, n_manq = factor(n_manq, levels = sort(unique(n_manq))))
  d <- d[order(-d$freq), ]
  d$patron <- factor(d$patron, levels = d$patron)
  nl <- length(levels(d$n_manq))
  pal <- if (nl == 1) COUL_OBS else
    grDevices::colorRampPalette(c(COUL_OBS, "#f0b27a", COUL_NA))(nl)
  ggplot(d, aes(x = patron, y = freq, fill = n_manq)) +
    geom_col() +
    geom_text(aes(label = freq), vjust = -0.3, size = 3.8) +
    scale_fill_manual(values = setNames(pal, levels(d$n_manq)),
                      name = "Variables manquantes dans le patron") +
    sy_fr() +
    labs(title = "Fréquence des patrons de manquants",
         subtitle = "Numéro du patron = ligne du tableau des patrons",
         x = "Patron", y = "Nombre d'observations") +
    theme_gg()
}

# ---- 2.3 % de manquants par variable --------------------------------------
graphe_pct_variable <- function(df) {
  pct <- vapply(df, function(x) 100 * mean(is.na(x)), numeric(1))
  d <- data.frame(variable = factor(names(pct), levels = names(pct)),
                  pct = as.numeric(pct))
  ggplot(d, aes(x = variable, y = pct)) +
    geom_col(fill = COUL_NA) +
    geom_text(aes(label = pct_fr(pct)), vjust = -0.3, size = 3.8) +
    sy_fr() +
    labs(title = "Pourcentage de valeurs manquantes par variable",
         x = "Variable", y = "% manquant") +
    theme_gg()
}

# ---- 2.4 Distribution observé vs imputé (illustre le rétrécissement) ------
graphe_imputation <- function(ampute, var, methode = "moyenne", graine = 2024L,
                              complet = NULL) {
  obs <- ampute[[var]][!is.na(ampute[[var]])]
  manq_idx <- is.na(ampute[[var]])
  if (methode == "moyenne") {
    imp_vals <- rep(mean(obs), sum(manq_idx))
    titre <- "Imputation par la moyenne"
  } else if (methode == "reg_det") {
    imp_vals <- imputer_regression(ampute, "deterministe", graine)[[var]][manq_idx]
    titre <- "Imputation par régression déterministe (valeur prédite)"
  } else if (methode == "regression") {
    imp_vals <- imputer_regression(ampute, "stochastique", graine)[[var]][manq_idx]
    titre <- "Imputation par régression stochastique (prédite + résidu)"
  } else {
    set.seed(graine)
    mi1 <- mice::mice(ampute, m = 1L, printFlag = FALSE)
    imp_vals <- mice::complete(mi1, 1)[[var]][manq_idx]
    titre <- "Imputation multiple : un des m jeux (appariement prédictif)"
  }
  sd_txt <- sprintf("Écart-type des valeurs imputées : %s  ·  des valeurs observées : %s",
                    fmt_fr(if (length(imp_vals) > 1) sd(imp_vals) else 0, 2), fmt_fr(sd(obs), 2))
  if (!is.null(complet))
    sd_txt <- paste0(sd_txt, sprintf("\ndes vraies valeurs effacées (connues en simulation) : %s",
                                     fmt_fr(sd(complet[[var]][manq_idx]), 2)))
  d <- data.frame(
    valeur = c(obs, imp_vals),
    type = factor(c(rep("Observé", length(obs)), rep("Imputé", length(imp_vals))),
                  levels = c("Observé", "Imputé"))
  )
  ggplot(d, aes(x = valeur, fill = type)) +
    geom_histogram(position = "identity", alpha = 0.55, bins = 30) +
    scale_fill_manual(values = c(`Observé` = COUL_OBS, `Imputé` = COUL_NA),
                      name = NULL) +
    sx_fr() +
    labs(title = paste0(titre, " — ", var),
         subtitle = sd_txt, x = var, y = "Effectif") +
    theme_gg()
}

graphe_manquants_bivarie <- function(complet, ampute, var_manq, var_x) {
  d <- data.frame(x = complet[[var_x]], y = complet[[var_manq]],
                  statut = factor(ifelse(is.na(ampute[[var_manq]]), "Manquant", "Observé"),
                                  levels = c("Observé", "Manquant")))
  ggplot(d, aes(x, y, colour = statut)) +
    geom_point(alpha = 0.75, size = 2.6) +
    scale_colour_manual(values = c(`Observé` = COUL_OBS, `Manquant` = COUL_NA),
                        name = NULL) +
    sx_fr() + sy_fr() +
    labs(title = sprintf("Présents et manquants de %s selon %s", var_manq, var_x),
         subtitle = sprintf("Chaque point est coloré selon que %s est observée ou manquante", var_manq),
         x = var_x, y = sprintf("%s (valeur vraie, connue en simulation)", var_manq)) +
    theme_gg()
}

# ---- 2.5 Tableau des patrons et diagrammes du biais -------------------------
table_patrons <- function(df) {
  mp <- mice::md.pattern(df, plot = FALSE)
  k <- nrow(mp); v <- ncol(mp)
  freq <- as.integer(rownames(mp)[-k])
  corps <- as.data.frame(mp[-k, -v, drop = FALSE])
  nb_manq <- as.integer(mp[-k, v])
  out <- cbind(Effectif = freq, corps, `Variables manquantes` = nb_manq)
  rownames(out) <- NULL
  out
}

graphe_biais <- function(tab) {
  d <- tab[tab$Methode != "Vérité (données complètes)", ]
  d$Methode <- factor(d$Methode, levels = rev(d$Methode))
  ggplot(d, aes(x = Methode, y = Biais, fill = Biais > 0)) +
    geom_col() +
    geom_hline(yintercept = 0, linewidth = 0.4) +
    coord_flip() +
    scale_fill_manual(values = c(`TRUE` = COUL_NA, `FALSE` = COUL_OBS),
                      guide = "none") +
    sy_fr() +
    labs(title = "Biais (estimation − vérité)",
         subtitle = "Un seul jeu simulé : un petit écart\npeut tenir au hasard",
         x = NULL, y = "Biais") +
    theme_gg()
}

graphe_biais_regression <- function(tab) {
  meths <- intersect(c("Cas_complets", "Imp_moyenne", "Imp_reg_det", "Imp_reg_sto",
                       "Imp_multiple", "FIML"), names(tab))
  d <- do.call(rbind, lapply(meths, function(m)
    data.frame(Terme = tab$Terme, Methode = m, Biais = tab[[m]] - tab$Verite)))
  lib <- c(Cas_complets = "Cas complets", Imp_moyenne = "Imp. moyenne",
           Imp_reg_det = "Régr. déterm.", Imp_reg_sto = "Régr. stoch.",
           Imp_multiple = "Imp. multiple", FIML = "FIML")
  d$Methode <- factor(lib[d$Methode], levels = lib[meths])
  d$Terme <- factor(d$Terme, levels = tab$Terme)
  ggplot(d, aes(x = Methode, y = Biais, fill = Biais > 0)) +
    geom_col() +
    geom_hline(yintercept = 0, linewidth = 0.4) +
    facet_wrap(~ Terme, scales = "free_y") +
    scale_fill_manual(values = c(`TRUE` = COUL_NA, `FALSE` = COUL_OBS), guide = "none") +
    sy_fr() +
    labs(title = "Biais de chaque méthode, par coefficient",
         x = NULL, y = "Biais (estimation − vérité)") +
    theme_gg() +
    theme(axis.text.x = element_text(angle = 35, hjust = 1))
}

# ============================================================================
#  MOTEUR MONTE CARLO — étude de simulation (biais, efficacité, couverture)
# ============================================================================

# Vérité de référence : paramètre estimé sur un très grand échantillon complet
verite_reference <- function(cible, params, p, rho,
                             N_ref = 50000L, graine_ref = 999983L) {
  ref <- simuler_complet(N_ref, p, rho, graine = graine_ref)
  estim_ponctuel(ref, cible, params)
}

# Estimation ponctuelle scalaire sur un jeu SANS manquant
estim_ponctuel <- function(d, cible, params) {
  if (cible == "moyenne") {
    mean(d[[params$var]])
  } else if (cible == "correlation") {
    cor(d[[params$v1]], d[[params$v2]])
  } else {
    form <- reformulate(params$x, params$y)
    coef(lm(form, data = d))[[params$terme]]
  }
}

# Estimation + erreur-type + IC 95 % pour un jeu AMPUTÉ, selon la méthode
# méthodes : "cc" (cas complets/disponibles), "moyenne", "mi"
estim_ic <- function(amp, cible, params, methode, m = 10L, graine = 1L) {
  m <- as.integer(m)
  if (cible == "moyenne") {
    v <- params$var
    if (methode == "mi") {
      set.seed(graine)
      mi <- mice::mice(amp, m = m, printFlag = FALSE)
      fits <- lapply(seq_len(m),
                     function(k) lm(reformulate("1", v), data = mice::complete(mi, k)))
      s <- summary(mice::pool(mice::as.mira(fits)), conf.int = TRUE)
      return(c(est = s$estimate[1], se = s$std.error[1],
               lo = s[["2.5 %"]][1], hi = s[["97.5 %"]][1]))
    }
    if (methode == "moyenne") {
      x <- amp[[v]]; x[is.na(x)] <- mean(x, na.rm = TRUE)
    } else {
      x <- amp[[v]][!is.na(amp[[v]])]
    }
    n <- length(x); est <- mean(x); se <- sd(x) / sqrt(n); tq <- qt(0.975, n - 1)
    return(c(est = est, se = se, lo = est - tq * se, hi = est + tq * se))

  } else if (cible == "correlation") {
    v1 <- params$v1; v2 <- params$v2
    if (methode == "mi") {
      set.seed(graine)
      mi <- mice::mice(amp, m = m, printFlag = FALSE)
      rs <- vapply(seq_len(m), function(k) {
        d <- mice::complete(mi, k); cor(d[[v1]], d[[v2]])
      }, numeric(1))
      n <- nrow(amp)
      zk <- atanh(rs); uk <- 1 / (n - 3)
      qbar <- mean(zk); B <- if (m > 1) var(zk) else 0
      Tvar <- mean(uk) + (1 + 1 / m) * B
      sez <- sqrt(Tvar)
      r_est <- tanh(qbar)
      # IC calculé sur l'échelle z; erreur-type ramenée à l'échelle r (méthode delta)
      return(c(est = r_est, se = (1 - r_est^2) * sez,
               lo = tanh(qbar - 1.96 * sez), hi = tanh(qbar + 1.96 * sez)))
    }
    if (methode == "moyenne") {
      d <- amp
      for (v in c(v1, v2)) d[[v]][is.na(d[[v]])] <- mean(d[[v]], na.rm = TRUE)
      r <- cor(d[[v1]], d[[v2]]); n <- nrow(d)
    } else {
      ok <- stats::complete.cases(amp[, c(v1, v2)])
      r <- cor(amp[[v1]][ok], amp[[v2]][ok]); n <- sum(ok)
    }
    z <- atanh(r); sez <- 1 / sqrt(n - 3)
    return(c(est = r, se = (1 - r^2) * sez,
             lo = tanh(z - 1.96 * sez), hi = tanh(z + 1.96 * sez)))

  } else { # regression
    form <- reformulate(params$x, params$y); trm <- params$terme
    if (methode == "mi") {
      set.seed(graine)
      mi <- mice::mice(amp, m = m, printFlag = FALSE)
      fits <- lapply(seq_len(m), function(k) lm(form, data = mice::complete(mi, k)))
      s <- summary(mice::pool(mice::as.mira(fits)), conf.int = TRUE)
      i <- match(trm, as.character(s$term))
      return(c(est = s$estimate[i], se = s$std.error[i],
               lo = s[["2.5 %"]][i], hi = s[["97.5 %"]][i]))
    }
    if (methode == "moyenne") {
      d <- amp
      for (v in names(d)) d[[v]][is.na(d[[v]])] <- mean(d[[v]], na.rm = TRUE)
    } else {
      d <- amp   # lm applique listwise
    }
    fit <- lm(form, data = d); co <- summary(fit)$coefficients
    est <- co[trm, 1]; se <- co[trm, 2]
    df <- fit$df.residual; tq <- qt(0.975, df)
    return(c(est = est, se = se, lo = est - tq * se, hi = est + tq * se))
  }
}

# Étude Monte Carlo complète : R réplications, graine variée par réplication
etude_monte_carlo <- function(R, N, p, rho, prop, mech, cible, params,
                              m = 10L, graine_base = 2024L,
                              methodes = c("cc", "moyenne", "mi"),
                              var_cible = NULL, var_pilote = NULL,
                              direction = "RIGHT", patron = "ponctuel") {
  verite <- verite_reference(cible, params, p, rho)
  lignes <- vector("list", R * length(methodes)); k <- 0L
  for (r in seq_len(R)) {
    g <- as.integer(graine_base + r)          # <-- graine variée : essentiel
    cplt <- simuler_complet(N, p, rho, graine = g)
    amp  <- amputer(cplt, prop = prop, mech = mech, graine = g,
                    var_cible = var_cible, var_pilote = var_pilote,
                    direction = direction, patron = patron)$amp
    for (meth in methodes) {
      ic <- estim_ic(amp, cible, params, meth, m = m, graine = g)
      k <- k + 1L
      lignes[[k]] <- data.frame(
        replication = r, methode = meth,
        est = ic[["est"]], se = ic[["se"]],
        lo = ic[["lo"]], hi = ic[["hi"]],
        couvre = (verite >= ic[["lo"]]) & (verite <= ic[["hi"]]),
        stringsAsFactors = FALSE)
    }
  }
  brut <- do.call(rbind, lignes)
  attr(brut, "verite") <- verite
  brut
}

# Tableau de synthèse : biais, écart-type empirique, erreur-type moyenne, couverture
synthese_mc <- function(brut) {
  verite <- attr(brut, "verite")
  etiq <- c(cc = "Cas complets / disponibles",
            moyenne = "Imputation par la moyenne",
            mi = "Imputation multiple (Rubin)")
  agg <- lapply(split(brut, brut$methode), function(d) {
    data.frame(
      Methode = etiq[d$methode[1]],
      Biais_moyen = mean(d$est) - verite,
      Ecart_type_empirique = sd(d$est),
      Erreur_type_moyenne = mean(d$se),
      Couverture_pct = 100 * mean(d$couvre),
      stringsAsFactors = FALSE)
  })
  out <- do.call(rbind, agg)
  out <- out[match(names(etiq), names(agg)), ]
  rownames(out) <- NULL
  out[!is.na(out$Methode), ]
}

# ---- Graphiques Monte Carlo -------------------------------------------------
etiq_meth <- c(cc = "Cas complets / disponibles",
               moyenne = "Imputation par la moyenne",
               mi = "Imputation multiple (Rubin)")
coul_meth <- c("Cas complets / disponibles" = "#2c7fb8",
               "Imputation par la moyenne"   = "#d95f0e",
               "Imputation multiple (Rubin)" = "#31a354")

graphe_mc_distribution <- function(brut) {
  verite <- attr(brut, "verite")
  d <- brut
  d$Methode <- factor(etiq_meth[d$methode], levels = etiq_meth)
  ggplot(d, aes(x = est, fill = Methode)) +
    geom_density(alpha = 0.45, colour = NA) +
    geom_vline(xintercept = verite, linetype = "dashed", linewidth = 0.6) +
    annotate("text", x = verite, y = 0, label = "  vérité", hjust = 0,
             vjust = -0.5, size = 3.4) +
    scale_fill_manual(values = coul_meth, name = NULL) +
    sx_fr() + sy_fr() +
    guides(fill = guide_legend(nrow = 2)) +
    labs(title = "Distribution des estimations",
         subtitle = sprintf("%d réplications\nDécalage par rapport à la vérité = biais; largeur = efficacité",
                            max(brut$replication)),
         x = "Estimation", y = "Densité") +
    theme_gg()
}

graphe_mc_couverture <- function(brut) {
  s <- synthese_mc(brut)
  s$Methode <- factor(s$Methode, levels = etiq_meth)
  ggplot(s, aes(x = Methode, y = Couverture_pct, fill = Methode)) +
    geom_col(width = 0.6) +
    geom_hline(yintercept = 95, linetype = "dashed", linewidth = 0.6) +
    geom_text(aes(label = pct_fr(Couverture_pct, 0)),
              vjust = -0.4, size = 4) +
    annotate("text", x = 3.45, y = 99, label = "cible : 95 %", hjust = 1, size = 3.6) +
    scale_fill_manual(values = coul_meth, guide = "none") +
    scale_x_discrete(labels = function(x) gsub(" (", "\n(", gsub(" / ", " /\n", gsub("Imputation ", "Imputation\n", x)), fixed = TRUE)) +
    coord_cartesian(ylim = c(0, 102)) +
    labs(title = "Couverture des IC à 95 %",
         subtitle = "Sous 95 % : intervalles trop étroits ou décalés",
         x = NULL, y = "Couverture (%)") +
    theme_gg()
}

graphe_caterpillar <- function(brut, n_montre = 25) {
  verite <- attr(brut, "verite")
  n_montre <- min(n_montre, max(brut$replication))
  d <- brut[brut$replication <= n_montre, ]
  d$Methode <- factor(etiq_meth[d$methode], levels = etiq_meth)
  ggplot(d, aes(y = replication)) +
    geom_vline(xintercept = verite, linewidth = 0.7, colour = "#333333") +
    geom_segment(aes(x = lo, xend = hi, yend = replication, colour = couvre),
                 linewidth = 0.5) +
    geom_point(aes(x = est, colour = couvre), size = 1) +
    facet_wrap(~ Methode) +
    scale_colour_manual(values = c(`TRUE` = COUL_OBS, `FALSE` = "#d7191c"),
                        labels = c(`TRUE` = "couvre la vérité", `FALSE` = "rate la vérité"),
                        name = NULL) +
    scale_y_reverse() + sx_fr() +
    labs(title = sprintf("Intervalles de confiance à 95 %% (%d premières réplications)", n_montre),
         subtitle = "Trait vertical = vraie valeur ; en rouge, les IC qui la ratent",
         x = "Estimation et intervalle de confiance", y = "Réplication") +
    theme_gg()
}

# ============================================================================
#  Test de Little (MCAR) — EM pour (mu, Sigma) + statistique du chi-carré
#  Little, R.J.A. (1988). JASA, 83(404), 1198–1202.  Sans dépendance externe.
# ============================================================================

# EM (multinormale, données manquantes) — E-step groupé par patron (rapide)
em_mvn <- function(X, max_iter = 300L, tol = 1e-6) {
  X <- as.matrix(X); n <- nrow(X); p <- ncol(X)
  M <- is.na(X)
  key <- apply(M, 1, function(r) paste(as.integer(!r), collapse = ""))
  groupes <- split(seq_len(n), key)
  mu <- colMeans(X, na.rm = TRUE)
  S <- stats::cov(X, use = "pairwise.complete.obs")
  S[is.na(S)] <- 0; d <- diag(S); d[is.na(d) | d <= 0] <- 1; diag(S) <- d
  rid <- 1e-8
  it <- 0L
  repeat {
    it <- it + 1L; mu_old <- mu; S_old <- S
    T1 <- numeric(p); T2 <- matrix(0, p, p)
    for (g in groupes) {
      obs <- !M[g[1], ]; mis <- !obs; nj <- length(g)
      Xg <- X[g, , drop = FALSE]
      if (all(obs)) {
        Xh <- Xg
      } else {
        So  <- S[obs, obs, drop = FALSE]
        Som <- S[obs, mis, drop = FALSE]
        Soi <- solve(So + diag(rid, nrow(So)))
        centre <- sweep(Xg[, obs, drop = FALSE], 2, mu[obs], "-")
        B <- sweep(centre %*% Soi %*% Som, 2, mu[mis], "+")   # nj x |mis|
        Xh <- matrix(0, nj, p); Xh[, obs] <- Xg[, obs]; Xh[, mis] <- B
        Cmm <- S[mis, mis, drop = FALSE] - t(Som) %*% Soi %*% Som
        corr <- matrix(0, p, p); corr[mis, mis] <- nj * Cmm
        T2 <- T2 + corr
      }
      T1 <- T1 + colSums(Xh)
      T2 <- T2 + t(Xh) %*% Xh
    }
    mu <- T1 / n
    S <- T2 / n - tcrossprod(mu); S <- (S + t(S)) / 2
    if (max(abs(mu - mu_old)) < tol && max(abs(S - S_old)) < tol) break
    if (it >= max_iter) break
  }
  list(mu = mu, sigma = S, iter = it)
}

# Test de Little : renvoie statistique, ddl, valeur-p, nb de patrons
little_mcar_test <- function(X) {
  X <- as.matrix(X); n <- nrow(X); p <- ncol(X)
  M <- is.na(X)
  if (!any(M)) {
    return(list(statistic = NA_real_, df = 0L, p.value = NA_real_,
                n_patterns = 1L, message = "Aucune donnée manquante."))
  }
  em <- em_mvn(X); mu <- em$mu; S <- em$sigma
  key <- apply(M, 1, function(r) paste(as.integer(!r), collapse = ""))
  pats <- unique(key)
  d2 <- 0; somme_pj <- 0L
  for (k in pats) {
    idx <- which(key == k); obs <- !M[idx[1], ]; pj <- sum(obs)
    if (pj == 0) next
    nj <- length(idx)
    xbar <- colMeans(X[idx, obs, drop = FALSE])
    diff <- xbar - mu[obs]
    Sji <- solve(S[obs, obs, drop = FALSE] + diag(1e-8, pj))
    d2 <- d2 + nj * as.numeric(t(diff) %*% Sji %*% diff)
    somme_pj <- somme_pj + pj
  }
  df <- somme_pj - p
  pval <- if (df > 0) stats::pchisq(d2, df, lower.tail = FALSE) else NA_real_
  list(statistic = d2, df = as.integer(df), p.value = pval,
       n_patterns = length(pats), iter = em$iter)
}

# Démonstration pédagogique : le test sous chaque mécanisme (cible X1)
demonstration_little <- function(N, prop, rho, graine, p = 4L) {
  scen <- list(
    list(nom = "MCAR",                               mech = "MCAR", rho = rho, pil = NULL),
    list(nom = "MAR (piloté par X2)",                mech = "MAR",  rho = rho, pil = "X2"),
    list(nom = sprintf("MNAR (ρ = %s)", fmt_fr(rho, 2)), mech = "MNAR", rho = rho, pil = NULL),
    list(nom = "MNAR (variables décorrélées, ρ = 0)", mech = "MNAR", rho = 0,   pil = NULL)
  )
  rows <- lapply(scen, function(s) {
    d <- simuler_complet(N, p, s$rho, graine = graine)
    a <- amputer(d, prop, s$mech, graine = graine,
                 var_cible = "X1", var_pilote = s$pil, direction = "RIGHT")$amp
    t <- little_mcar_test(a)
    data.frame(Scenario = s$nom, Statistique = t$statistic, ddl = t$df,
               p = t$p.value,
               Verdict = ifelse(is.na(t$p.value), "—",
                          ifelse(t$p.value < 0.05, "MCAR rejeté", "MCAR non rejeté")),
               stringsAsFactors = FALSE)
  })
  do.call(rbind, rows)
}

# ============================================================================
#  Diagnostics d'imputation (mice) — densités, convergence, effet de m
# ============================================================================

# Ajustement MICE partagé (calculé une seule fois, réutilisé par les 3 volets)
ajuster_mice <- function(ampute, m = 10L, maxit = 20L, graine = 2024L) {
  set.seed(graine)
  mice::mice(ampute, m = as.integer(m), maxit = as.integer(maxit),
             printFlag = FALSE)
}

vars_imputees <- function(mids) {
  names(mids$imp)[vapply(mids$imp,
    function(x) !is.null(x) && nrow(x) > 0, logical(1))]
}

# 1) Densités : valeurs observées vs imputées (une courbe par jeu imputé)
graphe_obs_imp <- function(mids, var) {
  obs <- mids$data[[var]][!is.na(mids$data[[var]])]
  impmat <- as.matrix(mids$imp[[var]])            # n_manquant x m
  m <- ncol(impmat)
  d_obs <- data.frame(valeur = obs)
  d_imp <- data.frame(valeur = as.vector(impmat),
                      serie = rep(paste0("imp", seq_len(m)), each = nrow(impmat)))
  ggplot() +
    geom_density(data = d_imp,
                 aes(x = valeur, group = serie, colour = "Imputé"),
                 linewidth = 0.4, key_glyph = "path") +
    geom_density(data = d_obs, aes(x = valeur, colour = "Observé"),
                 linewidth = 1.3, key_glyph = "path") +
    scale_colour_manual(values = c(`Observé` = COUL_OBS, `Imputé` = COUL_NA),
                        breaks = c("Observé", "Imputé"), name = NULL) +
    sx_fr() + sy_fr() +
    labs(title = paste0("Valeurs observées vs imputées — ", var),
         subtitle = sprintf("Épaisse : observé · fines : les %d jeux imputés", m),
         x = var, y = "Densité") +
    theme_gg()
}

# 1b) Nuage : la variable imputée selon une autre variable (tous les jeux imputés)
#     Les valeurs imputées doivent prolonger la relation observée.
graphe_obs_imp_nuage <- function(mids, var, x) {
  long <- mice::complete(mids, "long")
  manq <- is.na(mids$data[[var]])
  long$statut <- factor(ifelse(manq[long$.id], "Imputé", "Observé"),
                        levels = c("Observé", "Imputé"))
  obs <- long[long$.imp == 1 & long$statut == "Observé", ]
  imp <- long[long$statut == "Imputé", ]
  ggplot() +
    geom_point(data = obs, aes(.data[[x]], .data[[var]], colour = statut),
               alpha = 0.45, size = 1.6) +
    geom_point(data = imp, aes(.data[[x]], .data[[var]], colour = statut),
               alpha = 0.35, size = 1.6) +
    scale_colour_manual(values = c(`Observé` = COUL_OBS, `Imputé` = COUL_NA), name = NULL) +
    guides(colour = guide_legend(override.aes = list(alpha = 1, size = 3))) +
    sx_fr() + sy_fr() +
    labs(title = sprintf("%s selon %s", var, x),
         subtitle = sprintf("Bleu : observé · orange : imputé (%d jeux)", mids$m),
         x = x, y = var) +
    theme_gg()
}

# 2) Convergence : traces des moyennes imputées par itération et par chaîne
graphe_convergence <- function(mids) {
  vi <- vars_imputees(mids)
  cm <- mids$chainMean
  niter <- dim(cm)[2]; m <- dim(cm)[3]
  rows <- list(); k <- 0L
  for (v in vi) {
    for (ch in seq_len(m)) {
      k <- k + 1L
      rows[[k]] <- data.frame(variable = v, iteration = seq_len(niter),
                              chaine = factor(ch), moyenne = cm[v, , ch])
    }
  }
  d <- do.call(rbind, rows)
  ggplot(d, aes(x = iteration, y = moyenne, colour = chaine, group = chaine)) +
    geom_line(linewidth = 0.5) +
    facet_wrap(~ variable, scales = "free_y") +
    scale_colour_viridis_d(guide = "none") +
    sy_fr() +
    labs(title = "Convergence de l'algorithme MICE",
         subtitle = "Moyenne imputée par itération (une couleur par chaîne) : on cherche une stabilisation sans tendance",
         x = "Itération", y = "Moyenne imputée") +
    theme_gg()
}

# 3) Effet du nombre d'imputations sur l'estimation et son erreur-type (Rubin)
effet_m <- function(mids, var) {
  M <- mids$m; n <- nrow(mids$data)
  q <- numeric(M); u <- numeric(M)
  for (k in seq_len(M)) {
    dk <- mice::complete(mids, k)[[var]]
    q[k] <- mean(dk); u[k] <- stats::var(dk) / n
  }
  res <- lapply(2:M, function(mm) {
    qb <- mean(q[1:mm]); B <- stats::var(q[1:mm]); Ub <- mean(u[1:mm])
    Tv <- Ub + (1 + 1 / mm) * B
    data.frame(m = mm, estimation = qb, se = sqrt(Tv))
  })
  do.call(rbind, res)
}

graphe_effet_m <- function(tab, var) {
  d <- rbind(data.frame(m = tab$m, valeur = tab$estimation,
                        panneau = paste0("Estimation de la moyenne de ", var)),
             data.frame(m = tab$m, valeur = tab$se,
                        panneau = "Erreur-type (règles de Rubin)"))
  d$panneau <- factor(d$panneau, levels = unique(d$panneau))
  # Échelles comparables : l'estimation est montrée sur ± 1 erreur-type finale,
  # l'erreur-type sur ± 30 % de sa valeur finale (pour juger l'ampleur des variations)
  se_f <- tab$se[nrow(tab)]; est_f <- tab$estimation[nrow(tab)]
  bornes <- data.frame(m = min(tab$m),
                       valeur = c(est_f - se_f, est_f + se_f, 0.7 * se_f, 1.3 * se_f),
                       panneau = factor(rep(levels(d$panneau), each = 2), levels = levels(d$panneau)))
  ggplot(d, aes(x = m, y = valeur)) +
    geom_blank(data = bornes) +
    geom_line(colour = COUL_OBS, linewidth = 0.8) +
    geom_point(colour = COUL_OBS, size = 1.8) +
    facet_wrap(~ panneau, scales = "free_y") +
    sy_fr() +
    labs(title = paste0("Effet du nombre d'imputations — moyenne de ", var),
         subtitle = "Résultat mis en commun sur les m premiers jeux; on cherche le m à partir duquel les courbes cessent de bouger.\nÉchelles : ± 1 erreur-type pour l'estimation, ± 30 % pour l'erreur-type.",
         x = "Nombre d'imputations (m)", y = NULL) +
    theme_gg()
}

# ============================================================================
#  Jeux de données réels avec manquants authentiques (aucune vérité connue)
# ============================================================================

charger_reel <- function(nom) {
  if (nom == "airquality") {
    d <- datasets::airquality
  } else if (nom == "nhanes") {
    d <- mice::nhanes
  } else {                       # boys : sous-ensemble numérique (bmi retiré, colinéaire)
    d <- mice::boys[, c("age", "hgt", "wgt", "hc", "tv")]
  }
  d[, vapply(d, is.numeric, logical(1)), drop = FALSE]
}

desc_reel <- function(nom) {
  switch(nom,
    airquality = paste0(
      "<b>Qualité de l'air à New York (airquality).</b> Mesures quotidiennes de ",
      "mai à septembre 1973 (153 jours). L'ozone (Ozone, ppb) et le rayonnement ",
      "solaire (Solar.R) comportent des valeurs manquantes ; le vent, la ",
      "température, le mois et le jour sont complets. ",
      "<i>Source : jeu airquality de R (paquet datasets).</i>"),
    nhanes = paste0(
      "<b>NHANES (nhanes, paquet mice).</b> Extrait de 25 observations du National ",
      "Health and Nutrition Examination Survey. Variables : age (groupe d'âge, ",
      "complet), bmi (indice de masse corporelle), hyp (hypertension), chl ",
      "(cholestérol). Petit jeu canonique d'enseignement de l'imputation. ",
      "<i>Source : Schafer (1997), diffusé avec mice.</i>"),
    boys = paste0(
      "<b>Croissance de garçons néerlandais (boys, paquet mice).</b> Échantillon ",
      "transversal de 748 garçons. Sous-ensemble numérique retenu : âge (age), ",
      "taille (hgt), poids (wgt), périmètre crânien (hc) et volume testiculaire ",
      "(tv). Le volume testiculaire manque surtout chez les plus jeunes — un ",
      "manque lié à l'âge, variable observée (illustration typique de MAR). ",
      "<i>Source : Fredriks et al. (2000), diffusé avec mice.</i>")
  )
}

rd_lecon <- function(jeu, rejet, pval = NA_real_) {
  ptxt <- if (is.na(pval)) "" else if (pval < 0.001) " (p < 0,001)" else paste0(" (p = ", fmt_fr(pval, 3), ")")
  verdict <- if (is.na(rejet)) "n'est pas calculable ici"
             else if (rejet) "rejette MCAR" else "ne rejette pas MCAR"
  bloc <- function(titre, ...) tagList(
    tags$p(style = "margin-bottom:3px;", strong(titre)), tags$p(...))
  switch(jeu,
    airquality = tagList(
      bloc("Le mécanisme plausible.",
        "L'ozone manque certains jours de l'été 1973. Son absence peut tenir aux conditions ",
        "météo — elles-mêmes ", strong("observées"), " (vent, température, rayonnement). ",
        "Un manque MAR est donc crédible, sans qu'on puisse le prouver."),
      bloc("Ce que dit le test ici.",
        "Le test de Little ", strong(verdict), ptxt, " : le manque ", strong("n'est pas"),
        " complètement aléatoire. Mais il ne tranche pas entre MAR et MNAR — cette limite est absolue."),
      bloc("La leçon.",
        "Sur données réelles, MCAR est rarement crédible. Les variables météo complètes sont des ",
        "prédicteurs naturels pour imputer l'ozone : une imputation qui les exploite (MAR) est ",
        "raisonnable, tout en restant une ", em("hypothèse"), " qu'aucun test ne confirme.")),
    nhanes = tagList(
      bloc("Le mécanisme plausible.",
        "Sur ce petit extrait (25 observations), l'IMC, le cholestérol et l'hypertension manquent. ",
        "On ignore pourquoi — et surtout, on manque cruellement de données pour le savoir."),
      bloc("Ce que dit le test ici — attention au piège.",
        "Le test ", strong(verdict), ptxt, ". On serait tenté de conclure « MCAR ». ",
        strong("Ce serait une erreur"), " : avec 25 observations, le test n'a presque aucune ",
        strong("puissance"), " — il ne rejette pas faute de pouvoir détecter quoi que ce soit."),
      bloc("La leçon.",
        "Absence de preuve n'est pas preuve d'absence. Un test de Little non significatif sur un ",
        "petit échantillon ne justifie ", strong("jamais"), " de conclure à MCAR — c'est la leçon ",
        "que porte nhanes.")),
    boys = tagList(
      bloc("Le mécanisme plausible.",
        "Le volume testiculaire (tv) manque surtout chez les ", strong("plus jeunes"),
        " : le manque dépend de l'âge, une variable ", strong("observée"), ". C'est le MAR d'école."),
      bloc("Ce que dit le test ici.",
        "Le test ", strong(verdict), " très nettement", ptxt, " : le manque est ",
        "fortement structuré, cohérent avec une dépendance à l'âge — sans prouver MAR."),
      bloc("La leçon.",
        "Un MAR plausible et interprétable, gouverné par une variable de domaine — l'âge est ici le ",
        strong("prédicteur du manque"), ". L'inclure dans le modèle d'imputation est essentiel ; ",
        "l'oublier compromettrait l'imputation de tv, exactement comme au bac à sable. ",
        "Le nuage « tv selon age » permet de vérifier que les valeurs imputées suivent la croissance."))
  )
}

# ============================================================================
#  Analyse de sensibilité MNAR par delta-ajustement (van Buuren, 2018, ch. 9)
#  Impute sous MAR, décale les valeurs IMPUTÉES de var d'un δ, met en commun.
# ============================================================================
sensibilite_delta <- function(ampute, var, deltas_sd, m = 10L, graine = 2024L) {
  m <- as.integer(m)
  set.seed(graine)
  mids <- mice::mice(ampute, m = m, printFlag = FALSE)
  n <- nrow(ampute)
  miss <- is.na(ampute[[var]])
  sdv <- stats::sd(ampute[[var]], na.rm = TRUE)
  comps <- lapply(seq_len(m), function(k) mice::complete(mids, k)[[var]])
  res <- lapply(deltas_sd, function(dsd) {
    d_raw <- dsd * sdv
    q <- numeric(m); u <- numeric(m)
    for (k in seq_len(m)) {
      xk <- comps[[k]]; xk[miss] <- xk[miss] + d_raw
      q[k] <- mean(xk); u[k] <- stats::var(xk) / n
    }
    qb <- mean(q); B <- if (m > 1) stats::var(q) else 0
    se <- sqrt(mean(u) + (1 + 1 / m) * B)
    data.frame(delta = dsd, estimation = qb, se = se,
               lo = qb - 1.96 * se, hi = qb + 1.96 * se)
  })
  do.call(rbind, res)
}

# δ (en écarts-types) qui récupère la valeur vraie, par interpolation linéaire
delta_recuperateur <- function(tab, vrai) {
  if (vrai < min(tab$estimation) || vrai > max(tab$estimation)) return(NA_real_)
  stats::approx(x = tab$estimation, y = tab$delta, xout = vrai, ties = "ordered")$y
}

graphe_sensibilite <- function(tab, var, vrai = NULL) {
  g <- ggplot(tab, aes(x = delta, y = estimation)) +
    geom_hline(yintercept = 0, colour = "#9a9a9a", linewidth = 0.4) +
    geom_ribbon(aes(ymin = lo, ymax = hi), fill = COUL_OBS, alpha = 0.2) +
    geom_line(colour = COUL_OBS, linewidth = 0.9) +
    geom_vline(xintercept = 0, linetype = "dotted") +
    annotate("text", x = 0, y = min(tab$lo), label = "  δ = 0 (MAR)",
             hjust = 0, vjust = 0, size = 3.6) +
    sx_fr() + sy_fr() +
    labs(title = paste0("Analyse de sensibilité MNAR — moyenne de ", var),
         subtitle = "Estimation (bande : IC à 95 %) selon le décalage δ\nappliqué aux seules valeurs imputées",
         x = "δ (en écarts-types de la variable)",
         y = paste0("Moyenne estimée de ", var)) +
    theme_gg()
  if (!is.null(vrai)) {
    g <- g +
      geom_hline(yintercept = vrai, linetype = "dashed",
                 colour = COUL_ACC, linewidth = 0.8) +
      annotate("text", x = min(tab$delta), y = vrai,
               label = "vraie moyenne (connue en simulation)",
               hjust = 0, vjust = -0.4, colour = COUL_ACC, size = 3.6)
  }
  g
}

# Seuils de fragilité : plus petit |δ| (de chaque côté) où la décision
# « l'IC à 95 % exclut 0 » change par rapport à δ = 0
seuils_fragilite <- function(tab) {
  dec <- ifelse(tab$lo > 0, "positive", ifelse(tab$hi < 0, "négative", "non significative"))
  i0 <- which.min(abs(tab$delta)); d0 <- dec[i0]
  pos <- which(tab$delta > tab$delta[i0] & dec != d0)
  neg <- which(tab$delta < tab$delta[i0] & dec != d0)
  list(decision0 = d0,
       haut = if (length(pos)) tab$delta[min(pos)] else NA_real_,
       dec_haut = if (length(pos)) dec[min(pos)] else NA_character_,
       bas = if (length(neg)) tab$delta[max(neg)] else NA_real_,
       dec_bas = if (length(neg)) dec[max(neg)] else NA_character_)
}

# ============================================================================
#  3. BIBLIOGRAPHIE (références réelles, vérifiables)
# ============================================================================
# ============================================================================
#  BAC À SABLE (mini-laboratoire de l'Introduction)
# ============================================================================
# ============================================================================
#  MATRICE manque <-> variables (corrélation point-bisériale) + co-absence
# ============================================================================
matrice_manque <- function(complet, ampute) {
  vars <- names(ampute); miss <- vars[colSums(is.na(ampute)) > 0]
  if (length(miss) == 0) return(NULL)
  M <- matrix(NA_real_, length(miss), length(vars) + 1,
              dimnames = list(miss, c(vars, "vraie")))
  for (j in miss) {
    Rj <- as.numeric(is.na(ampute[[j]]))
    for (k in vars) {
      if (k == j) next
      xk <- ampute[[k]]; ok <- !is.na(xk)
      if (sd(Rj[ok]) > 0 && sd(xk[ok]) > 0) M[j, k] <- cor(Rj[ok], xk[ok])
    }
    M[j, "vraie"] <- cor(Rj, complet[[j]])
  }
  M
}

matrice_coabsence <- function(ampute) {
  miss <- names(ampute)[colSums(is.na(ampute)) > 0]
  if (length(miss) < 2) return(NULL)
  cor(sapply(miss, function(j) as.numeric(is.na(ampute[[j]]))))
}

cellule_cor <- function(v, seuil = 0.5, extra = "") {
  st <- "padding:5px 12px; text-align:center; border:1px solid #e3e3e3;"
  if (is.na(v)) return(tags$td("\u2014", style = paste0(st, "color:#bbb;")))
  a <- min(abs(v) / seuil, 1); bg <- sprintf("rgba(217,95,14,%.2f)", a)
  fg <- if (a > 0.55) "#ffffff" else "#333333"
  tags$td(fmt_fr(v, 2, signe = TRUE), style = paste0(st, "background:", bg, ";color:", fg, ";", extra))
}

table_matrice_manque <- function(complet, ampute) {
  M <- matrice_manque(complet, ampute)
  if (is.null(M)) return(div(class = "encadre", "Aucune variable manquante \u00e0 analyser."))
  cols <- colnames(M); vars <- setdiff(cols, "vraie")
  sth <- "padding:5px 12px; text-align:center;"
  vraie_st <- "border-left:2px solid #999; font-style:italic;"
  entete <- tags$tr(tags$th("Manque de\u2026", style = paste0(sth, "text-align:left;")),
    lapply(vars, function(v) tags$th(v, style = sth)),
    tags$th("sa valeur vraie", style = paste0(sth, vraie_st, " color:#666;")))
  corps <- lapply(rownames(M), function(j) tags$tr(
    tags$td(sprintf("M(%s)", j), style = "padding:5px 12px; font-weight:bold;"),
    lapply(vars, function(k) cellule_cor(M[j, k])),
    cellule_cor(M[j, "vraie"], extra = vraie_st)))
  tbl <- tags$table(style = "border-collapse:collapse; font-family:monospace; font-size:0.92rem;",
                    entete, corps)
  co <- matrice_coabsence(ampute); bloc_co <- NULL
  if (!is.null(co)) {
    mv <- rownames(co)
    ent2 <- tags$tr(tags$th("", style = sth),
      lapply(mv, function(v) tags$th(sprintf("M(%s)", v), style = sth)))
    cor2 <- lapply(mv, function(i) tags$tr(
      tags$td(sprintf("M(%s)", i), style = "padding:5px 12px; font-weight:bold;"),
      lapply(mv, function(k) if (i == k)
        tags$td("\u2014", style = "padding:5px 12px; text-align:center; color:#bbb;")
        else cellule_cor(co[i, k]))))
    bloc_co <- tagList(
      tags$p(style = "margin-top:14px;",
             tags$b("Co-absence (\u03c6) \u2014 quelles variables manquent ensemble")),
      tags$table(style = "border-collapse:collapse; font-family:monospace; font-size:0.92rem;",
                 ent2, cor2))
  }
  tagList(tbl, bloc_co)
}

bac_donnees <- local({
  set.seed(101)
  R <- matrix(c(1,   .78, .35, .35,          # X2 = prédicteur dominant de X1
                .78, 1,   .15, .15,          # X3, X4 : faiblement liés à X1 et à X2
                .35, .15, 1,   .30,
                .35, .15, .30, 1), 4, 4, byrow = TRUE)
  M <- MASS::mvrnorm(20, rep(0, 4), R)
  setNames(as.data.frame(round(50 + 10 * M)), paste0("X", 1:4))
})

# Sous MCAR, les lignes perdues sont tirées au sort; « tirage » fixe ce tirage
# (110 = tirage n° 1, représentatif; « Autre tirage » passe à 111, 112, ...).
BAC_TIRAGE_DEFAUT <- 110L
bac_lignes_manq <- function(mech, k, tirage = BAC_TIRAGE_DEFAUT) {
  d <- bac_donnees
  if (mech == "MCAR") { set.seed(tirage); sort(sample(seq_len(nrow(d)), k)) }
  else if (mech == "MAR")  order(d$X2, decreasing = TRUE)[seq_len(k)]
  else                     order(d$X1, decreasing = TRUE)[seq_len(k)]
}

bac_predicteurs <- function(modele) {
  if (identical(modele, "incomplet")) c("X3", "X4") else c("X2", "X3", "X4")
}

bac_valeurs_imputees <- function(mech, k, imput, modele = "complet",
                                 tirage = BAC_TIRAGE_DEFAUT) {
  d <- bac_donnees; n <- nrow(d)
  miss <- bac_lignes_manq(mech, k, tirage); obs <- setdiff(seq_len(n), miss)
  if (imput == "moyenne") {
    v <- rep(mean(d$X1[obs]), length(miss))
  } else {
    fit <- lm(stats::reformulate(bac_predicteurs(modele), "X1"), data = d[obs, ])
    v <- as.numeric(predict(fit, d[miss, ]))
    if (imput == "reg_bruit") { set.seed(2024); v <- v + rnorm(length(miss), 0, sigma(fit)) }
  }
  setNames(round(v), as.character(miss))
}

table_bac <- function(mech, k, vue, imput = "none", modele = "complet",
                      tirage = BAC_TIRAGE_DEFAUT) {
  d <- bac_donnees; n <- nrow(d); vars <- names(d)
  rows_miss <- bac_lignes_manq(mech, k, tirage)
  driver <- switch(mech, MAR = "X2", NMAR = "X1", MCAR = NA_character_)
  rows_red <- if (is.na(driver)) integer(0) else rows_miss
  imp_val <- NULL
  if (vue == "apres" && imput != "none")
    imp_val <- bac_valeurs_imputees(mech, k, imput, modele, tirage)
  cellule <- function(i, col) {
    val <- d[i, col]
    rouge  <- (!is.na(driver) && col == driver && i %in% rows_red)
    manque <- (vue == "apres" && col == "X1" && i %in% rows_miss)
    st <- "padding:4px 12px; text-align:center; border:1px solid #e3e3e3;"
    if (manque) {
      if (!is.null(imp_val)) return(tags$td(imp_val[[as.character(i)]], style = paste0(st,
        "background:#94c9f3; color:#0b5563; font-weight:bold;")))
      return(tags$td("NA", style = paste0(st,
        "background:#fde0e0; color:#b30000; font-style:italic; border:1px dashed #d95f0e;")))
    }
    if (rouge)  return(tags$td(val, style = paste0(st,
      "background:#f6a8a8; color:#6b0000; font-weight:bold;")))
    tags$td(val, style = st)
  }
  entete <- tags$tr(
    tags$th("#", style = "padding:4px 12px; color:#888;"),
    lapply(vars, function(v) tags$th(v, style = "padding:4px 12px;")))
  corps <- lapply(seq_len(n), function(i)
    tags$tr(tags$td(i, style = "padding:4px 12px; color:#aaa;"),
            lapply(vars, function(col) cellule(i, col))))
  tags$table(style = "border-collapse:collapse; font-family:monospace; font-size:0.92rem;",
             entete, corps)
}

bac_legende <- function(mech, vue, imput = "none", modele = "complet") {
  base <- switch(paste(mech, vue),
    "MCAR avant" = tagList(strong("MCAR. "),
      "Le manque n'est lié à aucune valeur : aucun motif. Passez à ",
      strong("Après"), " pour voir disparaître des cellules ", strong("au hasard"), "."),
    "MCAR apres" = tagList(strong("MCAR. "),
      "Les valeurs manquantes de X1 sont ", strong("dispersées au hasard"),
      " — indépendantes de toute variable."),
    "MAR avant"  = tagList(strong("MAR. "), "Les valeurs ", strong("élevées de X2"),
      " (en rouge) vont provoquer le manque de X1 dans ces mêmes lignes. Passez à ",
      strong("Après"), "."),
    "MAR apres"  = tagList(strong("MAR. "), "X1 a disparu ",
      strong("là où X2 est élevé"), " (encore en rouge). Le manque dépend d'une variable ",
      strong("observée"), " : on peut le corriger en s'appuyant sur X2."),
    "NMAR avant" = tagList(strong("NMAR. "), "Les valeurs ",
      strong("élevées de X1 elles-mêmes"),
      " (en rouge) vont provoquer leur propre manque. Passez à ", strong("Après"), "."),
    "NMAR apres" = tagList(strong("NMAR. "),
      "Les valeurs élevées de X1 se sont ", strong("effacées elles-mêmes"),
      " : les cellules rouges sont devenues NA. ", strong("La preuve a disparu"),
      " — voilà pourquoi NMAR est indétectable à partir des seules données observées.")
  )
  if (vue == "apres" && imput != "none") {
    preds <- if (identical(modele, "incomplet")) "X3 et X4 seulement (X2 est oublié)" else "X2, X3 et X4"
    note <- if (imput == "moyenne")
      tagList(tags$br(), tags$span(style = "color:#0b5563;",
        strong("Imputation par la moyenne : "),
        "chaque case comblée (en bleu) reçoit la MÊME valeur — la moyenne des X1 observés. La variabilité est écrasée."))
    else if (imput == "regression")
      tagList(tags$br(), tags$span(style = "color:#0b5563;",
        strong("Imputation par régression : "),
        "chaque case (en bleu) reçoit la valeur PRÉDITE à partir de ", preds,
        ". Elle suit la tendance, mais tombe exactement sur la droite : aucune dispersion autour d'elle."))
    else
      tagList(tags$br(), tags$span(style = "color:#0b5563;",
        strong("Imputation par régression + erreur : "),
        "valeur prédite à partir de ", preds, ", plus un résidu tiré au hasard. ",
        "La dispersion naturelle revient — au prix d'un tirage aléatoire qui fait bouger un peu les résultats."))
    return(tagList(base, note))
  }
  base
}

bac_analyse <- function(mech, k, imput, modele = "complet", tirage = BAC_TIRAGE_DEFAUT) {
  d <- bac_donnees; n <- nrow(d)
  miss <- bac_lignes_manq(mech, k, tirage); obs <- setdiff(seq_len(n), miss)
  est <- function(x1) {
    fit <- lm(x1 ~ d$X2)
    c(
      moyenne = mean(x1, na.rm = TRUE),
      ecart_type = sd(x1, na.rm = TRUE),
      pente = unname(coef(fit)[[2]]),
      et_pente = unname(summary(fit)$coefficients[2, "Std. Error"])
    )
  }
  vrai <- est(d$X1)
  x1cc <- d$X1; x1cc[miss] <- NA
  cc <- est(x1cc)
  imp <- NULL
  if (imput != "none") {
    x1i <- d$X1
    x1i[miss] <- bac_valeurs_imputees(mech, k, imput, modele, tirage)
    imp <- est(x1i)
  }
  list(vrai = vrai, cc = cc, imp = imp)
}

graphe_bac_nuage <- function(mech, k, imput, modele = "complet", tirage = BAC_TIRAGE_DEFAUT) {
  d <- bac_donnees; n <- nrow(d)
  miss <- bac_lignes_manq(mech, k, tirage); obs <- setdiff(seq_len(n), miss)
  fit2 <- lm(X1 ~ X2, data = d[obs, ])
  p <- ggplot() +
    geom_abline(intercept = coef(fit2)[1], slope = coef(fit2)[2],
                linetype = "dashed", colour = "#b0b0b0") +
    geom_point(data = d[miss, ], aes(X2, X1), shape = 1, colour = "#8a8a8a",
               size = 2.8, stroke = 1) +
    geom_point(data = d[obs, ], aes(X2, X1), colour = COUL_OBS, size = 2.4)
  if (imput != "none") {
    vimp <- as.numeric(bac_valeurs_imputees(mech, k, imput, modele, tirage))
    seg <- data.frame(X2 = d$X2[miss], y1 = d$X1[miss], y2 = vimp)
    p <- p +
      geom_segment(data = seg, aes(x = X2, xend = X2, y = y1, yend = y2),
                   colour = COUL_NA, alpha = 0.35, linewidth = 0.4) +
      geom_point(data = data.frame(X2 = d$X2[miss], X1 = vimp), aes(X2, X1),
                 colour = COUL_NA, size = 2.8)
  }
  p + labs(x = if (mech == "MAR") "X2 (prédicteur du manque)" else "X2",
           y = "X1",
           title = "Où tombent les valeurs imputées ?",
           subtitle = "Bleu : observé  ·  Cercle gris : vraie valeur (effacée)  ·  Orange : imputé\nTirets : droite de régression des cas complets") +
    theme_gg()
}

table_analyse <- function(mech, k, imput, modele = "complet", tirage = BAC_TIRAGE_DEFAUT) {
  a <- bac_analyse(mech, k, imput, modele, tirage)
  vert  <- "color:#177245; font-weight:bold;"
  rouge <- "color:#b30000; font-weight:bold;"
  ambre <- "color:#b36b00; font-weight:bold;"
  coul <- function(val, vrai, tol) if (abs(val - vrai) <= tol) vert else rouge
  # Erreur-type : trop petite = rouge (trop optimiste); trop grande = ambre (moins précis)
  coul_et <- function(val, vrai) {
    r <- val / vrai
    if (r < 0.8) rouge else if (r > 1.2) ambre else vert
  }
  st <- "padding:5px 10px; text-align:center; border:1px solid #e3e3e3;"
  ligne <- function(nom, idx, dec, tol = NULL) {
    cfun <- function(v) if (is.null(tol)) coul_et(v, a$vrai[idx]) else coul(v, a$vrai[idx], tol)
    cells <- list(
      tags$td(nom, style = paste0(st, "text-align:left;")),
      tags$td(fmt_fr(a$vrai[idx], dec), style = paste0(st, "font-weight:bold;")),
      tags$td(fmt_fr(a$cc[idx], dec), style = paste0(st, cfun(a$cc[idx]))))
    if (!is.null(a$imp)) cells <- c(cells, list(tags$td(fmt_fr(a$imp[idx], dec),
      style = paste0(st, cfun(a$imp[idx])))))
    do.call(tags$tr, cells)
  }
  entete <- do.call(tags$tr, c(
    list(tags$th("Analyse", style = paste0(st, "text-align:left;")),
         tags$th("Vérité", style = st), tags$th("Cas complets", style = st)),
    if (!is.null(a$imp)) list(tags$th(
      if (imput == "moyenne") "Après moyenne"
      else if (identical(modele, "incomplet") && imput == "reg_bruit") "Après régr. + erreur (sans X2)"
      else if (identical(modele, "incomplet")) "Après régr. (sans X2)"
      else if (imput == "reg_bruit") "Après régr. + erreur"
      else "Après régression", style = st))))
  tagList(
    tags$table(style = "border-collapse:collapse; font-size:0.92rem;",
      entete,
      ligne("Moyenne de X1", "moyenne", 1, 1.0),
      ligne("Écart-type de X1", "ecart_type", 1, 0.8),
      ligne("Pente de X1 ~ X2", "pente", 2, 0.10),
      ligne("Erreur-type naïve de la pente", "et_pente", 3)),
    tags$small(style = "display:block; margin-top:7px; color:#59636d; line-height:1.45;",
      tags$span(style = vert, "Vert"), " : proche de la vérité · ",
      tags$span(style = rouge, "rouge"), " : écart notable (pour l'erreur-type : trop petite, donc trop optimiste) · ",
      tags$span(style = ambre, "ambre"), " : erreur-type plus grande que la vraie (perte de précision).", tags$br(),
      "Après une imputation unique, l'erreur-type est calculée comme si les valeurs imputées avaient été observées : ",
      "elle omet l'incertitude d'imputation, que seule l'imputation multiple intègre.")
  )
}

bac_interpretation <- function(mech, imput = "none", modele = "complet") {
  base <- switch(mech,
    MCAR = tagList("Le manque est purement aléatoire : aucune estimation ne s'écarte de la vérité ",
      strong("de façon systématique"), ". Un écart isolé reste possible dans 20 lignes (essayez ",
      em("Autre tirage"), ") : il change de sens d'un tirage à l'autre — c'est une fluctuation, ",
      "pas un biais. Les cas complets perdent surtout de la ", strong("précision"),
      " (erreur-type plus grande)."),
    MAR  = tagList("Les cas complets ", strong("sous-estiment la moyenne"),
      " : les X1 effacés étaient justement les plus grands, puisque X2 et X1 vont de pair."),
    NMAR = tagList("Le manque dépend de X1 lui-même : ", strong("toutes les corrections restent biaisées"),
      ". La régression rapproche un peu la moyenne grâce au lien X1–X2, sans l'atteindre."))
  ajout <- NULL
  if (imput == "moyenne") {
    ajout <- tagList("L'imputation par la ", strong("moyenne"),
      " ne corrige pas la moyenne, ", strong("rétrécit l'écart-type"), " et ",
      strong("aplatit la pente"), (if (mech == "MCAR") " — même sous MCAR." else "."))
  } else if (imput %in% c("regression", "reg_bruit") && mech == "MAR") {
    ajout <- if (identical(modele, "incomplet")) {
      tagList("Sans X2, la régression ", strong("ne voit pas la sélection"),
              " : la moyenne reste aussi biaisée qu'avec les cas complets.")
    } else if (imput == "regression") {
      tagList("La régression sur X2 ", strong("récupère la moyenne et l'écart-type"),
              ". Mais son erreur-type naïve est ", strong("trop petite"),
              " : les valeurs prédites sont traitées comme certaines.")
    } else {
      tagList("Le résidu ajouté rend la ", strong("dispersion"), " naturelle; mais un ",
              strong("seul tirage"), " déplace un peu la moyenne et la pente, et l'erreur-type ",
              "reste trop petite. D'où l'idée de répéter l'imputation : l'imputation multiple.")
    }
  }
  if (is.null(ajout)) base else tagList(base, tags$br(), ajout)
}

bac_explication <- function(mech) {
  q <- function(txt) tags$p(style = "margin-bottom:4px;", strong(txt))
  switch(mech,
    MCAR = tagList(
      q("Pourquoi tout reste juste ?"),
      tags$p("Sous MCAR, les lignes qui perdent X1 sont tirées au sort : elles ne diffèrent ",
             "en rien des autres. La sous-population observée est une mini-réplique fidèle de ",
             "l'ensemble — mêmes moyennes, mêmes relations."),
      tags$p("C'est le ", strong("seul cas"), " où supprimer les lignes incomplètes ",
             "(cas complets) ne coûte que de la ", em("précision"), ", jamais du ",
             em("biais"), ".")),
    MAR = tagList(
      q("Pourquoi les cas complets échouent — et pourquoi la régression sauve la mise ?"),
      tags$p("X1 disparaît ", strong("là où X2 est élevé"), ". Comme X1 et X2 vont de pair, ",
             "les X1 manquants étaient justement les plus grands : les cas complets ne voient ",
             "que les petits, d'où une moyenne trop basse."),
      tags$p("Mais X2, lui, est ", strong("observé"), " : la régression s'en sert pour ",
             "reconstituer les grands X1 absents, et se rapproche de la vérité. L'imputation par la ",
             strong("moyenne"), " rebouche au contraire tous les trous avec une même valeur ",
             "tiède — la moyenne reste aussi basse, l'écart-type rétrécit et ",
             strong("la relation X1–X2 s'aplatit"), " (la pente s'effondre)."),
      tags$p("Encore faut-il ", strong("inclure X2 dans le modèle d'imputation"),
             " : passez le ", em("Modèle d'imputation"), " sur « incomplet » et regardez la ",
             "récupération s'effondrer. Le manque est toujours MAR — il dépend d'une variable ",
             "observée — mais une imputation qui ignore cette variable ne peut pas en tirer parti : ",
             "le biais reste alors ", strong("aussi entier que sous NMAR"), ".")),
    NMAR = tagList(
      q("Pourquoi rien ne suffit ?"),
      tags$p("Ici, X1 disparaît ", strong("là où X1 lui-même est élevé"), ". ",
             "L'information qui permettrait de deviner les valeurs manquantes est ",
             em("précisément celle qui manque"), "."),
      tags$p("Aucune variable observée ne porte toute la trace de la sélection : même la ",
             "régression ne récupère qu'en partie. On ne reconstitue pas ce dont il ne reste ",
             strong("aucune empreinte"), " — la vérité est hors de portée, et seule une ",
             strong("analyse de sensibilité"), " permet d'en cerner l'incertitude.")))
}

BIBLIO <- list(
  ouvrages = list(
    "Little, R. J. A., & Rubin, D. B. (2019). Statistical Analysis with Missing Data (3e éd.). Wiley.",
    "van Buuren, S. (2018). Flexible Imputation of Missing Data (2e éd.). CRC Press.",
    "Enders, C. K. (2022). Applied Missing Data Analysis (2e éd.). Guilford Press.",
    "Carpenter, J. R., & Kenward, M. G. (2013). Multiple Imputation and its Application. Wiley.",
    "Allison, P. D. (2001). Missing Data. Sage (Quantitative Applications in the Social Sciences, 136)."
  ),
  articles = list(
    "Rubin, D. B. (1976). Inference and missing data. Biometrika, 63(3), 581–592.",
    "Rubin, D. B. (1987). Multiple Imputation for Nonresponse in Surveys. Wiley.",
    "Little, R. J. A. (1988). A test of missing completely at random for multivariate data with missing values. Journal of the American Statistical Association, 83(404), 1198–1202.",
    "Schafer, J. L., & Graham, J. W. (2002). Missing data: Our view of the state of the art. Psychological Methods, 7(2), 147–177.",
    "Graham, J. W. (2009). Missing data analysis: Making it work in the real world. Annual Review of Psychology, 60, 549–576.",
    "van Buuren, S., & Groothuis-Oudshoorn, K. (2011). mice: Multivariate Imputation by Chained Equations in R. Journal of Statistical Software, 45(3), 1–67.",
    "Rosseel, Y. (2012). lavaan: An R Package for Structural Equation Modeling. Journal of Statistical Software, 48(2), 1–36."
  ),
  ressources = list(
    list(titre = "Flexible Imputation of Missing Data (S. van Buuren) — l'ouvrage de référence, en ligne et libre",
         url = "https://stefvanbuuren.name/fimd/"),
    list(titre = "mice — site officiel du paquet : vignettes, exemples, documentation",
         url = "https://amices.org/mice/"),
    list(titre = "naniar — exploration et visualisation des données manquantes en R",
         url = "https://naniar.njtierney.com/"),
    list(titre = "lavaan — équations structurelles et estimation FIML en R",
         url = "https://lavaan.ugent.be/")
  )
)

# ============================================================================
#  Banque de scénarios du quiz de classification (MCAR / MAR / NMAR)
# ============================================================================
QUIZ <- list(
  list(texte = "Un capteur météo tombe en panne quelques jours, de façon purement aléatoire et sans lien avec la température qu'il aurait mesurée.",
       rep = "MCAR",
       expl = "La panne est indépendante de toute donnée, observée ou non : les jours manquants forment un sous-échantillon représentatif."),
  list(texte = "Dans une enquête, les personnes âgées répondent moins souvent à la question sur le revenu. L'âge est enregistré pour tous les répondants.",
       rep = "MAR",
       expl = "Le manque dépend de l'âge, qui est observé — et non du revenu lui-même. C'est la définition de MAR."),
  list(texte = "Les personnes à revenu élevé refusent plus souvent de déclarer leur revenu.",
       rep = "NMAR",
       expl = "Le manque dépend de la valeur du revenu, précisément la donnée manquante : mécanisme non aléatoire (NMAR)."),
  list(texte = "La dernière page d'un questionnaire papier est arrachée par accident pour un lot d'exemplaires tirés au hasard.",
       rep = "MCAR",
       expl = "L'accident est sans aucun lien avec les réponses : perte complètement aléatoire."),
  list(texte = "Les patients les plus déprimés quittent l'étude avant de remplir l'échelle de dépression finale.",
       rep = "NMAR",
       expl = "L'abandon dépend du niveau de dépression finale — justement la valeur non observée. C'est le piège classique du décrochage lié à l'issue : NMAR."),
  list(texte = "Les hommes sautent plus souvent la question sur le poids que les femmes. Le sexe est connu pour tous.",
       rep = "MAR",
       expl = "Le manque dépend du sexe (observé), pas du poids. Conditionnellement au sexe, il redevient aléatoire : MAR."),
  list(texte = "Un appareil de mesure n'enregistre pas les concentrations situées sous son seuil de détection.",
       rep = "NMAR",
       expl = "La valeur manque parce qu'elle est trop basse : le manque dépend de la valeur elle-même (censure), donc NMAR."),
  list(texte = "Les élèves absents le jour d'un test ont, en moyenne, de moins bonnes notes aux évaluations antérieures, lesquelles sont enregistrées.",
       rep = "MAR",
       expl = "L'absence est liée aux notes antérieures observées, non à la note manquante du jour : conditionnellement à l'historique, MAR."),
  list(texte = "Une balance électronique n'affiche rien pour les objets dont le poids dépasse sa capacité maximale.",
       rep = "NMAR",
       expl = "Le manque dépend du poids réel (trop élevé), c'est-à-dire de la valeur manquante : NMAR."),
  list(texte = "Un laboratoire égare 5 % des tubes d'échantillon, au hasard, sans lien avec leur contenu.",
       rep = "MCAR",
       expl = "Perte purement aléatoire, indépendante du contenu : MCAR."),
  list(texte = "Sur une plateforme de commerce, les clients les plus insatisfaits ne laissent presque jamais d'avis.",
       rep = "NMAR",
       expl = "La note manque à cause du niveau de satisfaction qu'elle exprimerait — la valeur manquante elle-même : NMAR."),
  list(texte = "Les revenus manquent plus souvent dans les régions rurales ; la région est connue pour chaque ménage.",
       rep = "MAR",
       expl = "Le manque dépend de la région, qui est observée. Conditionnellement à la région, il redevient aléatoire : MAR."),
  list(texte = "Une coupure de courant efface les relevés d'une heure choisie au hasard dans la journée.",
       rep = "MCAR",
       expl = "L'effacement est sans lien avec les valeurs relevées : complètement aléatoire, MCAR."),
  list(texte = "Les patients les plus malades sont transférés vers un autre service et ne remplissent pas le questionnaire de sortie ; leur score de gravité à l'admission figure au dossier.",
       rep = "MAR",
       expl = "Piège classique : « les plus malades » évoque NMAR, mais la gravité est OBSERVÉE à l'admission. Conditionnellement à ce score, le manque est aléatoire : MAR."),
  list(texte = "Un thermomètre connecté cesse d'émettre quand sa batterie faiblit ; or la batterie se vide plus vite par grand froid.",
       rep = "NMAR",
       expl = "Le manque survient quand il fait froid, donc dépend de la température — précisément la valeur non enregistrée : NMAR."),
  list(texte = "Les utilisateurs les plus lents d'une appli de course désactivent le partage de leur fréquence cardiaque ; leur vitesse, elle, est toujours enregistrée.",
       rep = "MAR",
       expl = "La fréquence cardiaque manque selon la vitesse, qui est observée — pas selon la fréquence elle-même : MAR."),
  list(texte = "Les étudiants ayant obtenu une mauvaise note demandent plus souvent son retrait du tableau d'affichage public.",
       rep = "NMAR",
       expl = "Le retrait dépend de la valeur de la note — justement la donnée qui manque : NMAR."),
  list(texte = "Dans un panel, les ménages déménagent et sortent de l'étude ; l'adresse antérieure et les vagues précédentes sont enregistrées.",
       rep = "MAR",
       expl = "L'attrition est liée à des informations observées (historique, localisation), non à la valeur manquante de la vague courante : conditionnellement à l'historique, MAR.")
)

# ============================================================================
#  4. INTERFACE UTILISATEUR
# ============================================================================
ui <- page_navbar(
  title = "simManquant",
  theme = theme_app,
  fillable = FALSE,
  header = tags$head(tags$style(HTML("
    .encadre {background:#eef6fb;border-left:4px solid #2c7fb8;padding:12px 16px;border-radius:4px;margin:10px 0;}
    .alerte  {background:#fff3e6;border-left:4px solid #d95f0e;padding:12px 16px;border-radius:4px;margin:10px 0;}
    .cle     {background:#eafaf0;border-left:4px solid #31a354;padding:12px 16px;border-radius:4px;margin:10px 0;}
    h4 {color:#2c7fb8;margin-top:18px;}
  "))),

  # ---------------------------------------------------------------- INTRO ---
  nav_panel(
    title = "Introduction",
    icon = icon("book-open"),
    withMathJax(),
    navset_tab(
      nav_panel(
        "Concepts",
    card(
      card_header("Le problème des données manquantes"),
      card_body(
        layout_column_wrap(
          width = 1/2,
          p("Dans presque toute étude empirique, une partie des données est absente : ",
            "un répondant saute une question, un instrument tombe en panne, un ",
            "participant abandonne l'étude. Ignorer ce phénomène — ou le traiter ",
            "naïvement — peut ", strong("biaiser les estimations"), ", ",
            strong("gonfler ou réduire la variance"), " et ",
            strong("fausser les tests d'hypothèses"), "."),
          div(class = "encadre",
            "La question centrale n'est pas « ", strong("combien"), " » de données ",
            "manquent, mais « ", strong("pourquoi"), " » elles manquent. Ce ",
            em("mécanisme"), " détermine quelles méthodes sont valides — c'est le fil ",
            "conducteur de toute l'application.")
        )
      )
    ),
    h4("Les trois mécanismes de données manquantes (Rubin, 1976)"),
    p("Soit \\(Y\\) les données (partie observée \\(Y_{obs}\\), partie manquante ",
      "\\(Y_{mis}\\)) et \\(R\\) l'indicateur de réponse (\\(R=1\\) si observé). ",
      "Les mécanismes se distinguent par ce dont dépend la probabilité de manquer."),
    layout_column_wrap(
      width = 1/3,
      card(class = "h-100",
        card_header("MCAR", class = "bg-success text-white"),
        card_body(
          p(em("Missing Completely At Random")),
          p("Le manque ne dépend d'", strong("aucune"), " donnée :"),
          p("$$ P(R \\mid Y_{obs}, Y_{mis}) = P(R) $$"),
          div(class = "cle", strong("Exemple. "),
            "un tube d'échantillon se casse au hasard. Les cas complets restent ",
            em("représentatifs"), " : l'analyse perd surtout de la précision.")
        )),
      card(class = "h-100",
        card_header("MAR", class = "bg-primary text-white"),
        card_body(
          p(em("Missing At Random")),
          p("Le manque dépend des données ", strong("observées"), " :"),
          p("$$ P(R \\mid Y_{obs}, Y_{mis}) = P(R \\mid Y_{obs}) $$"),
          div(class = "cle", strong("Exemple. "),
            "les hommes déclarent moins souvent leur poids que les femmes. Le manque ",
            "dépend du sexe (observé) ; une imputation multiple ou un FIML qui ",
            "inclut les bonnes variables peut alors être valide.")
        )),
      card(class = "h-100",
        card_header("MNAR (NMAR)", class = "bg-warning"),
        card_body(
          p(em("Missing Not At Random")),
          p("Le manque dépend des valeurs ", strong("manquantes"), " elles-mêmes :"),
          p("$$ P(R \\mid Y_{obs}, Y_{mis}) \\ \\text{dépend de}\\ Y_{mis} $$"),
          div(class = "alerte", strong("Exemple. "),
            "les hauts revenus refusent de déclarer leur revenu. Aucune méthode ",
            "standard n'est garantie sans modéliser explicitement le manque.")
        ))
    ),
    card(
      card_header("Tableau de synthèse"),
      card_body(
        tags$table(class = "table table-sm table-striped mb-0",
          tags$thead(tags$tr(
            tags$th("Mécanisme"), tags$th("Le manque dépend de…"),
            tags$th("Analyse des cas complets"),
            tags$th("Imputation multiple / FIML"))),
          tags$tbody(
            tags$tr(tags$td(strong("MCAR")), tags$td("rien"),
              tags$td("valide; moins précis"),
              tags$td(tags$span(style = "color:#31a354;", "valide si le modèle est adéquat"))),
            tags$tr(tags$td(strong("MAR")), tags$td("valeurs observées"),
              tags$td(tags$span(style = "color:#d95f0e;", "souvent biaisé")),
              tags$td(tags$span(style = "color:#31a354;", "valide si MAR et modèle adéquat"))),
            tags$tr(tags$td(strong("MNAR")), tags$td("valeurs manquantes"),
              tags$td(tags$span(style = "color:#d95f0e;", "souvent biaisé")),
              tags$td(tags$span(style = "color:#d95f0e;", "analyse de sensibilité requise"))))
        )
      )
    ),
    layout_column_wrap(
      width = 1/2,
      card(class = "h-100",
        card_header("Impact sur les analyses"),
        card_body(
          tags$ul(
            tags$li(strong("Perte de puissance : "),
                    "supprimer les cas incomplets réduit la taille d'échantillon."),
            tags$li(strong("Biais : "),
                    "sous MAR/MNAR, l'analyse des cas complets peut biaiser moyennes, ",
                    "corrélations et coefficients de régression."),
            tags$li(strong("Sous-estimation de la variance : "),
                    "une valeur unique imputée (moyenne, régression déterministe) est ",
                    "traitée comme certaine et ", em("rétrécit"), " la dispersion et les ",
                    "erreurs-types."),
            tags$li(strong("Distorsion des relations : "),
                    "l'imputation par la moyenne atténue les corrélations vers zéro.")
          )
        )),
      card(class = "h-100",
        card_header("Approches de traitement"),
        card_body(
          tags$ul(
            tags$li(strong("Cas complets (listwise) : "),
                    "simple ; valide sous MCAR ; sinon perte d'information et biais."),
            tags$li(strong("Valeurs disponibles (pairwise) : "),
                    "exploite chaque paire ; matrices parfois non cohérentes."),
            tags$li(strong("Imputation par la moyenne : "),
                    "à éviter — rétrécit la variance et atténue les liens."),
            tags$li(strong("Imputation par régression : "),
                    "déterministe (trop lisse) ou stochastique (résidu aléatoire ajouté)."),
            tags$li(strong("Imputation multiple (Rubin) : "),
                    "plusieurs jeux complétés, analysés puis mis en commun ; souvent ",
                    "appropriée sous MAR si le modèle d'imputation est bien construit."),
            tags$li(strong("Maximum de vraisemblance (FIML) : "),
                    "estime directement sur les données incomplètes ; sous MAR et avec ",
                    "un modèle adéquat, il est une bonne alternative à l'imputation multiple.")
          )
        ))
    ),
    card(
      card_header("Comment utiliser ce laboratoire"),
      card_body(
        p("Ce laboratoire accompagne le ",
          tags$b("Guide de laboratoire \u2014 Analyse des données manquantes"),
          ", disponible \u00e0 ",
          tags$a("dcstats.org/simMISSING", href = "https://dcstats.org/simMISSING",
                 target = "_blank"),
          ". Lisez un concept dans le guide, puis manipulez-le ici : on retient ce qu'on manipule."),
        p("Commencez par le sous-onglet ", tags$b("Bac \u00e0 sable"), " (à droite de celui-ci) et le ",
          tags$b("Quiz"), ", puis suivez les onglets de la barre du haut, de gauche \u00e0 droite :"),
        tags$table(class = "table table-sm mb-0",
          tags$tbody(
            tags$tr(tags$td(strong("Données & visualisation")),
              tags$td("simuler un jeu, choisir un mécanisme, voir la carte, les patrons et la sélectivité du manque")),
            tags$tr(tags$td(strong("Traitements")),
              tags$td("comparer cas complets, imputations simples, imputation multiple et FIML face à la vérité")),
            tags$tr(tags$td(strong("Mécanisme du manque \u2192 Diagnostic")),
              tags$td("tester MCAR (test de Little) et comprendre pourquoi MAR et MNAR sont indiscernables")),
            tags$tr(tags$td(strong("Mécanisme du manque \u2192 Sensibilité MNAR")),
              tags$td("encadrer la conséquence d'un manque non aléatoire par delta-ajustement")),
            tags$tr(tags$td(strong("Évaluer et valider \u2192 Étude de simulation")),
              tags$td("répéter le cycle pour distinguer biais, efficacité et couverture des intervalles")),
            tags$tr(tags$td(strong("Évaluer et valider \u2192 Qualité de l'imputation")),
              tags$td("diagnostiquer une imputation sans vérité : densités, convergence, effet de m")),
            tags$tr(tags$td(strong("Données réelles")),
              tags$td("appliquer la démarche à airquality, NHANES et boys"))
          ))
      )
    )
      ),
      nav_panel(
        "Bac à sable",
        layout_sidebar(
          sidebar = sidebar(
            width = 320, title = "Réglages du bac à sable",
            radioButtons("bac_mech", "Mécanisme de manque",
                         c("MCAR (au hasard)" = "MCAR",
                           "MAR (dépend d'une variable observée)" = "MAR",
                           "NMAR (dépend de la variable manquante)" = "NMAR"),
                         selected = "MAR"),
            sliderInput("bac_k", "Nombre de valeurs manquantes", 1, 10, 5, step = 1),
            conditionalPanel(
              "input.bac_mech == 'MCAR'",
              div(style = "margin:-4px 0 12px;",
                actionButton("bac_retirer", "Autre tirage au hasard", icon = icon("dice"),
                             class = "btn-outline-secondary btn-sm"),
                actionButton("bac_tirage_defaut", "Tirage initial",
                             class = "btn-outline-secondary btn-sm"),
                uiOutput("bac_tirage_txt"))),
            radioButtons("bac_vue", "Vue",
                         c("Avant (les valeurs en cause)" = "avant",
                           "Après (le manque appliqué)" = "apres"),
                         selected = "avant"),
            radioButtons("bac_imput", "Imputation (vue « Après »)",
                         c("Aucune (montrer NA)" = "none",
                           "Par la moyenne" = "moyenne",
                           "Par régression" = "regression",
                           "Par régression + erreur" = "reg_bruit"),
                         selected = "none"),
            conditionalPanel(
              "input.bac_imput == 'regression' || input.bac_imput == 'reg_bruit'",
              radioButtons("bac_modele", "Modèle d'imputation",
                           c("Complet — inclut le prédicteur X2" = "complet",
                             "Incomplet — oublie X2" = "incomplet"),
                           selected = "complet")),
            helpText("Un même petit jeu de 20 lignes. Sous MAR et NMAR, le manque frappe ",
                     "les k valeurs les plus élevées du prédicteur, pour que le motif soit net; ",
                     "sous MCAR, les k lignes sont tirées au sort.")
          ),
          div(class = "encadre", uiOutput("bac_legende")),
          layout_columns(
            col_widths = c(4, 8),
            card(
              card_header("Un petit jeu de 20 observations (4 variables corrélées)"),
              card_body(
                div(style = "overflow-x:auto;", uiOutput("bac_table")),
                div(style = "margin-top:12px; font-size:0.85rem; color:#555; display:flex; flex-wrap:wrap; gap:6px;",
                  tags$span(style = "background:#f6a8a8; padding:2px 8px; border-radius:3px; white-space:nowrap;",
                            "Valeur en cause"),
                  tags$span(style = "background:#fde0e0; border:1px dashed #d95f0e; padding:2px 8px; border-radius:3px; white-space:nowrap;",
                            "NA (manquant)"),
                  tags$span(style = "background:#94c9f3; padding:2px 8px; border-radius:3px; white-space:nowrap;",
                            "Valeur imputée"))
              )
            ),
            div(
              card(
                card_header("Impact sur une analyse simple (ces mêmes 20 lignes)"),
                card_body(
                  div(style = "overflow-x:auto;", uiOutput("bac_resultats")),
                  div(class = "encadre", style = "margin-top:10px;",
                      uiOutput("bac_interpretation"))
                )
              ),
              card(
                card_header("Signature : de quoi dépend le manque ?"),
                card_body(
                  div(style = "overflow-x:auto;", uiOutput("bac_matrice")),
                  div(class = "encadre", style = "margin-top:8px; font-size:0.85rem;",
                    "M(X1) vaut 1 quand X1 manque, 0 sinon. Le tableau donne sa corrélation avec chaque variable. ",
                    tags$b("X2 s'allume"), " \u2192 MAR (prédicteur observé) ; ",
                    tags$b("sa valeur vraie s'allume"), " \u2192 NMAR (visible seulement en simulation).")
                )
              ),
              card(
                card_header("Où tombent les valeurs imputées ?"),
                plotOutput("bac_nuage", height = "300px")
              ),
              card(
                card_header("Pour bien comprendre"),
                card_body(uiOutput("bac_explication"))
              )
            )
          )
        )
      ),
        # -------------------------------------------------------------------- QUIZ
  nav_panel(
    title = "Quiz",
    icon = icon("circle-question"),
    layout_column_wrap(
      width = 1,
      card(
        card_header(textOutput("quiz_titre")),
        card_body(
          uiOutput("quiz_scenario"),
          radioButtons("quiz_reponse", "Ce mécanisme est :",
                       choices = c("MCAR", "MAR", "NMAR"),
                       selected = character(0), inline = TRUE),
          div(class = "mt-2",
            actionButton("quiz_valider", "Valider",
                         icon = icon("check"), class = "btn-primary"),
            actionButton("quiz_suivant", "Question suivante",
                         icon = icon("arrow-right"), class = "btn-outline-secondary ms-2")),
          uiOutput("quiz_feedback")
        )
      ),
      card(
        card_header("Score"),
        card_body(
          uiOutput("quiz_score"),
          actionButton("quiz_recommencer", "Recommencer (mélanger les questions)",
                       icon = icon("rotate"), class = "btn-outline-secondary")
        )
      )
    )
  ),
    # ---------------------------------------------------------- BIBLIOGRAPHIE -
  nav_panel(
    title = "Bibliographie",
    icon = icon("graduation-cap"),
    layout_column_wrap(
      width = 1,
      card(
        card_header("Ouvrages de référence"),
        card_body(
          tags$ul(lapply(BIBLIO$ouvrages,
                         function(r) tags$li(r, style = "margin-bottom:6px;")))
        )
      ),
      card(
        card_header("Articles fondateurs"),
        card_body(
          tags$ul(lapply(BIBLIO$articles,
                         function(r) tags$li(r, style = "margin-bottom:6px;")))
        )
      ),
      card(
        card_header("Ressources en ligne"),
        card_body(
          tags$ul(lapply(BIBLIO$ressources, function(r)
            tags$li(tags$a(r$titre, href = r$url, target = "_blank"),
                    style = "margin-bottom:6px;")))
        )
      )
    )
  )
)),

  # --------------------------------------------- DONNÉES & VISUALISATION ---
  nav_panel(
    title = "Données & visualisation",
    icon = icon("table-cells"),
    layout_sidebar(
      sidebar = sidebar(
        width = 330,
        title = "Simulation",
        accordion(
          open = c("Structure", "Manque"),
          accordion_panel(
            "Structure", icon = icon("cubes"),
            sliderInput("N", "Taille d'échantillon (N)", 50, 2000, 400, step = 50, sep = " "),
            sliderInput("p", "Nombre de variables", 2, 8, 4, step = 1),
            sliderInput("rho", "Corrélation entre variables (ρ)",
                        0, 0.9, 0.5, step = 0.05),
            numericInput("graine", "Graine aléatoire", value = sample.int(1e5, 1), min = 1, step = 1),
            checkboxInput("graine_fixe", "Fixer la graine (reproductible)", FALSE)
          ),
          accordion_panel(
            "Manque", icon = icon("eraser"),
            radioButtons("mech", "Mécanisme de manque",
                         c("MCAR" = "MCAR", "MAR" = "MAR", "NMAR (MNAR)" = "MNAR"),
                         selected = "MAR",
                         inline = TRUE),
            sliderInput("prop", "Proportion de cas incomplets",
                        0.05, 0.80, 0.35, step = 0.05),
            checkboxInput("pilote_on",
                          "Mode guidé : rendre le mécanisme visible", TRUE),
            conditionalPanel(
              "input.pilote_on",
              radioButtons("pilote_portee", "Portée du manque",
                           c("Une seule" = "une",
                             "Plusieurs" = "plusieurs"),
                           selected = "une",
                           inline = TRUE),
              conditionalPanel(
                "input.pilote_portee == 'une'",
                selectInput("pilote_cible", "Variable rendue manquante",
                            choices = NULL)),
              conditionalPanel(
                "input.pilote_portee == 'plusieurs'",
                checkboxGroupInput("pilote_cibles", "Variables rendues manquantes",
                                   choices = NULL, inline = TRUE),
                radioButtons("pilote_patron", "Patron de non-réponse",
                             c("Réponses ponctuelles réparties" = "ponctuel",
                               "Omission groupée d'un bloc" = "bloc"),
                             selected = "ponctuel"),
                helpText("Le premier patron répartit les omissions entre les variables; ",
                         "le second représente une section entière non répondue.")),
              conditionalPanel(
                "input.mech == 'MAR'",
                selectInput("pilote_var", "Prédicteur du manque",
                            choices = NULL)),
              conditionalPanel(
                "input.mech != 'MCAR'",
                radioButtons("pilote_dir", "Sens du mécanisme",
                             c("Valeurs élevées vers manque" = "RIGHT",
                               "Valeurs faibles vers manque" = "LEFT"),
                             selected = "RIGHT"))),
            helpText("Ce réglage vise la proportion de CAS touchés par au moins une ",
                     "valeur manquante. En mode multi-variables, le % de cellules ",
                     "manquantes dépend aussi du patron choisi et est affiché en haut.")
          )
        ),
        actionButton("simuler", "Simuler / rééchantillonner",
                     icon = icon("dice"), class = "btn-primary w-100")
      ),
      layout_column_wrap(
        width = 1/2, heights_equal = "row",
        value_box(title = "Cellules manquantes",
                  value = textOutput("vb_pct_cell"),
                  showcase = icon("percent"), theme = "warning"),
        value_box(title = "Cas incomplets",
                  value = textOutput("vb_pct_cas"),
                  showcase = icon("user-slash"), theme = "info")
      ),
      uiOutput("simulation_note"),
      navset_tab(
        # ---- Sous-onglet : Données ----
        nav_panel(
          "Données",
          card(
            full_screen = TRUE,
            card_header("Jeu de données amputé (les cellules vides sont manquantes)"),
            DTOutput("table_donnees")
          ),
          card(
            card_header("Résumé numérique : valeurs observées et vérité simulée"),
            DTOutput("resume_donnees")
          )
        ),
        # ---- Sous-onglet : Structure du manque ----
        nav_panel(
          "Structure du manque",
          layout_column_wrap(
            width = 1/2,
            card(card_header("Carte des valeurs manquantes"),
                 plotOutput("plot_carte", height = "420px")),
            card(
              card_header("Patrons de manquants"),
              card_body(
                p("Un ", strong("patron"), " est une combinaison de variables ",
                  "présentes et manquantes. Dans le tableau, ", strong("1 = observé"),
                  " et ", strong("0 = manquant"), " ; l'effectif indique combien ",
                  "d'observations suivent ce patron."),
                DTOutput("table_patrons")
              )
            )
          ),
          layout_column_wrap(
            width = 1/2,
            card(card_header("Fréquence des patrons"),
                 plotOutput("plot_patrons", height = "340px")),
            card(card_header("% manquant par variable"),
                 plotOutput("plot_pctvar", height = "340px"))
          )
        ),
        # ---- Sous-onglet : Mécanisme et traitement ----
        nav_panel(
          "Mécanisme et traitement",
          card(
            card_header("Manquants et présents selon une autre variable"),
            card_body(
              layout_column_wrap(
                width = 1/2,
                selectInput("viz_ymanq", "Variable dont on examine le manque",
                            choices = NULL),
                selectInput("viz_x", "Variable en abscisse", choices = NULL)
              ),
              plotOutput("viz_scatter", height = "380px"),
              div(class = "encadre",
                "Les points ", strong("orange"), " sont les cas où la variable examinée manque. ",
                "S'ils se regroupent le long de l'", strong("axe horizontal"), " (une plage de l'abscisse), ",
                "le manque dépend de cette variable observée : signature MAR. ",
                "S'ils se regroupent le long de l'", strong("axe vertical"), " (une plage de la valeur vraie), ",
                "le manque dépend de la valeur manquante elle-même : signature MNAR. ",
                "Cette seconde lecture n'est possible qu'en simulation, car l'axe vertical montre la ",
                strong("valeur vraie"), ", même là où elle a été effacée.")
            )
          ),
          card(
            card_header("Effet d'une imputation sur la distribution"),
            layout_sidebar(
              sidebar = sidebar(
                position = "right", width = 260,
                selectInput("viz_var", "Variable", choices = NULL),
                radioButtons("viz_meth", "Méthode d'imputation illustrée",
                             c("Par la moyenne" = "moyenne",
                               "Par régression (déterministe)" = "reg_det",
                               "Par régression (stochastique)" = "regression",
                               "Imputation multiple (1 jeu)" = "multiple"),
                             selected = "moyenne")
              ),
              plotOutput("plot_imput", height = "360px"),
              div(class = "encadre", textOutput("txt_imput"))
            )
          )
        )
      )
    )
  ),

  # ------------------------------------------------------------ TRAITEMENTS -
  nav_panel(
    title = "Traitements",
    icon = icon("wrench"),
    layout_sidebar(
      sidebar = sidebar(
        width = 320, title = "Analyse cible",
        radioButtons("cible", "Quantité à estimer",
                     c("Moyenne d'une variable" = "moyenne",
                       "Corrélation entre deux variables" = "correlation",
                       "Coefficients de régression" = "regression"),
                     selected = "regression"),
        conditionalPanel(
          "input.cible == 'moyenne'",
          selectInput("t_var", "Variable", choices = NULL)
        ),
        conditionalPanel(
          "input.cible == 'correlation'",
          selectInput("t_v1", "Variable 1", choices = NULL),
          selectInput("t_v2", "Variable 2", choices = NULL)
        ),
        conditionalPanel(
          "input.cible == 'regression'",
          selectInput("t_y", "Variable dépendante", choices = NULL),
          selectInput("t_x", "Prédicteurs", choices = NULL, multiple = TRUE)
        ),
        sliderInput("m_imp", "Nombre d'imputations (m)", 5, 50, 10, step = 5),
        actionButton("lancer", "Comparer les méthodes",
                     icon = icon("play"), class = "btn-primary w-100"),
        hr(),
        downloadButton("dl_rapport", "Télécharger le rapport HTML",
                       class = "btn-outline-secondary w-100")
      ),
      uiOutput("traitements_pret"),
      card(
        card_header("Comparaison des estimations"),
        DTOutput("table_comp")
      ),
      uiOutput("bloc_biais_interpretation")
    )
  ),

  # ------------------------------------------------------- MÉCANISME DU MANQUE
  nav_panel(
    title = "Mécanisme du manque",
    icon = icon("stethoscope"),
    navset_tab(
      # ---- Sous-onglet : Diagnostic ----
      nav_panel(
        "Diagnostic",
        layout_column_wrap(
          width = 1/3,
          card(
            card_header("Test de Little (MCAR) sur les données courantes"),
            card_body(uiOutput("little_resultat"))
          ),
          card(
            card_header("Démonstration : ce que le test peut — et ne peut pas — détecter"),
            card_body(
              p("Le même test appliqué à des jeux simulés (cible X1, N = 500, ",
                "proportion et ρ courants) sous chaque mécanisme :"),
              DTOutput("little_demo"),
              div(class = "alerte",
                strong("La leçon. "),
                "Le test distingue MCAR de « non-MCAR », mais il ne distingue ",
                strong("jamais"), " MAR de MNAR : les deux rejettent MCAR dès qu'une ",
                "structure est détectable. Et lorsque la variable manquante est ",
                "décorrélée des variables observées (dernière ligne), un mécanisme MNAR ",
                "devient totalement invisible — indiscernable de MCAR.")
            )
          ),
          card(
            card_header("La limite fondamentale"),
            card_body(
              p("Aucun test fondé sur les seules données observées ne peut confirmer MAR ",
                "plutôt que MNAR, car MNAR dépend précisément des valeurs qui manquent. ",
                "En pratique :"),
              tags$ul(
                tags$li(strong("MCAR rejeté : "), "le manque n'est pas complètement ",
                        "aléatoire — mais on ignore s'il est MAR ou MNAR."),
                tags$li(strong("MCAR non rejeté : "), "compatible avec MCAR, sans preuve ",
                        "(faible puissance, ou MNAR décorrélé indétectable)."),
                tags$li(strong("Décision : "), "le choix MAR vs MNAR relève d'un ",
                        em("argument de fond"), " sur le processus de manque, jamais d'un test.")
              )
            )
          )
        ),
        card(
          card_header("Matrice manque \u2194 variables \u2014 o\u00f9 sont les liens ?"),
          card_body(
            p("Chaque ligne suit l'indicateur de manque M d'une variable (M = 1 si la valeur manque, 0 sinon). ",
              "Chaque cellule donne la corr\u00e9lation (point-bis\u00e9riale) entre cet indicateur et les valeurs ",
              "observ\u00e9es d'une autre variable ; l'intensit\u00e9 de la couleur suit la force du lien. ",
              "Le second tableau (\u03c6) indique quelles variables manquent ensemble."),
            div(style = "overflow-x:auto;", uiOutput("diag_matrice")),
            div(class = "alerte", style = "margin-top:12px;",
              strong("Comment lire. "),
              "Une colonne qui s'allume d\u00e9signe un ", strong("pr\u00e9dicteur du manque"),
              " \u2014 \u00e0 inclure dans le mod\u00e8le d'imputation. La derni\u00e8re colonne, ",
              em("sa valeur vraie"), " (connue seulement en simulation), r\u00e9v\u00e8le la ",
              "signature : un lien fort sur une variable ", strong("observ\u00e9e"),
              " pointe vers MAR ; sur ", em("sa propre valeur vraie"), ", vers MNAR. Mais, comme ",
              "le test de Little, elle ne prouve jamais MAR contre MNAR \u00e0 partir des seules ",
              "donn\u00e9es observ\u00e9es.")
          )
        )
      ),
      # ---- Sous-onglet : Sensibilité MNAR ----
      nav_panel(
        "Sensibilité MNAR",
        layout_sidebar(
          sidebar = sidebar(
            width = 330, title = "Delta-ajustement",
            selectInput("sv_var", "Variable suspectée MNAR", choices = NULL),
            sliderInput("sv_range", "Amplitude de δ explorée (± écarts-types)",
                        0.5, 3, 1.5, step = 0.5),
            sliderInput("sv_m", "Nombre d'imputations (m)", 5, 30, 10, step = 5),
            actionButton("sv_lancer", "Lancer l'analyse",
                         icon = icon("play"), class = "btn-primary w-100"),
            helpText("On impute sous MAR, puis on décale les valeurs imputées de δ ",
                     "écarts-types. δ encode l'ampleur du manque non aléatoire redouté. ",
                     "δ = 0 correspond exactement à l'hypothèse MAR.")
          ),
          div(class = "encadre",
            "Le mécanisme MNAR étant indétectable (voir le sous-onglet Diagnostic), on ne ",
            "peut pas le ", strong("tester"), " — mais on peut en ",
            strong("encadrer la conséquence"), ". On fait varier δ sur une plage plausible ",
            "et l'on observe si les conclusions restent stables. En simulation, la courbe ",
            "traverse la vraie moyenne au δ correspondant au véritable mécanisme."),
          layout_columns(
            card(
              card_header("Courbe de sensibilité"),
              plotOutput("sv_plot", height = "420px")
              ),
            card(
              card_header("Interprétation"),
              card_body(uiOutput("sv_interpretation"))
            ),
            col_widths = c(6, 6)
          )
        )
      )
    )
  ),

  # -------------------------------------------------------- ÉVALUER ET VALIDER
  nav_panel(
    title = "Évaluer et valider",
    icon = icon("clipboard-check"),
    navset_tab(
      # ---- Sous-onglet : Étude de simulation ----
      nav_panel(
        "Étude de simulation",
        layout_sidebar(
          sidebar = sidebar(
            width = 340, title = "Réplications Monte Carlo",
            radioButtons("mc_cible", "Quantité étudiée",
                         c("Moyenne d'une variable" = "moyenne",
                           "Corrélation entre deux variables" = "correlation",
                           "Coefficient de régression" = "regression"),
                         selected = "correlation"),
            conditionalPanel(
              "input.mc_cible == 'moyenne'",
              selectInput("mc_var", "Variable", choices = NULL)),
            conditionalPanel(
              "input.mc_cible == 'correlation'",
              selectInput("mc_v1", "Variable 1", choices = NULL),
              selectInput("mc_v2", "Variable 2", choices = NULL)),
            conditionalPanel(
              "input.mc_cible == 'regression'",
              selectInput("mc_y", "Variable dépendante", choices = NULL),
              selectInput("mc_x", "Prédicteurs", choices = NULL, multiple = TRUE),
              selectInput("mc_terme", "Coefficient suivi", choices = NULL)),
            sliderInput("mc_R", "Nombre de réplications (R)", 50, 500, 200, step = 50),
            sliderInput("mc_m", "Imputations par réplication (m)", 5, 20, 5, step = 5),
            actionButton("mc_lancer", "Lancer l'étude",
                         icon = icon("play"), class = "btn-primary w-100"),
            helpText("Chaque réplication simule, ampute et traite un nouveau jeu ",
                     "(graine variée à chaque fois). Les réglages N, ρ, mécanisme et ",
                     "proportion proviennent de l'onglet Données & visualisation.")
          ),
          div(class = "encadre",
            "On répète ", strong("R"), " fois le cycle simuler → amputer → traiter. ",
            "La ", strong("dispersion"), " des estimations mesure l'efficacité, leur ",
            strong("décalage"), " par rapport à la vérité mesure le biais, et le ",
            strong("taux de couverture"), " révèle si les erreurs-types sont honnêtes."),
          layout_columns(
            card(
              card_header("Distribution d'échantillonnage des estimations"),
              plotOutput("mc_plot_dist", height = "360px")
              ),
            card(
              card_header("Couverture des intervalles de confiance à 95 %"),
              plotOutput("mc_plot_couv", height = "360px")
              ),
            col_widths = c(6, 6)
          ),
          card(
            card_header("Voir la couverture : chaque intervalle, réplication par réplication"),
            plotOutput("mc_plot_caterpillar", height = "440px")
          ),
          layout_columns(
            card(
              card_header("Synthèse"),
              DTOutput("mc_table")
              ),
            card(
              card_header("Interprétation"),
              uiOutput("mc_interpretation")
              ),
            col_widths = c(6, 6)
            )
        )
      ),

      # ---- Sous-onglet : Qualité de l'imputation ----
      nav_panel(
        "Qualité de l'imputation",
        layout_sidebar(
          sidebar = sidebar(
            width = 320, title = "Imputation à diagnostiquer",
            selectInput("qi_var", "Variable imputée à examiner", choices = NULL),
            sliderInput("qi_m", "Nombre d'imputations (m)", 5, 30, 10, step = 5),
            sliderInput("qi_maxit", "Itérations MICE (maxit)", 5, 40, 20, step = 5),
            actionButton("qi_lancer", "Lancer l'imputation",
                         icon = icon("play"), class = "btn-primary w-100"),
            helpText("L'imputation est calculée une fois puis diagnostiquée sous trois ",
                     "angles. Les variables listées sont celles qui comportent des ",
                     "valeurs manquantes.")
          ),
          div(class = "encadre",
            "En simulation, on connaît la vérité ; sur de vraies données, non. Ces trois ",
            "diagnostics permettent de ", strong("vérifier une imputation sans vérité de référence"),
            " : les valeurs imputées sont-elles plausibles, l'algorithme a-t-il convergé, ",
            "et le nombre d'imputations est-il suffisant ?"),
          card(card_header("1. Valeurs observées vs imputées"),
               plotOutput("qi_plot_obsimp", height = "360px"),
               div(class = "encadre", textOutput("qi_txt_obsimp"))),
          card(card_header("2. Convergence de l'algorithme"),
               plotOutput("qi_plot_conv", height = "380px")),
          card(card_header("3. Effet du nombre d'imputations"),
               plotOutput("qi_plot_effetm", height = "340px"))
        )
      )
    )
  ),

  # ------------------------------------------------------------ DONNÉES RÉELLES
  nav_panel(
    title = "Données réelles",
    icon = icon("database"),
    layout_sidebar(
      sidebar = sidebar(
        width = 330, title = "Jeu réel",
        selectInput("rd_jeu", "Jeu de données",
                    c("Qualité de l'air (airquality)" = "airquality",
                      "NHANES — santé (mice)" = "nhanes",
                      "Croissance de garçons (mice)" = "boys")),
        selectInput("rd_var", "Variable pour observé/imputé", choices = NULL),
        selectInput("rd_x", "Variable en abscisse (nuage)", choices = NULL),
        sliderInput("rd_m", "Nombre d'imputations (m)", 5, 20, 10, step = 5),
        actionButton("rd_lancer", "Imputer et diagnostiquer",
                     icon = icon("play"), class = "btn-primary w-100"),
        helpText("Sur données réelles, aucune vérité n'existe : impossible de mesurer ",
                 "le biais. On peut néanmoins visualiser le manque, tester MCAR et ",
                 "diagnostiquer l'imputation — la démarche exacte du terrain.")
      ),
      uiOutput("rd_desc"),
      card(
        card_header("Démarche terrain : ce que l'on peut vraiment conclure"),
        card_body(
          layout_column_wrap(
            width = 1/5,
            div(strong("1. Décrire"), br(),
                "Quantifier et cartographier les manquants."),
            div(strong("2. Expliquer"), br(),
                "Documenter le processus de collecte et les variables associées au manque."),
            div(strong("3. Traiter"), br(),
                "Inclure dans l'imputation l'issue, les prédicteurs et les auxiliaires pertinents."),
            div(strong("4. Vérifier"), br(),
                "Examiner plausibilité, convergence et stabilité de l'imputation."),
            div(strong("5. Sensibiliser"), br(),
                "Évaluer si une hypothèse MNAR raisonnable changerait la conclusion."))
        )
      ),
      layout_column_wrap(
        width = 1/2, heights_equal = "row",
        value_box(title = "Cellules manquantes", value = textOutput("rd_vb_cell"),
                  showcase = icon("percent"), theme = "warning"),
        value_box(title = "Cas incomplets", value = textOutput("rd_vb_cas"),
                  showcase = icon("user-slash"), theme = "info")
      ),
      card(full_screen = TRUE, card_header("Aperçu des données"),
           DTOutput("rd_table")),
      card(card_header("Carte des valeurs manquantes"),
           plotOutput("rd_carte", height = "400px")),
      layout_column_wrap(
        width = 1/2,
        card(card_header("Patrons de manquants"),
             plotOutput("rd_patrons", height = "340px")),
        card(card_header("% manquant par variable"),
             plotOutput("rd_pctvar", height = "340px"))
      ),
      card(card_header("Test de Little (MCAR)"),
           card_body(uiOutput("rd_little"))),
      layout_column_wrap(
        width = 1/2,
        card(card_header("Observé vs imputé : distributions"),
             plotOutput("rd_obsimp", height = "360px")),
        card(card_header("Observé vs imputé : selon une autre variable"),
             plotOutput("rd_nuage", height = "360px"))),
      div(class = "encadre",
        "Cliquez « Imputer et diagnostiquer ». Sans vérité, on juge la ",
        "plausibilité : les valeurs imputées doivent rester dans le registre ",
        "de l'observé et respecter les relations connues (par exemple, le volume ",
        "testiculaire croît avec l'âge dans ", em("boys"), ")."),
      card(card_header("Ce que ce jeu nous apprend"),
           card_body(uiOutput("rd_lecon")))
    )
  ),

  nav_spacer(),
  nav_item(tags$a(icon("globe"), " dcstats.org",
                  href = "https://dcstats.org", target = "_blank"))
)

# ============================================================================
#  5. SERVEUR
# ============================================================================
server <- function(input, output, session) {

  rv <- reactiveValues(complet = NULL, ampute = NULL,
                       pct_cell = NA, pct_cas = NA,
                       var_cible = NULL, var_pilote = NULL,
                       patron_manque = NULL, params = NULL,
                       simulation_id = 0L)
  bac_tirage <- reactiveVal(BAC_TIRAGE_DEFAUT)

  # Réglages qui définissent une simulation (comparés aux réglages courants)
  reglages_courants <- function() list(
    N = input$N, p = input$p, rho = input$rho, mech = input$mech, prop = input$prop,
    pilote_on = isTRUE(input$pilote_on), portee = input$pilote_portee,
    cible = input$pilote_cible, cibles = sort(as.character(input$pilote_cibles)),
    patron = input$pilote_patron, pilote = input$pilote_var,
    dir = if (is.null(input$pilote_dir)) "RIGHT" else input$pilote_dir)

  generer <- function() {
    graine <- if (isTRUE(input$graine_fixe)) as.integer(input$graine) else sample.int(1e5, 1)
    if (!isTRUE(input$graine_fixe)) updateNumericInput(session, "graine", value = graine)
    cplt <- simuler_complet(input$N, input$p, input$rho, graine = graine)
    vars_init <- names(cplt)
    portee <- if (identical(input$pilote_portee, "plusieurs")) "plusieurs" else "une"
    cible <- if (!isTRUE(input$pilote_on)) {
      NULL
    } else if (identical(portee, "plusieurs")) {
      unique(input$pilote_cibles)
    } else {
      input$pilote_cible
    }
    cible <- unique(cible[cible %in% vars_init])
    if (isTRUE(input$pilote_on) && !length(cible)) {
      cible <- if (identical(portee, "plusieurs")) vars_init[seq_len(min(2L, length(vars_init)))] else vars_init[1]
    }
    pilote <- if (isTRUE(input$pilote_on)) input$pilote_var else NULL
    if (!is.null(pilote) && (!(pilote %in% vars_init) || pilote %in% cible)) pilote <- NULL
    if (isTRUE(input$pilote_on) && identical(input$mech, "MAR") && is.null(pilote)) {
      pilote <- setdiff(vars_init, cible)[1]
    }
    if (isTRUE(input$pilote_on) && identical(input$mech, "MAR") &&
        (is.null(pilote) || is.na(pilote))) {
      showNotification("Sous MAR, laissez au moins une variable observée comme prédicteur du manque.",
                       type = "warning")
      return(invisible(NULL))
    }
    amp <- amputer(cplt, prop = input$prop, mech = input$mech,
                   graine = graine,
                   var_cible  = cible,
                   var_pilote = pilote,
                   direction  = if (is.null(input$pilote_dir)) "RIGHT" else input$pilote_dir,
                   patron = if (identical(portee, "plusieurs") && !is.null(input$pilote_patron))
                     input$pilote_patron else "ponctuel")
    rv$complet <- cplt
    rv$ampute  <- amp$amp
    rv$pct_cell <- amp$pct_cellules
    rv$pct_cas  <- amp$pct_cas
    rv$var_cible  <- amp$var_cible
    rv$var_pilote <- amp$var_pilote
    rv$patron_manque <- if (identical(portee, "plusieurs") && identical(input$pilote_patron, "bloc")) "bloc" else "ponctuel"
    params <- c(reglages_courants(), list(graine = graine))
    rv$simulation_id <- rv$simulation_id + 1L
    vars <- names(cplt)
    updateSelectInput(session, "viz_var", choices = vars, selected = vars[1])
    na_vars <- names(cplt)[colSums(is.na(amp$amp)) > 0]
    if (length(na_vars) == 0) na_vars <- names(cplt)
    updateSelectInput(session, "viz_ymanq", choices = na_vars, selected = na_vars[1])
    updateSelectInput(session, "viz_x", choices = vars,
                      selected = setdiff(vars, na_vars[1])[1])
    updateSelectInput(session, "t_var", choices = vars, selected = vars[1])
    updateSelectInput(session, "t_v1", choices = vars, selected = vars[1])
    updateSelectInput(session, "t_v2", choices = vars, selected = vars[2])
    updateSelectInput(session, "t_y", choices = vars, selected = vars[1])
    updateSelectInput(session, "t_x", choices = vars, selected = vars[-1])
    updateSelectInput(session, "mc_var", choices = vars, selected = vars[1])
    updateSelectInput(session, "mc_v1", choices = vars, selected = vars[1])
    updateSelectInput(session, "mc_v2", choices = vars, selected = vars[2])
    updateSelectInput(session, "mc_y", choices = vars, selected = vars[1])
    updateSelectInput(session, "mc_x", choices = vars, selected = vars[-1])
    updateSelectInput(session, "mc_terme",
                      choices = c("(Intercept)", vars[-1]), selected = vars[2])
    cible_unique <- if (length(cible) == 1L) cible else vars[1]
    cibles_multi <- if (length(cible) > 1L) cible else vars[seq_len(min(2L, length(vars)))]
    # Le prédicteur MAR doit rester observé : on exclut la ou les variables amputées
    pilotes_possibles <- setdiff(vars, if (identical(portee, "plusieurs")) cibles_multi else cible_unique)
    pilote_unique <- if (!is.null(pilote) && pilote %in% pilotes_possibles) pilote else pilotes_possibles[1]
    updateSelectInput(session, "pilote_cible", choices = vars, selected = cible_unique)
    updateCheckboxGroupInput(session, "pilote_cibles", choices = vars, selected = cibles_multi)
    updateSelectInput(session, "pilote_var", choices = pilotes_possibles, selected = pilote_unique)
    # On mémorise les valeurs effectivement utilisées (celles que les menus affichent désormais)
    params$cible <- cible_unique; params$cibles <- sort(cibles_multi); params$pilote <- pilote_unique
    rv$params <- params
    vars_na <- names(cplt)[colSums(is.na(amp$amp)) > 0]
    if (length(vars_na) == 0) vars_na <- names(cplt)
    updateSelectInput(session, "qi_var", choices = vars_na, selected = vars_na[1])
    updateSelectInput(session, "sv_var", choices = vars_na, selected = vars_na[1])
  }

  observeEvent(input$simuler, generer(), ignoreNULL = FALSE)

  # Les réglages ont-ils changé depuis la dernière simulation ?
  reglages_modifies <- reactive({
    req(rv$params)
    cur <- reglages_courants()
    old <- rv$params[names(cur)]
    cles <- c("N", "p", "rho", "mech", "prop", "pilote_on")
    if (isTRUE(cur$pilote_on)) cles <- c(cles, "portee", "dir",
      if (identical(cur$portee, "plusieurs")) c("cibles", "patron") else "cible",
      if (identical(cur$mech, "MAR")) "pilote")
    !isTRUE(all.equal(cur[cles], old[cles], check.attributes = FALSE))
  })
  bandeau_perime <- function() {
    if (!isTRUE(tryCatch(reglages_modifies(), error = function(e) FALSE))) return(NULL)
    div(class = "alerte", icon("triangle-exclamation"), strong(" Réglages modifiés. "),
        sprintf("Les résultats affichés correspondent encore à la simulation précédente (%s, N = %d, proportion de cas incomplets = %s). ",
                rv$params$mech, as.integer(rv$params$N), fmt_fr(rv$params$prop, 2)),
        "Cliquez « Simuler / rééchantillonner » pour les mettre à jour.")
  }

  # En mode MAR multi-variables, le prédicteur doit demeurer observé.
  observeEvent(list(input$pilote_cibles, input$pilote_portee, input$mech), {
    req(rv$complet)
    if (!identical(input$pilote_portee, "plusieurs") || !identical(input$mech, "MAR")) return()
    disponibles <- setdiff(names(rv$complet), input$pilote_cibles)
    selection <- if (!is.null(input$pilote_var) && input$pilote_var %in% disponibles) {
      input$pilote_var
    } else {
      disponibles[1]
    }
    updateSelectInput(session, "pilote_var", choices = disponibles, selected = selection)
  }, ignoreInit = TRUE)


  # ---- Données ----
  # Tirages successifs reproductibles : n° 1 (initial), 2, 3, ...
  observeEvent(input$bac_retirer, bac_tirage(bac_tirage() + 1L))
  observeEvent(input$bac_tirage_defaut, bac_tirage(BAC_TIRAGE_DEFAUT))
  output$bac_tirage_txt <- renderUI(tags$small(style = "color:#59636d;",
    if (bac_tirage() == BAC_TIRAGE_DEFAUT) "Tirage n° 1 (initial)." else
      sprintf("Tirage n° %d.", bac_tirage() - BAC_TIRAGE_DEFAUT + 1L)))
  bac_modele <- reactive(if (is.null(input$bac_modele)) "complet" else input$bac_modele)

  output$bac_table <- renderUI({
    table_bac(input$bac_mech, input$bac_k, input$bac_vue, input$bac_imput, bac_modele(), bac_tirage())
  })
  output$bac_legende <- renderUI({
    bac_legende(input$bac_mech, input$bac_vue, input$bac_imput, bac_modele())
  })

  output$bac_resultats <- renderUI({
    table_analyse(input$bac_mech, input$bac_k, input$bac_imput, bac_modele(), bac_tirage())
  })
  output$bac_interpretation <- renderUI({
    bac_interpretation(input$bac_mech, if (identical(input$bac_vue, "apres")) input$bac_imput else "none",
                       bac_modele())
  })

  output$bac_explication <- renderUI({
    bac_explication(input$bac_mech)
  })
  output$bac_matrice <- renderUI({
    amp <- bac_donnees
    amp$X1[bac_lignes_manq(input$bac_mech, input$bac_k, bac_tirage())] <- NA
    table_matrice_manque(bac_donnees, amp)
  })
  output$bac_nuage <- renderPlot({
    graphe_bac_nuage(input$bac_mech, input$bac_k,
                     if (identical(input$bac_vue, "apres")) input$bac_imput else "none",
                     bac_modele(), bac_tirage())
  })
  
  output$vb_pct_cell <- renderText(pct_fr(100 * rv$pct_cell))
  output$vb_pct_cas  <- renderText(pct_fr(100 * rv$pct_cas))

  output$simulation_note <- renderUI({
    req(rv$ampute, rv$params)
    tagList(bandeau_perime(), note_simulation())
  })
  note_simulation <- function() {
    cible <- rv$var_cible
    mech <- rv$params$mech
    if (is.null(cible)) {
      return(div(class = "encadre",
        strong("Mode diffus. "),
        "Le manque est réparti entre plusieurs variables. Ce mode est utile pour explorer des patrons complexes, mais le mécanisme est moins facile à lire au premier coup d'œil."))
    }
    if (length(cible) > 1L) {
      vars_txt <- paste(cible, collapse = ", ")
      patron_txt <- if (identical(rv$patron_manque, "bloc")) {
        "Les mêmes répondants omettent ce bloc de questions."
      } else {
        "Les omissions sont réparties entre ces variables : un même répondant ne perd pas nécessairement toutes ses réponses."
      }
      texte <- switch(mech,
        MCAR = paste0("Des valeurs de ", vars_txt, " sont retirées au hasard. ", patron_txt),
        MAR = paste0("Des valeurs de ", vars_txt, " deviennent manquantes selon ",
                     rv$var_pilote, ", qui demeure observée. ", patron_txt),
        MNAR = paste0("Les valeurs de ", vars_txt,
                      " deviennent manquantes selon leurs propres valeurs vraies. ", patron_txt,
                      " Cette dépendance est visible ici seulement parce que l'application connaît la vérité simulée.")
      )
      return(div(class = if (identical(mech, "MNAR")) "alerte" else "encadre",
        strong("Mode multi-variables. "), texte,
        " Changez les réglages, puis cliquez « Simuler / rééchantillonner »."))
    }
    sens <- if (identical(rv$params$dir, "LEFT")) "faibles" else "élevées"
    texte <- switch(mech,
      MCAR = paste0("Mode guidé : des valeurs de ", cible,
                    " sont retirées au hasard. Aucune variable ne pilote leur absence."),
      MAR = paste0("Mode guidé : des valeurs de ", cible, " deviennent manquantes lorsque ",
                   rv$var_pilote, " prend des valeurs ", sens,
                   ". Le prédicteur du manque reste observé."),
      MNAR = paste0("Mode guidé : des valeurs ", sens, " de ", cible,
                    " deviennent elles-mêmes manquantes. Cette information est visible ici seulement parce que l'application connaît la vérité simulée.")
    )
    div(class = if (identical(mech, "MNAR")) "alerte" else "encadre",
        strong("Ce que vous avez simulé. "), texte,
        " Changez les réglages, puis cliquez « Simuler / rééchantillonner ».")
  }

  output$table_donnees <- renderDT({
    req(rv$ampute)
    datatable(rv$ampute, options = list(pageLength = 8, scrollX = TRUE),
              rownames = TRUE) |>
      formatRound(names(rv$ampute), 3, dec.mark = ",") |>
      formatStyle(names(rv$ampute), backgroundColor = "white")
  })

  output$resume_donnees <- renderDT({
    req(rv$ampute)
    a <- rv$ampute; cp <- rv$complet
    tab <- data.frame(
      Variable = names(a),
      Manquantes = colSums(is.na(a)),
      `% manquant` = 100 * colMeans(is.na(a)),
      `Moyenne (observées)` = vapply(a, mean, numeric(1), na.rm = TRUE),
      `Moyenne (vérité)` = vapply(cp, mean, numeric(1)),
      `Écart-type (observées)` = vapply(a, sd, numeric(1), na.rm = TRUE),
      `Écart-type (vérité)` = vapply(cp, sd, numeric(1)),
      check.names = FALSE)
    datatable(tab, rownames = FALSE,
              caption = sprintf("Mécanisme %s · N = %d · graine %s. Les colonnes « vérité » ne sont connues qu'en simulation.",
                                rv$params$mech, nrow(a), rv$params$graine),
              options = list(dom = "t", ordering = FALSE)) |>
      formatRound("% manquant", 1, dec.mark = ",") |>
      formatRound(4:7, 3, dec.mark = ",")
  })

  # ---- Visualisation ----
  output$plot_carte  <- renderPlot(carte_manquants(req(rv$ampute)))
  output$plot_patrons <- renderPlot(graphe_patrons(req(rv$ampute)))
  output$plot_pctvar <- renderPlot(graphe_pct_variable(req(rv$ampute)))

  output$table_patrons <- renderDT({
    req(rv$ampute)
    datatable(table_patrons(rv$ampute), rownames = FALSE,
              options = list(dom = "t", ordering = FALSE, pageLength = 25))
  })

  output$viz_scatter <- renderPlot({
    req(rv$complet, rv$ampute, input$viz_ymanq, input$viz_x)
    validate(need(input$viz_ymanq != input$viz_x,
                  "Choisissez deux variables différentes."))
    graphe_manquants_bivarie(rv$complet, rv$ampute, input$viz_ymanq, input$viz_x)
  })

  output$plot_imput <- renderPlot({
    req(rv$ampute, input$viz_var)
    graphe_imputation(rv$ampute, input$viz_var, input$viz_meth,
                      graine = as.integer(rv$params$graine), complet = rv$complet)
  })
  output$txt_imput <- renderText({
    if (identical(input$viz_meth, "moyenne")) {
      "L'imputation par la moyenne concentre toutes les valeurs imputées sur un seul point : la distribution se déforme, la variance chute et les corrélations s'atténuent."
    } else if (identical(input$viz_meth, "reg_det")) {
      "La régression déterministe impute la valeur prédite : les valeurs imputées varient (elles suivent les autres variables), mais restent trop resserrées, car elles tombent toutes exactement sur la droite de régression. Comparez leur écart-type à celui des vraies valeurs effacées."
    } else if (identical(input$viz_meth, "regression")) {
      "La régression stochastique ajoute à la prédiction un résidu aléatoire : les valeurs imputées retrouvent un étalement comparable à celui des vraies valeurs effacées."
    } else {
      "Un jeu issu de l'imputation multiple (appariement prédictif) emprunte des valeurs réellement observées chez des cas voisins : il reproduit l'étalement. Répéter l'opération m fois permet ensuite de mesurer l'incertitude d'imputation."
    }
  })

  # ---- Diagnostic (test de Little) ----
  output$diag_matrice <- renderUI({
    req(rv$complet, rv$ampute)
    table_matrice_manque(rv$complet, rv$ampute)
  })

  output$little_resultat <- renderUI({
    req(rv$ampute)
    t <- little_mcar_test(rv$ampute)
    if (is.na(t$p.value)) {
      return(div(class = "encadre",
                 "Aucune donnée manquante : le test ne s'applique pas."))
    }
    rejet <- t$p.value < 0.05
    verdict <- tags$span(
      style = paste0("font-weight:bold;color:", if (rejet) "#d95f0e" else "#31a354", ";"),
      if (rejet) "MCAR rejeté" else "MCAR non rejeté")
    tags$div(
      bandeau_perime(),
      tags$p(sprintf("Statistique d\u00b2 = %s, ddl = %d (%d patrons de manquants).",
                     fmt_fr(t$statistic, 2), t$df, t$n_patterns)),
      tags$p("Valeur-p = ", tags$b(p_fr(t$p.value)), " \u2192 ", verdict),
      div(class = if (rejet) "alerte" else "encadre",
        if (rejet) {
          "L'hypothèse MCAR est rejetée : le manque n'est pas complètement aléatoire. Le test ne dit toutefois PAS s'il s'agit de MAR ou de MNAR."
        } else {
          "On ne rejette pas MCAR : les données sont compatibles avec un manque complètement aléatoire — sans que cela le prouve."
        })
    )
  })

  output$little_demo <- renderDT({
    req(rv$complet)
    req(rv$params)
    withProgress(message = "Test de Little sous chaque mécanisme…", value = 0.5, {
      tab <- demonstration_little(500L, rv$params$prop, rv$params$rho, as.integer(rv$params$graine))
    })
    tab$p <- p_fr(tab$p)
    names(tab) <- c("Scénario", "d\u00b2", "ddl", "valeur-p", "Verdict")
    datatable(tab, rownames = FALSE, options = list(dom = "t", ordering = FALSE)) |>
      formatRound("d\u00b2", 2, dec.mark = ",")
  })
  # ---- Traitements ----
  comparaison <- eventReactive(input$lancer, {
    req(rv$complet, rv$ampute)
    m <- as.integer(input$m_imp)
    g <- as.integer(rv$params$graine)
    res <- withProgress(message = "Comparaison des méthodes en cours…", value = 0.5, {
    if (input$cible == "moyenne") {
      req(input$t_var)
      tab <- comparer_moyenne(rv$complet, rv$ampute, input$t_var, m, g)
      list(type = "moyenne", tab = tab, var_disp = input$t_var,
           libelle = paste0("moyenne de ", input$t_var),
           disp = effet_dispersion(rv$complet, rv$ampute, input$t_var, attr(tab, "mids"), g))
    } else if (input$cible == "correlation") {
      req(input$t_v1, input$t_v2)
      validate(need(input$t_v1 != input$t_v2,
                    "Choisissez deux variables différentes."))
      tab <- comparer_correlation(rv$complet, rv$ampute, input$t_v1, input$t_v2, m, g)
      list(type = "correlation", tab = tab, var_disp = input$t_v1,
           libelle = paste0("corrélation entre ", input$t_v1, " et ", input$t_v2),
           disp = effet_dispersion(rv$complet, rv$ampute, input$t_v1, attr(tab, "mids"), g))
    } else {
      req(input$t_y, input$t_x)
      validate(need(length(input$t_x) >= 1, "Choisissez au moins un prédicteur."),
               need(!(input$t_y %in% input$t_x),
                    "La variable dépendante ne peut pas être un prédicteur."))
      tab <- comparer_regression(rv$complet, rv$ampute, input$t_y, input$t_x, m, g)
      list(type = "regression", tab = tab, var_disp = input$t_y,
           libelle = paste0("régression de ", input$t_y, " sur ", paste(input$t_x, collapse = ", ")),
           disp = effet_dispersion(rv$complet, rv$ampute, input$t_y, attr(tab, "mids"), g))
    }
    })
    res$simulation_id <- rv$simulation_id
    res$m <- m
    res
  })

  comparaison_courante <- reactive({
    res <- comparaison()
    if (is.null(res) || !identical(res$simulation_id, rv$simulation_id)) return(NULL)
    res
  })

  output$traitements_pret <- renderUI({
    tagList(bandeau_perime(), traitements_message())
  })
  traitements_message <- function() {
    res <- comparaison_courante()
    if (is.null(res)) {
      return(div(class = "encadre",
        strong("Étape suivante. "),
        "Choisissez une quantité à estimer, puis cliquez « Comparer les méthodes ». ",
        "Après toute nouvelle simulation, relancez cette comparaison : les résultats précédents ne correspondent plus au jeu courant."))
    }
    div(class = "cle",
      strong("Résultats du jeu courant. "),
      "La ligne « Vérité » est disponible uniquement parce que les données sont simulées. ",
      "En pratique, on compare plutôt des hypothèses et des diagnostics, jamais une estimation à une vérité observable.")
  }

  output$table_comp <- renderDT({
    res <- comparaison_courante()
    validate(need(!is.null(res), "Lancez une comparaison pour le jeu courant."))
    if (res$type == "regression") {
      tab <- res$tab
      names(tab) <- c("Terme", COLS_REG[setdiff(names(tab), "Terme")])
      datatable(tab, rownames = FALSE,
                options = list(dom = "t", ordering = FALSE, scrollX = TRUE)) |>
        formatRound(setdiff(names(tab), "Terme"), 3, dec.mark = ",") |>
        formatStyle("Vérité", fontWeight = "bold")
    } else {
      tab <- res$tab
      names(tab) <- c("Méthode", "Estimation", "Biais")
      datatable(tab, rownames = FALSE, options = list(dom = "t", ordering = FALSE)) |>
        formatRound(c("Estimation", "Biais"), 3, dec.mark = ",") |>
        formatStyle("Méthode", target = "row",
                    fontWeight = styleEqual("Vérité (données complètes)", "bold"))
    }
  })

  output$plot_comp_biais <- renderPlot({
    res <- comparaison_courante()
    req(res)
    req(res$type != "regression")
    graphe_biais(res$tab)
  })

  output$plot_comp_biais_reg <- renderPlot({
    res <- comparaison_courante()
    req(res)
    req(res$type == "regression")
    graphe_biais_regression(res$tab)
  })

  output$table_disp <- renderDT({
    res <- comparaison_courante()
    validate(need(!is.null(res), "Lancez une comparaison pour le jeu courant."))
    d <- res$disp; names(d) <- c(paste0("Situation (", res$var_disp, ")"), "Écart-type")
    datatable(d, rownames = FALSE, options = list(dom = "t", ordering = FALSE)) |>
      formatRound("Écart-type", 3, dec.mark = ",")
  })

output$bloc_biais_interpretation <- renderUI({
    res <- comparaison_courante()
    if (is.null(res)) return(NULL)
    graphe <- if (identical(res$type, "regression")) {
      card(card_header("Biais par méthode et par coefficient"),
           plotOutput("plot_comp_biais_reg", height = "380px"))
    } else {
      card(card_header("Biais par méthode"),
           plotOutput("plot_comp_biais", height = "340px"))
    }
    tagList(
      layout_columns(
        graphe,
        card(card_header("Effet sur la dispersion (écart-type)"), DTOutput("table_disp")),
        col_widths = c(7, 5), fill = FALSE),
      card(card_header("Interprétation"), uiOutput("interpretation"))
    )
  })
  output$interpretation <- renderUI({
    res <- comparaison_courante()
    req(res)
    mech <- rv$params$mech
    intro <- switch(mech,
      MCAR = "Sous MCAR, les cas complets restent représentatifs : on attend surtout une perte de précision, pas de biais systématique.",
      MAR  = "Sous MAR, supprimer les cas incomplets peut biaiser l'estimation. L'imputation multiple et le FIML peuvent corriger ce biais s'ils exploitent les variables observées dont dépend le manque.",
      MNAR = "Sous MNAR, aucune méthode standard n'est garantie : le manque dépend des valeurs non observées elles-mêmes."
    )
    # Méthodes les plus et les moins éloignées de la vérité
    if (res$type == "regression") {
      tab <- res$tab; cols <- setdiff(names(tab), c("Terme", "Verite"))
      ecart <- vapply(cols, function(cn) max(abs(tab[[cn]] - tab$Verite), na.rm = TRUE), numeric(1))
      noms <- COLS_REG[cols]
      mesure <- "plus grand écart absolu, tous coefficients confondus"
    } else {
      tab <- res$tab[-1, ]
      ecart <- abs(tab$Biais); noms <- tab$Methode
      mesure <- "écart absolu à la vérité"
    }
    ok <- is.finite(ecart)
    o <- order(ecart[ok]); noms <- noms[ok][o]; ecart <- ecart[ok][o]
    d <- res$disp
    sd_vrai <- d$Ecart_type[1]; sd_moy <- d$Ecart_type[3]
    tags$div(
      tags$p(strong(paste0("Mécanisme simulé : ", mech, ". ")), intro),
      tags$p(HTML(sprintf(
        "<b>Dans ce tirage</b> (%s; %s) : la méthode la plus proche est <em>%s</em> (%s); la plus éloignée est <em>%s</em> (%s).",
        res$libelle, mesure, noms[1], fmt_fr(ecart[1], 3), noms[length(noms)], fmt_fr(ecart[length(ecart)], 3)))),
      tags$p(HTML(sprintf(
        "<b>Dispersion.</b> L'imputation par la moyenne ramène l'écart-type de %s de %s à %s; la régression déterministe le réduit aussi (%s), tandis que la régression stochastique (%s) et l'imputation multiple (%s) le préservent.",
        res$var_disp, fmt_fr(sd_vrai, 2), fmt_fr(sd_moy, 2), fmt_fr(d$Ecart_type[4], 2),
        fmt_fr(d$Ecart_type[5], 2), fmt_fr(d$Ecart_type[6], 2)))),
      tags$p(style = "color:#59636d; font-size:0.9rem;",
             "Un seul jeu simulé ne classe pas les méthodes : un petit écart peut tenir au hasard. ",
             "Pour juger biais, efficacité et couverture, lancez l'", strong("Étude de simulation"), ".")
    )
  })

  # ---- Étude de simulation (Monte Carlo) ----
  observeEvent(list(input$mc_x, input$mc_y), {
    req(input$mc_x)
    updateSelectInput(session, "mc_terme",
                      choices = c("(Intercept)", input$mc_x),
                      selected = input$mc_x[1])
  })

  mc_resultat <- eventReactive(input$mc_lancer, {
    req(rv$complet)
    cible <- input$mc_cible
    if (cible == "moyenne") {
      req(input$mc_var)
      params <- list(var = input$mc_var)
    } else if (cible == "correlation") {
      req(input$mc_v1, input$mc_v2)
      validate(need(input$mc_v1 != input$mc_v2,
                    "Choisissez deux variables différentes."))
      params <- list(v1 = input$mc_v1, v2 = input$mc_v2)
    } else {
      req(input$mc_y, input$mc_x, input$mc_terme)
      validate(need(!(input$mc_y %in% input$mc_x),
                    "La variable dépendante ne peut pas être un prédicteur."))
      params <- list(y = input$mc_y, x = input$mc_x, terme = input$mc_terme)
    }
    res <- withProgress(message = "Réplications Monte Carlo en cours…", value = 0.5, {
      P <- rv$params
      etude_monte_carlo(R = as.integer(input$mc_R), N = P$N, p = P$p,
                        rho = P$rho, prop = P$prop, mech = P$mech,
                        cible = cible, params = params,
                        m = as.integer(input$mc_m),
                        graine_base = as.integer(P$graine),
                        var_cible  = if (isTRUE(P$pilote_on)) rv$var_cible else NULL,
                        var_pilote = if (isTRUE(P$pilote_on)) rv$var_pilote else NULL,
                        direction  = P$dir,
                        patron = if (isTRUE(P$pilote_on)) rv$patron_manque else "ponctuel")
    })
    attr(res, "simulation_id") <- rv$simulation_id
    res
  })

  mc_courant <- reactive({
    res <- mc_resultat()
    if (is.null(res) || !identical(attr(res, "simulation_id"), rv$simulation_id)) return(NULL)
    res
  })

  output$mc_plot_dist <- renderPlot({ res <- mc_courant(); req(res); graphe_mc_distribution(res) })
  output$mc_plot_couv <- renderPlot({ res <- mc_courant(); req(res); graphe_mc_couverture(res) })

  output$mc_plot_caterpillar <- renderPlot({ res <- mc_courant(); req(res); graphe_caterpillar(res) })

  output$mc_table <- renderDT({
    res <- mc_courant()
    validate(need(!is.null(res), "Lancez l'étude pour le jeu courant."))
    s <- synthese_mc(res)
    names(s) <- c("Méthode", "Biais moyen", "Écart-type empirique",
                  "Erreur-type moyenne", "Couverture (%)")
    datatable(s, rownames = FALSE, options = list(dom = "t", ordering = FALSE)) |>
      formatRound(c("Biais moyen", "Écart-type empirique",
                    "Erreur-type moyenne"), 3, dec.mark = ",") |>
      formatRound("Couverture (%)", 1, dec.mark = ",")
  })

  output$mc_interpretation <- renderUI({
    res <- mc_courant()
    req(res)
    s <- synthese_mc(res)
    R <- max(res$replication)
    # Un biais est « net » s'il dépasse 2 erreurs de Monte Carlo (ET empirique / √R)
    s$mcse <- s$Ecart_type_empirique / sqrt(R)
    s$biais_net <- abs(s$Biais_moyen) > 2 * s$mcse
    s$ratio <- s$Erreur_type_moyenne / s$Ecart_type_empirique
    ligne <- function(i) {
      r <- s[i, ]
      tags$li(HTML(sprintf(
        "<b>%s</b> : biais moyen %s %s; erreur-type moyenne = %s %% de la dispersion réelle; couverture %s.",
        r$Methode, fmt_fr(r$Biais_moyen, 3),
        if (r$biais_net) "(net : plus de 2 erreurs de Monte Carlo)" else "(non distinguable de 0)",
        fmt_fr(100 * r$ratio, 0), pct_fr(r$Couverture_pct, 1))))
    }
    naive <- s[grepl("moyenne", s$Methode), ]
    mi_row <- s[grepl("multiple", s$Methode), ]
    cv <- naive$Couverture_pct
    rythme <- if (cv >= 99.5) "ne rate pratiquement jamais la vérité" else
      sprintf("rate la vérité dans %s des réplications (5 %% visés)", pct_fr(100 - cv, 1))
    tags$div(
      tags$p(strong(sprintf("%d réplications, mécanisme %s. ", R, rv$params$mech)),
             "Pour chaque méthode :"),
      tags$ul(lapply(seq_len(nrow(s)), ligne)),
      tags$p(strong("Lecture. "),
        "Un biais net déplace la distribution; une erreur-type moyenne inférieure à la dispersion réelle ",
        "rend les intervalles trop étroits. Dans les deux cas, la couverture tombe sous 95 %. ",
        sprintf("Ici, l'imputation par la moyenne %s; l'imputation multiple couvre %s.",
                rythme, pct_fr(mi_row$Couverture_pct, 1))),
      tags$p(style = "color:#59636d; font-size:0.9rem;",
        sprintf("Avec R = %d, une couverture estimée est elle-même incertaine d'environ ± %s points (2 erreurs de Monte Carlo autour de 95 %%).",
                R, fmt_fr(200 * sqrt(0.95 * 0.05 / R), 1)))
    )
  })
  # ---- Qualité de l'imputation ----
  qi_mids <- eventReactive(input$qi_lancer, {
    req(rv$ampute)
    validate(need(any(is.na(rv$ampute)),
                  "Aucune valeur manquante à imputer dans le jeu courant."))
    withProgress(message = "Imputation MICE en cours…", value = 0.5, {
      mi <- ajuster_mice(rv$ampute, m = as.integer(input$qi_m),
                         maxit = as.integer(input$qi_maxit),
                         graine = as.integer(rv$params$graine))
      attr(mi, "simulation_id") <- rv$simulation_id
      mi
    })
  })

  qi_courant <- reactive({
    mi <- qi_mids()
    if (is.null(mi) || !identical(attr(mi, "simulation_id"), rv$simulation_id)) return(NULL)
    mi
  })

  output$qi_plot_obsimp <- renderPlot({
    mi <- qi_courant(); req(mi, input$qi_var %in% vars_imputees(mi))
    graphe_obs_imp(mi, input$qi_var)
  })
  output$qi_txt_obsimp <- renderText({
    paste0("Une imputation raisonnable produit des valeurs dans le même registre ",
           "que l'observé. Un décalage systématique est attendu sous MAR (l'imputation ",
           "corrige justement le biais en exploitant les variables observées) ; ",
           "en revanche, des valeurs aberrantes ou hors échelle signalent un problème.")
  })
  output$qi_plot_conv <- renderPlot({ mi <- qi_courant(); req(mi); graphe_convergence(mi) })
  output$qi_plot_effetm <- renderPlot({
    mi <- qi_courant(); req(mi, input$qi_var %in% vars_imputees(mi))
    graphe_effet_m(effet_m(mi, input$qi_var), input$qi_var)
  })
  # ---- Données réelles ----
  rd_data <- reactive({ req(input$rd_jeu); charger_reel(input$rd_jeu) })

  observeEvent(rd_data(), {
    d <- rd_data(); vna <- names(d)[colSums(is.na(d)) > 0]
    if (length(vna) == 0) vna <- names(d)
    sel <- if (input$rd_jeu == "boys" && "tv" %in% vna) "tv" else vna[1]
    updateSelectInput(session, "rd_var", choices = vna, selected = sel)
  })
  observeEvent(list(rd_data(), input$rd_var), {
    d <- rd_data(); autres <- setdiff(names(d), input$rd_var)
    defaut <- switch(input$rd_jeu, boys = "age", airquality = "Temp", nhanes = "age")
    sel <- if (!is.null(defaut) && defaut %in% autres) defaut else autres[1]
    updateSelectInput(session, "rd_x", choices = autres, selected = sel)
  })

  output$rd_desc  <- renderUI({ req(input$rd_jeu)
    div(class = "encadre", HTML(desc_reel(input$rd_jeu))) })
  output$rd_lecon <- renderUI({ req(input$rd_jeu)
    d <- rd_data(); t <- tryCatch(little_mcar_test(d), error = function(e) NULL)
    rejet <- if (is.null(t) || is.na(t$p.value)) NA else t$p.value < 0.05
    div(class = "encadre", rd_lecon(input$rd_jeu, rejet,
                                    if (is.null(t)) NA_real_ else t$p.value)) })
  output$rd_vb_cell <- renderText(pct_fr(100 * mean(is.na(rd_data()))))
  output$rd_vb_cas  <- renderText(pct_fr(100 * mean(!complete.cases(rd_data()))))
  output$rd_table <- renderDT(
    datatable(rd_data(), options = list(pageLength = 6, scrollX = TRUE)))
  output$rd_carte   <- renderPlot(carte_manquants(rd_data()))
  output$rd_patrons <- renderPlot(graphe_patrons(rd_data()))
  output$rd_pctvar  <- renderPlot(graphe_pct_variable(rd_data()))

  output$rd_little <- renderUI({
    d <- rd_data()
    t <- tryCatch(little_mcar_test(d), error = function(e) NULL)
    if (is.null(t) || is.na(t$p.value)) {
      return(div(class = "encadre", "Test non calculable sur ce jeu."))
    }
    rejet <- t$p.value < 0.05
    verdict <- tags$span(
      style = paste0("font-weight:bold;color:", if (rejet) "#d95f0e" else "#31a354", ";"),
      if (rejet) "MCAR rejeté" else "MCAR non rejeté")
    tags$div(
      tags$p(sprintf("Statistique d\u00b2 = %s, ddl = %d (%d patrons, N = %d).",
                     fmt_fr(t$statistic, 2), t$df, t$n_patterns, nrow(d))),
      tags$p("Valeur-p = ", tags$b(p_fr(t$p.value)), " \u2192 ", verdict),
      div(class = if (rejet) "alerte" else "encadre",
        if (rejet) {
          "MCAR est rejeté : le manque n'est pas complètement aléatoire. Le test ne dit pas s'il s'agit de MAR ou de MNAR."
        } else {
          "On ne rejette pas MCAR : compatible avec un manque complètement aléatoire (attention à la puissance sur les petits échantillons)."
        })
    )
  })

  rd_mi <- eventReactive(input$rd_lancer, {
    req(rd_data())
    validate(need(any(is.na(rd_data())), "Aucune valeur manquante dans ce jeu."))
    mi <- withProgress(message = "Imputation MICE en cours…", value = 0.5, {
      ajuster_mice(rd_data(), m = as.integer(input$rd_m), maxit = 15L, graine = 2024L)
    })
    attr(mi, "jeu") <- input$rd_jeu
    mi
  })
  rd_mi_courant <- reactive({
    validate(need(input$rd_lancer > 0, "Cliquez « Imputer et diagnostiquer »."))
    mi <- rd_mi()
    validate(need(identical(attr(mi, "jeu"), input$rd_jeu),
                  "Jeu modifié : relancez « Imputer et diagnostiquer »."))
    mi
  })
  output$rd_obsimp <- renderPlot({
    mi <- rd_mi_courant(); req(input$rd_var %in% vars_imputees(mi))
    graphe_obs_imp(mi, input$rd_var)
  })
  output$rd_nuage <- renderPlot({
    mi <- rd_mi_courant(); req(input$rd_var %in% vars_imputees(mi), input$rd_x %in% names(mi$data))
    graphe_obs_imp_nuage(mi, input$rd_var, input$rd_x)
  })
  # ---- Sensibilité MNAR (delta-ajustement) ----
  sv_res <- eventReactive(input$sv_lancer, {
    req(rv$ampute, rv$complet, input$sv_var)
    validate(need(any(is.na(rv$ampute[[input$sv_var]])),
                  "Cette variable ne comporte pas de valeurs manquantes."))
    grille <- seq(-input$sv_range, input$sv_range, length.out = 25L)
    withProgress(message = "Analyse de sensibilité en cours…", value = 0.5, {
      tab <- sensibilite_delta(rv$ampute, input$sv_var, grille,
                               m = as.integer(input$sv_m),
                               graine = as.integer(rv$params$graine))
    })
    list(simulation_id = rv$simulation_id, tab = tab,
         var = input$sv_var, vrai = mean(rv$complet[[input$sv_var]]))
  })

  sensibilite_courante <- reactive({
    res <- sv_res()
    if (is.null(res) || !identical(res$simulation_id, rv$simulation_id)) return(NULL)
    res
  })

  output$sv_plot <- renderPlot({
    r <- sensibilite_courante(); req(r); graphe_sensibilite(r$tab, r$var, r$vrai)
  })

  output$sv_interpretation <- renderUI({
    r <- sensibilite_courante(); req(r); tab <- r$tab
    i0 <- which.min(abs(tab$delta))
    mar <- tab$estimation[i0]
    biais <- mar - r$vrai
    dstar <- delta_recuperateur(tab, r$vrai)
    f <- seuils_fragilite(tab)
    dec_txt <- function(d) switch(d,
      "positive" = "l'IC à 95 % est entièrement au-dessus de 0 (moyenne significativement positive)",
      "négative" = "l'IC à 95 % est entièrement sous 0 (moyenne significativement négative)",
      "non significative" = "l'IC à 95 % contient 0 (moyenne non distinguable de 0)")
    cote <- function(delta, dec) if (is.na(delta)) NULL else
      tags$li(sprintf("dès δ = %s : %s.", fmt_fr(delta, 2, signe = TRUE), dec_txt(dec)))
    tags$div(
      tags$p(style = "color:#59636d; font-size:0.9rem;",
        sprintf("Données courantes : mécanisme simulé %s, %s %% de valeurs manquantes sur %s.",
                rv$params$mech, fmt_fr(100 * mean(is.na(rv$ampute[[r$var]])), 1), r$var)),
      tags$p(strong("Point de départ (MAR, δ = 0). "),
        sprintf("L'imputation sous MAR estime la moyenne de %s à %s [%s ; %s] : %s.",
                r$var, fmt_fr(mar, 3), fmt_fr(tab$lo[i0], 3), fmt_fr(tab$hi[i0], 3), dec_txt(f$decision0))),
      tags$p(strong("Seuils de fragilité. "),
        if (is.na(f$haut) && is.na(f$bas)) {
          sprintf("Sur toute la plage explorée (± %s ET), cette conclusion ne change pas : elle est robuste à un écart de cette ampleur par rapport à MAR.",
                  fmt_fr(max(tab$delta), 1))
        } else {
          tagList("La conclusion change :", tags$ul(cote(f$bas, f$dec_bas), cote(f$haut, f$dec_haut)),
                  "Si de tels décalages sont plausibles pour ce processus de non-réponse, la conclusion est fragile.")
        }),
      tags$p(strong("Ce que seule la simulation révèle. "),
        sprintf("La vraie moyenne vaut %s : le biais sous MAR est %s. ", fmt_fr(r$vrai, 3), fmt_fr(biais, 3, signe = TRUE)),
        if (is.na(dstar)) {
          "Elle se situe hors de la plage explorée : élargissez l'amplitude de δ pour la localiser."
        } else {
          sprintf("Un décalage δ ≈ %s écart-type l'aurait retrouvée. En pratique, ce δ est inconnu : on juge seulement si un décalage de cette taille est plausible.",
                  fmt_fr(dstar, 2, signe = TRUE))
        }),
      div(class = "encadre",
        "Démarche standard (van Buuren, 2018, ch. 9) : à défaut de tester MNAR, ",
        "on rapporte comment la conclusion varie selon δ, et l'on justifie la ",
        "plage de δ jugée plausible par un argument de fond sur le processus de manque.")
    )
  })
  # ---- Quiz de classification MCAR / MAR / NMAR ----
  quiz <- reactiveValues(ordre = NULL, idx = 1L, score = 0L,
                         repondu = FALSE, donnee = NULL, faites = 0L)

  init_quiz <- function() {
    quiz$ordre <- sample(seq_along(QUIZ))
    quiz$idx <- 1L; quiz$score <- 0L; quiz$repondu <- FALSE
    quiz$donnee <- NULL; quiz$faites <- 0L
    updateRadioButtons(session, "quiz_reponse", selected = character(0))
  }
  observe({ if (is.null(quiz$ordre)) init_quiz() })
  observeEvent(input$quiz_recommencer, init_quiz())

  question_courante <- reactive({ req(quiz$ordre); QUIZ[[ quiz$ordre[quiz$idx] ]] })

  output$quiz_titre <- renderText(sprintf("Question %d / %d", quiz$idx, length(QUIZ)))
  output$quiz_scenario <- renderUI(div(class = "encadre", question_courante()$texte))

  observeEvent(input$quiz_valider, {
    req(input$quiz_reponse)
    if (isTRUE(quiz$repondu)) return(invisible(NULL))
    q <- question_courante()
    quiz$donnee <- input$quiz_reponse
    quiz$repondu <- TRUE
    quiz$faites <- quiz$faites + 1L
    if (identical(quiz$donnee, q$rep)) quiz$score <- quiz$score + 1L
  })

  output$quiz_feedback <- renderUI({
    if (!isTRUE(quiz$repondu)) return(NULL)
    q <- question_courante()
    bon <- identical(quiz$donnee, q$rep)
    div(class = if (bon) "cle" else "alerte", style = "margin-top:12px;",
      strong(if (bon) "Correct ! " else paste0("Réponse attendue : ", q$rep, ". ")),
      q$expl)
  })

  observeEvent(input$quiz_suivant, {
    if (quiz$idx < length(QUIZ)) {
      quiz$idx <- quiz$idx + 1L
      quiz$repondu <- FALSE; quiz$donnee <- NULL
      updateRadioButtons(session, "quiz_reponse", selected = character(0))
    }
  })

  output$quiz_score <- renderUI({
    fini <- quiz$idx == length(QUIZ) && isTRUE(quiz$repondu)
    base <- sprintf("Bonnes réponses : %d sur %d %s (banque de %d).",
                    quiz$score, quiz$faites,
                    if (quiz$faites > 1) "questions traitées" else "question traitée", length(QUIZ))
    if (fini) {
      tags$div(tags$p(strong(base)),
        tags$p(sprintf("Quiz terminé — score final : %d / %d. Cliquez « Recommencer » pour un nouvel ordre.",
                       quiz$score, length(QUIZ))))
    } else {
      tags$p(strong(base))
    }
  })
  # ---- Rapport HTML ----
  output$dl_rapport <- downloadHandler(
    filename = function() paste0("rapport_manquants_", Sys.Date(), ".html"),
    content = function(file) {
      res <- comparaison_courante()
      validate(need(!is.null(res),
                    "Lancez une comparaison pour le jeu courant avant de télécharger le rapport."))
      tab_html <- function(df) {
        df[] <- lapply(df, function(x) if (is.numeric(x)) fmt_fr(x, 3) else x)
        if (all(names(COLS_REG) %in% names(df)))
          names(df) <- c("Terme", COLS_REG[setdiff(names(df), "Terme")])
        paste0("<table border='1' cellpadding='6' cellspacing='0' ",
               "style='border-collapse:collapse;'>",
               "<tr>", paste0("<th>", names(df), "</th>", collapse = ""), "</tr>",
               paste0(apply(df, 1, function(r)
                 paste0("<tr>", paste0("<td>", r, "</td>", collapse = ""), "</tr>")),
                 collapse = ""),
               "</table>")
      }
      html <- paste0(
        "<!DOCTYPE html><html lang='fr'><head><meta charset='utf-8'>",
        "<title>Rapport — données manquantes</title>",
        "<style>body{font-family:sans-serif;max-width:820px;margin:40px auto;color:#333;}",
        "h1,h2{color:#2c7fb8;} .box{background:#eef6fb;padding:12px 16px;",
        "border-left:4px solid #2c7fb8;border-radius:4px;}</style></head><body>",
        "<h1>Rapport — traitement des données manquantes</h1>",
        "<p>Généré le ", Sys.Date(), " avec <em>simManquant</em>.</p>",
        "<div class='box'><b>Mécanisme :</b> ", rv$params$mech,
        " &nbsp;|&nbsp; <b>Analyse :</b> ", res$libelle,
        " &nbsp;|&nbsp; <b>Graine :</b> ", rv$params$graine,
        " &nbsp;|&nbsp; <b>N :</b> ", nrow(rv$ampute),
        " &nbsp;|&nbsp; <b>Variables :</b> ", ncol(rv$ampute),
        " &nbsp;|&nbsp; <b>Cellules manquantes :</b> ",
        pct_fr(100 * rv$pct_cell), "</div>",
        "<h2>Comparaison des méthodes</h2>", tab_html(res$tab),
        "<h2>Effet sur la dispersion</h2>", tab_html(res$disp),
        "<h2>Note de lecture</h2><p>La ligne « Vérité » n'existe qu'en simulation. ",
        "Un biais calculé sur un seul jeu peut tenir au hasard : seule une étude de simulation ",
        "(onglet Évaluer et valider) juge le biais, l'efficacité et la couverture d'une méthode. ",
        "L'imputation par la moyenne et la régression déterministe rétrécissent l'écart-type.</p>",
        "</body></html>"
      )
      writeLines(html, file, useBytes = TRUE)
    }
  )
}

shinyApp(ui, server)
