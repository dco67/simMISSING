# _setup.R — configuration commune, à sourcer UNE fois en tête de chaque chapitre :
#   ```{r}
#   #| include: false
#   source("_setup.R")
#   ```
# Le code R n'est jamais affiché : seules les figures apparaissent.

suppressPackageStartupMessages({
  library(mice)
  library(ggplot2)
  library(MASS)
})

knitr::opts_chunk$set(
  echo       = FALSE,
  message    = FALSE,
  warning    = FALSE,
  fig.width  = 7,
  fig.height = 4.2,
  fig.align  = "center",
  dpi        = 150,
  out.width  = "90%",
  dev        = "ragg_png"
)

# Palette commune (mêmes teintes que l'application)
COUL_OBS <- "#2c7fb8"   # observé
COUL_NA  <- "#d95f0e"   # manquant / imputé
COUL_ACC <- "#31a354"   # vérité / accent
GRIS     <- "#7f7f7f"

theme_set(theme_minimal(base_size = 12) +
            theme(panel.grid.minor = element_blank(),
                  legend.position = "bottom",
                  strip.text = element_text(face = "bold")))

# Nombres à la française
frn <- function(x, d = 2) sub("^-", "\u2212", formatC(x, format = "f", digits = d, decimal.mark = ","))
axe_fr <- function(x) sub("^-", "\u2212", format(x, decimal.mark = ",", big.mark = " ",
                             scientific = FALSE, trim = TRUE, drop0trailing = TRUE))
sx_fr <- function(...) scale_x_continuous(labels = axe_fr, ...)
sy_fr <- function(...) scale_y_continuous(labels = axe_fr, ...)
pct_axe <- function(x) paste0(axe_fr(x), " %")

# Données multinormales équicorrélées (comme dans l'application)
sim_mvn <- function(N, p = 4, rho = 0.5, graine = 2416) {
  set.seed(graine)
  R <- matrix(rho, p, p); diag(R) <- 1
  d <- as.data.frame(MASS::mvrnorm(N, rep(0, p), R))
  names(d) <- paste0("X", seq_len(p)); d
}

# Probabilité logistique de manquer, centrée sur la médiane du prédicteur
p_manque <- function(z, force = 2.5, cible = 0.35) {
  a <- uniroot(function(a) mean(plogis(a + force * z)) - cible, c(-20, 20))$root
  plogis(a + force * z)
}

set.seed(2416)
