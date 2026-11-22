# _moteur.R — fonctions de calcul reprises telles quelles de l'application
# simManquant (simulation, amputation, test de Little, imputations simples, FIML),
# pour que les chiffres du guide soient identiques à ceux de l'écran.
# Ce fichier est sourcé par les chapitres qui en ont besoin; il n'est jamais affiché.

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

