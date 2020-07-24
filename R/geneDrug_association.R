mCI <- function(predictions,
                observations,
                delta.pred = 0,
                delta.obs = 0,
                alpha = 0.05,
                outx = TRUE,
                alternative = c("two.sided", "less", "greater"))
{
  alternative <- match.arg(alternative)
  predictions[which(is.nan(predictions))] <- NA
  observations[which(is.nan(observations))] <- NA
  cc.ix <- complete.cases(predictions, observations)
  predictions <- predictions[which(cc.ix)]
  observations <- observations[which(cc.ix)]
  N <- length(which(cc.ix))
  c <- d <- u <- matrix(0, nrow = 1, ncol = N)
  c.d.seq <- NULL
  for (i in seq(from = 1, to = N - 1)) {
    for (j in seq(from = i + 1, to = N)) {
      pair <- c(i, j)
      if (abs(predictions[i] - predictions[j]) > delta.pred &
          abs(observations[i] - observations[j]) > delta.obs)
      {
        pp <- (predictions[i] < predictions[j])
        oo <- (observations[i] < observations[j])
        if (pp == oo) {
          c[pair] <- c[pair] + 1
          c.d.seq <- c(c.d.seq, TRUE)
        } else {
          d[pair] <- d[pair] + 1
          c.d.seq <- c(c.d.seq, FALSE)
        }
      } else {
        if (outx) {
          u[pair] <- u[pair] + 1
        } else{
          d[pair] <- d[pair] + 1
          c.d.seq <- c(c.d.seq, FALSE)
        }
      }
    }
  }
  C <- sum(c)
  D <- sum(d)
  if (N < 3 || (C == 0 && D == 0)) {
    return(list(
      "cindex" = NA,
      "p.value" = NA,
      "lower" = NA,
      "upper" = NA,
      "relevant.pairs.no" = 0
    ))
  }

  cindex <- C / (C + D)
  CC <- sum(c * (c - 1))
  DD <- sum(d * (d - 1))
  CD <- sum(c * d)
  varp <-
    4 * ((D ^ 2 * CC) - (2 * C * D * CD) + (C ^ 2 * DD)) / (C + D) ^ 4 * N * (N - 1) / (N - 2)

  if (varp >= 0) {
    sterr <- sqrt(varp / N)
    ci <- qnorm(p = alpha / 2, lower.tail = FALSE) * sterr
    p <- pnorm((cindex - 0.5) / sterr)
  } else {
    return(
      list(
        "cindex" = cindex,
        "p.value" = 1,
        "lower" = 0,
        "upper" = 0,
        "relevant.pairs.no" = (C + D) / 2,
        "concordant.pairs" = c.d.seq
      )
    )
  }
  return(
    list(
      "cindex" = cindex,
      "p.value" = switch(
        alternative,
        less = p,
        greater = 1 - p,
        two.sided = 2 * min(p, 1 - p)
      ),
      "lower" = max(cindex - ci, 0),
      "upper" = min(cindex + ci, 1),
      "relevant.pairs.no" = (C + D) / 2,
      "concordant.pairs" = c.d.seq
    )
  )
}


#' @import scales
getCI <- function(pred,
                  obs,
                  delta.pred = 0,
                  delta.obs = 0)
{
  ciz <- tryCatch({
    mCI(pred, obs, delta.pred = 0, delta.obs = 0)
  },
  error = function(cond) {
    list(
      "cindex" = NA,
      "p.value" = NA,
      "lower" = NA,
      "upper" = NA,
      "relevant.pairs.no" = 0
    )
  })
  ciz$dci <- NA
  if (!is.na(ciz$cindex))
  {
    ciz$dci <-
      scales::rescale(ciz$cindex, to = c(-1, 1), from = c(0, 1))
  }

  return(ciz)
}


getCor <- function(p, q, method = c("pearson", "spearman"))
{
  tryCatch({
    v <- stats::cor.test(p, q, method = method[1])
    list(cor = as.numeric(v$estimate), p.value = v$p.value)
  }, error = function(e)
  {
    list(cor = NA, p.value = NA)
  })
}

#' @import doParallel
#' @import parallel
#' @import PharmacoGx
compute_association <-
  function(x,
           y,
           fit = c("lm", "CI", "pearson", "spearman"),
           nthread = 1,
           type = NULL,
           standardize = 'SD',
           verbose = TRUE)
  {
    fit = fit[1]
    if (is(x, "matrix") == FALSE)
    {
      stop("x must be a matrix")
    }

    if (is.null(colnames(x)))
    {
      colnames(x) <- seq_len(ncol(x))
    }

    if (standardize == "SD") {
      x <- scale(x)[, ]
    }
    if (standardize == "rescale") {
      x <- as.matrix(apply(x, 2, .normalize01))
    }

    if (fit == "lm")
    {
      rr <- PharmacoGx:::rankGeneDrugSensitivity(data = x, drugpheno = y,
                                                 type = type,
        batch = rep("batch", length(y)), single.type = FALSE,
        standardize = standardize, nthread = nthread, verbose = verbose )

      rr <- data.frame(rr[[1]], stringsAsFactors = FALSE)
      rr$feature <- colnames(x)
      rr <- .reorderCol(rr, "feature", 1)
      return(rr)
    }

    ##------ for other fits --------
    if (nthread != 1) {
      availcore <- parallel::detectCores()
      if (missing(nthread) || nthread < 1 || nthread > availcore) {
        nthread <- availcore
      }
    }

    splitix <- parallel::splitIndices(nx = ncol(x), ncl = nthread)
    splitix <-
      splitix[vapply(splitix, length, FUN.VALUE = numeric(1)) > 0]

    if (fit == "CI")
    {
      mcres <- parallel::mclapply(splitix, function(i, x, y) {
        res <- apply(x[, i, drop = FALSE], 2, getCI, obs = y)
        rtx <- data.frame()
        for (gn in names(res))
        {
          rtx <- rbind(
            rtx,
            data.frame(
              feature = gn,
              ci = res[[gn]]$cindex,
              p.value = res[[gn]]$p.value,
              stringsAsFactors = FALSE
            )
          )
        }
        return(rtx)

      }, x = x, y = y, mc.cores = nthread)

      rest <- do.call(rbind, mcres)
      return(rest)
    }

    if (fit %in% c("pearson", "spearman"))
    {
      mcres <- parallel::mclapply(splitix, function(i, x, y) {
        res <- apply(x[, i, drop = FALSE], 2, getCor, q = y, method = fit)
        rtx <- data.frame()
        for (gn in names(res))
        {
          rtx <- rbind(
            rtx,
            data.frame(
              feature = gn,
              cor = res[[gn]]$cor,
              p.value = res[[gn]]$p.value,
              stringsAsFactors = FALSE
            )
          )
        }
        return(rtx)
      }, x = x, y = y, mc.cores = nthread)

      rest <- do.call(rbind, mcres)
      return(rest)
    }
  }
