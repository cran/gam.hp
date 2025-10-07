#' Hierarchical Partitioning of Adjusted R2 and Explained Deviance for Generalized Additive Models (GAM and BAM)
#'
#' @description
#' This function conducts hierarchical partitioning to calculate the individual contributions of each predictor
#' towards total adjusted R^2 and explained deviance for Generalized Additive Models fitted by either
#' \code{\link[mgcv]{gam}} or \code{\link[mgcv]{bam}} in the \pkg{mgcv} package.
#'
#' @param mod Fitted \code{"gam"} or \code{"bam"} model object from the \pkg{mgcv} package.
#' @param iv Optional. A list specifying groups of predictor variables for assessing group-wise relative importance.
#' Each element of the list should contain the names of variables belonging to a specific group, corresponding
#' to the predictor names defined in the model (\code{mod}).
#' @param type Character. The type of R-square for GAM/BAM models, either \code{"dev"} or \code{"adjR2"}.
#' \code{"dev"} represents the explained deviance, and \code{"adjR2"} represents the adjusted R-square.
#' The default is \code{"dev"}.
#' @param commonality Logical; if \code{TRUE}, the function returns the results of commonality analysis
#' (i.e., 2^N - 1 fractions for N predictors). Default is \code{FALSE}.
#' @param data Optional. The dataset used to fit the model. If not provided, the function will attempt to extract
#' the data directly from the fitted \code{mod} object. This argument is mainly useful for \code{\link[mgcv]{bam}}
#' models, where the model object may not always store the full dataset (especially when using large datasets or
#' parallel fitting options).
#'
#' @details
#' The function supports both \code{\link[mgcv]{gam}} and \code{\link[mgcv]{bam}} model objects.
#' It decomposes the total explained deviance or adjusted R^2 into unique and shared contributions of
#' individual predictors or groups of predictors using hierarchical partitioning.
#' The adjusted R^2 and explained deviance values are extracted from \code{summary.gam()} or \code{summary.bam()}.
#'
#' @return
#' \item{dev}{The R2 (either explained deviance or adjusted R^2) for the full model.}
#' \item{hierarchical.partitioning}{A matrix containing the individual effects and the percentage
#' contribution of each predictor to total explained deviance or adjusted R^2.}
#'
#' @author
#' {Jiangshan Lai} \email{lai@njfu.edu.cn}
#'
#' @references
#' \itemize{
#' \item Lai J., Tang J., Li T., Zhang A., Mao L. (2024). Evaluating the relative importance of predictors in Generalized Additive Models using the gam.hp R package. *Plant Diversity*, 46(4): 542-546. <DOI:10.1016/j.pld.2024.06.002>
#' \item Lai J., Zhu W., Cui D., Mao L. (2023). Extension of the glmm.hp package to zero-inflated generalized linear mixed models and multiple regression. *Journal of Plant Ecology*, 16(6): rtad038. <DOI:10.1093/jpe/rtad038>
#' \item Lai J., Zou Y., Zhang S., Zhang X., Mao L. (2022). glmm.hp: an R package for computing individual effect of predictors in generalized linear mixed models. *Journal of Plant Ecology*, 15(6): 1302-1307. <DOI:10.1093/jpe/rtac096>
#' \item Lai J., Zou Y., Zhang J., Peres-Neto P. (2022). Generalizing hierarchical and variation partitioning in multiple regression and canonical analyses using the rdacca.hp R package. *Methods in Ecology and Evolution*, 13(4): 782-788. <DOI:10.1111/2041-210X.13800>
#' \item Chevan, A. & Sutherland, M. (1991). Hierarchical partitioning. *The American Statistician*, 45, 90-96. <DOI:10.1080/00031305.1991.10475776>
#' \item Nimon, K., Oswald, F. L. & Roberts, J. K. (2013). Yhat: Interpreting regression effects. R package version 2.0.0.
#' }
#'
#' @export
#' @examples
#' library(mgcv)
#'
#' ## Example 1: Using gam()
#' mod1 <- gam(Sepal.Length ~ s(Petal.Length) + s(Petal.Width) + Sepal.Width,
#'             data = iris)
#' summary(mod1)
#' gam.hp(mod1)
#' gam.hp(mod1, type = "adjR2")
#' gam.hp(mod1, commonality = TRUE)
#' iv <- list(env1 = c("s(Petal.Length)", "s(Petal.Width)"), env2 = "Sepal.Width")
#' gam.hp(mod1, iv, type = "adjR2")
#' gam.hp(mod1, iv, commonality = TRUE)
#'
#' ## Example 2: Using bam()
#' mod2 <- bam(Sepal.Length ~ s(Petal.Length) + s(Petal.Width) + Sepal.Width,
#'             data = iris)
#' summary(mod2)
#' gam.hp(mod2)
#' gam.hp(mod2, type = "adjR2")
#' gam.hp(mod2, commonality = TRUE)
#' ## Explicitly specifying data (useful for bam)
#' gam.hp(mod2, data = iris)
#'
gam.hp <- function(mod, iv = NULL, type = "dev", commonality = FALSE, data = NULL)
{
  ## Support both gam and bam objects
  if (!(inherits(mod, "gam") || inherits(mod, "bam")))
    stop("gam.hp currently supports 'gam' and 'bam' objects from the mgcv package.")

  ## Check RHS for '*' (encourage using ':' for interactions)
  rhs_char <- tryCatch(as.character(formula(mod))[3], error = function(e) "")
  if (nzchar(rhs_char) && grepl("\\*", rhs_char, fixed = FALSE))
    stop("Please put the interaction term as a new variable (use colon ':' rather than asterisk '*') in the original model formula.")

  ## Extract term labels robustly
  ivname <- attr(terms(mod), "term.labels")
  if (is.null(ivname)) ivname <- character(0)

  ## Check for offset term (if any)
  offsetterm <- NULL
  offpos <- attr(terms(mod), "offset")
  if (!is.null(offpos) && length(offpos) > 0 && offpos > 0) {
    offsetterm <- as.character(attr(mod$terms, "variables")[[offpos + 1]])
  }

  iv.name <- ivname

  ## Extract R^2 or explained deviance from summary
  smod <- summary(mod)
  if (type == "adjR2") outr2 <- smod$r.sq else if (type == "dev") outr2 <- smod$dev.expl else stop("type must be 'dev' or 'adjR2'")

  r2type <- row.names(outr2)
  nr2type <- length(r2type)
  if (nr2type == 0) {
    nr2type <- 1
    if (commonality) r2type <- "commonality.analysis" else r2type <- "hierarchical.partitioning"
  }

  ## Try to retrieve the original data:
  dat <- NULL
  if (!is.null(data)) {
    dat <- data
  } else {
    if (!is.null(mod$model) && is.data.frame(mod$model)) {
      dat <- mod$model
    } else {
      if (!is.null(mod$call$data)) {
        dat_try <- tryCatch(eval(mod$call$data, envir = environment(formula(mod))), error = function(e) NULL)
        if (is.null(dat_try)) dat_try <- tryCatch(eval(mod$call$data, envir = parent.frame()), error = function(e) NULL)
        if (!is.null(dat_try) && is.data.frame(dat_try)) dat <- dat_try
      }
      if (is.null(dat)) {
        dat_try2 <- tryCatch(model.frame(formula(mod), data = environment(formula(mod))), error = function(e) NULL)
        if (is.null(dat_try2)) dat_try2 <- tryCatch(model.frame(mod), error = function(e) NULL)
        if (!is.null(dat_try2) && is.data.frame(dat_try2)) dat <- dat_try2
      }
    }
  }

  if (is.null(dat) || !is.data.frame(dat)) {
    stop("Cannot locate the data.frame used to fit the model. Refit the model with a named 'data' argument or call gam.hp(mod, data = your_data).")
  }
  dat <- na.omit(dat)

  ## Build null model (remove all predictor variables)
  to_del <- paste(paste("-", iv.name, sep = ""), collapse = " ")
  modnull <- stats::update(stats::formula(mod), paste(". ~ . ", to_del, sep = ""))
  mod_null <- stats::update(object = mod, formula. = modnull, data = dat)

  ## Main hierarchical partitioning logic (supports iv == NULL or grouped iv)
  if (is.null(iv)) {

    nvar <- length(iv.name)
    if (nvar < 2)
      stop("Analysis not conducted. Insufficient number of predictors.")

    totalN <- 2^nvar - 1
    binarymx <- matrix(0, nvar, totalN)
    for (i in 1:totalN) {
      binarymx <- creatbin(i, binarymx)
    }

    outputList <- list()
    outputList[[1]] <- outr2
    for (k in 1:nr2type) {
      commonM <- matrix(nrow = totalN, ncol = 3)
      for (i in 1:totalN) {
        tmp.name <- iv.name[as.logical(binarymx[, i])]

        to_add_terms <- tmp.name
        if (!is.null(offsetterm)) to_add_terms <- c(to_add_terms, offsetterm)

        to_add <- paste("~", paste(to_add_terms, collapse = " + "))
        modnew <- stats::update(object = mod_null, data = dat, formula = as.formula(to_add))

        if (type == "dev") commonM[i, 2] <- summary(modnew)$dev.expl
        if (type == "adjR2") commonM[i, 2] <- summary(modnew)$r.sq
      }

      commonlist <- vector("list", totalN)

      seqID <- vector()
      for (i in 1:nvar) {
        seqID[i] <- 2^(i - 1)
      }

      for (i in 1:totalN) {
        bit <- binarymx[1, i]
        if (bit == 1)
          ivname_work <- c(0, -seqID[1])
        else ivname_work <- seqID[1]
        for (j in 2:nvar) {
          bit <- binarymx[j, i]
          if (bit == 1) {
            alist <- ivname_work
            blist <- genList(ivname_work, -seqID[j])
            ivname_work <- c(alist, blist)
          }
          else ivname_work <- genList(ivname_work, seqID[j])
        }
        ivname_work <- ivname_work * -1
        commonlist[[i]] <- ivname_work
      }

      for (i in 1:totalN) {
        r2list <- unlist(commonlist[i])
        numlist <- length(r2list)
        ccsum <- 0
        for (j in 1:numlist) {
          indexs <- r2list[[j]]
          indexu <- abs(indexs)
          if (indexu != 0) {
            ccvalue <- commonM[indexu, 2]
            if (indexs < 0)
              ccvalue <- ccvalue * -1
            ccsum <- ccsum + ccvalue
          }
        }
        commonM[i, 3] <- ccsum
      }

      orderList <- vector("list", totalN)
      index <- 0
      for (i in 1:nvar) {
        for (j in 1:totalN) {
          nbits <- sum(binarymx[, j])
          if (nbits == i) {
            index <- index + 1
            commonM[index, 1] <- j
          }
        }
      }

      outputcommonM <- matrix(nrow = totalN + 1, ncol = 2)
      totalRSquare <- sum(commonM[, 3])
      for (i in 1:totalN) {
        outputcommonM[i, 1] <- round(commonM[commonM[i, 1], 3], digits = 4)
        outputcommonM[i, 2] <- round((commonM[commonM[i, 1], 3] / totalRSquare) * 100, digits = 2)
      }
      outputcommonM[totalN + 1, 1] <- round(totalRSquare, digits = 4)
      outputcommonM[totalN + 1, 2] <- round(100, digits = 4)
      rowNames <- NULL
      for (i in 1:totalN) {
        ii <- commonM[i, 1]
        nbits <- sum(binarymx[, ii])
        cbits <- 0
        if (nbits == 1)
          rowName <- "Unique to "
        else rowName <- "Common to "
        for (j in 1:nvar) {
          if (binarymx[j, ii] == 1) {
            if (nbits == 1)
              rowName <- paste(rowName, iv.name[j], sep = "")
            else {
              cbits <- cbits + 1
              if (cbits == nbits) {
                rowName <- paste(rowName, "and ", sep = "")
                rowName <- paste(rowName, iv.name[j], sep = "")
              } else {
                rowName <- paste(rowName, iv.name[j], sep = "")
                rowName <- paste(rowName, ", ", sep = "")
              }
            }
          }
        }
        rowNames <- c(rowNames, rowName)
      }
      rowNames <- c(rowNames, "Total")
      rowNames <- format.default(rowNames, justify = "left")
      colNames <- format.default(c("Fractions", " % Total"), justify = "right")
      dimnames(outputcommonM) <- list(rowNames, colNames)

      VariableImportance <- matrix(nrow = nvar, ncol = 4)
      for (i in 1:nvar) {
        VariableImportance[i, 3] <- round(sum(binarymx[i, ] * (commonM[, 3] / apply(binarymx, 2, sum))), digits = 4)
      }

      VariableImportance[, 1] <- outputcommonM[1:nvar, 1]
      VariableImportance[, 2] <- VariableImportance[, 3] - VariableImportance[, 1]

      total <- round(sum(VariableImportance[, 3]), digits = 4)
      VariableImportance[, 4] <- round(100 * VariableImportance[, 3] / total, 2)
      dimnames(VariableImportance) <- list(iv.name, c("Unique", "Average.share", "Individual", "I.perc(%)"))

      if (commonality) {
        outputList[[k + 1]] <- outputcommonM
      } else {
        outputList[[k + 1]] <- VariableImportance
      }
    }
  } else {
    ## Grouped iv branch
    nvar <- length(iv)
    ilist <- names(iv)
    if (is.null(ilist)) {
      names(iv) <- paste("X", 1:nvar, sep = "")
    } else {
      whichnoname <- which(ilist == "")
      if (length(whichnoname) > 0) names(iv)[whichnoname] <- paste("X", whichnoname, sep = "")
    }

    ilist <- names(iv)

    ivlist <- ilist
    iv.name <- ilist

    ivID <- matrix(nrow = nvar, ncol = 1)
    for (i in 0:nvar - 1) {
      ivID[i + 1] <- 2^i
    }

    totalN <- 2^nvar - 1

    binarymx <- matrix(0, nvar, totalN)
    for (i in 1:totalN) {
      binarymx <- creatbin(i, binarymx)
    }

    outputList <- list()
    outputList[[1]] <- outr2
    for (k in 1:nr2type) {
      commonM <- matrix(nrow = totalN, ncol = 3)
      for (i in 1:totalN) {
        tmpname <- iv.name[as.logical(binarymx[, i])]
        tmp.name <- unlist(iv[names(iv) %in% tmpname])
        to_add_terms <- tmp.name
        if (!is.null(offsetterm)) to_add_terms <- c(to_add_terms, offsetterm)
        to_add <- paste("~", paste(to_add_terms, collapse = " + "))
        modnew <- stats::update(object = mod_null, data = dat, formula = as.formula(to_add))

        if (type == "dev") commonM[i, 2] <- summary(modnew)$dev.expl
        if (type == "adjR2") commonM[i, 2] <- summary(modnew)$r.sq
      }

      commonlist <- vector("list", totalN)

      seqID <- vector()
      for (i in 1:nvar) {
        seqID[i] <- 2^(i - 1)
      }

      for (i in 1:totalN) {
        bit <- binarymx[1, i]
        if (bit == 1)
          ivname_work <- c(0, -seqID[1])
        else ivname_work <- seqID[1]
        for (j in 2:nvar) {
          bit <- binarymx[j, i]
          if (bit == 1) {
            alist <- ivname_work
            blist <- genList(ivname_work, -seqID[j])
            ivname_work <- c(alist, blist)
          } else ivname_work <- genList(ivname_work, seqID[j])
        }
        ivname_work <- ivname_work * -1
        commonlist[[i]] <- ivname_work
      }

      for (i in 1:totalN) {
        r2list <- unlist(commonlist[i])
        numlist <- length(r2list)
        ccsum <- 0
        for (j in 1:numlist) {
          indexs <- r2list[[j]]
          indexu <- abs(indexs)
          if (indexu != 0) {
            ccvalue <- commonM[indexu, 2]
            if (indexs < 0)
              ccvalue <- ccvalue * -1
            ccsum <- ccsum + ccvalue
          }
        }
        commonM[i, 3] <- ccsum
      }

      orderList <- vector("list", totalN)
      index <- 0
      for (i in 1:nvar) {
        for (j in 1:totalN) {
          nbits <- sum(binarymx[, j])
          if (nbits == i) {
            index <- index + 1
            commonM[index, 1] <- j
          }
        }
      }

      outputcommonM <- matrix(nrow = totalN + 1, ncol = 2)
      totalRSquare <- sum(commonM[, 3])
      for (i in 1:totalN) {
        outputcommonM[i, 1] <- round(commonM[commonM[i, 1], 3], digits = 4)
        outputcommonM[i, 2] <- round((commonM[commonM[i, 1], 3] / totalRSquare) * 100, digits = 2)
      }
      outputcommonM[totalN + 1, 1] <- round(totalRSquare, digits = 4)
      outputcommonM[totalN + 1, 2] <- round(100, digits = 4)
      rowNames <- NULL
      for (i in 1:totalN) {
        ii <- commonM[i, 1]
        nbits <- sum(binarymx[, ii])
        cbits <- 0
        if (nbits == 1)
          rowName <- "Unique to "
        else rowName <- "Common to "
        for (j in 1:nvar) {
          if (binarymx[j, ii] == 1) {
            if (nbits == 1)
              rowName <- paste(rowName, iv.name[j], sep = "")
            else {
              cbits <- cbits + 1
              if (cbits == nbits) {
                rowName <- paste(rowName, "and ", sep = "")
                rowName <- paste(rowName, iv.name[j], sep = "")
              } else {
                rowName <- paste(rowName, iv.name[j], sep = "")
                rowName <- paste(rowName, ", ", sep = "")
              }
            }
          }
        }
        rowNames <- c(rowNames, rowName)
      }
      rowNames <- c(rowNames, "Total")
      rowNames <- format.default(rowNames, justify = "left")
      colNames <- format.default(c("Fractions", " % Total"), justify = "right")
      dimnames(outputcommonM) <- list(rowNames, colNames)

      VariableImportance <- matrix(nrow = nvar, ncol = 4)
      for (i in 1:nvar) {
        VariableImportance[i, 3] <- round(sum(binarymx[i, ] * (commonM[, 3] / apply(binarymx, 2, sum))), digits = 4)
      }

      VariableImportance[, 1] <- outputcommonM[1:nvar, 1]
      VariableImportance[, 2] <- VariableImportance[, 3] - VariableImportance[, 1]

      total <- round(sum(VariableImportance[, 3]), digits = 4)
      VariableImportance[, 4] <- round(100 * VariableImportance[, 3] / total, 2)
      dimnames(VariableImportance) <- list(iv.name, c("Unique", "Average.share", "Individual", "I.perc(%)"))

      if (commonality) {
        outputList[[k + 1]] <- outputcommonM
      } else {
        outputList[[k + 1]] <- VariableImportance
      }
    }
  }

  if (type == "adjR2") {
    names(outputList) <- c("adjusted.R2", r2type)
  }
  if (type == "dev") {
    names(outputList) <- c("Explained.deviance", r2type)
  }
  outputList$variables <- iv.name
  if (commonality) {
    outputList$type <- "commonality.analysis"
  } else {
    outputList$type <- "hierarchical.partitioning"
  }
  class(outputList) <- "gamhp"
  outputList
}
