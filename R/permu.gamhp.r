#' Permutation Test of Hierarchical Partitioning for GAM Analysis
#' @param  iv  optional The relative importance of predictor groups will be assessed. The input for iv should be a list, where each element contains the names of variables belonging to a specific group. These variable names must correspond to the predictor variables defined in the model (mod).
#' @param  mod gam model generated by mgcv::gam()
#' @param  type The type of total explained variation, either "dev" or "adjR2", in which "dev" is deviance explained and "adjR2" is adjusted R-square, the default is "adjR2".

#' @param  permutations An integer; Number of permutations for computing p value of individual contribution for the randomized dataset.

#' @details This function is a permutation test of hierarchical partitioning for gam analysis. It returns a matrix of I values (the individual contribution towards total explained variation) for all values from permutations randomizations. For each permutation, the values in each variable (i.e each column of iv) are randomized independently, and gam.hp is run on the randomized iv. As well as the randomized I matrix, the function returns a summary table listing the observed I values, the p value of I for the randomized dataset.

#' @return a data.frame containing a summary table listing the observed individual contribution, the p value of individual contribution for the randomized dataset

#' @author {Jiangshan Lai} \email{lai@njfu.edu.cn}


#'@export
#'@examples
#'library(mgcv)
#'mod1 <- gam(Sepal.Length ~ s(Petal.Length) + s(Petal.Width) + Sepal.Width,data = iris)
#'permu.gamhp(mod=mod1,type="dev",permutations=10)
#'iv <- list(env1=c("s(Petal.Length)","s(Petal.Width)"),env2="Sepal.Width")
#'permu.gamhp(mod=mod1,iv,type="dev",permutations=10)


permu.gamhp <- function (mod = NULL,iv = NULL, type = "dev", permutations = 10)
{
    y_name = as.character(mod$formula)[2]
    x_name = names(mod$var.summary)
    d_y = as.data.frame(mod$y)
    colnames(d_y) <- y_name
    Iv <- as.data.frame(mod$model[, x_name])
    obs <- gam.hp(mod = mod, iv=iv, type = type, commonality = FALSE)
    r2q_p <- obs$hierarchical.partitioning[, 3]
    n <- dim(Iv)[1]
    nvar <- dim(Iv)[2]
    for (i in seq(1, permutations - 1, 1)) {
        newiv <- Iv
        for (j in 1:nvar) {
            perms <- sample(1:n, n)
            newiv[, j] <- Iv[, j][perms]
        }

        dat1 <- as.data.frame(cbind(d_y, newiv))
        mod2 <- stats::update(object = mod, formula. = . ~ ., data = dat1)
        simu_p = gam.hp(mod = mod2,iv=iv, type = type, commonality = FALSE)
        r2q_p = cbind(r2q_p, simu_p$hierarchical.partitioning[,3])
    }
    Signi <- function(x) {
        pval = round(1 - ecdf(x)(x[1]) + 1/(permutations + 1),
            nchar(permutations))
        if (pval <= 0.001) {
            return(noquote(paste(pval, "***", sep = " ")))
        }
        if (pval <= 0.01) {
            return(noquote(paste(pval, " **", sep = " ")))
        }
        if (pval <= 0.05) {
            return(noquote(paste(pval, "  *", sep = " ")))
        }
        else {
            return(noquote(paste(pval, "   ", sep = " ")))
        }
    }
    p.R2_p = apply(r2q_p, 1, Signi)
    result = data.frame(obs$hierarchical.partitioning[, 3], Pr = p.R2_p)
    colnames(result) <- c("Individual", "Pr(>I)")
    return(result)
}
