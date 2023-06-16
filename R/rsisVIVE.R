#' @title Robustified Some Valid Some Invalid Instrumental Variable Estimator
#'
#' @description
#' This function provides the Robustified Some Valid Some Invalid Instrumental Variable Estimator (rsisVIVE) coefficients for step one (alpha) and step two (beta) of the instrumental variable estimation. The algorithm closely follows Kang et. al.'s some invalid some valid instrumental variable estimator (2016) however, the rsisVIVE can tolerate large proportions of contamination both in the exposure and/or outcome variables in addition to tolerating some invalid instrumental variables.
#'
#' @param Y matrix array nx1 specifying the outcome values
#' @param D matrix array nx1 specifying the exposure values
#' @param Z matrix array nxL specifying the instruments where L is the number of instruments
#' @param method character string specifying the type of method to obtain the step one (alpha) coefficients. If `BIC` is chosen, it uses the robust BIC as specified in the Sparse Least Trimmed Squares method (Alfons et. al., 2013). If `PE_MIN` or `PE_SE` is chosen, the function uses a 10-fold cross-validation minimizing the Root Trimmed Mean Squared Prediction Error and selecting the minimum or the one standard error, correspondingly. Note that the selection procedure follows the majority rule as outlined in the
#' @param ncores numeric indicating the number of cores used for the cross-validation
#' @return description
#' @examples
#' n <- 50
#' L <- 10
#' s <- 3
#' contam <- 0.2
#' Sigma <- diag(1, nrow = L)
#'
#' a0 <- 1
#' a1 <- rep(0.05, L)
#'
#' b0 <- 1
#' b1 <- c(rep(1,s), rep(0, L-s)) # alpha
#' b2 <- 1                        # beta
#'
#' error <- MASS::mvrnorm(n, c(0, 0), matrix(c(1, 0.9, 0.9, 1), 2, 2))
#'
#' Z <- MASS::mvrnorm(n, mu = rep(0, L), Sigma)
#' D <- a0 + Z %*% a1 + error[, 1]
#' Y <- b0 + Z %*% b1 + D * b2 + error[, 2]
#' Y[sample(1:n, m)] <- rnorm(m, 70, 1) # contamination
#'
#' rsisVIVE(Y, D, Z, method = 'PE_SE', ncores = 1)
#' @import dplyr
#' @import MASS
#' @import robustbase
#' @import robustHD
#' @import perry
#' @export
rsisVIVE <- function(
        Y,
        D,
        Z,
        method = c("BIC", "PE_MIN", "PE_SE"),
        ncores = 1) {
    # robust standardize
    Y <- robustHD::robStandardize(Y, centerFun = median, scaleFun = mad)
    D <- robustHD::robStandardize(D, centerFun = median, scaleFun = mad)
    Z <- robustHD::robStandardize(Z, centerFun = median, scaleFun = mad)

    # projections
    QR <- qr(Z)
    Yhat <- qr.fitted(QR, Y)

    # use robust regression to get robust Dhat values
    # patchwork for the rare case where robust regression doesn't converge
    tryCatch(
        {
            # Code block where an error might occur
            Dhat <- predict(MASS::rlm(D ~ Z, method = "MM"), data.frame(Z))
        },
        error = function(e) {
            # Error handling code
            Dhat <- qr.fitted(QR, D)
        }
    )

    # transform Y and Z according to paper
    Z_transformed <- Z - Dhat %*% t(Dhat) %*% Z / sum(Dhat^2)
    Y_transformed <- Yhat - Dhat * (sum(Dhat * Yhat) / sum(Dhat^2))

    # STEP 1: SPARSE LTS
    # LAMBDA SEQUENCE
    lambda_0_hat <- robustHD::lambda0(Z_transformed, Y_transformed) * 2
    lambda_seq <- seq(0, lambda_0_hat, by = lambda_0_hat * 0.025)

    if (method == "BIC") {
        fit <- robustHD::sparseLTS(
            x = Z_transformed,
            y = as.vector(Y_transformed),
            intercept = FALSE,
            normalize = TRUE,
            crit = "BIC",
            lambda = lambda_seq,
            mode = "lambda",
            ncores = ncores
        )
        coefs <- fit$coefficients
        idx <- which(colSums(coefs != 0) <= floor(ncol(Z_transformed) / 2))
        idx_min <- which.min(fit$crit$values[, "reweighted"][idx])
        alpha <- fit$coefficients[, idx_min]
    } else if (method %in% c("PE_MIN", "PE_SE")) {
        fits <- list()
        for (i in 1:length(lambda_seq)) {
            lambda <- lambda_seq[i]
            fit <- robustHD::sparseLTS(
                x = Z_transformed,
                y = as.vector(Y_transformed),
                intercept = FALSE,
                normalize = TRUE,
                lambda = lambda,
                crit = "PE",
                ncores = ncores
            )
            fits[[i]] <- fit
        }

        # data frame of lambda, scale, and number of nonzero coefficients
        fits_properties <- lapply(fits, function(fit) {
            lambda <- fit$lambda
            scale <- fit$scale
            num_nonzero <- sum(fit$coefficients != 0)
            data.frame(lambda, scale, num_nonzero)
        }) %>% bind_rows()

        # get new restricted lambda sequence
        tryCatch(
            {
                new_lambdas <- fits_properties %>%
                    filter(num_nonzero <= floor(ncol(Z_transformed) / 2)) %>%
                    summarise(
                        max_lambda = max(lambda),
                        min_lambda = min(lambda),
                    )
            },
            error = function(e) {
                new_lambdas <- data.frame(
                    min_lambda = 1e-5,
                    max_lambda = lambda_0_hat / 4)
            }
        )
        lambda_seq_pe <- seq(
            new_lambdas$min_lambda,
            new_lambdas$max_lambda,
            length.out = 30
        )
        fit <- robustHD::sparseLTS(
            x = Z_transformed,
            y = as.vector(Y_transformed),
            intercept = FALSE,
            normalize = TRUE,
            lambda = lambda_seq_pe,
            mode = "lambda",
            splits = foldControl(K = 10, R = 1, type = "random"),
            crit = "PE",
            ncores = ncores
        )
        if (method == "PE_MIN") {
            alpha <- coef(fit)
        } else if (method == "PE_SE") {
            pe_values <- round((fit$pe$reweighted), 3)
            lambda_values <- unlist((fit$tuning))
            min_error <- min(pe_values)
            se <- sd(pe_values) / sqrt(length(pe_values))
            idx <- which(pe_values < (min_error + se))[1]
            lambda_pe_se <- lambda_values[idx]
            fit_pe_se <- robustHD::sparseLTS(
                x = Z_transformed,
                y = as.vector(Y_transformed),
                intercept = FALSE,
                normalize = TRUE,
                lambda = lambda_pe_se,
                crit = "PE",
                ncores = ncores
            )
            alpha <- coef(fit_pe_se)
        }
    } else {
        stop("method must be one of 'BIC', 'PE_MIN', 'PE_SE'")
    }

    # STEP 2: LTS
    Y_new <- (Y - (Z %*% alpha))
    fit_2 <- robustbase::ltsReg(Y_new ~ Dhat)
    beta <- as.numeric(coef(fit_2)[2])

    # packaging Object
    out <- list(
        alpha = alpha,
        beta = beta
    )
    return(out)
}
