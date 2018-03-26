#' sim_z_cov
#'
#' Simulate survival data
#' @param nInd Number of individuals added to population each year (if \code{staggered = FALSE}, nInd * nYears individuals are added in year 1)
#' @param nYears Total number of years in simulation (can be shortened in createCH)
#' @param phi.sum Monthly summer survival probability
#' @param phi.win Monthly winter survival probability
#' @param mu.phi.aut Mean monthly autumn survival probability
#' @param lsigma.phi.aut SD of variation in autumn survival
#' @param rel.diff Relative difference between spring and fall migration survival (i.e, mu.phi.spr = mu.phi.aut * rel.diff)
#' @param cor.spr Covariance between spring and fall migration
#' @param staggered Should nInd new individuals be added to simulation each year (if FALSE, all individuals added in year 1)
#' @export

sim_z_cov <- function(nInd = 50, nYears = 12, phi.sum = 0.97, phi.win = 0.98,
                      mu.phi.aut = 0.9, r.var = 0.1, beta.aut, beta.spr, rel.diff, staggered = TRUE){
  ## TODO: generate separate z matrices for summer and winter
  if(staggered){
    ## Logit(mean) spring and fall migration survival
    mean.phi.aut <- log(mu.phi.aut) - log(1 - mu.phi.aut)
    mu.phi.spr <- mu.phi.aut * rel.diff
    mean.phi.spr <- log(mu.phi.spr) - log(1 - mu.phi.spr)

    ## Annual deviation from mean
    X.aut <- rnorm(nYears)
    X.spr <- rnorm(nYears)
    
    ## Realized fall survival
    lphi.aut <- mean.phi.aut + beta.aut * X.aut + rnorm(nYears, 0, r.var ^ 0.5)
    phi.aut <- 1 / (1 + exp(-lphi.aut))
    PHI.aut <- phi.aut ^ 2

    ## Realized spring survival
    lphi.spr <- mean.phi.spr + beta.spr * X.spr + rnorm(nYears, 0, r.var ^ 0.5)
    PHI.spr <- 1 / (1 + exp(-lphi.spr))

    ## Stationary period survival
    PHI.sum <- phi.sum ^ 4
    PHI.win <- phi.win ^ 5
    
    ## Realized annual survival
    PHI.annual <- PHI.sum * PHI.aut * PHI.win * PHI.spr

    ## Realized migration survival
    PHI.mig <- PHI.spr * PHI.aut

    ## Phi matrix
    phi.mat <- numeric(length = nYears * 4)
    year.end <- 4 * (seq(from = 1, to = nYears) - 1)
    phi.mat[seq(1, nYears * 4, by = 4)] <- PHI.sum
    phi.mat[seq(2, nYears * 4, by = 4)] <- PHI.aut
    phi.mat[seq(3, nYears * 4, by = 4)] <- PHI.win
    phi.mat[seq(4, nYears * 4, by = 4)] <- PHI.spr
    

    ## Simulate survival data
    z.sum <- matrix(0, nrow = nInd * nYears, ncol = (nYears * 4 + 1))
    z.win <- matrix(0, nrow = nInd * nYears, ncol = ((nYears - 1) * 4 + 2))

    f.sum <- 1 + 4 * rep(seq(from = 0, to = nYears - 1), each = nInd)
    f.win <- 1 + 4 * rep(seq(from = 0, to = nYears - 1), each = nInd)

    for(i in 1:(nInd * nYears)){
      z.sum[i, f.sum[i]] <- 1
      z.win[i, f.win[i]] <- 1
      for(j in (f.sum[i] + 1):ncol(z.sum)){
        z.sum[i, j] <- rbinom(n = 1, 1, prob = phi.mat[j - 1] * z.sum[i, j - 1])
      }

      for(j in (f.win[i] + 1):ncol(z.win)){
        z.win[i, j] <- rbinom(n = 1, 1, prob = phi.mat[j + 1] * z.win[i, j - 1])
      }
    }

    ll <- list(z.sum = z.sum, z.win = z.win, PHI.annual = PHI.annual, PHI.mig = PHI.mig,
               PHI.spr = PHI.spr, PHI.aut = PHI.aut, PHI.sum = PHI.sum, PHI.win = PHI.win,
               mean.phi.sum = phi.sum, mean.phi.win = phi.win,
               nYears = nYears, nInd.sum = nInd,  nInd.win = nInd, X.aut = X.aut, X.spr = X.spr)
  }else{
    mean.phi.aut <- log(mu.phi.aut) - log(1 - mu.phi.aut)
    mu.phi.spr <- mu.phi.aut * rel.diff
    mean.phi.spr <- log(mu.phi.spr) - log(1 - mu.phi.spr)
    
    ## Annual deviation from mean
    X.aut <- rnorm(nYears)
    X.spr <- rnorm(nYears)
    
    ## Realized fall survival
    lphi.aut <- mean.phi.aut + beta.aut * X.aut + rnorm(nYears, 0, r.var ^ 0.5)
    phi.aut <- 1 / (1 + exp(-lphi.aut))
    PHI.aut <- phi.aut ^ 2
    
    ## Realized spring survival
    lphi.spr <- mean.phi.spr + beta.spr * X.spr + rnorm(nYears, 0, r.var ^ 0.5)
    PHI.spr <- 1 / (1 + exp(-lphi.spr))

    ## Stationary period survival
    PHI.sum <- phi.sum ^ 4
    PHI.win <- phi.win ^ 5
    
    ## Realized annual survival
    PHI.annual <- PHI.sum * PHI.aut * PHI.win * PHI.spr
    
    ## Realized migration survival
    PHI.mig <- PHI.spr * PHI.aut

    ## Phi matrix
    phi.mat <- numeric(length = nYears * 12)
    year.end <- 12 * (seq(from = 1, to = nYears) - 1)
    phi.mat[sort(unlist(lapply(seq(from = 1, to = 4), function(x) x + year.end)))] <- phi.sum
    phi.mat[sort(unlist(lapply(seq(from = 7, to = 11), function(x) x + year.end)))] <- phi.win

    for(i in 1:nYears){
      phi.mat[5:6 + 12 * (i - 1)] <- phi.aut[i]
      phi.mat[12 + 12 * (i - 1)] <- phi.spr[i]
    }

    ## Simulate survival data
    z.sum <- matrix(0, nrow = nInd * nYears, ncol = (nYears * 12) - 8)
    z.win <- matrix(0, nrow = nInd * nYears, ncol = (nYears * 12) - 6)

    z.sum[, 1] <- 1
    z.win[, 1] <- 1

    for(i in 1:(nInd * nYears)){
      for(j in 2:ncol(z.sum)){
        z.sum[i, j] <- rbinom(n = 1, 1, prob = phi.mat[j] * z.sum[i, j - 1])
      }

      for(j in 2:ncol(z.win)){
        z.win[i, j] <- rbinom(n = 1, 1, prob = phi.mat[j + 5] * z.win[i, j - 1])
      }
    }

    ll <- list(z.sum = z.sum, z.win = z.win, phi.annual = phi.annual, phi.mig = phi.mig,
               phi.spr = phi.spr, phi.aut = phi.aut, phi.sum = phi.sum, phi.win = phi.win,
               nYears = nYears, nInd = nInd,  X.aut = X.aut, X.spr = X.spr)
  }

  ll
}

