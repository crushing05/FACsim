#' sim_z
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

sim_z <- function(nInd = 150, omega = 65, nYears = 12, 
                  phi.sum = 0.97, phi.win = 0.98,
                  mu.phi.aut = 0.9, lsigma.phi.aut, 
                  rel.diff, cor.spr){
  
    nImm <- rpois(nYears - 1, lambda = omega)
    
    ## Logit(mean) spring and fall migration survival
    mean.phi.aut <- log(mu.phi.aut) - log(1 - mu.phi.aut)
    mu.phi.spr <- mu.phi.aut * rel.diff
    mean.phi.spr <- log(mu.phi.spr) - log(1 - mu.phi.spr)

    ## Annual deviation from mean
    lsigma.phi.spr <- (1 - cor.spr ^ 2) * ((mu.phi.aut * (1 - mu.phi.aut) * lsigma.phi.aut) / (mu.phi.spr * (1 - mu.phi.spr)))^2
    
    Sigma <- matrix(c(lsigma.phi.aut ^ 2, cor.spr * sqrt(lsigma.phi.aut ^ 2 * lsigma.phi.spr), 
                      cor.spr * sqrt(lsigma.phi.aut ^ 2 * lsigma.phi.spr), lsigma.phi.spr), nrow = 2, byrow = TRUE)
    epsilon <- MASS::mvrnorm(n = nYears, mu = c(0, 0), Sigma = Sigma, empirical = TRUE)

    ## Realized fall survival
    lphi.aut <- mean.phi.aut + epsilon[, 1]
    phi.aut <- 1 / (1 + exp(-lphi.aut))
    PHI.aut <- phi.aut ^ 2

    ## Realized spring survival
    lphi.spr <- mean.phi.spr + epsilon[, 2]
    PHI.spr <- 1 / (1 + exp(-lphi.spr))

    ## Stationary period survival
    PHI.sum <- phi.sum ^ 4
    PHI.win <- phi.win ^ 5
    
    ## Realized annual survival
    PHI.annual <- PHI.sum * PHI.aut * PHI.win * PHI.spr

    ## Realized migration survival
    PHI.mig <- PHI.spr * PHI.aut

    ## Phi matrices
    year.end <- 4 * (seq(from = 1, to = nYears - 1))
    phi.mat.sum <- numeric(length = nYears * 4)
    phi.mat.sum[seq(1, nYears * 4, by = 4)] <- PHI.sum
    phi.mat.sum[seq(2, nYears * 4, by = 4)] <- PHI.aut
    phi.mat.sum[seq(3, nYears * 4, by = 4)] <- PHI.win
    phi.mat.sum[seq(4, nYears * 4, by = 4)] <- PHI.spr
    
    phi.mat.win <- numeric(length = (nYears - 1) * 4)
    phi.mat.win[seq(1, (nYears - 1) * 4, by = 4)] <- PHI.win
    phi.mat.win[seq(2, (nYears - 1) * 4, by = 4)] <- PHI.spr[1:(nYears - 1)]
    phi.mat.win[seq(3, (nYears - 1) * 4, by = 4)] <- PHI.sum
    phi.mat.win[seq(4, (nYears - 1) * 4, by = 4)] <- PHI.aut[2:nYears]

    ## Simulate survival data
    z.sum <- matrix(0, nrow = sum(nInd, nImm), ncol = length(phi.mat.sum) + 1)
    z.win <- matrix(0, nrow = sum(nInd, nImm[-length(nImm)]), ncol = length(phi.mat.win) + 1)

    f.sum <- 1 + 4 * rep(seq(from = 0, to = nYears - 1), c(nInd, nImm))
    f.win <- 1 + 4 * rep(seq(from = 0, to = nYears - 2), c(nInd, nImm[-length(nImm)]))

    ## True state of summer individuals
    for(i in 1:sum(nInd, nImm)){
      z.sum[i, f.sum[i]] <- 1
      for(j in (f.sum[i] + 1):ncol(z.sum)){
        z.sum[i, j] <- rbinom(n = 1, 1, prob = phi.mat.sum[j - 1] * z.sum[i, j - 1])
      }
    }
    
    ## Add immigrants that arrive in final breeding season
    z.last <- matrix(0, nrow = rpois(1, omega), ncol = ncol(z.sum))
    z.last[,ncol(z.last)] <- 1
    z.sum <- rbind(z.sum, z.last)
    
    ## True state of winter individuals
    for(i in 1:sum(nInd, nImm[-length(nImm)])){
      z.win[i, f.win[i]] <- 1
      
      for(j in (f.win[i] + 1):ncol(z.win)){
        z.win[i, j] <- rbinom(n = 1, 1, prob = phi.mat.win[j - 1] * z.win[i, j - 1])
      }
    }
    
    ## Add immigrants that arrive in final winter
    z.last <- matrix(0, nrow = nImm[length(nImm)], ncol = ncol(z.win))
    z.last[,ncol(z.last)] <- 1
    z.win <- rbind(z.win, z.last)
    
    N.sum <- apply(z.sum, 2, sum)[seq(1, ncol(z.sum), 4)]
    N.win <- apply(z.win, 2, sum)[seq(1, ncol(z.win), 4)]

    
    ll <- list(z.sum = z.sum, z.win = z.win, 
               N.sum = N.sum, 
               N.win = N.win,
               PHI.annual = PHI.annual,
               PHI.mig = PHI.mig,
               PHI.spr = PHI.spr, 
               PHI.aut = PHI.aut, 
               PHI.sum = PHI.sum, 
               PHI.win = PHI.win,
               mean.phi.sum = phi.sum, mean.phi.win = phi.win,
               nYears = nYears, nInd.sum = sum(nInd, nImm), 
               nInd.win = sum(nInd, nImm[-length(nImm)]))
  
    return(ll) 
}

