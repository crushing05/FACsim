#' createCH
#'
#' Create capture history from z matrix
#' @param z list containing the output from sim_data()
#' @param p.sum Summer detection probability
#' @param p.win Winter detection probability
#' @param sigma.y Observation error for count data
#' @export

createCH <- function(z, p.sum = 0.6, p.win = 0.4, sigma.y = 2, monthly = FALSE, cov = FALSE){

  nYears <- z$nYears

  year.end <- 4 * (seq(from = 1, to = nYears) - 1)

    ## first and last summer/winter occasions in z matrix
    z.summer.occ <- c(sort(unlist(lapply(c(1, 2), function(x) x + year.end))), nYears * 4 + 1)
    z.winter.occ <- sort(unlist(lapply(c(1, 2), function(x) x + year.end)))

    ## Subset z matrix, including only first and last summer/winter occasions
    z.summer <- z$z.sum[, z.summer.occ]
    z.winter <- z$z.win[, z.winter.occ[-length(z.winter.occ)]]

    ## First capture of each individual
    f.summer <- apply(z.summer, 1, function(x) min(which(x == 1)))
    f.winter <- apply(z.winter, 1, function(x) min(which(x == 1)))


    ## Create empty CHs
    CH.summer <- matrix(0, nrow = z$nInd.sum, ncol = dim(z.summer)[2])
    CH.winter <- matrix(0, nrow = z$nInd.win, ncol = dim(z.winter)[2])

    ## Generate CH from z matrix and detection probabilities
    for(i in 1:z$nInd.sum){
      CH.summer[i, f.summer[i]] <- 1

      ## Simulate summer captures based on true state
      for(j in (f.summer[i] + 1):ncol(CH.summer)){
        CH.summer[i, j] <- rbinom(n = 1, size = 1, prob = z.summer[i, j] * p.sum)
      }
    }
  
    CH.summer <- rbind(CH.summer, z.summer[(z$nInd.sum + 1):nrow(z.summer),])
    
    for(i in 1:z$nInd.win){
      CH.winter[i, f.winter[i]] <- 1
      
      ## Simulate winter captures based on true state
      for(j in (f.winter[i] + 1):ncol(CH.winter)){
        CH.winter[i, j] <- rbinom(n = 1, size = 1, prob = z.winter[i, j] * p.win)
      }
    }
    
    CH.winter <- rbind(CH.winter, z.winter[(z$nInd.win + 1):nrow(z.winter),])
  
  ## Convert capture histories to m-arrays
  marr.sum <- FAC.sim::marray(CH.summer)
  r.sum <- rowSums(marr.sum)

  marr.win <- FAC.sim::marray(CH.winter)
  r.win <- rowSums(marr.win)

  if(!cov){
    ## Sample breeding and winter populations
    C.sum <- round(rnorm(n = length(z$N.sum), mean = z$N.sum, sd = sigma.y))
    C.win <- round(rnorm(n = length(z$N.win), mean = z$N.win, sd = sigma.y))
    
    ll <- list(marr.sum = marr.sum, r.sum = r.sum, z.summer = z.summer,
               marr.win = marr.win, r.win = r.win, z.winter = z.winter,
               C.sum= C.sum, C.win = C.win)
  }else{
    ll <- list(marr.sum = marr.sum, r.sum = r.sum, z.summer = z.summer,
               marr.win = marr.win, r.win = r.win, z.winter = z.winter)
  }


  return(ll)
}


#' marray
#'
#' Convert capture history to m-array format
#' @param CH nind x noccasions capture history
#' @export

marray <- function(CH){
  nind <- dim(CH)[1]
  n.occasions <- dim(CH)[2]
  m.array <- matrix(data = 0, ncol = n.occasions+1, nrow = n.occasions)
  # Calculate the number of released individuals at each time period
  for (t in 1:n.occasions){
    m.array[t,1] <- sum(CH[,t])
  }
  for (i in 1:nind){
    pos <- which(CH[i,]!=0)
    g <- length(pos)
    for (z in 1:(g-1)){
      m.array[pos[z],pos[z+1]] <- m.array[pos[z],pos[z+1]] + 1
    } #z
  } #i
  # Calculate the number of individuals that were never recaptured
  for (t in 1:n.occasions){
    m.array[t,n.occasions+1] <- m.array[t,1] - sum(m.array[t,2:n.occasions])
  }
  out <- m.array[1:(n.occasions-1),2:(n.occasions+1)]
  return(out)
}
