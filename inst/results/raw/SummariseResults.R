### 'Basic' model ----
## Bias, RMSE, and power
dat <- read.csv(here::here("inst/results/FAC_results.csv"))
param_df <- expand.grid(rel_diff = c(1, 0.875, 0.75),
                        lsigma_phi_aut = c(0.02, 0.25, 0.5),
                        cov_spr = c(0, 0.4, 0.8))
param_df$sim <- seq(1:nrow(param_df))

dat1 <- dplyr::filter(dat, ts == "mean" & Season %in% c("Spring", "Autumn"))

dat2 <- dplyr::left_join(dat1, param_df)
dat2 <- dplyr::group_by(dat2, sim, Season)
dat3 <- dplyr::summarise(dat2, rel.bias = mean(bias),
                         rmse = sqrt(sum((est - true)^2)/100))
dat4 <- dplyr::left_join(dat3, param_df)

dat2 <- dplyr::ungroup(dat2)
dat5 <- dplyr::select(dat2, sim, it, Season, est)
dat6 <- tidyr::spread(dat5, key = Season, value = est)
dat7 <- dplyr:::mutate(dat6, direction = Autumn > Spring)
dat8 <- dplyr::group_by(dat7, sim)
dat9 <- dplyr::summarise(dat8, power = mean(direction))
dat10 <- dplyr::left_join(dat4, dat9)
write.csv(dat10, "inst/results/basic_mu.csv")

## Annual correlation
f <- function(true, est) cor(true, est)

dat1 <- dplyr::filter(dat, ts == "annual" & Season %in% c("Spring", "Autumn"))

dat2 <- dplyr::group_by(dat1, Season, sim, it, year)
dat3 <- dplyr::summarise(dat2, true = mean(true), est = mean(est))
dat3 <- dplyr::ungroup(dat3)
dat4 <- dplyr::group_by(dat3, sim, Season, it)
dat5 <- dplyr::summarise(dat4, rho = f(true, est))
dat5 <- dplyr::ungroup(dat5)
dat5 <- dplyr::group_by(dat5, sim, Season)
dat6 <- dplyr::summarise(dat5, LCI = quantile(rho, probs = 0.025),
                         UCI = quantile(rho, probs = 0.975), rho = mean(rho))

dat7 <- dplyr::left_join(dat6, param_df)
write.csv(dat7, "inst/results/basic_corr.csv")



### 'Covariate' model ----
## Bias, RMSE, and power
dat <- read.csv(here::here("inst/results/FAC_results_cov.csv"))
beta_df <- expand.grid(beta.spr = c(0.0, 0.5, 1.0),
                       beta.aut = c(0.0, 0.5, 1.0),
                       r.var = c(0.02, 0.25, 0.50))
beta_df$sim <- seq(1:nrow(beta_df))

dat1 <- dplyr::filter(dat, ts == "mean" & Season %in% c("Spring", "Autumn"))

dat2 <- dplyr::left_join(dat1, beta_df)
dat2 <- dplyr::group_by(dat2, sim, Season)
dat3 <- dplyr::summarise(dat2, rel.bias = mean(bias),
                         rmse = sqrt(sum((est - true)^2)/100))
dat4 <- dplyr::left_join(dat3, beta_df)

dat2 <- dplyr::ungroup(dat2)
dat5 <- dplyr::select(dat2, sim, it, Season, est)
dat6 <- tidyr::spread(dat5, key = Season, value = est)
dat7 <- dplyr:::mutate(dat6, direction = Autumn > Spring)
dat8 <- dplyr::group_by(dat7, sim)
dat9 <- dplyr::summarise(dat8, power = mean(direction))
dat10 <- dplyr::left_join(dat4, dat9)
write.csv(dat10, "inst/results/cov_mu.csv")


## Annual correlation
f <- function(true, est) cor(true, est)

dat1 <- dplyr::filter(dat, ts == "annual" & Season %in% c("Spring", "Autumn"))

dat2 <- dplyr::group_by(dat1, Season, sim, it, year)
dat3 <- dplyr::summarise(dat2, true = mean(true), est = mean(est))
dat4 <- dplyr::group_by(dat3, sim, Season, it)
dat5 <- dplyr::summarise(dat4, rho = f(true, est))
dat6 <- dplyr::summarise(dat5, LCI = quantile(rho, probs = 0.025),
                         UCI = quantile(rho, probs = 0.975), rho = mean(rho))

dat6 <- dplyr::left_join(dat6, beta_df)

write.csv(dat6, "inst/results/cov_corr.csv")
