#' plot_joint1
#'
#' Plot bias & RMSE of basic integrated CJS model
#' @export

plot_joint1 <- function(){
  dat <- read.csv(here::here("inst/results/raw/FAC_results.csv"))
  param_df <- expand.grid(rel_diff = c(1, 0.875, 0.75),
                          lsigma_phi_aut = c(0.02, 0.25, 0.5),
                          cov_spr = c(0, 0.4, 0.8))
  param_df$sim <- seq(1:nrow(param_df))

  dat1 <- dplyr::filter(dat, ts == "mean" & Season %in% c("Spring", "Autumn"))

  dat2 <- suppressMessages(dplyr::left_join(dat1, param_df))
  dat2 <- dplyr::group_by(dat2, sim, Season)
  dat3 <- dplyr::summarise(dat2, rel.bias = mean(bias),
                           rmse = sqrt(sum((est - true)^2)/250))
  dat4 <- suppressMessages(dplyr::left_join(dat3, param_df))

  dat4$sigma <- paste("sigma^2 == ", dat4$lsigma_phi_aut)

  dat4 <- dplyr::filter(dat4, cov_spr == 0)

  p <- ggplot(dat4, aes(x = rel_diff, y = rmse, color = Season, linetype = Season, group = Season))
  p <- p + facet_grid( ~ sigma, labeller = label_parsed)
  p <- p + scale_y_continuous("Root mean square error", limits = c(0, 0.2),
                              breaks = seq(from = 0, to = 0.2, by = 0.025),
                              labels = c(0, "", 0.05, "", 0.10, "", 0.15, "", 0.2))
  p <- p + scale_x_continuous(expression(paste("Relative difference (", Delta, ")")))
  p <- p + geom_path()
  p <- p + geom_point(color = "white", size = 5)
  p <- p + geom_point(size = 4)
  p <- p + scale_color_manual(values = c('#cb4b16', '#e5ca28'))
  p <- p + theme(panel.grid.major = element_line(color = "grey90", size = 0.5),
                 panel.grid.major.x = element_blank(),
                 strip.text.x = element_blank(),
                 axis.line.x = element_line(color = "black"))
  p <- p + guides(color = FALSE, linetype = FALSE)
  p <- p + annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf, size = 1.5)
  p <- p + annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf, size = 1.5)
  p <- p + geom_hline(yintercept = 0, color = "#777777")

  q <- ggplot(dat4, aes(x = rel_diff, y = rel.bias, color = Season, linetype = Season, group = Season))
  q <- q + facet_grid( ~ sigma, labeller = label_parsed)
  q <- q + scale_y_continuous("Relative bias", limits = c(-0.1, 0.20),
                              breaks = seq(from = -0.2, to = 0.2, by = 0.05),
                              labels = c(-0.2, "", -0.1, "", 0, "", 0.1, "", 0.2))
  q <- q + scale_x_continuous("")
  q <- q + geom_hline(yintercept = 0, color = "#777777")
  q <- q + geom_path()
  q <- q + geom_point(color = "white", size = 5)
  q <- q + geom_point(size = 4)
  q <- q + scale_color_manual(values = c('#cb4b16', '#e5ca28'))
  q <- q + theme(panel.grid.major = element_line(color = "grey90", size = 0.5),
                 panel.grid.major.x = element_blank(),
                 axis.line.x = element_line(color = "black"),
                 legend.position = c(.94, .75), legend.background = element_rect(fill = NA),
                 legend.text = element_text(size = 12), legend.key = element_rect(color = "grey70", size = 0.5))
  q <- q + annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf, size = 1.5)
  q <- q + annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf, size = 1.5)
  cowplot::plot_grid(q, p, align = "h", nrow = 2)
}


#' plot_cor1
#'
#' Plot correlation between true and estimated annual survival estimates from basic integrated CJS model
#' @export

plot_cor1 <- function(){
  dat <- read.csv(here::here("inst/results/raw/FAC_results.csv"))
  param_df <- expand.grid(rel_diff = c(1, 0.875, 0.75),
                          lsigma_phi_aut = c(0.02, 0.25, 0.5),
                          cov_spr = c(0, 0.4, 0.8))
  param_df$sim <- seq(1:nrow(param_df))
  f <- function(true, est) cor(true, est)

  dat1 <- dplyr::filter(dat, ts == "annual" & Season %in% c("Spring", "Autumn"))

  dat2 <- dplyr::group_by(dat1, Season, sim, it, year)
  dat3 <- dplyr::summarise(dat2, true = mean(true), est = mean(est))
  dat4 <- dplyr::group_by(dat3, sim, Season, it)
  dat5 <- dplyr::summarise(dat4, rho = f(true, est))
  dat6 <- dplyr::summarise(dat5, LCI = quantile(rho, probs = 0.025),
                           UCI = quantile(rho, probs = 0.975), rho = mean(rho))

  dat6 <- dplyr::mutate(dat6, label = paste("~italic(r)[",substr(Season, 1, 3),"] ==",
                                            round(rho, digits = 2), " ~~(",
                                            round(LCI, digits = 2), ":",
                                            round(UCI, digits = 2), ")", sep = ""))
  dat4 <- suppressMessages(dplyr::left_join(dat4, param_df))
  dat4 <- dplyr::filter(dat4, rel_diff == 0.75 & cov_spr == 0)
  dat4$sigma <- paste("sigma^2 == ", dat4$lsigma_phi_aut)
  dat6$x <- max(dat4$true) - 0.25
  dat6$y <- ifelse(dat6$Season == "Autumn", min(dat4$est) + 0.05, 0.27)
  dat6 <- dplyr::select(dat6, -Season)
  dat6 <- suppressMessages(dplyr::left_join(dat6, param_df))
  dat6 <- dplyr::filter(dat6, rel_diff == 0.75 & cov_spr == 0)
  dat6$sigma <- paste("sigma^2 == ", dat6$lsigma_phi_aut)

  dat7 <- data.frame(x = rep(2, 2), y = rep(2, 2), Season = c("Spring", "Autumn"))

  p <- ggplot(dat4, aes(x = true, y = est))
  p <- p + geom_abline(slope = 1, intercept = 0, color = "grey60")
  p <- p + geom_point(data = dat7, aes(x = x, y = y, color = Season, group = Season))
  p <- p + facet_wrap(~ sigma, nrow = 1, labeller = label_parsed)
  p <- p + scale_y_continuous(expression(paste("Estimated ", phi)), limits = c(0.25, 1))
  p <- p + scale_x_continuous(expression(paste("True ", phi)), limits = c(0.25, 1))
  p <- p + geom_point(size = 3, alpha = 0.2, aes(color = Season))
  p <- p + stat_smooth(method = "lm", aes(linetype = Season, group = Season), color = "black", se = FALSE)
  p <- p + scale_color_manual(values = c('#cb4b16', '#e5ca28'))
  p <- p + geom_text(data = dat6, aes(x = x, y = y, label = label, group = label), parse = TRUE, size = 6)
  p <- p + theme(legend.key = element_rect(fill = NA))
  p
}


#' plot_joint2
#'
#' Plot bias & RMSE of integrated CJS model w/ covariates
#' @export

plot_joint2 <- function(){
  dat <- read.csv(here::here("inst/results/raw/FAC_results_cov.csv"))

  beta_df <- expand.grid(beta.aut = c(0.0, 0.5, 1.0),
                         beta.spr = c(0.0, 0.5, 1.0),
                         r.var = c(0.02, 0.25, 0.50))
  beta_df$sim <- seq(1:nrow(beta_df))
  dat1 <- dplyr::filter(dat, ts == "mean" & Season %in% c("Spring", "Autumn"))

  dat2 <- suppressMessages(dplyr::left_join(dat1, beta_df))
  dat2 <- dplyr::group_by(dat2, sim, Season)
  dat3 <- dplyr::summarise(dat2, rel.bias = mean(bias),
                           rmse = sqrt(sum((est - true)^2)/250))
  dat4 <- suppressMessages(dplyr::left_join(dat3, beta_df))

  dat4$beta1 <- paste("beta[Autumn] == ", dat4$beta.aut)
  dat4$beta2 <- paste("beta[Spring] == ", dat4$beta.spr)
  dat4$sigma <- paste("sigma^2 == ", dat4$r.var)

  dat5 <- dplyr::filter(dat4, sim %in% c(1, 5, 9, 10, 14, 18, 19, 23, 27))
  dat5$beta <- dat5$beta.aut

  p <- ggplot(dat5, aes(x = beta, y = rmse, color = Season, linetype = Season, group = Season))
  p <- p + facet_grid( ~ sigma, labeller = label_parsed)
  p <- p + scale_y_continuous("Root mean square error", limits = c(0, 0.2),
                              breaks = seq(from = 0, to = 0.2, by = 0.025),
                              labels = c(0, "", 0.05, "", 0.10, "", 0.15, "", 0.2))
  p <- p + scale_x_continuous(expression(paste("Strength of covariate (", beta, ")")))
  p <- p + geom_path()
  p <- p + geom_point(color = "white", size = 5)
  p <- p + geom_point(size = 4)
  p <- p + scale_color_manual(values = c('#cb4b16', '#e5ca28'))
  p <- p + theme(panel.grid.major = element_line(color = "grey90", size = 0.5),
                 panel.grid.major.x = element_blank(),
                 strip.text.x = element_blank(),
                 axis.line.x = element_line(color = "black"))
  p <- p + guides(color = FALSE, linetype = FALSE)
  p <- p + annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf, size = 1.5)
  p <- p + annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf, size = 1.5)
  p <- p + geom_hline(yintercept = 0, color = "#777777")

  q <- ggplot(dat5, aes(x = beta, y = rel.bias, color = Season, linetype = Season, group = Season))
  q <- q + facet_grid( ~ sigma, labeller = label_parsed)
  q <- q + scale_y_continuous("Relative bias", limits = c(-0.1, 0.20),
                              breaks = seq(from = -0.2, to = 0.2, by = 0.05),
                              labels = c(-0.2, "", -0.1, "", 0, "", 0.1, "", 0.2))
  q <- q + scale_x_continuous("")
  q <- q + geom_hline(yintercept = 0, color = "#777777")
  q <- q + geom_path()
  q <- q + geom_point(color = "white", size = 5)
  q <- q + geom_point(size = 4)
  q <- q + scale_color_manual(values = c('#cb4b16', '#e5ca28'))
  q <- q + theme(panel.grid.major = element_line(color = "grey90", size = 0.5),
                 panel.grid.major.x = element_blank(),
                 axis.line.x = element_line(color = "black"),
                 legend.position = c(.94, .75), legend.background = element_rect(fill = NA),
                 legend.text = element_text(size = 10), legend.key = element_rect(color = "grey70", size = 0.5))
  q <- q + annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf, size = 1.5)
  q <- q + annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf, size = 1.5)
  cowplot::plot_grid(q, p, align = "h", nrow = 2)
}

#' plot_cor2
#'
#' Plot correlation between true and estimated annual survival estimates from integrated CJS model with covariates
#' @export

plot_cor2 <- function(){
  dat <- read.csv(here::here("inst/results/raw/FAC_results_cov.csv"))

  beta_df <- expand.grid(beta.aut = c(0.0, 0.5, 1.0),
                         beta.spr = c(0.0, 0.5, 1.0),
                         r.var = c(0.02, 0.25, 0.50))
  beta_df$sim <- seq(1:nrow(beta_df))
  f <- function(true, est) cor(true, est)

  dat1 <- dplyr::filter(dat, ts == "annual" & Season %in% c("Spring", "Autumn"))

  dat2 <- dplyr::group_by(dat1, Season, sim, it, year)
  dat3 <- dplyr::summarise(dat2, true = mean(true), est = mean(est))
  dat4 <- dplyr::group_by(dat3, sim, Season, it)
  dat5 <- dplyr::summarise(dat4, rho = f(true, est))
  dat6 <- dplyr::summarise(dat5, LCI = quantile(rho, probs = 0.025),
                           UCI = quantile(rho, probs = 0.975), rho = mean(rho))

  dat6 <- dplyr::mutate(dat6, label = paste("~italic(r)[",substr(Season, 1, 3),"] ==",
                                            round(rho, digits = 2), " ~~(",
                                            round(LCI, digits = 2), ":",
                                            round(UCI, digits = 2), ")", sep = ""))
  dat4 <- suppressMessages(dplyr::left_join(dat4, beta_df))
  dat4 <- dplyr::filter(dat4, r.var == 0.25)
  dat4$sigma <- paste("sigma^2 == ", dat4$r.var)
  dat4$beta1 <- paste("beta[Autumn] == ", dat4$beta.aut)
  dat4$beta2 <- paste("beta[Spring] == ", dat4$beta.spr)

  dat6$x <- max(dat4$true) - 0.35
  dat6$y <- ifelse(dat6$Season == "Autumn", min(dat4$est) + 0.1, min(dat4$est) + 0.01)
  dat6 <- dplyr::select(dat6, -Season)
  dat6 <- suppressMessages(dplyr::left_join(dat6, beta_df))
  dat6 <- dplyr::filter(dat6, r.var == 0.25)
  dat6$sigma <- paste("sigma^2 == ", dat6$r.var)
  dat6$beta1 <- paste("beta[Autumn] == ", dat6$beta.aut)
  dat6$beta2 <- paste("beta[Spring] == ", dat6$beta.spr)


  dat7 <- dplyr::filter(dat4, sim %in% c(1, 5, 9, 10, 14, 18, 19, 23, 27))
  dat7$beta <- paste("beta == ", dat7$beta.aut)

  dat8 <- dplyr::filter(dat6, sim %in% c(1, 5, 9, 10, 14, 18, 19, 23, 27))
  dat8$beta <- paste("beta == ", dat8$beta.aut)

  dat9 <- data.frame(x = rep(2, 2), y = rep(2, 2), Season = c("Spring", "Autumn"))

  p <- ggplot(dat7, aes(x = true, y = est, group = Season))
  p <- p + geom_abline(slope = 1, intercept = 0, color = "grey60")
  p <- p + geom_point(data = dat9, aes(x = x, y = y, color = Season))
  p <- p + facet_grid( ~ beta, labeller = label_parsed)
  p <- p + scale_y_continuous(expression(paste("Estimated ", phi)), limits = c(0, 1))
  p <- p + scale_x_continuous(expression(paste("True ", phi)), limits = c(0, 1))
  p <- p + geom_point(size = 3, alpha = 0.2, aes(color = Season))
  p <- p + stat_smooth(method = "lm", aes(linetype = Season, group = Season), color = "black", se = FALSE)
  p <- p + scale_color_manual(values = c('#cb4b16', '#e5ca28'))
  p <- p + geom_text(data = dat8, aes(x = x, y = y, label = label, group = label), parse = TRUE, size = 5)
  p
}


#' plot_joint3
#'
#' Plot bias & RMSE of nYear models
#' @export

plot_joint3 <- function(){
  dat <- read.csv(here::here("inst/results/raw/FAC_results_nYears.csv"))


  dat1 <- dplyr::filter(dat, ts == "mean" & Season %in% c("Spring", "Autumn"))

  dat2 <- dplyr::group_by(dat1, sim, Season, nYears)
  dat3 <- dplyr::summarise(dat2, rel.bias = mean(bias),
                           rmse = sqrt(sum((est - true)^2)/250))
  dat3$rel.bias[dat3$Season=="Spring"] <- abs(dat3$rel.bias[dat3$Season=="Spring"])

  p <- ggplot(dat3, aes(x = nYears, y = rmse, color = Season, linetype = Season, group = Season))
  p <- p + scale_y_continuous("Root mean square error", limits = c(0, 0.2),
                              breaks = seq(from = 0, to = 0.2, by = 0.025),
                              labels = c(0, "", 0.05, "", 0.10, "", 0.15, "", 0.2))
  p <- p + scale_x_continuous("Number of Years")
  p <- p + geom_path()
  p <- p + geom_point(color = "white", size = 5)
  p <- p + geom_point(size = 4)
  p <- p + scale_color_manual(values = c('#cb4b16', '#e5ca28'))
  p <- p + theme(panel.grid.major = element_line(color = "grey90", size = 0.5),
                 panel.grid.major.x = element_blank(),
                 strip.text.x = element_blank(),
                 axis.line.x = element_line(color = "black"))
  p <- p + guides(color = FALSE, linetype = FALSE)
  p <- p + annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf, size = 1.5)
  p <- p + annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf, size = 1.5)
  p <- p + geom_hline(yintercept = 0, color = "#777777")

  q <- ggplot(dat3, aes(x = nYears, y = rel.bias, color = Season, linetype = Season, group = Season))
  q <- q + scale_y_continuous("Relative bias", limits = c(-0.1, 0.20),
                              breaks = seq(from = -0.2, to = 0.2, by = 0.05),
                              labels = c(-0.2, "", -0.1, "", 0, "", 0.1, "", 0.2))
  q <- q + scale_x_continuous("")
  q <- q + geom_hline(yintercept = 0, color = "#777777")
  q <- q + geom_path()
  q <- q + geom_point(color = "white", size = 5)
  q <- q + geom_point(size = 4)
  q <- q + scale_color_manual(values = c('#cb4b16', '#e5ca28'))
  q <- q + theme(panel.grid.major = element_line(color = "grey90", size = 0.5),
                 panel.grid.major.x = element_blank(),
                 axis.line.x = element_line(color = "black"),
                 legend.position = c(0.9, 0.8), legend.background = element_rect(fill = NA),
                 legend.text = element_text(size = 12), legend.key = element_rect(color = "grey70", size = 0.5))
  q <- q + annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf, size = 1.5)
  q <- q + annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf, size = 1.5)
  cowplot::plot_grid(q, p, align = "h", nrow = 2)
}

#' plot_cor3
#'
#' Plot correlation between true and estimated annual survival estimates from integrated CJS model with 1 season count data
#' @export

plot_cor3 <- function(){
  dat <- read.csv(here::here("inst/results/raw/FAC_results_nYears.csv"))

  f <- function(true, est) cor(true, est)

  dat1 <- dplyr::filter(dat, ts == "annual" & Season %in% c("Spring", "Autumn"))

  dat2 <- dplyr::group_by(dat1, Season, nYears, it, year)
  dat3 <- dplyr::summarise(dat2, true = mean(true), est = mean(est))
  dat4 <- dplyr::group_by(dat3, nYears, Season, it)
  dat4$label <- factor(paste("Number of years =", dat4$nYears))
  levels(dat4$label) <- unique(dat4$label)
  dat5 <- dplyr::summarise(dat4, rho = f(true, est))
  dat6 <- dplyr::summarise(dat5, LCI = quantile(rho, probs = 0.025),
                           UCI = quantile(rho, probs = 0.975), rho = mean(rho))

  dat6$label <- factor(paste("Number of years =", dat6$nYears))
  levels(dat6$label) <- unique(dat6$label)
  dat6 <- dplyr::mutate(dat6, label2 = paste("~italic(r)[",substr(Season, 1, 3),"] ==",
                                            round(rho, digits = 2), " ~~(",
                                            round(LCI, digits = 2), ":",
                                            round(UCI, digits = 2), ")", sep = ""))

  dat6$x <- max(dat4$true) - 0.28
  dat6$y <- ifelse(dat6$Season == "Autumn", min(dat4$est) + 0.03, 0.265)
  dat6 <- dplyr::select(dat6, -Season)
  dat7 <- data.frame(x = rep(2, 2), y = rep(2, 2), Season = c("Spring", "Autumn"))

  p <- ggplot(dat4, aes(x = true, y = est, group = Season))
  p <- p + geom_abline(slope = 1, intercept = 0, color = "grey60")
  p <- p + geom_point(data = dat7, aes(x = x, y = y, color = Season))
  p <- p + facet_wrap(~ label, nrow = 2)#, labeller = label_parsed)
  p <- p + scale_y_continuous(expression(paste("Estimated ", phi)), limits = c(0.25, 1))
  p <- p + scale_x_continuous(expression(paste("True ", phi)), limits = c(0.25, 1))
  p <- p + geom_point(size = 3, alpha = 0.2, aes(color = Season))
  p <- p + stat_smooth(method = "lm", aes(linetype = Season, group = Season), color = "black", se = FALSE)
  p <- p + scale_color_manual(values = c('#cb4b16', '#e5ca28'))
  p <- p + geom_text(data = dat6, aes(x = x, y = y, label = label2, group = label2), parse = TRUE, size = 5)
  p <- p + theme(strip.text.x = element_text(size = 12))
  p
}



#' plot_S1
#'
#' Plot relative bias of basic integrated CJS model
#' @export

plot_S1 <- function(){
  dat <- read.csv(here::here("inst/results/raw/FAC_results.csv"))
  param_df <- expand.grid(rel_diff = c(1, 0.875, 0.75),
                          lsigma_phi_aut = c(0.02, 0.25, 0.5),
                          cov_spr = c(0, 0.4, 0.8))
  param_df$sim <- seq(1:nrow(param_df))

  dat1 <- suppressWarnings(dplyr::filter(dat, ts == "mean" & Season %in% c("Spring", "Autumn")))

  dat2 <- suppressMessages(dplyr::left_join(dat1, param_df))
  dat2 <- dplyr::group_by(dat2, sim, Season)
  dat3 <- dplyr::summarise(dat2, rel.bias = mean(bias))
  dat4 <- suppressMessages(dplyr::left_join(dat3, param_df))

  dat4$rho <- paste("rho == ", dat4$cov_spr)
  dat4$delta <- paste("Delta == ", dat4$rel_diff)

  p <- ggplot(dat4, aes(x = lsigma_phi_aut, y = rel.bias, color = Season, group = Season))
  p <- p + facet_grid(rho ~ delta, labeller = label_parsed)
  p <- p + scale_y_continuous("Relative bias", limits = c(-0.2, 0.250),
                              breaks = seq(from = -0.2, to = 0.25, by = 0.05),
                              labels = c(-0.2, "", -0.1, "", 0, "", 0.1, "", 0.2, ""))
  p <- p + scale_x_continuous(expression(paste("Annual variation (", sigma^2, ")")))
  p <- p + geom_hline(yintercept = 0, color = "#777777")
  p <- p + geom_path()
  p <- p + geom_point(color = "white", size = 4)
  p <- p + geom_point()
  p <- p + scale_color_manual(values = c('#cb4b16', '#e5ca28'))
  p <- p + theme(panel.grid.major = element_line(color = "grey90", size = 0.25),
                 panel.grid.major.x = element_blank())
  p
}


#' plot_S2
#'
#' Plot RMSE of basic integrated CJS model
#' @export

plot_S2 <- function(){
  dat <- read.csv(here::here("inst/results/raw/FAC_results.csv"))
  param_df <- expand.grid(rel_diff = c(1, 0.875, 0.75),
                          lsigma_phi_aut = c(0.02, 0.25, 0.5),
                          cov_spr = c(0, 0.4, 0.8))
  param_df$sim <- seq(1:nrow(param_df))

  dat1 <- dplyr::filter(dat, ts == "mean" & Season %in% c("Spring", "Autumn"))

  dat2 <- suppressMessages(dplyr::left_join(dat1, param_df))
  dat2 <- dplyr::group_by(dat2, sim, Season)
  dat3 <- dplyr::summarise(dat2, rel.bias = mean(bias),
                           rmse = sqrt(sum((est - true)^2)/250))
  dat4 <- suppressMessages(dplyr::left_join(dat3, param_df))

  dat4$rho <- paste("rho == ", dat4$cov_spr)
  dat4$delta <- paste("Delta == ", dat4$rel_diff)


  p <- ggplot(dat4, aes(x = lsigma_phi_aut, y = rmse, color = Season, group = Season))
  p <- p + facet_grid(rho ~ delta, labeller = label_parsed)
  p <- p + scale_y_continuous("Root mean square error", limits = c(0, 0.2))
  p <- p + scale_x_continuous(expression(paste("Relative difference (", Delta, ")")))
  p <- p + geom_path()
  p <- p + geom_point(color = "white", size = 4)
  p <- p + geom_point()
  p <- p + scale_color_manual(values = c('#cb4b16', '#e5ca28'))
  p <- p + theme(panel.grid.major = element_line(color = "grey90"),
                 panel.grid.major.x = element_blank())
  p <- p + geom_hline(yintercept = 0)
  p
}

#' plot_S3
#'
#' Plot correlation between true and estimated annual survival estimates from basic integrated CJS model
#' @export

plot_S3 <- function(){
  dat <- read.csv(here::here("inst/results/raw/FAC_results.csv"))
  param_df <- expand.grid(rel_diff = c(1, 0.875, 0.75),
                          lsigma_phi_aut = c(0.02, 0.25, 0.5),
                          cov_spr = c(0, 0.4, 0.8))
  param_df$sim <- seq(1:nrow(param_df))
  f <- function(true, est) cor(true, est)

  dat1 <- dplyr::filter(dat, ts == "annual" & Season %in% c("Spring", "Autumn"))

  dat2 <- dplyr::group_by(dat1, Season, sim, it, year)
  dat3 <- dplyr::summarise(dat2, true = mean(true), est = mean(est))
  dat4 <- dplyr::group_by(dat3, sim, Season, it)
  dat5 <- dplyr::summarise(dat4, rho = f(true, est))
  dat6 <- dplyr::summarise(dat5, LCI = quantile(rho, probs = 0.025),
                           UCI = quantile(rho, probs = 0.975), rho = mean(rho))

  dat6 <- dplyr::mutate(dat6, label = paste("~italic(r)[",substr(Season, 1, 3),"] ==",
                                            round(rho, digits = 2), " ~~(",
                                            round(LCI, digits = 2), ":",
                                            round(UCI, digits = 2), ")", sep = ""))
  dat4 <- suppressMessages(dplyr::left_join(dat4, param_df))
  dat4 <- dplyr::filter(dat4, cov_spr == 0)
  dat4$sigma <- paste("sigma^2 == ", dat4$lsigma_phi_aut)
  dat4$delta <- paste("Delta == ", dat4$rel_diff)
  dat6$x <- max(dat4$true) - 0.25
  dat6$y <- ifelse(dat6$Season == "Autumn", min(dat4$est) + 0.06, 0.26)
  dat6 <- dplyr::select(dat6, -Season)
  dat6 <- suppressMessages(dplyr::left_join(dat6, param_df))
  dat6 <- dplyr::filter(dat6, cov_spr == 0)
  dat6$sigma <- paste("sigma^2 == ", dat6$lsigma_phi_aut)
  dat6$delta <- paste("Delta == ", dat6$rel_diff)

  dat7 <- data.frame(x = rep(2, 2), y = rep(2, 2), Season = c("Spring", "Autumn"))

  p <- ggplot(dat4, aes(x = true, y = est, group = Season))
  p <- p + geom_abline(slope = 1, intercept = 0, color = "grey60")
  p <- p + geom_point(data = dat7, aes(x = x, y = y, color = Season))
  p <- p + facet_grid(delta ~ sigma, labeller = label_parsed)
  p <- p + scale_y_continuous(expression(paste("Estimated ", phi)), limits = c(0.25, 1))
  p <- p + scale_x_continuous(expression(paste("True ", phi)), limits = c(0.25, 1))
  p <- p + geom_point(size = 3, alpha = 0.2, aes(color = Season))
  p <- p + stat_smooth(method = "lm", aes(linetype = Season, group = Season), color = "black", se = FALSE)
  p <- p + scale_color_manual(values = c('#cb4b16', '#e5ca28'))
  p <- p + geom_text(data = dat6, aes(x = x, y = y, label = label, group = label), parse = TRUE, size = 5)
  p <- p + theme(legend.key = element_rect(fill = NA))
  p
}

#' plot_S4
#'
#' Plot relative bias of integrated CJS model with covariates
#' @export

plot_S4 <- function(){
  dat <- read.csv(here::here("inst/results/raw/FAC_results_cov.csv"))

  beta_df <- expand.grid(beta.spr = c(0.0, 0.5, 1.0),
                         beta.aut = c(0.0, 0.5, 1.0),
                         r.var = c(0.02, 0.25, 0.50))
  beta_df$sim <- seq(1:nrow(beta_df))

  dat1 <- dplyr::filter(dat, ts == "mean" & Season %in% c("Spring", "Autumn"))

  dat2 <- dplyr::group_by(dat1, sim, Season)
  dat3 <- dplyr::summarise(dat2, rel.bias = mean(bias))
  dat4 <- suppressMessages(dplyr::left_join(dat3, beta_df))

  dat4$beta1 <- paste("beta[Aut] == ", dat4$beta.aut)
  dat4$beta2 <- paste("beta[Spr] == ", dat4$beta.spr)
  dat4$sigma <- paste("sigma^2 == ", dat4$r.var)


  p <- ggplot(dat4, aes(x = r.var, y = rel.bias, color = Season, linetype = Season, group = Season, label = sim))
  p <- p + facet_grid(beta1 ~ beta2, labeller = label_parsed)
  p <- p + scale_y_continuous("Relative bias", limits = c(-0.1, 0.2),
                              breaks = seq(from = -0.1, to = 0.2, by = 0.05),
                              labels = c(-0.1, "", 0, "", 0.1, "", 0.2))
  p <- p + scale_x_continuous(expression(paste("Residual variation (", sigma^2, ")")))
  p <- p + geom_hline(yintercept = 0, color = "#777777")
  p <- p + geom_path()
  p <- p + geom_point(color = "white", size = 6)
  p <- p + geom_point(size = 4.5)
  p <- p + scale_color_manual(values = c('#cb4b16', '#e5ca28'))
  p <- p + theme(panel.grid.major = element_line(size = 0.25, color = "grey90"),
                 panel.grid.major.x = element_blank())
  p
}


#' plot_S5
#'
#' Plot RMSE of integrated CJS model with covariates
#' @export

plot_S5 <- function(){
  dat <- read.csv(here::here("inst/results/raw/FAC_results_cov.csv"))

  beta_df <- expand.grid(beta.spr = c(0.0, 0.5, 1.0),
                         beta.aut = c(0.0, 0.5, 1.0),
                         r.var = c(0.02, 0.25, 0.50))
  beta_df$sim <- seq(1:nrow(beta_df))

  dat1 <- dplyr::filter(dat, ts == "mean" & Season %in% c("Spring", "Autumn"))

  dat2 <- dplyr::group_by(dat1, sim, Season)
  dat3 <- dplyr::summarise(dat2, rmse = sqrt(sum((est - true)^2)/250))
  dat4 <- suppressMessages(dplyr::left_join(dat3, beta_df))

  dat4$beta1 <- paste("beta[Aut] == ", dat4$beta.aut)
  dat4$beta2 <- paste("beta[Spr] == ", dat4$beta.spr)
  dat4$sigma <- paste("sigma^2 == ", dat4$r.var)

  p <- ggplot(dat4, aes(x = r.var, y = rmse, color = Season, group = Season))
  p <- p + facet_grid(beta1 ~ beta2, labeller = label_parsed)
  p <- p + scale_y_continuous("Root mean square error", limits = c(0, 0.20))
  p <- p + scale_x_continuous(expression(paste("Residual variation (", sigma^2, ")")))
  p <- p + geom_path()
  p <- p + geom_point(color = "white", size = 6)
  p <- p + geom_point(size = 4.5)
  p <- p + scale_color_manual(values = c('#cb4b16', '#e5ca28'))
  p <- p + theme(panel.grid.major = element_line(color = "grey80"),
                 panel.grid.major.x = element_blank())
  p
}

#' plot_S6
#'
#' Plot correlation between true and estimated annual survival estimates from integrated CJS model with covariates
#' @export

plot_S6 <- function(){
  dat <- read.csv(here::here("inst/results/raw/FAC_results_cov.csv"))

  beta_df <- expand.grid(beta.spr = c(0.0, 0.5, 1.0),
                         beta.aut = c(0.0, 0.5, 1.0),
                         r.var = c(0.02, 0.25, 0.50))
  beta_df$sim <- seq(1:nrow(beta_df))
  f <- function(true, est) cor(true, est)

  dat1 <- dplyr::filter(dat, ts == "annual" & Season %in% c("Spring", "Autumn"))

  dat2 <- dplyr::group_by(dat1, Season, sim, it, year)
  dat3 <- dplyr::summarise(dat2, true = mean(true), est = mean(est))
  dat4 <- dplyr::group_by(dat3, sim, Season, it)
  dat5 <- dplyr::summarise(dat4, rho = f(true, est))
  dat6 <- dplyr::summarise(dat5, LCI = quantile(rho, probs = 0.025),
                           UCI = quantile(rho, probs = 0.975), rho = mean(rho))

  dat6 <- dplyr::mutate(dat6, label = paste("~italic(r)[",substr(Season, 1, 3),"] ==",
                                            round(rho, digits = 2), " ~~(",
                                            round(LCI, digits = 2), ":",
                                            round(UCI, digits = 2), ")", sep = ""))
  dat4 <- suppressMessages(dplyr::left_join(dat4, beta_df))
  dat4 <- dplyr::filter(dat4, r.var == 0.25)
  dat4$beta1 <- paste("beta[Aut] == ", dat4$beta.aut)
  dat4$beta2 <- paste("beta[Spr] == ", dat4$beta.spr)

  dat6$x <- max(dat4$true) - 0.35
  dat6$y <- ifelse(dat6$Season == "Autumn", min(dat4$est) + 0.1, min(dat4$est) + 0.01)
  dat6 <- dplyr::select(dat6, -Season)
  dat6 <- suppressMessages(dplyr::left_join(dat6, beta_df))
  dat6 <- dplyr::filter(dat6, r.var == 0.25)
  dat6$beta1 <- paste("beta[Aut] == ", dat6$beta.aut)
  dat6$beta2 <- paste("beta[Spr] == ", dat6$beta.spr)


  dat9 <- data.frame(x = rep(2, 2), y = rep(2, 2), Season = c("Spring", "Autumn"))

  p <- ggplot(dat4, aes(x = true, y = est, group = Season))
  p <- p + geom_abline(slope = 1, intercept = 0, color = "grey60")
  p <- p + geom_point(data = dat9, aes(x = x, y = y, color = Season))
  p <- p + facet_grid(beta1 ~ beta2, labeller = label_parsed)
  p <- p + scale_y_continuous(expression(paste("Estimated ", phi)), limits = c(0, 1))
  p <- p + scale_x_continuous(expression(paste("True ", phi)), limits = c(0, 1))
  p <- p + geom_point(size = 3, alpha = 0.2, aes(color = Season))
  p <- p + stat_smooth(method = "lm", aes(linetype = Season, group = Season), color = "black", se = FALSE)
  p <- p + scale_color_manual(values = c('#cb4b16', '#e5ca28'))
  p <- p + geom_text(data = dat6, aes(x = x, y = y, label = label, group = label), parse = TRUE, size = 5)
  p
}
