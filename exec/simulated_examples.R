## illustrative simulated example of OAT with the VI model
library(gridExtra)
library(ggplot2)


set.seed(1308)
cov_par <- list("sigma" = 2, "l" = 1, "tau" = 1)
xy <- matrix(ncol = 1, data = sort(runif(n = 300, min = 0, max = 10)))
mu <- rep(0, times = nrow(xy))
cov_fun <- "sqexp"

Sigma <- make_cov_matC(x = xy, x_pred = matrix(), cov_fun = cov_fun, cov_par = cov_par, delta = 1e-6)

cov_par2 <- cov_par
cov_par2$tau <- 0
f <- as.numeric(mvtnorm::rmvnorm(n = 1, mean = mu, 
                                 sigma = make_cov_matC(x = xy, x_pred = matrix(), cov_par = cov_par2, cov_fun = "sqexp", delta = 1e-6)))
y <- f + rnorm(n = length(f), mean = 0, sd = cov_par$tau)

plot(xy, y)
points(x = xy, y = f, col = "red", type = "l")

dcov_fun_dtheta <- list("sigma" = dsqexp_dsigma, "l" = dsqexp_dl, "tau" = dsqexp_dtau)
dcov_fun_dknot <- dsqexp_dx2

## showing spikes of VI objective function
cov_par_start <- cov_par
set.seed(1308)
# fp_spike <- paste("TestCode_and_Results/vi_results/simulated_experiments/vi_spike_model.rds")
# system.time(spike_model <- optimize_gp(cov_par_start = cov_par_start,
#                                                 family = "gaussian",
#                                                 cov_fun = "sqexp",
#                                                 nugget = TRUE, 
#                                                 sparse = TRUE, 
#                                                 xu_opt = "fixed", 
#                                                 xu = c(1,3,5,7,9),
#                                                 xy = xy, vi = TRUE,
#                                                 y = y,
#                                                 mu = rep(0, times = length(y)),
#                                                 muu = rep(0, times = 5),
#                                                 opt = list("optim_method" = "adadelta", 
#                                                            "maxknot" = maxknot,
#                                                            "grad_tol" = Inf, 
#                                                            tol_knot = 1e-3, 
#                                                            obj_tol = 1e-3,
#                                                            "epsilon" = 1e-4,
#                                                            "TTmax" = maxknot, 
#                                                            maxit = 1000,
#                                                            "TTmin" = 10),
#                                                 verbose = TRUE, file_path = fp_spike)
# )
spike_model <- readRDS(file = "TestCode_and_Results/vi_results/simulated_experiments/vi_spike_model.rds")

xpred <- c(xy, seq(from = 0, to = 10, by = 0.01))
pred_spike <- predict_gp(mod = spike_model, 
                        x_pred = xpred, 
                        mu_pred = rep(0, times = length(xpred)), 
                        full_cov = FALSE, vi = TRUE)


plot_spike_pred <- pred_plot_fun(y = c(y, rep(NA, times = length(xpred) - length(y))), 
                                 xy = xpred, 
                              xu = spike_model$results$xu, 
                              xu_init = matrix(c(1,3,5,7,9), ncol = 1), 
                              pred_results = pred_spike$pred, 
                              ci_level = 0.95, 
                              family = "gaussian", 
                              alpha = 0.2, 
                              title = "", sparse = TRUE)
plot_spike_pred


## spike plot of the ELBO
elbo_vals <- numeric()
xu_seq <- seq(from = 0, to = 10, by = 0.01)

for(i in 1:length(xu_seq))
{
  elbo_vals[i] <- optimize_gp(cov_par_start = spike_model$results$cov_par,
                                                         family = "gaussian",
                                                         cov_fun = "sqexp",
                                                         nugget = TRUE, 
                                                         sparse = TRUE, 
                                                         xu_opt = "fixed", 
                                                         xu = c(1,3,5,7,9, xu_seq[i]),
                                                         xy = xy, vi = TRUE,
                                                         y = y,
                                                         mu = rep(0, times = length(y)),
                                                         muu = rep(0, times = 5 + 1),
                                                         opt = list("optim_method" = "adadelta", 
                                                                    # "maxknot" = 6,
                                                                    "grad_tol" = Inf, 
                                                                    tol_knot = 1e-3, 
                                                                    obj_tol = 1e-3,
                                                                    "epsilon" = 1e-4,
                                                                    # "TTmax" = 6, 
                                                                    maxit = 0,
                                                                    "TTmin" = 10),
                                                         verbose = FALSE, file_path = NULL)$results$obj_fun
  print(i)
}

hist(elbo_vals)

plot(xu_seq, elbo_vals, type = "l")

elbo_vals_plot <- ggplot() +
  geom_line(mapping = aes(x = xu_seq, y = elbo_vals)) +
  geom_hline(mapping = aes(yintercept = spike_model$results$obj_fun[length(spike_model$results$obj_fun)]), colour = "grey") +
  xlab("6th Knot") +
  ylab("ELBO") +
  theme_bw() 
  # theme(text = element_text(size = 12))

elbo_vals_plot


grid.arrange(plot_spike_pred, elbo_vals_plot, nrow = 2)

gA <- ggplotGrob(plot_spike_pred)
gB <- ggplotGrob(elbo_vals_plot)
maxWidth = grid::unit.pmax(gA$widths[2:5], gB$widths[2:5])
gA$widths[2:5] <- as.list(maxWidth)
gB$widths[2:5] <- as.list(maxWidth)

pdf("C://Users/Nate/Documents/Bayes_Research/GPs/GPrepo/Papers/GPpaper_knot_selection_vi/figures/vi_spike_plot.pdf", 
    width = 5, height = 3.5)
grid.arrange(gA, gB, ncol=1)
dev.off()

# ggsave(filename = "C://Users/Nate/Documents/Bayes_Research/GPs/GPrepo/Papers/GPpaper_knot_selection_vi/figures/vi_spike_plot.eps",
#        plot = arrangeGrob(gA, gB, ncol=1), device = , width = 5, height = 3.5)

cairo_ps(filename = "C://Users/Nate/Documents/Bayes_Research/GPs/GPrepo/Papers/GPpaper_knot_selection_vi/figures/vi_spike_plot.eps", 
         width = 5, height = 3.5, pointsize = 12)
grid.draw(arrangeGrob(gA, gB, ncol=1))
dev.off()

## try and see what happens when a knot is added to a bad spot
set.seed(1111)
bad_knot <- optimize_gp(cov_par_start = spike_model$results$cov_par,
                        family = "gaussian",
                        cov_fun = "sqexp",
                        nugget = TRUE, 
                        sparse = TRUE, 
                        xu_opt = "random", 
                        xu = c(1,3,5,7,9),
                        xy = xy, vi = TRUE,
                        y = y,
                        mu = rep(0, times = length(y)),
                        muu = rep(0, times = 5),
                        opt = list("optim_method" = "adadelta", 
                                   "maxknot" = 6,
                                   "grad_tol" = Inf, 
                                   tol_knot = 1e-3, 
                                   obj_tol = 1e-3,
                                   "epsilon" = 1e-4,
                                   "TTmax" = 1, 
                                   maxit = 1000,
                                   "TTmin" = 1),
                        verbose = TRUE, file_path = NULL)

bad_knot$results$xu

set.seed(1111)
good_knot <- optimize_gp(cov_par_start = spike_model$results$cov_par,
                        family = "gaussian",
                        cov_fun = "sqexp",
                        nugget = TRUE, 
                        sparse = TRUE, 
                        xu_opt = "random", 
                        xu = c(1,3,5,7,9),
                        xy = xy, vi = TRUE,
                        y = y,
                        mu = rep(0, times = length(y)),
                        muu = rep(0, times = 5),
                        opt = list("optim_method" = "adadelta", 
                                   "maxknot" = 6,
                                   "grad_tol" = Inf, 
                                   tol_knot = 1e-3, 
                                   obj_tol = 1e-3,
                                   "epsilon" = 1e-4,
                                   "TTmax" = 40, 
                                   maxit = 1000,
                                   "TTmin" = 10),
                        verbose = TRUE, file_path = NULL)

good_knot$results$xu

## fit the models for the first figure
maxknot <- 30
k_init <- 5
set.seed(1308)
xu_rand1 <- runif(n = k_init, min = 0, max = 10)
set.seed(1309)
xu_rand2 <- runif(n = k_init, min = 0, max = 10)
set.seed(1310)
xu_rand3 <- runif(n = k_init, min = 0, max = 10)

muu <- rep(0, times = length(xu_rand1))

## OAT VI 1
cov_par_start <- cov_par
set.seed(1308)
fp_oat_vi1 <- paste("TestCode_and_Results/vi_results/simulated_experiments/gaussian_sim_xu5_rand_oat1.rda")
# t_oat_vi1 <- system.time(oat_vi1 <- optimize_gp(cov_par_start = cov_par_start,
#                                                    family = "gaussian",
#                                                    cov_fun = "sqexp",
#                                                    nugget = TRUE, 
#                                                    sparse = TRUE, 
#                                                    xu_opt = "oat", 
#                                                    xu = xu_rand1,
#                                                    xy = xy, vi = TRUE,
#                                                    y = y,
#                                                    mu = mu,
#                                                    muu = muu,
#                                                    opt = list("optim_method" = "adadelta", 
#                                                               "maxknot" = maxknot,
#                                                               "grad_tol" = Inf, 
#                                                               tol_knot = 1e-3, 
#                                                               obj_tol = 1e-3,
#                                                               "epsilon" = 1e-4,
#                                                               "TTmax" = maxknot, 
#                                                               "TTmin" = 10),
#                                                    verbose = TRUE, file_path = fp_oat_vi1)
# )
oat_vi1 <- readRDS(fp_oat_vi1)


## OAT VI 2
cov_par_start <- cov_par
set.seed(1308)
fp_oat_vi2 <- paste("TestCode_and_Results/vi_results/simulated_experiments/gaussian_sim_xu5_rand_oat2.rda")
# t_oat_vi2 <- system.time(oat_vi2 <- optimize_gp(cov_par_start = cov_par_start,
#                                                 family = "gaussian",
#                                                 cov_fun = "sqexp",
#                                                 nugget = TRUE, 
#                                                 sparse = TRUE, 
#                                                 xu_opt = "oat", 
#                                                 xu = xu_rand2,
#                                                 xy = xy, vi = TRUE,
#                                                 y = y,
#                                                 mu = mu,
#                                                 muu = muu,
#                                                 opt = list("optim_method" = "adadelta", 
#                                                            "maxknot" = maxknot,
#                                                            "grad_tol" = Inf, 
#                                                            tol_knot = 1e-3, 
#                                                            obj_tol = 1e-3,
#                                                            "epsilon" = 1e-4,
#                                                            "TTmax" = maxknot, 
#                                                            "TTmin" = 10),
#                                                 verbose = TRUE, file_path = fp_oat_vi2)
# )
oat_vi2 <- readRDS(fp_oat_vi2)

## OAT VI 3
cov_par_start <- cov_par
set.seed(1308)
fp_oat_vi3 <- paste("TestCode_and_Results/vi_results/simulated_experiments/gaussian_sim_xu5_rand_oat3.rda")
# t_oat_vi3 <- system.time(oat_vi3 <- optimize_gp(cov_par_start = cov_par_start,
#                                                 family = "gaussian",
#                                                 cov_fun = "sqexp",
#                                                 nugget = TRUE, 
#                                                 sparse = TRUE, 
#                                                 xu_opt = "oat", 
#                                                 xu = xu_rand3,
#                                                 xy = xy, vi = TRUE,
#                                                 y = y,
#                                                 mu = mu,
#                                                 muu = muu,
#                                                 opt = list("optim_method" = "adadelta", 
#                                                            "maxknot" = maxknot,
#                                                            "grad_tol" = Inf, 
#                                                            tol_knot = 1e-3, 
#                                                            obj_tol = 1e-3,
#                                                            "TTmax" = maxknot, 
#                                                            "TTmin" = 10),
#                                                 verbose = TRUE, file_path = fp_oat_vi3)
# )
oat_vi3 <- readRDS(fp_oat_vi3)

## Simultaneous optimizations using the knot locations found by OAT
xu_simult1 <- oat_vi1$results$xu
fp_simult_vi1 <- paste("TestCode_and_Results/vi_results/simulated_experiments/gaussian_sim_xu5_rand_simult1.rda")
# t_simult_vi1 <- system.time(simult_vi1 <- optimize_gp(y = y, xy = xy, 
#                                                         xu = xu_simult1, 
#                                                         sparse = TRUE,
#                                                         cov_fun = "sqexp", 
#                                                         cov_par_start = oat_vi1$results$cov_par, 
#                                                         mu = mu, muu = muu, vi = TRUE,
#                                                         family = "gaussian", 
#                                                         xu_opt = "simultaneous", nugget = TRUE, 
#                                                         opt = list("optim_method" = "adadelta",
#                                                                    "grad_tol" = Inf,
#                                                                    tol_knot = 1e-3, 
#                                                                    "epsilon" = 1e-4,
#                                                                    obj_tol = 1e-3), 
#                                                         verbose = TRUE, file_path = fp_simult_vi1))
simult_vi1 <- readRDS(fp_simult_vi1)

## Simultaneous optimizations using the knot locations found by OAT
xu_simult2 <- oat_vi2$results$xu
fp_simult_vi2 <- paste("TestCode_and_Results/vi_results/simulated_experiments/gaussian_sim_xu5_rand_simult2.rda")
# t_simult_vi2 <- system.time(simult_vi2 <- optimize_gp(y = y, xy = xy, 
#                                                       xu = xu_simult2, 
#                                                       sparse = TRUE,
#                                                       cov_fun = "sqexp", 
#                                                       cov_par_start = oat_vi2$results$cov_par, 
#                                                       mu = mu, muu = muu, vi = TRUE,
#                                                       family = "gaussian", 
#                                                       xu_opt = "simultaneous", nugget = TRUE, 
#                                                       opt = list("optim_method" = "adadelta",
#                                                                  "grad_tol" = Inf,
#                                                                  tol_knot = 1e-3,
#                                                                  "epsilon" = 1e-4,
#                                                                  obj_tol = 1e-3), 
#                                                       verbose = TRUE, file_path = fp_simult_vi2))
simult_vi2 <- readRDS(fp_simult_vi2)


## Simultaneous optimizations using the knot locations found by OAT
xu_simult3 <- oat_vi3$results$xu
fp_simult_vi3 <- paste("TestCode_and_Results/vi_results/simulated_experiments/gaussian_sim_xu5_rand_simult3.rda")
# t_simult_vi3 <- system.time(simult_vi3 <- optimize_gp(y = y, xy = xy, 
#                                                       xu = xu_simult3, 
#                                                       sparse = TRUE,
#                                                       cov_fun = "sqexp", 
#                                                       cov_par_start = oat_vi3$results$cov_par, 
#                                                       mu = mu, muu = muu, vi = TRUE,
#                                                       family = "gaussian", 
#                                                       xu_opt = "simultaneous", nugget = TRUE, 
#                                                       opt = list("optim_method" = "adadelta",
#                                                                  "grad_tol" = Inf,
#                                                                  tol_knot = 1e-3, 
#                                                                  "epsilon" = 1e-4,
#                                                                  obj_tol = 1e-3), 
#                                                       verbose = TRUE, file_path = fp_simult_vi3))
simult_vi3 <- readRDS(fp_simult_vi3)

## get predictions
pred_oat1 <- predict_gp(mod = oat_vi1, 
                        x_pred = xy, 
                        mu_pred = mu, 
                        full_cov = FALSE, vi = TRUE)
pred_oat2 <- predict_gp(mod = oat_vi2, 
                        x_pred = xy, 
                        mu_pred = mu, 
                        full_cov = FALSE, vi = TRUE)
pred_oat3 <- predict_gp(mod = oat_vi3, 
                        x_pred = xy, 
                        mu_pred = mu, 
                        full_cov = FALSE, vi = TRUE)

pred_simult1 <- predict_gp(mod = simult_vi1, 
                        x_pred = xy, 
                        mu_pred = mu, 
                        full_cov = FALSE, vi = TRUE)
pred_simult2 <- predict_gp(mod = simult_vi2, 
                        x_pred = xy, 
                        mu_pred = mu, 
                        full_cov = FALSE, vi = TRUE)
pred_simult3 <- predict_gp(mod = simult_vi3, 
                        x_pred = xy, 
                        mu_pred = mu, 
                        full_cov = FALSE, vi = TRUE)

## plots of predictions
plot_oat1 <- pred_plot_fun(y = y, xy = xy, 
                           xu = oat_vi1$results$xu, 
                           xu_init = xu_rand1, 
                           pred_results = pred_oat1$pred, 
                           ci_level = 0.95, 
                           family = "gaussian", 
                           alpha = 0.2, 
                           title = "", sparse = TRUE)

plot_oat2 <- pred_plot_fun(y = y, xy = xy, 
                           xu = oat_vi2$results$xu, 
                           xu_init = xu_rand2, 
                           pred_results = pred_oat2$pred, 
                           ci_level = 0.95, 
                           family = "gaussian", 
                           alpha = 0.2, 
                           title = "", sparse = TRUE)

plot_oat3 <- pred_plot_fun(y = y, xy = xy, 
                           xu = oat_vi3$results$xu, 
                           xu_init = xu_rand3, 
                           pred_results = pred_oat3$pred, 
                           ci_level = 0.95, 
                           family = "gaussian", 
                           alpha = 0.2, 
                           title = "", sparse = TRUE)

plot_simult1 <- pred_plot_fun(y = y, xy = xy, 
                           xu = simult_vi1$results$xu, 
                           xu_init = xu_simult1, 
                           pred_results = pred_simult1$pred, 
                           ci_level = 0.95, 
                           family = "gaussian", 
                           alpha = 0.2, 
                           title = "", sparse = TRUE)

plot_simult2 <- pred_plot_fun(y = y, xy = xy, 
                           xu = simult_vi2$results$xu, 
                           xu_init = xu_simult2, 
                           pred_results = pred_simult2$pred, 
                           ci_level = 0.95, 
                           family = "gaussian", 
                           alpha = 0.2, 
                           title = "", sparse = TRUE)

plot_simult3 <- pred_plot_fun(y = y, xy = xy, 
                           xu = simult_vi3$results$xu, 
                           xu_init = xu_simult3, 
                           pred_results = pred_simult3$pred, 
                           ci_level = 0.95, 
                           family = "gaussian", 
                           alpha = 0.2, 
                           title = "", sparse = TRUE)

grid.arrange(plot_oat1, plot_oat2, plot_oat3,
             plot_simult1, plot_simult2, plot_simult3, nrow = 2)

pdf(file = "C://Users/Nate/Documents/Bayes_Research/GPs/GPrepo/Code/TestCode_and_Results/vi_results/simulated_experiments/oat_vs_simult_gaussian_example.pdf")
results_plot <- grid.arrange(plot_oat1, plot_oat2, plot_oat3,
                             plot_simult1, plot_simult2, plot_simult3, nrow = 2)
dev.off()

pdf("C://Users/Nate/Documents/Bayes_Research/GPs/GPrepo/Papers/GPpaper_knot_selection_vi/figures/oat_vs_simult_gaussian_example.pdf", width = 8, height = 4)
results_plot <- grid.arrange(plot_oat1, plot_oat2, plot_oat3,
                             plot_simult1, plot_simult2, plot_simult3, nrow = 2)
dev.off()

cairo_ps(filename = "C://Users/Nate/Documents/Bayes_Research/GPs/GPrepo/Papers/GPpaper_knot_selection_vi/figures/oat_vs_simult_gaussian_example.eps", 
         width = 8, height = 4, pointsize = 12)
results_plot <- grid.arrange(plot_oat1, plot_oat2, plot_oat3,
                             plot_simult1, plot_simult2, plot_simult3, nrow = 2)
dev.off()
