# The Airfoil data is available for download at
# https://archive.ics.uci.edu/ml/datasets/Airfoil+Self-Noise

## NOTE: Results should closely, but not exactly match those in the paper.
##  This is due to slight modifications of optimization settings for
##  different models. These changes were either to (a) try to speed up
##  optimization slightly in sparse models or (b) because delta = 1e-6
##  results in numerical problems for some models and thus was changed to
##  delta = 1e-4. These changes should be very minimal.

## read in raw Airfoil data
# airfoil_raw <- read.table(file = "Data/airfoil_self_noise.dat", header = FALSE)
airfoil <- airfoil_raw

colnames(airfoil) <- c(paste("X",1:5, sep = ""), "y")

## function to center and scale the data
my_scale <- function(x, mean_vec, sd_vec)
{
  mean_mat <- matrix(rep(x = mean_vec, times = nrow(x)),
                     ncol = ncol(x), nrow = nrow(x), byrow = TRUE)
  sd_mat <- matrix(rep(x = sd_vec, times = nrow(x)),
                   ncol = ncol(x), nrow = nrow(x), byrow = TRUE)

  return((x - mean_mat) / sd_mat)
}

ptrain <- 0.8
runs <- 5
cnames <- c("full", "oat-bo-vi", "oat-rs-vi", "oat-bo", "simult-vi", "adjust-vi")

xu_opt <- c("NA", "oat", "random", "oat", "simultaneous", "simultaneous")
vi <- c(FALSE, TRUE, TRUE, FALSE, TRUE, TRUE)
sparse <- c(FALSE, TRUE, TRUE, TRUE, TRUE, TRUE)
adjust <- c(FALSE, FALSE, FALSE, FALSE, FALSE, TRUE)
train_time <- numeric()
K <- numeric()
nlp <- numeric()
kl <- numeric()
srmse <- numeric()

## optimization parameters
maxit = 1000
obj_tol = 1e-3
grad_tol = Inf
delta = 1e-4
epsilon <- 1e-4
k_init <- 5
maxknot <- 80
TTmax <- maxknot / 2
TTmin <- 5

## directory to which to save models
# dir <- "./"

results_table <- data.frame(xu_opt = "NA", vi = FALSE, sparse = FALSE, adjust = FALSE, run = 1,
                            "MNLP" = 0,
                            "AUKL" = 0,
                            "SRMSE" = 0,
                            train_time = 1, K = 1)
counter <- 0
for(i in 1:runs)
{
  set.seed(1307 + i)
  train_rows <- sample.int(n = nrow(airfoil), size = round(ptrain * nrow(airfoil)), replace = FALSE)
  airfoil$train <- (1:nrow(airfoil)) %in% train_rows

  X_train <- as.matrix(airfoil[airfoil$train == TRUE,1:5])
  train_mean <- apply(X_train, MARGIN = 2, FUN = mean)
  train_sd <- apply(X_train, MARGIN = 2, FUN = sd)


  X_train_sc <- my_scale(x = X_train, mean_vec = train_mean, sd_vec = train_sd)
  y_train <- airfoil[airfoil$train == TRUE,]$y

  # y_train_sc <- (y_train - mean(y_train)) / sd(y_train)

  X_test <- as.matrix(airfoil[airfoil$train == FALSE,1:5])
  X_test_sc <- my_scale(x = X_test, mean_vec = train_mean, sd_vec = train_sd)
  y_test <- airfoil[airfoil$train == FALSE,]$y

  for(j in 1:length(xu_opt))
  {
    counter <- counter + 1
    d <- ncol(X_train_sc)
    cov_par_start <- list("sigma" = sqrt(var(y_train) / 2), "l" = 1, "tau" = sqrt(var(y_train) / 2))

    set.seed(1308)
    xu_init <- kmeans(x = X_train_sc, centers = k_init)$centers
    mu <- rep(mean(y_train), times = length(y_train))
    muu_init <- rep(mean(y_train), times = k_init)

    if(xu_opt[j] == "simultaneous" & adjust[j] == FALSE)
    {
      # m_vi_oat <- readRDS(file = paste(dir, "_xu_opt", "oat",
      #                                  "_vi", TRUE,
      #                                  "_sparse", TRUE,
      #                                  "_adjust", FALSE,
      #                                  "_", i,
      #                                  ".rds", sep = ""))
      xu_init <- kmeans(x = X_train_sc, centers = nrow(m_vi_oat$results$xu))$centers
      mu <- rep(mean(y_train), times = length(y_train))
      muu_init <- rep(mean(y_train), times = nrow(m_vi_oat$results$xu))
    }

    if(xu_opt[j] == "simultaneous" & adjust[j] == TRUE)
    {
      # m_vi_oat <- readRDS(file = paste(dir, "_xu_opt", "oat",
      #                                  "_vi", TRUE,
      #                                  "_sparse", TRUE,
      #                                  "_adjust", FALSE,
      #                                  "_", i,
      #                                  ".rds", sep = ""))
      xu_init <- m_vi_oat$results$xu
      mu <- rep(mean(y_train), times = length(y_train))
      muu_init <- rep(mean(y_train), times = nrow(m_vi_oat$results$xu))
    }


    # filename <- paste(dir, "_xu_opt", xu_opt[j],
    #                   "_vi", vi[j],
    #                   "_sparse", sparse[j],
    #                   "_adjust", adjust[j],
    #                   "_", i,
    #                   ".rds", sep = "")

    temp_time <- system.time(m <- try(optimize_gp(y = y_train,
                                                  xy = X_train_sc,
                                                  cov_fun = "sqexp",
                                                  cov_par_start = cov_par_start,
                                                  mu = mu,
                                                  family = "gaussian",
                                                  nugget = TRUE, sparse = sparse[j],
                                                  xu_opt = xu_opt[j],
                                                  xu = xu_init,
                                                  muu = muu_init, vi = vi[j],
                                                  opt = list(maxit = maxit,
                                                             obj_tol = obj_tol,
                                                             grad_tol = grad_tol,
                                                             delta = delta,
                                                             epsilon = epsilon,
                                                             TTmin = TTmin,
                                                             TTmax = TTmax,
                                                             maxknot = maxknot),
                                                  verbose = FALSE,
                                                  file_path = NULL)))

    if(class(m) == "try-error")
    {
      print(paste("run = ", i, " / mod = ", j, sep = ""))
      next
    }
    if(j == 1)
    {
      m_full <- m
    }
    if(j == 2)
    {
      m_vi_oat <- m
    }

    preds_full <- predict_gp(mod = m_full,
                             x_pred = X_test_sc,
                             mu_pred = rep(mean(y_train), times = length(y_test)),
                             full_cov = FALSE, vi = FALSE)

    preds <- predict_gp(mod = m,
                        x_pred = X_test_sc,
                        mu_pred = rep(mean(y_train), times = length(y_test)),
                        full_cov = FALSE, vi = vi[j])

    ## calculate SRMSE, MNLP, and AUKL from full GP
    nlp[counter] <- median(my_nlp(y = y_test, pred_mean = preds$pred$pred_mean,
                                  pred_var = preds$pred$pred_var, family = "gaussian",
                                  par = list("tau" = m$results$cov_par$tau),
                                  mc = FALSE, mv = FALSE)$nlp)

    kl[counter] <- my_kl(mean1 = preds_full$pred$pred_mean,
                         mean2 = preds$pred$pred_mean,
                         sigma1 = preds_full$pred$pred_var,
                         sigma2 = preds$pred$pred_var, mv = FALSE)

    srmse[counter] <- my_rmse(pred = preds$pred$pred_mean, actual = y_test) / sd(y_test)

    print(paste("run = ", i, " / mod = ", j, sep = ""))

    ## record times and knots
    temp_df <- data.frame("xu_opt" = xu_opt[j],
                          "vi" = vi[j],
                          "sparse" = sparse[j],
                          "adjust" = adjust[j],
                          "train_time" = temp_time[1],
                          "MNLP" = nlp[counter],
                          "AUKL" = kl[counter],
                          "SRMSE" = srmse[counter],
                          "K" = ifelse(is.null(nrow(m$results$xu)),
                                       yes = "NA",
                                       no = nrow(m$results$xu)),
                          "run" = i)
    results_table <- rbind(results_table,
                           temp_df)
  }
}

## print the results table
results <- results_table[-1,]

## create the Boston results figure from the paper
## create model index
results$xu_opt <- as.character(results$xu_opt)
results$xu_opt[(results$xu_opt) == "NA"] <- "-"
results$model_number <- 1*(results$xu_opt == "-") +
  2 * 1 * (results$xu_opt == "oat" & results$vi == TRUE) +
  3 * 1 * (results$xu_opt == "random" & results$vi == TRUE) +
  4 * 1 * (results$xu_opt == "oat" & results$vi == FALSE) +
  5 * 1 * (results$xu_opt == "simultaneous" & results$adjust == FALSE) +
  6 * 1 * (results$xu_opt == "simultaneous" & results$adjust == TRUE)

results$K <- as.numeric(results$K)
results$Model <- ifelse(test = results$model_number == 1, yes = "FGP",
                        no = ifelse(test = results$model_number == 2, yes = "OBVk",
                                    no = ifelse(test = results$model_number == 3, yes = "ORVk",
                                                no = ifelse(test = results$model_number == 4, yes = "OBFk",
                                                            no = ifelse(test = results$model_number == 5, yes = "SVk",
                                                                        no = "SVO")))))

results


## plots of results
## put runs on x-axis and connect horizontally with lines
colnames(results)[9] <- "Train Time"
colnames(results)[10] <- "# Knots"
# head(results)
results$`Knot Opt.` <- ifelse(test = results$xu_opt == "-", yes = "N/A",
                              no = ifelse(test = results$xu_opt %in% c("oat", "random"), yes = "OAT",
                                          no = "Simult."))
results_long <- reshape2::melt(data = results, id.vars = c("Model", "run", "model_number", "Knot Opt."),
                               measure.vars = c("MNLP","AUKL","SRMSE","Train Time", "# Knots"), variable.name = "Metric")
results_long$Model <- as.character(results_long$Model)

## fix the fact that I did not average AUKL before
for(i in 1:nrow(results_long))
{
  if(results_long$Metric[i] == "AUKL")
  {
    results_long$value[i] <- results_long$value[i] / 301
  }
}

for(i in 1:nrow(results_long))
{
  if(results_long$Metric[i] == "AUKL")
  {
    results_long$value[i] <- log(x = results_long$value[i], base = 10)
  }
}

for(i in 1:nrow(results_long))
{
  if(results_long$Metric[i] == "Train Time")
  {
    results_long$value[i] <- log(x = results_long$value[i], base = 10)
  }
}

levels(results_long$Metric) <- c("MNLP", "log(AUKL)","SRMSE","log(Train Time)", "# Knots")

results_long <- results_long[is.finite(results_long$value),]

cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
# The palette with black:
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#0072B2", "#D55E00", "#CC79A7")

results_plots <- ggplot(data = results_long) +
  geom_point(mapping = aes(x = run, y = value, colour = Model, shape = Model), size = 2) +
  facet_grid(Metric ~ ., scales = "free_y") +
  theme_bw() +
  # theme(text = element_text(size = 14)) +
  geom_line(mapping = aes(x = run, y = value, colour = Model, linetype = `Knot Opt.`), size = 0.5) +
  # scale_y_continuous(, trans = "log10")
  scale_colour_manual(values = cbbPalette) +
  scale_shape_manual(values = c(1,2,3,4,7,8)) +
  scale_linetype_manual(values = c(1:3))

results_plots

