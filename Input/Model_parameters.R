# Author: Sandor Piross - sandor.piross@gmail.com 

# The script estimates empirical and defines non-estimated parameters
# for the model

rm(list = ls())
#.rs.restartR()

library(fitdistrplus)
library(tidyverse)
library(readxl)


# Importing empirical data
parameter_estimation_edges <-
  readRDS("./Input/Empirical_data/empirical_edges.rds") |> 
  filter(dataset=="Estimation")

parameter_estimation_nodes <-
  readRDS("./Input/Empirical_data/empirical_nodes.rds") |> 
  filter(dataset=="Estimation")

# Data frame for estimated parameters -------------------------------------

estimated_parameters <- data.frame(
  parameter_name = NA,
  # Parameter for the simulation model
  distribution = NA,
  # Fitted distribution
  distr_param = NA,
  # Distribution parameter
  distr_param_est = NA,
  # Distribution parameter
  distr_param_est_se = NA # Distribution parameter standard error
)


# NETWORK STRUCTURE ------------------------------------------------------
# Number of nests per colony ----------------------------------------------
nests_per_colony <-
  parameter_estimation_nodes |>
  filter(species == "nest") |>
  group_by(colony, year) |>
  summarise(N = n()) |>
  dplyr::select(N) |>
  pull()

# Geometric distribution fit
nests_per_colony_fit_geom <-
  fitdist(nests_per_colony, distr = "geom", method = "mle")
plot(nests_per_colony_fit_geom)
summary(nests_per_colony_fit_geom)

# Poisson distribution fit
nests_per_colony_fit_pois <-
  fitdist(nests_per_colony, distr = "pois", method = "mle")
plot(nests_per_colony_fit_pois)
summary(nests_per_colony_fit_pois)

# Negative binomial distribution fit
nests_per_colony_fit_nb <-
  fitdist(nests_per_colony, distr = "nbinom", method = "mle")
plot(nests_per_colony_fit_nb)
summary(nests_per_colony_fit_nb)

# Saving results
estimated_parameters <-
  rbind(
    estimated_parameters,
    data.frame(
      # Parameter for the simulation model
      parameter_name = "Nests per colony",
      # Fitted distribution
      distribution = "nbinom",
      # Distribution parameter
      distr_param = "mu",
      # Distribution parameter
      distr_param_est = nests_per_colony_fit_nb$estimate["mu"],
      # Distribution parameter standard error
      distr_param_est_se = nests_per_colony_fit_nb$sd["mu"]
    )
  )
estimated_parameters <- estimated_parameters[-1, ]

estimated_parameters <-
  rbind(
    estimated_parameters,
    data.frame(
      # Parameter for the simulation model
      parameter_name = "Nests per colony",
      # Fitted distribution
      distribution = "nbinom",
      # Distribution parameter
      distr_param = "size",
      # Distribution parameter
      distr_param_est = nests_per_colony_fit_nb$estimate["size"],
      # Distribution parameter standard error
      distr_param_est_se = nests_per_colony_fit_nb$sd["size"]
    )
  )


# Internest trail length --------------------------------------------------
internest_length <-
  parameter_estimation_edges |>
  filter(trail_type == "internest") |>
  dplyr::select(distance) |>
  pull()

plot(density(internest_length))

# Gamma distribution fit
internest_length_fit_gamma <- fitdist(internest_length,
                                      distr = "gamma", method = "mle")
plot(internest_length_fit_gamma)
summary(internest_length_fit_gamma)

# Log-normal distribution fit
internest_length_fit_lnorm <- fitdist(internest_length,
                                      distr = "lnorm", method = "mle")
plot(internest_length_fit_lnorm)
summary(internest_length_fit_lnorm)

# Exponential distribution fit
internest_length_fit_exp <- fitdist(internest_length,
                                    distr = "exp", method = "mle")
plot(internest_length_fit_exp)
summary(internest_length_fit_exp)

# Saving results
estimated_parameters <-
  rbind(
    estimated_parameters,
    data.frame(
      # Parameter for the simulation model
      parameter_name = "Internest trail length",
      # Fitted distribution
      distribution = "gamma",
      # Distribution parameter
      distr_param = "shape",
      # Distribution parameter
      distr_param_est = internest_length_fit_gamma$estimate["shape"],
      # Distribution parameter standard error
      distr_param_est_se = internest_length_fit_gamma$sd["shape"]
    )
  )

estimated_parameters <-
  rbind(
    estimated_parameters,
    data.frame(
      # Parameter for the simulation model
      parameter_name = "Internest trail length",
      # Fitted distribution
      distribution = "gamma",
      # Distribution parameter
      distr_param = "rate",
      # Distribution parameter
      distr_param_est = internest_length_fit_gamma$estimate["rate"],
      # Distribution parameter standard error
      distr_param_est_se = internest_length_fit_gamma$sd["rate"]
    )
  )

# Foraging trail length --------------------------------------------------
foraging_length <-
  parameter_estimation_edges |>
  filter(trail_type == "foraging") |>
  dplyr::select(distance) |>
  pull()

plot(density(foraging_length))

# Gamma distribution fit
foraging_length_fit_gamma <- fitdist(foraging_length,
                                     distr = "gamma", method = "mle")
plot(foraging_length_fit_gamma)
summary(foraging_length_fit_gamma)

# Log-normal distribution fit
foraging_length_fit_lnorm <- fitdist(foraging_length,
                                     distr = "lnorm", method = "mle")
plot(foraging_length_fit_lnorm)
summary(foraging_length_fit_lnorm)

# Exponential distribution fit
foraging_length_fit_exp <- fitdist(foraging_length,
                                   distr = "exp", method = "mle")
plot(foraging_length_fit_exp)
summary(foraging_length_fit_exp)

# Saving results
estimated_parameters <-
  rbind(
    estimated_parameters,
    data.frame(
      # Parameter for the simulation model
      parameter_name = "Foraging trail length",
      # Fitted distribution
      distribution = "gamma",
      # Distribution parameter
      distr_param = "shape",
      # Distribution parameter
      distr_param_est = foraging_length_fit_gamma$estimate["shape"],
      # Distribution parameter standard error
      distr_param_est_se = foraging_length_fit_gamma$sd["shape"]
    )
  )

estimated_parameters <-
  rbind(
    estimated_parameters,
    data.frame(
      # Parameter for the simulation model
      parameter_name = "Foraging trail length",
      # Fitted distribution
      distribution = "gamma",
      # Distribution parameter
      distr_param = "rate",
      # Distribution parameter
      distr_param_est = foraging_length_fit_gamma$estimate["rate"],
      # Distribution parameter standard error
      distr_param_est_se = foraging_length_fit_gamma$sd["rate"]
    )
  )

# Nest sizes --------------------------------------------------
nest_size <-
  parameter_estimation_nodes |>
  filter(species == "nest") |>
  dplyr::select(nest.size) |>
  pull()

plot(density(nest_size))
hist(nest_size, breaks = 30)
summary(nest_size)
#nest_size <- nest_size[nest_size<quantile(nest_size,0.975)]

# Gamma distribution fit
nest_size_fit_gamma <- fitdist(
  nest_size,
  distr = "gamma",
  method = "mle",
  #lower = c(0, 0),
  start = list(rate = 0.1, shape = 0.5)
)
plot(nest_size_fit_gamma)
summary(nest_size_fit_gamma)

# Log-normal distribution fit
nest_size_fit_lnorm <- fitdist(nest_size,
                               distr = "lnorm", method = "mle")
plot(nest_size_fit_lnorm)
summary(nest_size_fit_lnorm)


# Exponential distribution fit
nest_size_fit_exp <- fitdist(nest_size,
                             distr = "exp", method = "mle")

plot(nest_size_fit_exp)
summary(nest_size_fit_exp)

estimated_parameters <-
  rbind(
    estimated_parameters,
    data.frame(
      # Parameter for the simulation model
      parameter_name = "Nest size",
      # Fitted distribution
      distribution = "gamma",
      # Distribution parameter
      distr_param = "shape",
      # Distribution parameter
      distr_param_est = nest_size_fit_gamma$estimate["shape"],
      # Distribution parameter standard error
      distr_param_est_se = nest_size_fit_gamma$sd["shape"]
    )
  )

estimated_parameters <-
  rbind(
    estimated_parameters,
    data.frame(
      # Parameter for the simulation model
      parameter_name = "Nest size",
      # Fitted distribution
      distribution = "gamma",
      # Distribution parameter
      distr_param = "rate",
      # Distribution parameter
      distr_param_est = nest_size_fit_gamma$estimate["rate"],
      # Distribution parameter standard error
      distr_param_est_se = nest_size_fit_gamma$sd["rate"],
      row.names = ""
    )
  )


# Number of foraging trails per nest --------------------------------------

foraging_trail_num <-
  parameter_estimation_edges |>
  group_by(date, from_ID, trail_type) |>
  summarise(N = n()) |>
  pivot_wider(names_from = "trail_type",
              values_from = "N") |>
  rename("foraging_num" = "foraging",
         "internest_num" = "internest") |>
  mutate(date_from_ID = paste(date, from_ID, sep = "_")) |>
  ungroup() |>
  dplyr::select(date_from_ID, foraging_num, internest_num) |>
  right_join(
    parameter_estimation_nodes |>
      mutate(date_node_ID = paste(date, node_ID, sep = "_")),
    by = c("date_from_ID" = "date_node_ID")
  ) |>
  mutate(
    foraging_num = replace_na(foraging_num, 0),
    internest_num = replace_na(internest_num, 0)
  ) |>
  filter(species == "nest") |>
  dplyr::select(foraging_num) |>
  pull()

table(foraging_trail_num)

# Poisson distribution fit
foraging_trail_num_fit_pois <- fitdist(foraging_trail_num,
                                       distr = "pois", method = "mle")
plot(foraging_trail_num_fit_pois)
summary(foraging_trail_num_fit_pois)


# Negative binomial distribution fit
foraging_trail_num_fit_nbinom <- fitdist(foraging_trail_num,
                                         distr = "nbinom", method = "mle")
plot(foraging_trail_num_fit_nbinom)
summary(foraging_trail_num_fit_nbinom)

# Geometric distribution
foraging_trail_num_fit_geom <- fitdist(foraging_trail_num,
                                       distr = "geom", method = "mle")
plot(foraging_trail_num_fit_geom)
summary(foraging_trail_num_fit_geom)

estimated_parameters <-
  rbind(
    estimated_parameters,
    data.frame(
      # Parameter for the simulation model
      parameter_name = "Number of foraging trails per nest",
      # Fitted distribution
      distribution = "pois",
      # Distribution parameter
      distr_param = "lambda",
      # Distribution parameter
      distr_param_est = foraging_trail_num_fit_pois$estimate["lambda"],
      # Distribution parameter standard error
      distr_param_est_se = foraging_trail_num_fit_pois$sd["lambda"]
    )
  )

# Number of trails per nest --------------------------------------

trail_num <-
  parameter_estimation_edges |>
  group_by(date, from_ID, trail_type) |>
  summarise(N = n()) |>
  pivot_wider(names_from = "trail_type",
              values_from = "N") |>
  rename("foraging_num" = "foraging",
         "internest_num" = "internest") |>
  mutate(date_from_ID = paste(date, from_ID, sep = "_")) |>
  ungroup() |>
  dplyr::select(date_from_ID, foraging_num, internest_num) |>
  right_join(
    parameter_estimation_nodes |>
      mutate(date_node_ID = paste(date, node_ID, sep = "_")),
    by = c("date_from_ID" = "date_node_ID")
  ) |>
  mutate(
    foraging_num = replace_na(foraging_num, 0),
    internest_num = replace_na(internest_num, 0)
  ) |>
  filter(species == "nest") |>
  dplyr::select(foraging_num, internest_num) |>
  mutate(trail_num = foraging_num + internest_num) |>
  pull()

table(trail_num)

# Poisson distribution fit
trail_num_fit_pois <- fitdist(trail_num,
                              distr = "pois", method = "mle")
plot(trail_num_fit_pois)
summary(trail_num_fit_pois)


# Negative binomial distribution fit
trail_num_fit_nbinom <- fitdist(trail_num,
                                distr = "nbinom", method = "mle")
plot(trail_num_fit_nbinom)
summary(trail_num_fit_nbinom)

# Geometric distribution
trail_num_fit_geom <- fitdist(trail_num,
                              distr = "geom", method = "mle")
plot(trail_num_fit_geom)
summary(trail_num_fit_geom)

estimated_parameters <-
  rbind(
    estimated_parameters,
    data.frame(
      # Parameter for the simulation model
      parameter_name = "Number of trails per nest",
      # Fitted distribution
      distribution = "pois",
      # Distribution parameter
      distr_param = "lambda",
      # Distribution parameter
      distr_param_est = trail_num_fit_pois$estimate["lambda"],
      # Distribution parameter standard error
      distr_param_est_se = trail_num_fit_pois$sd["lambda"]
    )
  )

# NETWORK DYNAMICS --------------------------------------------------------


# Gravity model parameters - Power distance decay  ----------------------------
gravity_pow_error <- function(p, data) {
  G <-  p[1] # "Gravitational" constant
  alpha <-  p[2] # Nest size exponent
  beta <- p[3] # Distance multiplier
  S <- data$strength # Observed strength
  N1 <- data$from_nest_size # From nest size
  N2 <- data$to_nest_size # To nest size (1 for trees)
  d <- data$distance # Distance
  
  # Exponential distance decay
  #MSE
  MSE <- mean((S - (G * (((N1 * N2) ^ (alpha)
  ) /
    (
      d ^ beta
    )))) ^ 2)
  
  return(MSE)
  
}

gravity_pow_par_foraging <-
  optim(
    par = c(0.0001, 0.0001, 0.0001),
    fn = gravity_pow_error,
    data = parameter_estimation_edges |> filter(trail_type == "foraging")
  )

gravity_pow_par_internest <-
  optim(
    par = c(0.0001, 0.0001, 0.0001),
    fn = gravity_pow_error,
    data = parameter_estimation_edges |> filter(trail_type == "internest")
  )

gravity_pow_parameters <-
  tibble(trail_type = c("foraging", "internest")) |>
  bind_cols(rbind(gravity_pow_par_foraging$par,
                  gravity_pow_par_internest$par)) |>
  rename(G_pow = ...2,
         alpha_pow = ...3,
         beta_pow = ...4)



# NOT ESTIMATED PARAMETERS ------------------------------------------------

# K - Nest carrying capacity
K <-
  qgamma(0.999,
         shape = estimated_parameters[estimated_parameters$parameter_name ==
                                        "Nest size" &
                                        estimated_parameters$distr_param ==
                                        "shape",
                                      "distr_param_est"],
         rate = estimated_parameters[estimated_parameters$parameter_name ==
                                       "Nest size" &
                                       estimated_parameters$distr_param ==
                                       "rate",
                                     "distr_param_est"])
# Maximum foraging trail length
foraging_max_length_m <-
  qgamma(0.99,
         shape = estimated_parameters[estimated_parameters$parameter_name ==
                                        "Foraging trail length" &
                                        estimated_parameters$distr_param ==
                                        "shape",
                                      "distr_param_est"],
         rate = estimated_parameters[estimated_parameters$parameter_name ==
                                       "Foraging trail length" &
                                       estimated_parameters$distr_param ==
                                       "rate",
                                     "distr_param_est"])


# Nest abandonment threshold (1000 ants)
nest_abandonment_threshold <-
  qgamma(0.05,
         shape = estimated_parameters[estimated_parameters$parameter_name ==
                                        "Nest size" &
                                        estimated_parameters$distr_param ==
                                        "shape",
                                      "distr_param_est"],
         rate = estimated_parameters[estimated_parameters$parameter_name ==
                                       "Nest size" &
                                       estimated_parameters$distr_param ==
                                       "rate",
                                     "distr_param_est"])


not_estimated_parameters <-
  list(
    K = K,
    foraging_max_length_m = foraging_max_length_m,
    nest_abandonment_threshold = nest_abandonment_threshold,
    map_width_m = 60,
    # Width of the map (m)
    rSSI_inhibit_dist_m = 2,
    # Tree inhibition distance (used in rSSI())
    rSSI_density_tree_per_sq = 0.01,
    # Density of trees (/m^2)
    internest_min_length_m = 1.5,
    # Minimum distance for newly placed nests (m)
    internest_redraw_max_length_m = 2,
    # Upper-bound used in new nest position generation distance redraw (m)
    initial_nest_size_cap = 0.95,
    # Maximum initial size for nests as a portion of K (carrying capacity)
    internest_rejection.nb = 100,
    # Number of times a potential new internest trail is tried to be placed
    colony_loss_rate = 0.075,
    # Constant nest size loss rate
    max_growth_rate = 0.125,
    # Maximum growth rate of nests
    max_budding_prob = 0.04,
    # Maximum probability of a nest budding in each iteration
    max_budding_tries = 100,
    # Number of times a potential new nest is tried to be placed
    new_nest_initial_size = 10,
    # Starting nest size of newly budded nests (1000 ants)
    new_foraging_trail_prob = 0.0015,
    # Probability of a nest forming a new foraging trail
    new_internest_trail_prob = 0.04,
    # Probability of a nest forming a new internest trail
    trail_abandonment_threshold = 0.0025,  # Abandoning trails weaker than 10 ant / 4 m (ants/mm)
    # Burn-in length (number of iterations, weeks) 
    burn_in_length = 8
  )


# Saving results ----------------------------------------------------------
gravity_parameters <-
  gravity_pow_parameters

write.csv(gravity_parameters, file = "./Output/Model_parameters/gravity_parameters.csv", row.names = FALSE)
saveRDS(gravity_parameters, file = "./Output/Model_parameters/gravity_parameters.rds")

write.csv(estimated_parameters, file = "./Output/Model_parameters/estimated_parameters.csv", row.names = FALSE)
saveRDS(estimated_parameters, file = "./Output/Model_parameters/estimated_parameters.rds")

write.csv(not_estimated_parameters,
          file = "./Output/Model_parameters/not_estimated_parameters.csv",
          row.names = FALSE)
saveRDS(not_estimated_parameters, file = "./Output/Model_parameters/not_estimated_parameters.rds")

# Exporting session information
sink(file = "./Output/Model_parameters/sessionInfo.txt")
sessionInfo()
sink()
write_lines(print(sessionInfo()), file = "./Output/Model_parameters/sessionInfo_all.txt")

