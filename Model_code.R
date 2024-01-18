# Author: Sandor Piross - sandor.piross@gmail.com 
# Contains code from:
# Lecheval, V., Larson, H., Burns, D.D., Ellis, S., Powell, S., Donaldson-Matasci, M.C.
# and Robinson, E.J., 2021. From foraging trails to transport networks: how the
# quality-distance trade-off shapes network structure. Proceedings of the Royal
# Society B, 288(1949), p.20210430. https://doi.org/10.1098/rspb.2021.0430

# The script contains the codebase for running simulation experiments including:
  # Initalisation of networks
    # generate_network()
      # my_rSSI() - generation of trees
      # generate_nests()  - generation of nests
      # generate_internest()  - generation of internest trails
      # generate_foraging()  - generation of foraging trails
      # calculate_strength() -  calculating trail strengths
  # Updating state variables (dynamics)
    # update_network()
  # Definition of the experiments
    # run_exclusion_experiment()
  # Calculating network statistics
    # calculate_experiment_network_stats()
      # calculate_experiment_network_stats()

library(spatstat)
library(tidyverse)
library(igraph)
library(brainGraph)
library(lazyeval)

options(dplyr.summarise.inform = FALSE)
# rm(list = ls())

# Model parameters ----------------------------------------------------
estimated_parameters <- readRDS("./Output/Model_parameters/estimated_parameters.rds")
gravity_parameters <- readRDS("./Output/Model_parameters/gravity_parameters.rds")
not_estimated_parameters <- readRDS("./Output/Model_parameters/not_estimated_parameters.rds")

# Utility functions -------------------------------------------------------
calculate_descriptive = function(data, y, ...) {
  # Calculates descriptive statistics on a numeric variable by groups
  y <- enquo(y)
  data %>%
    group_by(.dots = lazyeval::lazy_dots(...)) %>% 
    summarize(Mean = mean(!!y, na.rm=T),
              N=n(),
              SD=sd(!!y, na.rm=T),
              Var_to_Mean=var(!!y, na.rm=T)/mean(!!y, na.rm=T),
              CV=sd(!!y, na.rm=T)/mean(!!y, na.rm=T),
              Min=quantile(!!y,0, na.rm=T),
              LQ=quantile(!!y,0.25, na.rm=T),
              Median=quantile(!!y,0.5, na.rm=T),
              UQ=quantile(!!y,0.75, na.rm=T),
              Max=quantile(!!y,1, na.rm=T),
              IQR=quantile(!!y,0.75, na.rm=T)-quantile(!!y,0.25, na.rm=T)
    )
}

calculate_robustness <- function(network_graph) {
  # Calculates network robustness statistic from igraph objects
  
  # Robustness: Average (over all edges) of the proportion of network efficiency
  # retained, when each edge is removed Latora, V. and Marchiori, M. (2003)
  # Economic small-world behavior in weighted networks. The European Physical
  # Journal B - Condensed Matter and Complex Systems. 32 (2), 249–263.
  # https://doi.org/10.1140/epjb/e2003-00095-5
  
  return(mean(
    # Efficiencies after each edge removed one-by-one
    unlist(lapply(E(network_graph), function(x) {
      efficiency(delete.edges(network_graph, x), type = "global")
    })) /
      # Efficiency
      efficiency(network_graph, type = "global"),
    na.rm = T
    
  ))
}

is.intersecting <- function(x, y) {
  # Checks whether segments intersect
  
  ## given coordinates of two segments [AB] and [CD], returns whether
  ## segments intersect or not.
  ## input: xA = x[1], xB = x[2], xC = x[3], xD = x[4] — same goes for y
  ## algo https://stackoverflow.com/a/7069702
  if (length(x) != 4 &
      length(y) != 4)
    stop("Input vectors should have four elements!")
  
  if (any(is.na(c(x,y))))
    stop("NA-s in coordinates!")
  
  # If segments are the same, they intersect
  if (x[1] == x[3] & x[2] == x[4] & y[1] == y[3] & y[2] == y[4]) {
    return(TRUE)
  } else {
    s1 <-
      sign((x[4] - x[3]) * (y[1] - y[4]) - (y[4] - y[3]) * (x[1] - x[4]))
    s2 <-
      sign((x[4] - x[3]) * (y[2] - y[4]) - (y[4] - y[3]) * (x[2] - x[4]))
    s3 <-
      sign((x[2] - x[1]) * (y[3] - y[2]) - (y[2] - y[1]) * (x[3] - x[2]))
    s4 <-
      sign((x[2] - x[1]) * (y[4] - y[2]) - (y[2] - y[1]) * (x[4] - x[2]))
    ans <- s1 * s2 < 0 & s3 * s4 < 0
    return(ans)
  }
}


# Generating trees --------------------------------------------------------
my_rSSI <- function(width_m = not_estimated_parameters$map_width_m,
                    inhibit_dist_m = not_estimated_parameters$rSSI_inhibit_dist_m,
                    density_tree_per_sq = not_estimated_parameters$rSSI_density_tree_per_sq) {
  # Generates treescape using spatstat::rSSI() as part of the network initialisation
  res <- rSSI(r = inhibit_dist_m,
              n = ((density_tree_per_sq) * ((width_m ^ 2))),
              win = square(width_m))
  
  tree_list <- data.frame(
    ID = paste0("T", 1:length(res$x)),
    node_type = "tree",
    nest_size = 1,
    # Only added for illustrative purposes
    x = round(res$x - mean(res$x), 3),
    y = round(res$y - mean(res$y), 3)
  ) |>
    mutate(status = "active")
  
  return(tree_list)
}


# Generating nests --------------------------------------------------------
generate_nests <- function() {
  # Generates nests as part of the network initialisation

  ## Number of nest in the colony
  n.nests <- 0
  while (n.nests <= 4 | n.nests >= 20) {
    # Number of nests must be at least two
    n.nests <- rnbinom(1,
                       mu = estimated_parameters[estimated_parameters$parameter_name ==
                                                   "Nests per colony" &
                                                   estimated_parameters$distr_param ==
                                                   "mu",
                                                 "distr_param_est"],
                       size = estimated_parameters[estimated_parameters$parameter_name ==
                                                     "Nests per colony" &
                                                     estimated_parameters$distr_param ==
                                                     "size",
                                                   "distr_param_est"])
  }
  

  
  
  ## Distance between nests
  ## Exponential distribution
  r <- rgamma(n.nests - 1,
              shape = estimated_parameters[estimated_parameters$parameter_name ==
                                             "Internest trail length" &
                                             estimated_parameters$distr_param ==
                                             "shape",
                                           "distr_param_est"],
              rate = estimated_parameters[estimated_parameters$parameter_name ==
                                            "Internest trail length" &
                                            estimated_parameters$distr_param ==
                                            "rate",
                                          "distr_param_est"]) # Distances between nests
  
  r.tooshort <- r < not_estimated_parameters$internest_min_length_m # Are distances below the empirical minimum
  
  r[r.tooshort] <-
    runif(n = sum(r.tooshort),
          min = not_estimated_parameters$internest_min_length_m,
          max = not_estimated_parameters$internest_redraw_max_length_m) # Distances too short a redrawn from a uniform distribution betwee, 1.5 and 2 m
  phi <-
    runif(n = n.nests - 1, min = -pi, max = pi) # Drawing angles [rad] from a uniform distribution
  dx <- r * cos(phi) # relative x-coordinates
  dy <- r * sin(phi) # relative y-coordinates
  nest.ref <-
    sapply(1:(n.nests - 1), function (n) {
      sample(seq(n), size = 1)
    }) # Draw a number for between 1 and number of distances for each distance
  
  # New nodes are placed relative to this node's position
  x = rep(0, n.nests)
  y = rep(0, n.nests)
  for (inest in 2:n.nests) {
    iref <- inest - 1
    x[inest] <- x[nest.ref[iref]] + dx[iref]
    y[inest] <- y[nest.ref[iref]] + dy[iref]
  }
  
  # Generating nest sizes
  nest_size <- rgamma(n.nests,
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
  
  # Initial nest sizes cannot be grater than the initial nest size cap
  nest_size <-
    ifelse(
      nest_size > not_estimated_parameters$initial_nest_size_cap * not_estimated_parameters$K,
      not_estimated_parameters$initial_nest_size_cap * not_estimated_parameters$K,
      nest_size
    )
  
  nest_list <- data.frame(
    ID = paste0("N", 1:n.nests),
    node_type = "nest",
    nest_size = nest_size,
    x = round(x, 3),
    y = round(y, 3)
  ) |> 
    mutate(status="active")
  
  return(nest_list)
}



# Generating internest trails ---------------------------------------------
generate_internest <-
  function(nest_list,
           rejection.nb = not_estimated_parameters$internest_rejection.nb,
           av.rejection = 0) {
    # Creates internest trails as part of the network initialisation
    n.nests <- nrow(nest_list)
    x <- nest_list$x
    y <- nest_list$y
    X <- x
    Y <- y
    edges <- c(NULL, NULL)
    all.pairs <- combn(1:n.nests, 2) # All nest pairs
    all.dists <- sqrt((x[all.pairs[1,]] - x[all.pairs[2,]]) ** 2 +
                        (y[all.pairs[1,]] - y[all.pairs[2,]]) ** 2)  # All distances between nests

    for (inest in 1:n.nests) {
      # For all nests
      poss.pairs.cols <-
        (all.pairs[1,] == inest) |
        (all.pairs[2,] == inest) # Where is the focal nest in the combinations
      possible.dists <-
        all.dists[poss.pairs.cols] # Distances of possible connections of the focal nest
      dist.sort <-
        sort(possible.dists, index.return = TRUE) # Sorted possible distances
      failed.connections <- 0
      go <- TRUE
      while (go) {
        edge.intersection <- FALSE
        # Choosing a distance (node at random)
        rand.dist <-
          rgamma(1,
                 shape = estimated_parameters[estimated_parameters$parameter_name ==
                                                "Internest trail length" &
                                                estimated_parameters$distr_param ==
                                                "shape",
                                              "distr_param_est"],
                 rate = estimated_parameters[estimated_parameters$parameter_name ==
                                               "Internest trail length" &
                                               estimated_parameters$distr_param ==
                                               "rate",
                                             "distr_param_est"])
        
        rand.node <-
          which(abs(dist.sort$x - rand.dist) == min(abs(dist.sort$x - rand.dist))) # Which nest is closest to the randomly drawn distance
        if (rand.node <= n.nests - 1) {
          # If random node exists
          node.co <-
            all.pairs[, which(poss.pairs.cols)[dist.sort$ix[rand.node]]] # Finding the corresponding connection
          nb.duplicates <-
            sum(edges[, 1] == node.co[2] &
                  edges[, 2] == node.co[1]) +
            sum(edges[, 1] == node.co[1] & edges[, 2] == node.co[2])
          if (nb.duplicates == 0 &
              !is.null(nrow(edges))) {
            # If we have no duplicates and have edges
            for (iedge in 1:nrow(edges)) {
              # For all edges
              nests.indices <- c(node.co, edges[iedge,])
              if (is.intersecting(X[nests.indices], Y[nests.indices])) {
                # Examining that the trails intersect
                edge.intersection <- TRUE
                failed.connections <- failed.connections + 1
                av.rejection <- av.rejection + 1
                node.co <- c(NA, NA)
                break
                
              }
            }
            if (failed.connections < rejection.nb) {
              go <- edge.intersection # Repeat till edges don't intersect
            } else {
              go <- FALSE
            }
          } else if (failed.connections < rejection.nb) {
            failed.connections <- failed.connections + 1
            av.rejection <- av.rejection + 1
            node.co <- c(NA, NA)
          } else {
            go <- FALSE
          }
        }
      }
      edges <- na.omit(rbind(edges, node.co))
    }
    
    internest_list <- data.frame(
      from_ID = paste0("N", edges[, 1]),
      to_ID = paste0("N", edges[, 2]),
      trail_type = "internest"
    )
    
    internest_list <-
      internest_list |>
      left_join(nest_list |>
                  dplyr::select(ID, x, y), by = c("from_ID" = "ID")) |>
      rename(from_x = x,
             from_y = y) |>
      relocate(from_x, from_y, .after = last_col()) |>
      left_join(nest_list |>
                  dplyr::select(ID, x, y), by = c("to_ID" = "ID")) |>
      rename(to_x = x,
             to_y = y) |>
      relocate(to_x, to_y, .after = last_col())
    
    internest_list <- 
      internest_list |> 
      mutate(status="active")
    
    return(internest_list)
  }


# Generating foraging trails ----------------------------------------------
generate_foraging <-
  function(nest_list, tree_list, internest_list) {
    # Creates foraging trails as part of the network initialisation
    
    # Creating dataframe from foraging trails
    foraging_list <- internest_list[0, ]
    
    # Drawing the number of foraging trails for each nest
    # Ther must be at least one foraing trial in the colony
    n_foraging_sum <- 0
    while(n_foraging_sum<2){
    
    nest_list$n_foraging <- rpois(nrow(nest_list),
                                  estimated_parameters[estimated_parameters$parameter_name ==
                                                         "Number of foraging trails per nest" &
                                                         estimated_parameters$distr_param == "lambda",
                                                       "distr_param_est"])
    n_foraging_sum <- sum(nest_list$n_foraging)}
    
    # Which nests need trees? (Larger nests come first)
    trees_needed <-
      nest_list |>
      arrange(desc(nest_size)) |>
      mutate(row_ID = ID) |>
      column_to_rownames("row_ID")
    
    # Duplicating rows by the number of their foraging trails
    trees_needed <-
      trees_needed[rep(trees_needed$ID, trees_needed$n_foraging), c("ID", "x", "y")]
    
    
    for (i in 1:nrow(trees_needed)) {
      # Which trees are reachable for that nest without intersecting other trails?
      segments_to_all_trees <-
        tree_list |>
        mutate(nest.x = trees_needed[i, "x"],
               nest.y = trees_needed[i, "y"]) |>
        filter(!(ID %in% foraging_list[foraging_list$from_ID == trees_needed[i, "ID"], "to_ID"]))
      
      
      # All trails
      trail_list <- bind_rows(na.omit(foraging_list), internest_list)
      
      any_intersecting <- logical(nrow(segments_to_all_trees))
      for (j in 1:nrow(segments_to_all_trees)) {
        intersecting <- logical(nrow(trail_list))
        for (k in 1:nrow(trail_list)) {
          intersecting[k] <-
            is.intersecting(
              x = c(segments_to_all_trees[j, c("x", "nest.x")],
                    trail_list[k, c("from_x", "to_x")], recursive = T),
              y = c(segments_to_all_trees[j, c("y", "nest.y")],
                    trail_list[k, c("from_y", "to_y")], recursive = T)
            )
        }
        any_intersecting[j] <- any(intersecting)
      }
      
      segments_to_all_trees$any_intersecting <- any_intersecting
      
      # Trees available without intersecting other trails
      available_trees <- segments_to_all_trees[!any_intersecting, ]
      
      # Distances to available trees
      available_trees <-
        available_trees |>
        mutate(distance = sqrt((nest.x - x) ^ 2 + (nest.y - y) ^ 2)) |>
        arrange(distance) # Sorting trees based on dstance
      
      
      # Excluding trees from available list that are further than allowed max distance
      available_trees <- 
      available_trees |> 
        filter(distance<not_estimated_parameters$foraging_max_length_m)
      
      # If no trees are available
      if (nrow(available_trees) == 0) {
        # Adding to existing trails
        foraging_list[i, ] <-
          data.frame(
            from_ID = trees_needed[i, "ID"],
            to_ID = "No available tree",
            trail_type = "foraging",
            from_x = trees_needed[i, "x"],
            from_y = trees_needed[i, "y"],
            to_x = NA,
            to_y = NA
          )
      } else {

        # Choosing the closest available tree
        to_tree <- available_trees[1, "ID"]
        
        # Adding to existing trails
        foraging_list[i, ] <-
          data.frame(
            from_ID = trees_needed[i, "ID"],
            to_ID = to_tree,
            trail_type = "foraging",
            from_x = trees_needed[i, "x"],
            from_y = trees_needed[i, "y"],
            to_x = tree_list[tree_list$ID ==
                               to_tree, "x"],
            to_y = tree_list[tree_list$ID ==
                               to_tree, "y"]
          )
      }
      
    }
    
    foraging_list <- 
      foraging_list |> 
      mutate(status="active") |> 
      filter(to_ID != "No available tree")
    
    return(foraging_list)
  }


# Calculating strengths ---------------------------------------------------

calculate_strength <-
  function(nest_list,
           tree_list,
           internest_list,
           foraging_list) {
    # Calculates trail strengths using the gravity formulae
    
    trail_list <-
      bind_rows(internest_list,
                foraging_list) |>
      #Adding from nest sizes
      left_join(nest_list |>
                  dplyr::select(ID, nest_size),
                by = c("from_ID" = "ID")) |>
      rename(from_nest_size = nest_size) |>
      # Adding to nest sizes
      left_join(nest_list |>
                  dplyr::select(ID, nest_size),
                by = c("to_ID" = "ID")) |>
      rename(to_nest_size = nest_size) |>
      # Setting nest size of trees to 1
      mutate(to_nest_size = ifelse(trail_type == "foraging", 1, to_nest_size)) |>
      # Calculating distances
      mutate(distance = sqrt(((to_x - from_x) ^ 2) +
                               ((to_y - from_y) ^ 2))) |>
      # Adding gravity model parameters
      left_join(gravity_parameters, by="trail_type") |>
      # Calculating strength based on gravity model
      mutate(
             strength = (G_pow * (((from_nest_size * to_nest_size) ^ (alpha_pow)) / (distance ^ beta_pow)))
      ) |> 
      dplyr::select(-c(G_pow,alpha_pow,beta_pow))
    
    return(trail_list)
    
  }


# Generating complete networks --------------------------------------------

## Generating the network
generate_network <- function(tree_generating_function = my_rSSI(),
                             burn_in_length=not_estimated_parameters$burn_in_length) {
  
  # Generates the complete networks using the previous functions
  
  burn_in <- (-burn_in_length)
  
  nest_list <-  generate_nests() |> 
    mutate(t=burn_in) |> 
    relocate(t)
  #tree_generating_function <-  my_rSSI()
  tree_list <-  tree_generating_function |> 
    mutate(t=burn_in) |> 
    relocate(t)
  internest_list <-  generate_internest(nest_list)
  foraging_list <- generate_foraging(nest_list,
                                     tree_list,
                                     internest_list)
  trail_list <- calculate_strength(nest_list,
                                     tree_list,
                                     internest_list,
                                   foraging_list) |> 
    mutate(t=burn_in) |> 
    relocate(t)
  
  return(
    list(
      nest_list = nest_list,
      tree_list = tree_list,
      trail_list = trail_list
    ))

}


# Update network ----------------------------------------------------------

update_network <- function(network){
  # Updates state variables in each iteration
  
  # Getting lists from last time step
  network_last_t <- list(
    nest_list = network$nest_list |>
      filter(t == max(t)),
    tree_list = network$tree_list |>
      filter(t == max(t)),
    trail_list = network$trail_list |>
      filter(t == max(t))
  )
  
  
  # Updating nest sizes -----------------------------------------------------
  
  # Calculating foraging activities
  activities <-
    network_last_t$trail_list |>
    dplyr::select(trail_type,
                  from_ID,
                  to_ID,
                  from_nest_size,
                  to_nest_size,
                  from_x,
                  from_y,
                  to_x,
                  to_y,
                  distance,
                  strength,
                  t,
                  status) |>
    # Selecting internest trails
    filter(trail_type == "internest") |>
    # Reversing internest trails
    rename(
      to_ID = from_ID,
      from_ID = to_ID,
      from_nest_size = to_nest_size,
      to_nest_size = from_nest_size,
      from_x = to_x,
      from_y = to_y,
      to_x = from_x,
      to_y = from_y
    ) |>
    # Add reversed internest trails to trail list
    bind_rows(network_last_t$trail_list) |>
    # Calculating incoming (to) and outgoing (from) strengths
    mutate(
    strength_split = ifelse((from_nest_size + to_nest_size)==0,0,(strength * (
      from_nest_size / (from_nest_size + to_nest_size))))
    )

  # Calculating cumulative activities for each nest 
  node_activities <-
    full_join(
      activities |>
        group_by(from_ID) |>
        summarise(strength_split_from = sum(strength_split)),
      activities |>
        filter(trail_type == "internest") |>
        group_by(to_ID) |>
        summarise(strength_split_to = sum(strength_split)),
      by = c("from_ID" = "to_ID")
    ) |>
    rename(ID = from_ID)
  
  # Updating nest sizes with discrete logistic growth
  max_growth_rate <- not_estimated_parameters$max_growth_rate
  network_last_t$nest_list <- 
    network_last_t$nest_list |> 
    left_join(node_activities, by="ID") |> 
    mutate(growth_rate = 
             max_growth_rate*(ifelse((strength_split_from + strength_split_to)==0,0,(strength_split_from - strength_split_to)/(strength_split_from + strength_split_to)))-not_estimated_parameters$colony_loss_rate,
           new_nest_size = nest_size+
             (growth_rate*nest_size*
                ((not_estimated_parameters$K-nest_size)/not_estimated_parameters$K)),
           new_nest_size = ifelse(new_nest_size<0,0,new_nest_size),
           new_nest_size = ifelse(new_nest_size>not_estimated_parameters$K,not_estimated_parameters$K,new_nest_size)
    ) |> 
    dplyr::select(t, ID, node_type, new_nest_size, x, y, status) |> 
    rename(nest_size=new_nest_size)
  
  # Updating nest sizes in trail list
  network_last_t$trail_list <- 
    network_last_t$nest_list |> 
    dplyr::select(ID,nest_size) |> 
    right_join(network_last_t$trail_list, by=c("ID"="from_ID")) |>
    rename(from_ID=ID) |> 
    mutate(from_nest_size=nest_size) |> 
    dplyr::select(-nest_size)
  
  network_last_t$trail_list <- 
    network_last_t$nest_list |> 
    dplyr::select(ID,nest_size) |> 
    right_join(network_last_t$trail_list, by=c("ID"="to_ID")) |> 
    rename(to_ID=ID) |> 
    relocate(from_ID) |> 
    mutate(to_nest_size=nest_size) |>
    mutate(to_nest_size=ifelse(trail_type=="foraging",1,to_nest_size)) |> 
    dplyr::select(-nest_size)
  
  
  # Abandoning nests based on nest size  --------------------
  
  # Nests below threshold get abandoned 
  network_last_t$nest_list <- 
    network_last_t$nest_list |> 
    mutate(nest_size=ifelse(nest_size<=not_estimated_parameters$nest_abandonment_threshold,0,nest_size),
           status=ifelse(nest_size==0,"abandoned",status))
  
  # Trails connected to abandoned nests also get abandoned
  abandoned_trails <- (network_last_t$trail_list$from_ID %in% (network_last_t$nest_list |> 
                                                                 filter(status=="abandoned") |> 
                                                                 dplyr::select(ID) |> 
                                                                 pull()))|
    (network_last_t$trail_list$to_ID %in% (network_last_t$nest_list |> 
                                             filter(status=="abandoned") |> 
                                             dplyr::select(ID) |> 
                                             pull()))
  network_last_t$trail_list <- 
    network_last_t$trail_list |> 
    mutate(status=ifelse(abandoned_trails,"abandoned",status ),
           strength=ifelse(abandoned_trails,0,strength))
  
  
  # Budding new nests -------------------------------------------------------
  network_last_t$nest_list <- 
    network_last_t$nest_list |> 
    mutate(new = FALSE)
  
  # Finding active nests
  active_nests <- network_last_t$nest_list |>
    filter(status == "active") 
  
  # Finding active trails
  active_trails <- network_last_t$trail_list |>
    filter(status == "active")
  
  # If there are no active nests and trails left, skip trail formation
  if(nrow(active_nests>0)){
    if(nrow(active_trails>0)){
      
  # Number of new nests budded by each active nest
  budded_numbers <- numeric(nrow(active_nests)) 
  for (i in 1:nrow(active_nests)) {
    budded_numbers[i] <- 
      rbinom(
        1, size = 1,
        prob = (not_estimated_parameters$max_budding_prob * (
          pgamma(active_nests$nest_size[i],
                 shape = estimated_parameters[estimated_parameters$parameter_name ==
                                                "Nest size" &
                                                estimated_parameters$distr_param ==
                                                "shape",
                                              "distr_param_est"],
                 rate = estimated_parameters[estimated_parameters$parameter_name ==
                                               "Nest size" &
                                               estimated_parameters$distr_param ==
                                               "rate",
                                             "distr_param_est"]))))
    
  }
  
  # If there are nests to bud
  if (sum(budded_numbers) > 0) {
    # List of budder (parent nests)
    budders <-
      factor(rep(active_nests$ID, times = budded_numbers))
    
    # Sorting list by nest size
    budders <- factor(
      budders,
      levels = network_last_t$nest_list |>
        filter(ID %in% budders) |>
        arrange(desc(nest_size)) |>
        dplyr::select(ID) |> pull()
    )
    
    budders <- budders[order(budders)]
    
    
    # For each newly budded nest
    for (i in 1:length(budders)) {
      # Trying to find non-intersecting trail
      for (j in 1:not_estimated_parameters$max_budding_tries) {
        # Generating coordinates for new nests with non-intersecting internest trails to their parents
        r <- rgamma(1,
                    shape = estimated_parameters[estimated_parameters$parameter_name ==
                                                   "Internest trail length" &
                                                   estimated_parameters$distr_param ==
                                                   "shape",
                                                 "distr_param_est"],
                    rate = estimated_parameters[estimated_parameters$parameter_name ==
                                                  "Internest trail length" &
                                                  estimated_parameters$distr_param ==
                                                  "rate",
                                                "distr_param_est"]) # Distances between nests
        
        r.tooshort <- r < not_estimated_parameters$internest_min_length_m # Are distances below the empirical minimum
        
        r[r.tooshort] <-
          runif(n = sum(r.tooshort),
                min = not_estimated_parameters$internest_min_length_m,
                max = not_estimated_parameters$internest_redraw_max_length_m) # Distances too short a redrawn from a uniform distribution between, 1.5 and 2 m
        phi <-
          runif(n =  1,
                min = -pi,
                max = pi) # Drawing angles [rad] from a uniform distribution
        dx <- r * cos(phi) # relative x-coordinates
        dy <- r * sin(phi) # relative y-coordinates
        
        # New coordinates
        new_trail <-
          c(active_nests[active_nests$ID == budders[1], c("x", "y")][1,],
            active_nests[active_nests$ID == budders[1], c("x", "y")][1,] +
              c(dx, dy),
            recursive = TRUE)
        
        # Does the new trail intersect any active trails?
        if(nrow(active_trails)==0){
          intersects <- FALSE
        } else {
        
        intersects <- rep(NA,nrow(active_trails))
        
        if(nrow(active_trails)>0){
          intersects <- any(sapply(1:nrow(active_trails),
                                   function(z) {
                                     is.intersecting(
                                       x = c(new_trail[c(1, 3)], active_trails[z, c("from_x", "to_x")] , recursive =
                                               TRUE),
                                       y = c(new_trail[c(2, 4)], active_trails[z, c("from_y", "to_y")] , recursive =
                                               TRUE)
                                     )
                                   }))
        }

        }
        
        if (intersects == FALSE) {
          # Add to trail list
          new_trail_list_item <-
            data.frame(
              from_ID = budders[i],
              # Next nest ID (max ID number +1)
              to_ID = paste0("N", max(as.numeric(
                substr(network_last_t$nest_list$ID, 2, nchar(network_last_t$nest_list$ID))
              )) + 1),
              trail_type = "internest",
              from_x = new_trail[1],
              from_y = new_trail[2],
              to_x = new_trail[3],
              to_y = new_trail[4],
              status = "active",
              from_nest_size = active_nests[active_nests$ID == budders[i], "nest_size"],
              to_nest_size = not_estimated_parameters$new_nest_initial_size,
              distance = r,
              strength = 0,
              t = network_last_t$trail_list$t[1]
            )
          
          network_last_t$trail_list <- bind_rows(new_trail_list_item,network_last_t$trail_list)
          
          # Adding to nest list
          network_last_t$nest_list <-
            data.frame(
              ID = new_trail_list_item$to_ID,
              node_type = "nest",
              nest_size = new_trail_list_item$to_nest_size,
              x = new_trail_list_item$to_x,
              y = new_trail_list_item$to_y,
              status = "active",
              t = network_last_t$nest_list$t[1],
              new = TRUE
            ) |>
            bind_rows(network_last_t$nest_list)
          
          break
        } 
      }
      
    }
  }
  
  

# Adding foraging trails to newly budded nests and stochastic formation --------

  # Which nests need trees? (Larger nests come first)
  trees_needed <-
    network_last_t$nest_list |>
    filter(new==TRUE) |> 
    arrange(desc(nest_size)) |>
    mutate(row_ID = ID) |>
    column_to_rownames("row_ID")
  
  # Drawing foraging trail numbers for new nests
  trees_needed <- 
    trees_needed |>
    mutate(n_foraging = sign(rpois(nrow(trees_needed),
                              estimated_parameters[estimated_parameters$parameter_name ==
                                                     "Number of foraging trails per nest" &
                                                     estimated_parameters$distr_param == "lambda",
                                                   "distr_param_est"])))
  
  
  # Duplicating rows by the number of their foraging trails
  trees_needed <-
    trees_needed[rep(trees_needed$ID, trees_needed$n_foraging), c("ID", "x", "y")]
  
  # Stochastic foraging trail formation: each old (not newly budded nest has n%
  # chance of forming a new foraging trail)
  new_foraging_trails <- 
  network_last_t$nest_list |>
    filter(new==FALSE,
           status=="active") |>
    arrange(desc(nest_size)) |>
    mutate(new_foraging = as.logical(rbinom(
      network_last_t$nest_list |>
        filter(new==FALSE,
               status=="active") |>
        nrow()
      , 1, not_estimated_parameters$new_foraging_trail_prob
    ))) |> 
    filter(new_foraging==TRUE) |> 
    dplyr::select(ID,x,y) 
  
  trees_needed <-
    trees_needed |> 
    bind_rows(new_foraging_trails)
  
  # Adding new foraging trails
  if(nrow(trees_needed)>0){
    for (i in 1:nrow(trees_needed)) {
      foraging_list <- 
        network_last_t$trail_list |>
        filter(status == "active",
               trail_type == "foraging")
      
      # Which trees are reachable for that nest without intersecting other trails?
      segments_to_all_trees <-
        network_last_t$tree_list |>
        filter(status=="active") |> 
        mutate(nest.x = trees_needed[i, "x"],
               nest.y = trees_needed[i, "y"]) |>
        filter(!(ID %in% foraging_list[foraging_list$from_ID == trees_needed[i, "ID"], "to_ID"]))
      
      
      # All trails
      trail_list <- network_last_t$trail_list |> filter(status=="active")
      
      any_intersecting <- logical(nrow(segments_to_all_trees))
      for (j in 1:nrow(segments_to_all_trees)) {
        intersecting <- logical(nrow(trail_list))
        for (k in 1:nrow(trail_list)) {
          intersecting[k] <-
            is.intersecting(
              x = c(segments_to_all_trees[j, c("x", "nest.x")],
                    trail_list[k, c("from_x", "to_x")], recursive = T),
              y = c(segments_to_all_trees[j, c("y", "nest.y")],
                    trail_list[k, c("from_y", "to_y")], recursive = T)
            )
        }
        any_intersecting[j] <- any(intersecting)
      }
      
      segments_to_all_trees$any_intersecting <- any_intersecting
      
      # Trees available without intersecting other trails
      available_trees <- segments_to_all_trees[!any_intersecting, ]
      
      # Distances to available trees
      available_trees <-
        available_trees |>
        mutate(distance = sqrt((nest.x - x) ^ 2 + (nest.y - y) ^ 2)) |>
        arrange(distance) # Sorting trees based on distance
      
      # Excluding trees from available list that are further than allowed max distance
      available_trees <- 
        available_trees |> 
        filter(distance<not_estimated_parameters$foraging_max_length_m)
      
      # If no trees are available
      if (nrow(available_trees) == 0) {
        # # Adding to existing trails
        # network_last_t$trail_list[i, ] <-
        #   data.frame(
        #     from_ID = trees_needed[i, "ID"],
        #     to_ID = "No available tree",
        #     trail_type = "foraging",
        #     from_x = trees_needed[i, "x"],
        #     from_y = trees_needed[i, "y"],
        #     to_x = NA,
        #     to_y = NA
        #   )
      } else {
        
        
        # Choosing the closest available tree
        to_tree <- available_trees[1, "ID"]
        
        # Adding to existing trails
        network_last_t$trail_list <-
          data.frame(
            from_ID = trees_needed[i, "ID"],
            to_ID = to_tree,
            trail_type = "foraging",
            from_x = trees_needed[i, "x"],
            from_y = trees_needed[i, "y"],
            to_x = network_last_t$tree_list[network_last_t$tree_list$ID ==
                                              to_tree, "x"],
            to_y = network_last_t$tree_list[network_last_t$tree_list$ID ==
                                              to_tree, "y"],
            status = "active",
            t = network_last_t$trail_list$t[1]
          ) |> 
          bind_rows(network_last_t$trail_list)
      }
      
    }
  }
  

# Stochastic internest trail formation -------------------------------------
  # Stochastic internest trail formation: each old (not newly budded nest has n%
  # chance of forming a new internest trail)
  nests_needed <- 
    network_last_t$nest_list |>
    filter(new==FALSE,
           status=="active") |>
    arrange(desc(nest_size)) |>
    mutate(new_internest = as.logical(rbinom(
      network_last_t$nest_list |>
        filter(new==FALSE,
               status=="active") |>
        nrow()
      , 1, not_estimated_parameters$new_internest_trail_prob
    ))) |> 
    filter(new_internest==TRUE) |> 
    dplyr::select(ID,x,y)
  
  
  # Adding new internest trails
  if((nrow(nests_needed)>0)&(nrow(active_nests)>1)){
    for (i in 1:nrow(nests_needed)) {

      # All trails
      trail_list <- network_last_t$trail_list |> filter(status=="active")
      
      # Which nests are reachable for that nest without intersecting other trails?
      segments_to_all_nests <-
        network_last_t$nest_list |>
        mutate(nest.x = nests_needed[i, "x"],
               nest.y = nests_needed[i, "y"]) |>
        filter(!(ID %in% trail_list[trail_list$from_ID == nests_needed[i, "ID"], "to_ID"]))
      
      any_intersecting <- as.logical(rep(NA,nrow(segments_to_all_nests)))
      if((nrow(segments_to_all_nests)>0)&(nrow(trail_list)>0)){
        for (j in 1:nrow(segments_to_all_nests)) {
          intersecting <- logical(nrow(trail_list))
          for (k in 1:nrow(trail_list)) {
            intersecting[k] <-
              is.intersecting(
                x = c(segments_to_all_nests[j, c("x", "nest.x")],
                      trail_list[k, c("from_x", "to_x")], recursive = T),
                y = c(segments_to_all_nests[j, c("y", "nest.y")],
                      trail_list[k, c("from_y", "to_y")], recursive = T)
              )
          }
          any_intersecting[j] <- any(intersecting)
        }
        
        segments_to_all_nests$any_intersecting <- any_intersecting
        
        # Nests available without intersecting other trails
        available_nests <- segments_to_all_nests[!any_intersecting, ] |> 
          filter(ID != nests_needed[i,"ID"])
      }
      
      # If no nests are available
      if (nrow(available_nests) == 0) {
        # # Adding to existing trails
        # network_last_t$trail_list[i, ] <-
        #   data.frame(
        #     from_ID = nests_needed[i, "ID"],
        #     to_ID = "No available nest",
        #     trail_type = "internest",
        #     from_x = nests_needed[i, "x"],
        #     from_y = nests_needed[i, "y"],
        #     to_x = NA,
        #     to_y = NA
        #   )
      } else {
        # Distances to available trees
        available_nests <-
          available_nests |>
          mutate(distance = sqrt((nest.x - x) ^ 2 + (nest.y - y) ^ 2)) |>
          arrange(distance) # Sorting trees based on dstance
        
        # Choosing the closest available nest
        to_nest <- available_nests[1, "ID"]
        
        # Adding to existing trails
        network_last_t$trail_list <-
          data.frame(
            from_ID = nests_needed[i, "ID"],
            to_ID = to_nest,
            trail_type = "internest",
            from_x = nests_needed[i, "x"],
            from_y = nests_needed[i, "y"],
            to_x = network_last_t$nest_list[network_last_t$nest_list$ID ==
                                              to_nest, "x"],
            to_y = network_last_t$nest_list[network_last_t$nest_list$ID ==
                                              to_nest, "y"],
            status = "active",
            t = network_last_t$trail_list$t[1]
          ) |> 
          bind_rows(network_last_t$trail_list)
      }
      
    }
  }
  
  # Removing "new" variable from nest list
  network_last_t$nest_list <- 
    network_last_t$nest_list |> 
    dplyr::select(-new)
  
  # Filling out missing variables in trail list
  network_last_t$trail_list <- 
    network_last_t$trail_list |> 
    mutate(strength = ifelse(is.na(strength),0,strength),
           distance = sqrt(((to_x - from_x) ^ 2) +
                             ((to_y - from_y) ^ 2)))
  
  # Updating nest sizes in trail list
  network_last_t$trail_list <- 
    network_last_t$nest_list |> 
    dplyr::select(ID,nest_size) |> 
    right_join(network_last_t$trail_list, by=c("ID"="from_ID")) |>
    rename(from_ID=ID) |> 
    mutate(from_nest_size=nest_size) |> 
    dplyr::select(-nest_size)
  
  network_last_t$trail_list <- 
    network_last_t$nest_list |> 
    dplyr::select(ID,nest_size) |> 
    right_join(network_last_t$trail_list, by=c("ID"="to_ID")) |> 
    rename(to_ID=ID) |> 
    relocate(from_ID) |> 
    mutate(to_nest_size=nest_size) |>
    mutate(to_nest_size=ifelse(trail_type=="foraging",1,to_nest_size)) |> 
    dplyr::select(-nest_size)
  
    }
  }
  # Updating trail strengths ------------------------------------------------
  
  network_last_t$trail_list <- 
    network_last_t$trail_list |>
    left_join(gravity_parameters, by="trail_type") |>
    # Calculating strength based on gravity model
    mutate(
      strength = ifelse(from_nest_size * to_nest_size==0,0,
                        (G_pow * (((from_nest_size * to_nest_size) ^ (alpha_pow)) / (distance ^ beta_pow))))
    ) |> 
    dplyr::select(-c(G_pow,alpha_pow,beta_pow)) |> 
    mutate(strength=ifelse(status=="abandoned",0,strength)) |> 
  # Abandoning trails weaker than 10 ant / 4 m
  mutate(status=ifelse(strength<not_estimated_parameters$trail_abandonment_threshold,"abandoned",status),
         strength=ifelse(strength<not_estimated_parameters$trail_abandonment_threshold,0,strength)) 

  

  
  # Adding new timestep to network object -----------------------------------
  
  network$nest_list <- 
    network$nest_list |> 
    bind_rows(network_last_t$nest_list |> 
                mutate(t=max(network$nest_list$t)+1))
  
  network$trail_list <- 
    network$trail_list |> 
    bind_rows(network_last_t$trail_list |> 
                mutate(t=max(network$trail_list$t)+1))
  
  network$tree_list <- 
    network$tree_list |> 
    bind_rows(network_last_t$tree_list |> 
                mutate(t=max(network$tree_list$t)+1))

  
  return(network)
}


# Running simulations -----------------------------------------------------

run_exclusion_experiment <- function() {
  # Runs exclusion experiments
  
  # Generating a network
  my_network_control <- generate_network()
  
  # Running burn-in
  for (i in -(not_estimated_parameters$burn_in_length-1):0) {
    my_network_control <- update_network(my_network_control)
  }
  
  # Exclusion treatments on the network
  my_network_strongest <- my_network_control
  my_network_weakest <- my_network_control
  my_network_random <- my_network_control
  
  # Finding trees to exclude in the experiments
  tree_sum_strengths <- 
  my_network_control$trail_list |>
    filter(trail_type=="foraging",
           status=="active",
           t==0) |> 
    group_by(to_ID) |> 
  # Finding the trees with the strongest and weakest cumulative trail strengths
  # and a random tree
    summarise(sum_strength=sum(strength)) |> 
    mutate(exclusion_category="not_excluded") 
  
  strongest <- tree_sum_strengths[tree_sum_strengths$sum_strength==max(tree_sum_strengths$sum_strength),]
  strongest$exclusion_category <- "strongest"
  
  weakest <- tree_sum_strengths[tree_sum_strengths$sum_strength==min(tree_sum_strengths$sum_strength),]
  weakest$exclusion_category <- "weakest"
  
  random <- tree_sum_strengths[sample(1:nrow(tree_sum_strengths),1),]
  random$exclusion_category <- "random"
  
  excluded_trees <- 
    bind_rows(strongest, random, weakest)
  
  # Affected trails
  affected_trails <- 
  my_network_control$trail_list |>
    filter(t==0) |> 
    left_join(excluded_trees |> 
                select(-sum_strength), by="to_ID") |> 
    filter(!is.na(exclusion_category))
  
  # Removing trails to that tree
  # Strongest exclusion 
  my_network_strongest$trail_list <-
    my_network_strongest$trail_list |>
    mutate(
      status = ifelse((to_ID == excluded_trees |> 
                        filter(exclusion_category=="strongest") |> 
                        pull(to_ID)) & (t==0) ,
                      "abandoned", status),
      strength = ifelse((to_ID == excluded_trees |> 
                           filter(exclusion_category=="strongest") |> 
                           pull(to_ID)) & (t==0),
                        0, strength)
    )
  
  # Weakest exclusion 
  my_network_weakest$trail_list <-
    my_network_weakest$trail_list |>
    mutate(
      status = ifelse((to_ID == excluded_trees |> 
                         filter(exclusion_category=="weakest") |> 
                         pull(to_ID)) & (t==0),
                      "abandoned", status),
      strength = ifelse((to_ID == excluded_trees |> 
                           filter(exclusion_category=="weakest") |> 
                           pull(to_ID)) & (t==0),
                        0, strength)
    )
  
  # Random exclusion 
  my_network_random$trail_list <-
    my_network_random$trail_list |>
    mutate(
      status = ifelse((to_ID == excluded_trees |> 
                         filter(exclusion_category=="random") |> 
                         pull(to_ID)) & (t==0),
                      "abandoned", status),
      strength = ifelse((to_ID == excluded_trees |> 
                           filter(exclusion_category=="random") |> 
                           pull(to_ID)) & (t==0),
                        0, strength)
    )
  
  # Excluding tree from active tree list
  # Strongest
  my_network_strongest$tree_list <-
    my_network_strongest$tree_list |>
    mutate(status = ifelse((ID == excluded_trees |> 
                             filter(exclusion_category=="strongest") |> 
                             pull(to_ID)) & (t==0),
                           "abandoned", status))
  
  # Weakest
  my_network_weakest$tree_list <-
    my_network_weakest$tree_list |>
    mutate(status = ifelse((ID == excluded_trees |> 
                              filter(exclusion_category=="weakest") |> 
                              pull(to_ID)) & (t==0),
                           "abandoned", status))
  
  # Random
  my_network_random$tree_list <-
    my_network_random$tree_list |>
    mutate(status = ifelse((ID == excluded_trees |> 
                              filter(exclusion_category=="random") |> 
                              pull(to_ID)) & (t==0),
                           "abandoned", status))

  ## Running simulation: first year
  for (i in 1:28) {
    my_network_control <- update_network(my_network_control)
    my_network_strongest <- update_network(my_network_strongest)
    my_network_weakest <- update_network(my_network_weakest)
    my_network_random <- update_network(my_network_random)
  }
  
  # Splitting exclusion networks to reintroduction and sustained exclusion groups
  my_network_strongest_re <- my_network_strongest
  my_network_strongest_sus <- my_network_strongest
  
  my_network_weakest_re <- my_network_weakest
  my_network_weakest_sus <- my_network_weakest
  
  my_network_random_re <- my_network_random
  my_network_random_sus <- my_network_random
  
  # Reintroducing excluded tree
  my_network_strongest_re$tree_list <-
    my_network_strongest_re$tree_list |>
    mutate(status = ifelse(t == max(t) &
                             status == "abandoned", "active", status))
  
  my_network_weakest_re$tree_list <-
    my_network_weakest_re$tree_list |>
    mutate(status = ifelse(t == max(t) &
                             status == "abandoned", "active", status))
  
  my_network_random_re$tree_list <-
    my_network_random_re$tree_list |>
    mutate(status = ifelse(t == max(t) &
                             status == "abandoned", "active", status))
  
  # Running simulation: second and third years
  for (i in 1:56) {
    my_network_control <- update_network(my_network_control)
    
    my_network_strongest_re <- update_network(my_network_strongest_re)
    my_network_strongest_sus <- update_network(my_network_strongest_sus)
    
    my_network_weakest_re <- update_network(my_network_weakest_re)
    my_network_weakest_sus <- update_network(my_network_weakest_sus)
    
    my_network_random_re <- update_network(my_network_random_re)
    my_network_random_sus <- update_network(my_network_random_sus)
  }
  
  # Did excluded trees rejoin the network in the reintroduction groups?
  rejoin <-
    bind_rows(
      my_network_strongest_re$trail_list |>
        filter(
          t > 28,
          status == "active",
          to_ID == excluded_trees |>
            filter(exclusion_category == "strongest") |>
            pull(to_ID)
        ) |>
        group_by(from_ID) |>
        filter(t == min(t)) |>
        mutate(treatment = "strongest"),
      
      my_network_weakest_re$trail_list |>
        filter(
          t > 28,
          status == "active",
          to_ID == excluded_trees |>
            filter(exclusion_category == "weakest") |>
            pull(to_ID)
        ) |>
        group_by(from_ID) |>
        filter(t == min(t)) |>
        mutate(treatment = "weakest"),
      
      my_network_random_re$trail_list |>
        filter(
          t > 28,
          status == "active",
          to_ID == excluded_trees |>
            filter(exclusion_category == "random") |>
            pull(to_ID)
        ) |>
        group_by(from_ID) |>
        filter(t == min(t)) |>
        mutate(treatment = "random")
    ) |> suppressWarnings()
  
  return(list(
    networks = list(
      control = my_network_control,
      
      strongest_re = my_network_strongest_re,
      strongest_sus = my_network_strongest_sus,
      
      weakest_re = my_network_weakest_re,
      weakest_sus = my_network_weakest_sus,
      
      random_re = my_network_random_re,
      random_sus = my_network_random_sus
    ),
    exclusion_info = list(
      affected_trails = affected_trails,
      rejoin = rejoin
    )
  ))
  
}




# Network statistics ------------------------------------------------------

calculate_network_stats <- function(network, timestep) {
  # Calculates network measures for a given timestep
  
  nest_list <- network$nest_list |> 
    filter(t==timestep) |>
    relocate(t, .after = ID)
  
  tree_list <- network$tree_list |> 
    filter(t==timestep)|>
    relocate(t, .after = ID)
  
  trail_list <- network$trail_list |>
    filter(t==timestep)|>
    relocate(t, .after = to_ID)
  
  ## Calculating network stats
  # Creating igraph object for stats
  network_graph <-
    graph_from_data_frame(
      d = trail_list |>
        filter(to_ID !=
                 "No available tree",
               status == "active"),
      vertices = bind_rows(
        nest_list |> 
          filter(status == "active"),
        tree_list |> 
          filter(ID %in% (trail_list |>
                            filter(trail_type == "foraging",
                                   status == "active") |> 
                            pull(to_ID)))
      )
    )
  
  
  
  
  # Creating igraph object for stats - Nests only
  network_graph_nests <-
    graph_from_data_frame(d = trail_list |>
                            filter(trail_type == "internest",
                                   status == "active"),
                          vertices = nest_list |> 
                            filter(status == "active"))
  
  # Number of nests
  num_nests = nest_list |> 
    filter(status == "active") |> 
    nrow()
  
  # Number of trees
  num_trees = trail_list |>
    filter(trail_type == "foraging",
           status == "active",
           to_ID != "No available tree") |>
    dplyr::select(to_ID) |>
    unique() |>
    nrow()
  # Number of foraging trails
  num_foraging = trail_list |>
    filter(trail_type == "foraging",
           status == "active",
           to_ID != "No available tree") |>
    nrow()
  # Number of internest trails
  num_internest = trail_list |>
    filter(trail_type == "internest",
           status == "active") |>
    nrow()
  # Trees to nests ratio
  trees_to_nests_rat = num_trees / num_nests
  # Internest trails to nests ratio
  internest_to_nests_rat = num_internest / num_nests
  # Foraging trails to nests ratio
  foraging_to_nests_rat = num_foraging / num_nests
  # Network components - nests only
  num_components_nests =
    components(network_graph_nests)$no
  # Efficiency
  network_efficiency <-
    network_graph |> efficiency(
      type = "global",
      weights = bind_rows(trail_list |>
                            filter(to_ID !=
                                     "No available tree",
                                   status=="active")) |>
        mutate(distance = sqrt(((from_x - to_x) ^ 2
        ) + ((from_y - to_y) ^ 2
        ))) |>
        dplyr::select(distance) |>
        pull()
    )
  # Efficiency - nests only
  network_efficiency_nests <-
    network_graph_nests |> efficiency(
      type = "global",
      weights = trail_list |>
        filter(trail_type == "internest", status=="active") |>
        mutate(distance = sqrt(((from_x - to_x) ^ 2
        ) + ((from_y - to_y) ^ 2
        ))) |>
        dplyr::select(distance) |>
        pull()
    )
  # Robustness
  network_robustness <-
    calculate_robustness(network_graph)
  # Robustness - nests only
  network_robustness_nests <-
    calculate_robustness(network_graph_nests)
  # Network cost
  network_cost <- trail_list |>
    filter(to_ID !=
             "No available tree",
           status == "active") |>
    mutate(distance = sqrt(((from_x - to_x) ^ 2) + ((from_y - to_y) ^
                                                      2))) |>
    summarise(network_cost = sum(distance)) |>
    pull()
  # Network cost - nests only
  network_cost_nests <- trail_list |>
    filter(trail_type == "internest",
           status == "active") |>
    mutate(distance = sqrt(((from_x - to_x) ^ 2) + ((from_y - to_y) ^
                                                      2))) |>
    summarise(network_cost = sum(distance)) |>
    pull()
  
  
  
  # Gathering network-level stats
  networklevel_stats <-
    data.frame(
      t=timestep,
      num_nests,
      num_trees,
      num_foraging,
      num_internest,
      trees_to_nests_rat,
      internest_to_nests_rat,
      foraging_to_nests_rat,
      num_components_nests,
      network_efficiency,
      network_efficiency_nests,
      network_robustness,
      network_robustness_nests,
      network_cost,
      network_cost_nests
    )
  
  return(
    list(
      networklevel_stats = networklevel_stats
    )
  )
  
}

calculate_experiment_network_stats <-
  function(experiment,
           time_steps = (-not_estimated_parameters$burn_in_length):84 #c(0,28,56,84) 
           ) {
    # Calls the function calculating network measures for specified timesteps
    network_stats <- list()
    
    for (i in 1:length(experiment$networks)) {
      network_stats[[i]] <-
        sapply(time_steps, function(x) {
          calculate_network_stats(experiment$networks[[i]], timestep = x)
        }, simplify = TRUE) |>
        bind_rows() |>
        mutate(treatment = names(experiment$networks)[i]) |>
        relocate(treatment, .after = t)
    }
    network_stats <- bind_rows(network_stats)
    return(network_stats)
  }

