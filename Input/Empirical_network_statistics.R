# Author: Sandor Piross - sandor.piross@gmail.com 

# The script calculates network measures on the empirical data

source("Model_code.R")
options(dplyr.summarise.inform = FALSE)

# Importing empirical data
empirical_edges <- 
  read.csv("./Input/Empirical_data/empirical_edges.csv")

empirical_nodes <-
  read.csv("./Input/Empirical_data/empirical_nodes.csv") 


calculate_empirical_stats <- 
  function(focal_network_ID, edge_data, node_data) {

  focal_edges <- edge_data |>
    filter(network_ID == focal_network_ID) |>
    relocate(from,to)
  
  focal_nodes <- node_data |>
    filter(network_ID == focal_network_ID) |>
    relocate(ID)
  
  ## Calculating network stats
  # Creating igraph object for stats
  network_graph <-
    graph_from_data_frame(d = focal_edges,
                          vertices = focal_nodes)
  
  # Creating igraph object for stats - Nests only
  network_graph_nests <-
    graph_from_data_frame(
      d = focal_edges |>
        filter(trail_type == "internest"),
      vertices = focal_nodes |>
        filter(species == "nest")
    )
  
  # Number of nests
  num_nests = nrow(focal_nodes |>
                     filter(species == "nest"))
  # Number of trees
  num_trees = nrow(focal_nodes |>
                     filter(species != "nest"))
  # Number of foraging trails
  num_foraging = focal_edges |>
    filter(trail_type == "foraging") |>
    nrow()
  # Number of internest trails
  num_internest = focal_edges |>
    filter(trail_type == "internest") |>
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
  #Edge density
  edge_density <- network_graph |> edge_density()
  #Edge density - nests only
  edge_density_nests <- network_graph_nests |> edge_density()
  # Efficiency
  network_efficiency <-
    network_graph |> efficiency(type = "global",
                                weights = focal_edges$distance)
  # Efficiency - nests only
  network_efficiency_nests <-
    network_graph_nests |> efficiency(
      type = "global",
      weights = focal_edges |>
        filter(trail_type == "internest") |>
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
  network_cost <- sum(focal_edges$distance, na.rm = T)
  # Network cost - nests only
  network_cost_nests <- focal_edges |>
    filter(trail_type == "internest") |>
    summarise(cost = sum(distance, na.rm = T)) |>
    pull()
  
  # Gathering network-level stats
  networklevel_stats <-
    data.frame(
      dataset = focal_nodes$dataset[1],
      colony = focal_nodes$colony[1],
      date = focal_nodes$date[1],
      num_nests,
      num_trees,
      num_foraging,
      num_internest,
      trees_to_nests_rat,
      internest_to_nests_rat,
      foraging_to_nests_rat,
      num_components_nests,
      edge_density,
      edge_density_nests,
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

# Empirical data stats -----------------------------------------

# Network list
network_ID_list_empirical <-
  empirical_edges |>
  dplyr::select(network_ID) |>
  unique() |>
  arrange(network_ID) |>
  pull()

# Calculating data set statistics
empirical_stats_all <- lapply(network_ID_list_empirical,function(x){
  calculate_empirical_stats(x,
                             edge_data=empirical_edges,
                             node_data=empirical_nodes)})

# Getting network stats
empirical_stats <-
  as.data.frame(do.call(rbind, lapply(empirical_stats_all, function(x) {
    x$networklevel_stats
  })))


# Exporting empirical statistics
write.csv(empirical_stats, file = "./Output/Empirical_statistics/empirical_stats.csv",
          row.names = FALSE)
saveRDS(empirical_stats, file = "./Output/Empirical_statistics/empirical_stats.rds")



