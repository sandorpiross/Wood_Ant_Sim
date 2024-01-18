# Author: Sandor Piross - sandor.piross@gmail.com 

# The script runs the exclusion experiments prepares the validation and result
# reports and generates individual report for the experiments

rm(list = ls())

library(doParallel)
library(parallel)
library(rmarkdown)
library(spatstat)
library(tidyverse)
library(igraph)
library(brainGraph)
library(lazyeval)
#library(fitdistrplus)
library(readxl)
library(knitr)
library(patchwork)
library(ggrepel)
library(RColorBrewer)
library(flexdashboard)

# Loading simulation code
source("Model_code.R")

# Number of simulations to run
number_of_simulations <- 8

# Creating a directory for the outputs
new_dir <- paste("Simulation run", gsub(":", "_", Sys.time()))
dir.create(paste0("./Output/", new_dir), showWarnings = FALSE)
dir.create(paste0("./Output/", new_dir, "/Network objects"),
           showWarnings = FALSE)
dir.create(paste0("./Output/", new_dir, "/Network statistics"),
           showWarnings = FALSE)
dir.create(paste0("./Output/", new_dir, "/Network reports"),
           showWarnings = FALSE)

# Exporting session information
sink(file = paste0("./Output/", new_dir, "/sessionInfo.txt"))
Sys.time()
sessionInfo()
sink()
write_lines(print(sessionInfo()),
            file = paste0("./Output/", new_dir, "/sessionInfo_all.txt"))

# Exporting parameter list
render(
  input = "./List_of_parameters.Rmd",
  output_format = "pdf_document",
  output_file = paste0("List_of_parameters.pdf"),
  output_dir = paste0("./Output/", new_dir),
  params = list(
    estimated_parameters = estimated_parameters,
    not_estimated_parameters = not_estimated_parameters,
    gravity_parameters = gravity_parameters,
    new_dir = new_dir
  ),
  envir = parent.frame()
)

# Parallel computing settings
no_cores <- detectCores(logical = TRUE)
cl <- makeCluster(no_cores,
                  outfile = paste0("./Output/", new_dir, "/Cluster_outfile.txt"))
registerDoParallel(cl)

# Running simulations
run_exclusion_experiment_parallel <- function(x, new_dir = new_dir) {
  # Runs exclusion experiments on parallel cores
  print(paste("Experiment", x))
  
  options(dplyr.summarise.inform = FALSE)
  
  # Running an experiment
  experiment <- run_exclusion_experiment()
  
  # Saving experiment and statistics
  saveRDS(experiment,
          file = paste0("./Output/", new_dir,
                        "/Network objects/experiment_", x, ".rds"
          ))

}

calculate_experiment_network_stats_parallel <- function(x, new_dir = new_dir) {
  # Calculates network measures on parallel cores
  print(paste("Experiment", x))
  
  options(dplyr.summarise.inform = FALSE)
  
  # Reading experiment
  experiment <- readRDS(
          file = paste0("./Output/", new_dir,
                        "/Network objects/experiment_", x, ".rds"
          ))
  
  # Calculating network statistics
  experiment_network_stats <-
    calculate_experiment_network_stats(experiment)
  
  
  saveRDS(experiment_network_stats,
          file = paste0("./Output/", new_dir,
                        "/Network statistics/network_stats_", x, ".rds"
          ))
  
}


# Loading packages to cluster
clusterExport(varlist = ls(), cl = cl)

clusterEvalQ(cl, library(tidyverse))
clusterEvalQ(cl, library(ggrepel))
clusterEvalQ(cl, library(patchwork))
clusterEvalQ(cl, library(spatstat))
clusterEvalQ(cl, library(igraph))
clusterEvalQ(cl, library(brainGraph))
clusterEvalQ(cl, library(RColorBrewer))
clusterEvalQ(cl, library(rlang))
clusterEvalQ(cl, library(Matrix))

# Starting simulations
start_from <- 1
(start <- Sys.time())
clusterMap(
  cl = cl,

  fun = function(x, new_dir)
    tryCatch(
      run_exclusion_experiment_parallel(x, new_dir),
      error = function(e)
        e
    ),
  x = start_from:number_of_simulations,
  new_dir = new_dir
)

# Calculating network statistics
start_from <- 1
(start <- Sys.time())
clusterMap(
  cl = cl,

  fun = function(x, new_dir)
    tryCatch(
      calculate_experiment_network_stats_parallel(x, new_dir),
      error = function(e)
        e
    ),
  x = start_from:number_of_simulations,
  new_dir = new_dir
)

stopCluster(cl = cl)

# Started at
start
# End of running simulations
(end <- Sys.time())

# Gathering simulation statistics
simulation_stats <-
  lapply(list.files(paste0(
    "./Output/", new_dir, "/Network statistics"
  )),
  function(x) {
    readRDS(paste0("./Output/", new_dir, "/Network statistics/",
                   x)) |>
      mutate(dataset = "Simulation",
             colony = x) |>
      relocate(dataset, colony)
  }) |> bind_rows()

# Generating validation report
render(
  input = "./Validation_report_template.Rmd",
  output_format = "html_document",
  output_file = paste0("Validation_report.html"),
  output_dir = paste0("./Output/", new_dir),
  params = list(simulation_stats = simulation_stats),
  envir = parent.frame()
)

# Generating results report
render(
  input = "./Simulation_results_template.Rmd",
  output_format = "html_document",
  output_file = paste0("Simulation_results.html"),
  output_dir = paste0("./Output/", new_dir),
  params = list(simulation_stats = simulation_stats),
  envir = parent.frame()
)

# Generating individual reports

experiments <- c(1:2) # Experiment indices

lapply(experiments,
       function(x){
         focal_experiment <- 
           readRDS(
             file = paste0("./Output/", new_dir,
                           "/Network objects/experiment_", x, ".rds"
             ))
         
         focal_experiment_network_stats <- 
           readRDS(file = 
                     paste0("./Output/",new_dir,"/Network statistics/network_stats_",x,".rds"
                     ))
         
         rmarkdown::render(
           input = "./Experiment_report_template.Rmd",
           output_format = "html_document",
           output_file = paste0("Report_", x, ".html"),
           output_dir = paste0("./Output/", new_dir, "/Network reports"),
           params = list(
             experiment = focal_experiment,
             experiment_network_stats = focal_experiment_network_stats 
           ),
           envir = parent.frame()
         )
       })

