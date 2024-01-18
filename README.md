# Strong and weak environmental perturbations cause contrasting restructure of ant transportation networks

This project contains the code base and data used in Piross et al. "Strong and weak environmental perturbations cause contrasting restructure of ant transportation networks".

The study aims to improve our understanding of the structural changes happening in multi-sink multi-source transport networks in response to environmental perturbations. We designed the model to accurately describe the structure, dynamics and perturbation responses of wood ant (_Formica lugubris_) colony networks.

We provide a complete, detailed model description, following the ODD (Overview, Design concepts, Details) protocol [1–3] in the Supplementary Material document.

The project is implemented in R 4.2.3. See Supplementary Methods S1.4 and sessionInfo.txt for more details. Execute Run_simulations.R to run simulation experiments. The results will be saved in "./Output/" in a folder named by the date and time when the simulations started. 

# List of files

- Input
	- Empirical data
		- empirical_edges.csv - Empirical data: network edges
		- empirical_nodes - Empirical data: network nodes
	- Empirical_network_statistics.R - Calculation of network measures on the empirical networks
	- Model_parameters.R - Estimates empirical and defines non-estimated parameters

- Output
	- Empirical_statistics - The calculated empirical network measures
		- empirical_stats.csv - Empirical network measures
		- empirical_stats.rds - Empirical network measures
	- Model_parameters
		- estimated_parameters.csv - List of parameters estimated from the empirical data
		- estimated_parameters.rds - List of parameters estimated from the empirical data
		- gravity_parameters.csv - Estimated parameters for the gravity model
		- gravity_parameters.rds - Estimated parameters for the gravity model
		- not_estimated_parameters.csv - List of not estimated parameters
		- not_estimated_parameters.rds - List of not estimated parameters

- Supplementary_Methods.pdf - Detailed model description
- Supplementary_Methods.Rmd - Rmarkdown files preparing Supplementary Methods S1
- Supplementary_Results.Rmd - Rmarkdown files preparing Supplementary Methods S2
- Supplementary_references.bib - BibLatex file with the references for the Supplementary Materials
- Validation_report_template.Rmd - Rmarkdown template called by Run_simulations.R to prepare the validation report
- Experiment_report_template.Rmd - Rmarkdown template called by Run_simulations.R to prepare individual reports on the experiments
- List_of_parameters.Rmd - Rmarkdown template called by Run_simulations.R to prepare a list of the used parameters
- Model_code.R - The codebase for running simulation experiments
- Run_simulations.R - Runs the exclusion experiments, prepares the validation and result reports, and generates individual reports for the experiments
- sessionInfo.txt - Information on R and package versions
- sessionInfo_all.txt - Detailed information on packages
- Simulation_results_template.Rmd - Rmarkdown template called by Run_simulations.R to prepare the results report

The data of the 10000 simulation experiments used for the publication and to produce Supplementary_Methods.pdf and Supplementary_Results.pdf are available upon request (3.53 GB).

# References

1. Grimm V _et al._ 2006 A standard protocol for describing individual-based and agent-based models. _Ecol. Model._ **198**, 115–126. (doi:10.1016/j.ecolmodel.2006.04.023)

2. Grimm V, Berger U, DeAngelis DL, Polhill JG, Giske J, Railsback SF. 2010 The ODD protocol: A review and first update. _Ecol. Model._ **221**, 2760–2768. (doi:10.1016/j.ecolmodel.2010.08.019)

3. Grimm V _et al._ 2020 The ODD protocol for describing agent-based and other simulation models: A second update to improve clarity, replication, and structural realism. _J. Artif. Soc. Soc. Simul._ **23**. (doi:10.18564/jasss.4259)
