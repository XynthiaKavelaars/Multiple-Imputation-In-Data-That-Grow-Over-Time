# Multiple-Imputation-In-Data-That-Grow-Over-Time

This folder contains all necessary files to replicate the simulation study for the manuscript ”Multiple imputation in data that grow over time: A comparison of three strategies." by Xynthia Kavelaars, Joost van Ginkel, and Stef van Buuren

Description of files and folders:
########################################################################################################
#					Simulation - Linear regression				       #
########################################################################################################
\Simulation - Linear regression		Folder containing subfolders with scripts, functions, simulation conditions, and the final simulated data R-workspace for the simulation of linear regression (section 4.1).

####					Simulation 1 (20% missing values)			    ####
\\Simulation 1 (20% missing values):	Folder containing subfolders with scripts, functions, simulation conditions, and the final simulated data R-workspace for 20% missingness.
0. Simulation_setup.R			R-script to run a pilot study to obtain the required sample size and number of samples
1. Execute_simulation.R			R-script to run the simulation
2. Make_table_simulation.R		R-script to extract the tabulated and in-text data as presented in the manuscript
3. Make_figure_simulation.R		R-script to create plots of simulation results presented in manuscript

\\\Functions:				Folder containing functions used by "1. Execute_simulation.R"
Condition.R				Function to write files for each condition
Chol.R					Function to 
Df.R					Function to calculate degrees of freedom according to pooling rules
GetVariableDefinitions.R 		R-script containing variable definitions
Impute.R				Functions to impute incomplete data using three strategies
Resample.R				Function to repeat the simulation and extract results
Sigma.R					Function to check and correct positive definiteness of the correlation matrix
Simulate.R				Function to simulate data
	
\\\Conditions:				Folder for saving simulation conditions as obtained by running "1. Execute_simulation.R"

\\\Plots:				Folder for saving plots as obtained by running "3. Make_figure_simulation.R"
	
\\\Workspaces:				Folder containing final simulated data R-workspace as obtained by running "1. Execute_simulation.R"
Simulation_design.Rdata			Workspace containing output of "0. Simulation_setup.R"
Simulation.Rdata			Workspace containing output of "1. Execute_simulation.R"


####					Simulation 2 (50% missing values)			    ####
\\Simulation 2 (50% missing values):	Folder containing subfolders with scripts, functions, simulation conditions, and the final simulated data R-workspace for 50% missing values.
0. Simulation_setup.R			R-script to run a pilot study to obtain the required sample size and number of samples
1. Execute_simulation.R			R-script to run the simulation
2. Make_table_simulation.R		R-script to extract the tabulated and in-text data as presented in the manuscript
3. Make_figure_simulation.R		R-script to create plots of simulation results presented in manuscript

\\\Functions:				Folder containing functions used by "1. Execute_simulation.R"
Condition.R				Function to write files for each condition
Df.R					Function to calculate degrees of freedom according to pooling rules
Impute.R				Functions to impute incomplete data using three strategies
GetVariableDefinitions.R 		R-script containing variable definitions
Resample.R				Function to repeat the simulation and extract results
Sigma.R					Function to check and correct positive definiteness of the correlation matrix
Simulate.R				Function to simulate data

\\\Conditions:				Folder for saving simulation conditions as obtained by running "1. Execute_simulation.R"

\\\Plots:				Folder for saving plots as obtained by running "3. Make_figure_simulation.R"
		
\\\Workspaces:				Folder containing final simulated data R-workspace as obtained by running "1. Execute_simulation.R"
Simulation.Rdata			Workspace containing output of "1. Execute_simulation.R"

########################################################################################################
#					Simulation - Growth model				       #
########################################################################################################
\Simulation - Growth model		Folder containing subfolders with scripts, functions, and the final simulated data R-workspace of the growth model (Section 4.2).
1. Simulate.R				R-script to run the simulation
2. Evaluate.R				R-script to extract the in-text data and figure as presented in the manuscript

\\Functions:				Folder containing functions used by "1. Execute_simulation.R"
Functions.R				Functions used for simulation
Functions_evaluate.R			Functions used for evaluation of simulation output
VariableDefinitions.R 			R-script containing variable definitions

\\Workspaces:				Folder containing final simulated data R-workspace as obtained by running "1. Simulate.R"


########################################################################################################
#					Application 						       #
########################################################################################################
\Application:				Folder containing subfolders with scripts, functions, and the final simulated data R-workspace for application to the POPS dataset
0. Problem illustration.R		R-script to extract results in problem illustration
1. Execute_application.R		R-script to run the application
2. Make_table_application.R		R-script to extract the tabulated and in-text data as presented in the manuscript
POPS_data.csv				POPS dataset [[NOT INCLUDED DUE TO PRIVACY CONSIDERATIONS]]

\\Functions:				Folder containing functions used by "1. Execute_application.R"
Df.R					Function to calculate degrees of freedom according to pooling rules

\\Workspaces:				Folder containing final simulated data R-workspace as obtained by running "1. Execute_application.R"
Workspace_POPS.Rdata			Workspace containing output of 1. Execute_application.R [[NOT INCLUDED DUE TO PRIVACY CONSIDERATIONS]]

########################################################################################################
#					Important 						       #
########################################################################################################
Make sure the working directory is set to folder "Simulation 1 20% missing values" or "Simulation 2 50% missing values" when running scripts from the "Simulation - Linear regression" 
folder. Make sure the working directory is set to folder "Simulation - Growth model" when running scripts from the "Simulation - Growth model" 
folder. The working directory must be set to folder "Application" when running scripts from the
"Application" folder.

For any help with the files in this archive, please contact Xynthia Kavelaars (xynthia.kavelaars@ou.nl). 
