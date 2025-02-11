Sexual Selection and Population Structure Effects Analysis
This project analyzes how sexual selection and population structure affect network metrics derived from Ethovision recording data. The workflow includes data processing, network metric calculation, mixed-effects model fitting, and a randomization procedure for extracting model estimates.

Table of Contents
Overview
Directory Structure
Dependencies
How to Run
Functions Overview
Output Files
Additional Notes
License
Contact

Overview
This project uses custom R functions to process raw data, calculate network metrics, and fit statistical models. Key analyses include:

Data Processing:
Reading and merging raw data from text files, cleaning and reshaping the data, and performing pairwise comparisons of network weights.

Model Fitting:
Calculating individual and population-level network metrics, fitting mixed-effects models (using lmer or glmmTMB), and extracting beta estimates.

Randomization Tests:
Shuffling key grouping variables (e.g., sex selection and population subgroup) in a parallelized environment to generate beta estimates under random conditions. A progress bar provides live monitoring of the iterations.

Directory Structure
A suggested directory structure for the project is:

/<project_root>
├── functions/
│   ├── start_progress.R                    # Contains the start_progress function for data preparation
│   ├── load_or_install.R                   # Contains the load_or_install function for package management
│   └── random_netcal_model_betas_pro.R     # Contains the random_netcal_model_betas_pro function
├── text/
│   ├── manual_correction96.txt             # Raw network weight data from Ethovision recordings
│   └── recorddataALLCG.txt                 # Raw recording details data
├── analysis_script.R                       # Main R script that calls all functions and performs the analysis
└── README.md                               # This file

Dependencies
The analysis requires the following R packages:

Data Manipulation & Reshaping: dplyr, tidyr, purrr
Network Analysis: tnet, igraph
Statistical Modeling: lmerTest, lme4, glmmTMB
Parallel Processing & Monitoring: doParallel, foreach, pbapply, parallel

The load_or_install function automatically installs any missing packages from the CRAN mirror (set to "https://cloud.r-project.org").

How to Run
Set Up Your Working Directory:

The main script automatically sets the working directory based on your system’s home directory. The custom set_working_directory function (embedded in start_progress) creates the path using your specified folder structure (e.g., Dropbox/david_quevedo_SexNet/thesis/text for data files).

Source the Functions:

The main script uses a helper function (e.g., source_all_functions) to load all custom functions from the functions directory. Ensure that the working directory is set correctly so that these functions can be found.

Run the Main Analysis Script:

Execute analysis_script.R in R or RStudio. This script:

- Loads and processes raw data via start_progress.
- Calculates network metrics and fits mixed-effects models.
- Performs randomization iterations using random_netcal_model_betas_pro.
- Extracts and summarizes beta estimates.
- Exports the results as CSV files to the working directory.

Review the Output:

Several CSV files are generated (e.g., betas_mean_sd_indv.csv, betas_combined_indv.csv, results_df_individual.csv, etc.). These files contain the beta estimates and model fitting results.

Functions Overview
1. start_progress
Purpose:
Prepares the data for analysis by setting the working directory, reading and merging data files, cleaning and reshaping the data, and performing pairwise comparisons on network weights (both before and after matrix reshaping).

Key Steps:

- Set the working directory.
- Read data from manual_correction96.txt and recorddataALLCG.txt.
- Merge datasets and exclude specified videos.
- Reshape data for network analysis (using gather and spread).
- Compare directional pairs of network weights and report consistency.
- Return a list of processed data objects.

2. load_or_install
Purpose:
Checks if a required package is installed and loads it. If not installed, it installs the package from the specified CRAN mirror and then loads it.

3. random_netcal_model_betas_pro
Purpose:
Performs multiple iterations of a randomization procedure to shuffle key grouping variables, recalculate network metrics, fit models, and extract beta estimates.
Uses parallel processing with error handling and a progress bar.

Key Steps:

- Set up a parallel processing cluster.
- Export necessary data and functions to the worker nodes.
- Shuffle data, calculate network metrics, and fit models for each iteration.
- Extract beta estimates from individual and population-level analyses.
- Combine and return the beta estimates.

Output Files
After running the analysis, the following output files are generated in your working directory:

Beta Estimate Summaries:

- betas_mean_sd_indv.csv
- betas_mean_sd_pop.csv

Combined Beta Estimates:

- betas_combined_indv.csv
- betas_combined_pop.csv

Model Results:

- results_df_individual.csv
- results_df_populational.csv

These files provide detailed summaries and raw outputs from the network analysis and model fitting.

Additional Notes
Error Handling:
The randomization function uses tryCatch to ensure that errors in individual iterations do not halt the entire process.

Performance Considerations:
The script utilizes parallel processing to speed up the randomization iterations. Adjust the number of cores used if necessary (currently set to detectCores() - 2).

Customization:
You can adjust parameters such as the threshold for network weights (e.g., setting bcd values below 11 to 0) or modify the data shuffling procedure by editing the corresponding functions.

License
All this repository is owned by the Spanish Research Council and their affiliations. 

Contact
For questions or further assistance, please contact:

David P. Quevedo
Email: davquecol@outlook.es
