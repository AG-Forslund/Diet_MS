# Diet_MS
Base code for this publication: "Functional microbiome shifts are associated with neuroinflammatory outcomes under dietary interventions in multiple sclerosis"

## Repository Structure
This repository contains two folders, each with its own `renv.lock` file, since the scripts were run on two separate machines with different environments.

- **Metabolomics_Metadata_Creation**: Contains the scripts used to prepare and analyze metadata (e.g., age, medication) and to preprocess the metabolomics files.
- **Main_Analysis**: Contains the core analysis scripts, including most of the statistical analyses presented in the publication.

Both folders follow the same structure: the `renv.lock` file sits at the project root, alongside the `scripts`, `figures`, and `data` subfolders. The `scripts` subfolder contains all the scripts provided here.


| Script | Folder | Purpose |
|---|---|---|
| `Analysis_MetabolomicsRatios.Rmd` | Metabolomics_Metadata_Creation | Prepares the metabolomics data, incl. potential outlier filtering and option to group metabolites, e.g. along a pathway. Calculates differences over time, IDO activation, and serum-stool correlation. Formats the result for later use. Allows for basic plots. |
| `Baseline_Table_Hardcoded.Rmd` | Metabolomics_Metadata_Creation | Calcualtes correlations for both baseline tables. |
| `Metadata_Preparation.Rmd` | Metabolomics_Metadata_Creation | Prepares metadata, including filtering. Translates the medication data, and calcuates differences where applicable. Plots the correlation between BHB values from MS-lab and commercial lab. |
| `Quality_Control_Metabolites.Rmd` | Metabolomics_Metadata_Creation | Metabolomics quality control, including linearity filter, batch correction and dimensionality reduction analysis. |
| `Functions_Metabolites.R` | Metabolomics_Metadata_Creation | All analysis functions needed for the .Rmd scripts. |
| `Merge_dfs_serum.R` | Metabolomics_Metadata_Creation | Merges serum datasets after their preparation and quality control. |
| `Merge_dfs_stool.R` | Metabolomics_Metadata_Creation | Merges stool datasets after their preparation and quality control. |
| `Merge_TrpSCFA.R` | Metabolomics_Metadata_Creation | Merges tryptophan and SCFA datasets per matrix after their preparation and quality control. |
| `Packages_Metabolites.R` | Metabolomics_Metadata_Creation | All packages needed for the .Rmd scripts. |
| `Plots_Metabolites.R` | Metabolomics_Metadata_Creation | All plot functions needed for the .Rmd scripts. |
||||
| `Integrated_Analysis_NAMS.Rmd` | Main_Analysis | Correlates the microbiome with the metabolome in various ways, including interaction models and mediation, and plots the results. |
| `NAMS_Difference_Model_Analysis.Rmd` | Main_Analysis | Analyses baseline values and their relation to changes over time, including plots. |
| `NAMS_Microbiome_Analysis.Rmd` | Main_Analysis | Prepares Microbiome data, incl. assigning a taxonomic level and filtering. Analyses diversity metrics and performs outlier analysis, including imputation of diversity values. |
| `NAMS_Model_Analysis.Rmd` | Main_Analysis | Analyses level changes over time of metabolites, metabolite-over-tryptophan-ratios and microbial features, also with regard to diet differences and MS outcomes. Plots results and post hoc tests.|
| `1_Taxa_diversity_2.R` | Main_Analysis | Performs rarefaction and plots the results. Analyses diversity metrics (not used in publication). |
| `0_Taxa_diversity_2.R` | Main_Analysis | Processes sequencing results. |
| `Functions_NAMS_Difference_Analysis.R` | Main_Analysis | All analysis functions needed for the NAMS_Difference_Model_Analysis.Rmd scripts. |
| `Functions.R` | Main_Analysis | All analysis functions needed for .Rmd scripts except NAMS_Difference_Model_Analysis.Rmd. |
| `Packages.R` | Main_Analysis | All packages needed for the .Rmd scripts. |
| `Plots.R` | Main_Analysis | All packages needed for the .Rmd scripts. |


### Restoring the environments
Each folder uses [`renv`](https://rstudio.github.io/renv/) to manage its own R package versions. To reproduce the exact environment used for a given folder run `renv::restore()`. This will install the exact package versions recorded in that folder's `renv.lock` file.

Repeat this step separately for each folder, as they rely on different `renv.lock` files.

## Data availability
Individual participant data are not publicly available due to privacy and ethical considerations. However, they will be shared after de-identification and subject consent with researchers who provide a methodologically sound proposal on a related research question.

## Contact
Please contact friederike.gutmann@mdc-berlin.de for requests and further questions.
