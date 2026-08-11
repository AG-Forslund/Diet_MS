# Diet_MS
Base code for this publication: "Functional microbiome shifts are associated with neuroinflammatory outcomes under dietary interventions in multiple sclerosis"

## Repository Structure
This repository contains two folders, each with its own `renv.lock` file, since the scripts were run on two separate machines with different environments.

- **Metabolomics_Metadata_Creation**: Contains the scripts used to prepare and analyze metadata (e.g., age, medication) and to preprocess the metabolomics files.
- **Main_Analysis**: Contains the core analysis scripts, including most of the statistical analyses presented in the publication.

Both folders follow the same structure: the `renv.lock` file sits at the project root, alongside the `scripts`, `figures`, and `data` subfolders. The `scripts` subfolder contains all the scripts provided here.


| Script | Folder | Purpose |
|---|---|---|
| `Analysis_MetabolomicsRatios.Rmd` | Metabolomics_Metadata_Creation | XY |
| `Baseline_Table_Hardcoded.Rmd` | Metabolomics_Metadata_Creation | XY |
| `Metadata_Preparation.Rmd` | Metabolomics_Metadata_Creation | XY |
| `Quality_Control_Metabolites.Rmd` | Metabolomics_Metadata_Creation | Runs primary statistical models |
| `Functions_Metabolites.R` | Metabolomics_Metadata_Creation | All analysis functions needed for the .Rmd scripts |
| `Merge_dfs_serum.R` | Metabolomics_Metadata_Creation | Merges serum datasets after their preparation and quality control |
| `Merge_dfs_stool.R` | Metabolomics_Metadata_Creation | Merges stool datasets after their preparation and quality control |
| `Merge_TrpSCFA.R` | Metabolomics_Metadata_Creation | Merges tryptophan and SCFA datasets after their preparation and quality control |
| `Packages_Metabolites.R` | Metabolomics_Metadata_Creation | All packages needed for the .Rmd scripts |
| `Plots_Metabolites.R` | Metabolomics_Metadata_Creation | All plot functions needed for the .Rmd scripts |
| `QCRLSC_only.R` | Metabolomics_Metadata_Creation | Performs QCRLSC |
||||
| `Integrated_Analysis_NAMS.Rmd` | Main_Analysis | XY |


### Restoring the environments
Each folder uses [`renv`](https://rstudio.github.io/renv/) to manage its own R package versions. To reproduce the exact environment used for a given folder run `renv::restore()`. This will install the exact package versions recorded in that folder's `renv.lock` file.

Repeat this step separately for each folder, as they rely on different `renv.lock` files.

## Contact
The underlying data can be accessed upson reasonable request to friederike.gutmann@mdc-berlin.de.
