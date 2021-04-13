Supporting analyses for the paper "Characterization and detection of  bacterial meningitis localized epidemics in the African meningitis belt using high spatial resolution  surveillance data"
==============

#+date: April 2021

Here are the scripts used to derive the figures of the paper "Characterization and detection of  bacterial meningitis localized epidemics in the African meningitis belt using high spatial resolution  surveillance data" (manuscript in preparation). The aim is to make sure that the analyses are fully reproducible, that is:
- The full dataset will be available upon request. They remain the property of the Direction de la Lutte contre la Maladie, Burkina Faso
- The preprocessing and analysis scripts are made public here

This document is organized as follows:

1. Organization of the repository: where to find the analyses, the scripts, etc.
2. Processed data: we introduce the intermediate analysis files that were derived and then used for the analysis
3. Figures: we introduce the rationale and the key steps of the analysis behind each figure (this do not prevent you from reading the code)
4. Additional analyses: additional analysis that might be of interest, even if they didn't end up in the paper.


These scripts are released under the GNU General Public License v3 or newer (GPLv3+). This means that you have the right to use, to study and modify the code and to share the modifications. Please keep a mention of the initial author.

Also, the scripts are written in `R`. Sometimes the code might be hard to read/hard to follow, feel free to contact the authors for clarification. 

## Organization of the repository
### File architecture
*Disclaimer:* this repository do not contain the datasets. These datasets are available upon request at the DLM or from the corresponding author.

The repository is organized in several folders and subfolders:
1. `scripts`: this folder contains only the scripts to generate the data. These are the same scripts as the ones presented in this [other Github repository](https://github.com/MaximeMaW/MeningitisDustDynamics), from Woringer et al, 2018.
2. `figures` this folder contains the figures that were included in the manuscript. The output might not be exactly the same as the ones you obtain when running the script. One of the explanations is that sometimes, the figures have been edited with the Inkscape software (to adjust the layout, add arrows, change colours, etc).
3. `README.md`: this instructions file

The analyses are split by "categories":

1. `cluster_analysis`: about the derivation of a fine characterization of single epidemic clusters
2. `cluster_example`: about the computation of the the spatio-temporal K-Ripley functions, and the resamplings to seek for a null and a generative model.
3. `data_completeness`: provides various statistics about the quality of the dataset, including the validation with respect to other sources of data.
4. `data_statistics`: some statistics about the datasets

In the `scripts/` folder are the `.R` scripts and in the `figures/` folder the `.pdf` exports, and if needed, the `.svg` that creates the final figure.

| Fig. | Description                           | Folder              | Script                    |
|------+---------------------------------------+---------------------+---------------------------|
|    S1 | Completeness of the dataset           | `data_completeness` | `fig_data_completeness.R` |
|    1 | Temporal and spatial dynamics of meningitis clusters| `cluster_analysis`   | `fig_cluster_analysis.R` |
|    2 | K-Ripley vs. null model               | `cluster_example`   | `fig_clustering.R`        |
|    3 | K-Ripley vs. generative model         | `cluster_example`   | `fig_clustering.R`        |
|   S1 | Geolocalization of the HC             | `data_completeness` | `map_hc.R`                |


For some figures, I haven't taken the time to format a proper working code. Please contact us if you need it.
Some extra numbers mentioned in the text are generated in the `data_statistics/` folder.

## Running the scripts
All the scripts depend on a toolbox called `utilitaires.R` (included in the `scripts` folder). Make sure that `R` can access it. In many cases, one script generates several figures. It might be possible to run the chunks separately if you are only interested in one figure or in one panel (but be careful with side-effects). Also, sometimes you have the choice to either recompute a heavy analysis or to directly load a preprocessed version. The choice is usually performed by switching a flag to `TRUE` or `FALSE`.

The scripts have various dependencies. Make sure that they are installed and that you have placed the datasets at the right locations (which is usually NOT next to the analysis script) before you start. Dependencies include (but not only):
- `reshape`
- `zoo`
- `RColorBrewer`
- `classInt`
- `rgdal`
- `rgeos`
- `maptools`
- `viridis`
- `fields`
- `scales`
- `spatstat`
- `stpp` # Should be in version 1.0
- `ISOweek`
- `ncdf4`
- `reshape2`
- `parallel`
- `stringr`

On the system, the following softwares might be needed: BWidget, GDAL, OpenSSL: `apt install bwidget libgdal-dev libssl-dev` (this command-line applies to Debian-like systems)

### Installation of `stpp`, version < 2.0
To work, the scripts require the library `stpp` in version 1.0. At the time of this update, version 2.0 has been released and is the default version installed. To install an older version, you need to proceed as follow (inspired from [this page](https://support.rstudio.com/hc/en-us/articles/219949047-Installing-older-versions-of-packages) :

```
install.packages('devtools')
require(devtools)
install_version("stpp", version = "1.0-5") # 1.0-5 is the latest compatible version
```

## Processed datasets
These scripts depend on various sources of data that are described below:

## Figures

## Additional analyses
