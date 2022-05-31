# Sector-coupled Euro-Calliope Model with 98 European regions

This snakemake workflow builds upon Euro-Calliope v1.0 with the following:

1. Representation of transport, heat, and industry sectors. These sectors lead to the addition of unique carriers and technologies, and associated myriad of input datasets.
2. Bespoke clustering of NUTS3 regions to create an intermediate resolution model of Europe, consisting of 98 regions.
3. Inclusion of Iceland in the energy system model.
4. Use of grid transfer capacities (GTCs) from the E-Highways 2050 Euroepan project to place limits on line capacity between regions.
5. Inclusion of proposed/planned new AC and DC lines connecting regions, according to ENTSOE.

The resulting energy system transmission network is:

<img src="map.png" alt="98-region European energy system model, including all possible AC and DC lines connecting regions" width="800" height="800">

The key rules of the workflow shown graphically are:

<img src="rulegraph.png" alt="Rulegraph of workflow to build a sector-coupled, 98-region European Calliope energy system model" width="800">

This worflow currently relies on a custom version of the [euro-calliope](https://github.com/brynpickering/euro-calliope/tree/sector-coupled) workflow. The correct version of this subworkflows are given as git submodules in this repository.

## Building and running the model
To run the Sector-Coupled Euro-Calliope workflow, you need to undertake the following steps:

1. Clone this repository, including the submodules: `git clone --recurse-submodules git@github.com:calliope-project/sector-coupled-euro-calliope.git`
2. Populate the different subworkflow data directories with manual data.
4. Install the base conda environment: `conda env create -f environment.yaml` + `conda activate sector-coupled-euro-calliope`.
5. Build the Calliope model definition by runnning the workflow from the top level (same place as this file can be found): `snakemake --use-conda --cores 1` (You can change the number of cores, or run using the cluster profile by adding the argument `--profile config/euler`).
If you need to rebuild something from the `euro-calliope` subworkflow, then you should `cd` into the directory and run snakemake from there: `snakemake --configfile ../config/euro-calliope-2050.yaml --use-conda --cores 1 [path/to/file]`.
1. Optimise the built model. This can be done in one of two ways:
    1. use snakemake directly (e.g. `snakemake --use-conda --profile config/euler "build/ehighways/outputs/run_2018_2H_neutral.nc"`)
    2. use separate run scripts found in the directory `run_scripts`. These scripts are decoupled from the workflow since they can take a very long time to run and can clash with snakemake's file watching processes. The run scripts are suitable for use on a high performance computer using the LSF software.


## Troubleshooting

There are a number of issues you may come across when trying to run this workflow.
Many of these issues will be solved in upcoming releases.

1. Not enough disk space/memory to run scripts.
This full workflow is designed to be run on a high performance computing cluster.
If you want to create specific files that you know should have a low resource requirement, you can ask snakemake to create that file specifically: `snakemake [arguments] "build/path/to/file"`

2. Snakemake rules not being triggered in subworkflows.
Due to an enresolved bug, you may find that from the top-level workflow, snakemake doesn't think there is anything to do in the subworkflows, even if no data exists in them.
You can solve this by going into the individual workflow directories and running them directly from there.
See the `Building and running the model` above for more information.

3. Choosing between 2030 and 2050 projection year.
The workflow has two sets of assumptions for demand and supply technology characteristics, depending on the "projection" year.
You can simply change the `projection_year` option in the main workflow configuration file to switch between the two.
If running from the main workflow, this will pass the correct information down to the `euro-calliope` sub-workflow.
If running from within the `euro-calliope` sub-workflow, you should make sure to point to the correct configuration file (`../config/euro-calliope-2030.yaml` or `../config/euro-calliope-2050.yaml`).

4. Running for a subset of countries in Europe.
This is not exactly possible in the current formulation of the workflow, since all countries are needed to enable gap filling (where some countries are used to fill in for others).
Your best bet is to instead run for all countries and then exclude certain regions when building and running the Calliope model.

5. Changing the composition of regions.
The file `statistical_units_to_ehighways_regions.csv` maps NUTS3 regions to user-defined Euro-Calliope regions.
If you would like to update the grouping of regions, you can do so here.

6. Running the exact workflow to produce results as in the accompanying publication.
There have been many updates to the workflow since the version used to produce results in [DOI:10.1016/j.joule.2022.05.009](https://doi.org/10.1016/j.joule.2022.05.009).
These updates rarely change the data, but rather the structure of the workflow and the method of running the Calliope models.
The version of the workflow that was used can be found [here](https://github.com/calliope-project/sector-coupled-euro-calliope/tree/2021-11-29).
A pre-built version of the same model can be be found on [Zenodo](https://zenodo.org/record/5774988).
To rebuild the original version of the model, we recommend you use this current version of the workflow to build the model, making changes to YAML template files to match the data used in the older version.
Major changes to data and final Calliope model structure are:
    1. The way in which synthetic fuels can be transported between regions.
    The model has been updated to allow all synthetic fuels to transported between regions freely (using the `synfuel_transmission` Calliope override).
    Previously only industry synethtic fuel demands in one region could be met by other regions (using the `industry_fuel_shared` Calliope override).
    1.  The cost of the Hydrogen-to-Methanol technology has been corrected to not double count the electricity requirement for electrolysis to produce hydrogen.
    2.  The efficiency of the Biofuel-to-Methanol has been corrected upwards based on the type of biofuel being consumed.
    3.  A Hydrogen-fuelled combined heat and power (CHP) plant is now an available technology.
    4.  The model is now designed for optimisation in version Calliope 0.6.8, with an intermediate Python script to add custom constraints.
    Previously, a fork of Calliope version 0.6.6 was required to run the model such that custom constraints could be included.
    1. Most manual data download steps have been automated, by making data available on Zenodo.
    You can emulate this by downloading any missing data from the Zenodo respoitories given in this workflow version.

## Citation

If you use this workflow or the pre-built Sector-Coupled Euro-Calliope model, please cite:

Pickering, B., Lombardi, F., Pfenninger, S., 2022. Diversity of options to eliminate fossil fuels and reach carbon neutrality across the entire European energy system. Joule. [DOI:10.1016/j.joule.2022.05.009](https://doi.org/10.1016/j.joule.2022.05.009)
