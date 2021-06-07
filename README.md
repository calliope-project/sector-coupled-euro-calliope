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

This worflow currently relies on updated versions of the [euro-calliope](https://github.com/brynpickering/euro-calliope/tree/sector-coupled) and [solar-and-wind-potential](https://github.com/brynpickering/possibility-for-electricity-autarky/tree/custom-regions) workflows, as well as a bespoke version of [Calliope v0.6.6](https://github.com/calliope-project/calliope/tree/euro-calliope)
