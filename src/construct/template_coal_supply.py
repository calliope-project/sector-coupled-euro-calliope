import jinja2
import pandas as pd
import geopandas as gpd
import shapely.geometry

from eurocalliopelib import filters

TEMPLATE = """
overrides:
    coal_supply:
        techs:
            coal_power_plant:
                essentials:
                    name: Coal-fired power plant
                    carrier_in: coal
                    carrier_out: electricity
                    parent: conversion
                constraints:
                    lifetime: 15  # [Danish energy agency, 02 LTE existing plant, 2030]
                    energy_eff: 0.4  # https://www.ge.com/power/transform/article.transform.articles.2018.mar.come-hele-or-high-water
                costs:
                    monetary:
                        energy_cap: {{ 0.24 * 1e6 * scaling_factors.monetary / scaling_factors.power }}  # {{ (scaling_factors.power / scaling_factors.monetary) | unit("EUR/MW") }} | [Danish energy agency, 02 LTE existing plant, 2030]

        locations:
            {% for region, capacity in capacity_per_region.iteritems() %}
            {% if capacity > 0 %}
            {{ region }}:
                techs.coal_power_plant.constraints.energy_cap_max: {{ capacity * scaling_factors.power }} # {{ (1 / scaling_factors.power) | unit("MW") }}
            {% endif %}
            {% endfor %}
"""


def generate_emissions_scenarios(
    path_to_power_plants, path_to_units, scaling_factors, path_to_result

):
    """Generate a file that represents links in Calliope."""

    scaling_factors["specific_costs"] = scaling_factors["monetary"] / scaling_factors["power"]

    power_plant_db = pd.read_csv(path_to_power_plants, header=0)
    units = gpd.read_file(path_to_units)
    capacity_current = power_plant_db[
        (power_plant_db.type_g.str.lower().str.find("coal") > -1) &
        power_plant_db.year_decommissioned.isna()
    ].dropna(subset=["lat", "lon"])

    capacity_current_points = [
        shapely.geometry.Point(xy)
        for xy in zip(capacity_current.lon, capacity_current.lat)
    ]
    capacity_current_gdf = gpd.GeoDataFrame(
        capacity_current, geometry=capacity_current_points, crs='epsg:4326'
    )
    capacity_per_region = (
        gpd.overlay(capacity_current_gdf, units)
        .groupby('id').sum()
        .capacity_g
    )

    env = jinja2.Environment(lstrip_blocks=True, trim_blocks=True)
    env.filters["unit"] = filters.unit

    scenarios = env.from_string(TEMPLATE).render(
        capacity_per_region=capacity_per_region,
        scaling_factors=scaling_factors
    )
    with open(path_to_result, "w") as result_file:
        result_file.write(scenarios)


if __name__ == "__main__":
    generate_emissions_scenarios(
        path_to_power_plants=snakemake.input.power_plants,
        path_to_units=snakemake.input.units,
        scaling_factors=snakemake.params.scaling_factors,
        path_to_result=snakemake.output[0],
    )
