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
                    name: Coal supply
                    carrier: electricity
                    parent: supply
                constraints:
                    lifetime: 15  # [Danish energy agency, 02 LTE existing plant, 2030]
                    energy_eff: 0.4  # https://www.ge.com/power/transform/article.transform.articles.2018.mar.come-hele-or-high-water
                costs:
                    monetary:
                        om_con: {{ fuel_cost.Coal * scaling_factors.specific_costs }}  # {{ (1 / scaling_factors.specific_costs) | unit("EUR2015/MWh") }}
                        energy_cap: {{ 0.24 * 1e6 * scaling_factors.monetary / scaling_factors.power }}  # {{ (scaling_factors.power / scaling_factors.monetary) | unit("EUR/MW") }} | [Danish energy agency, 02 LTE existing plant, 2030]
                    co2.om_con: {{ 0.34 * scaling_factors.co2_cost / scaling_factors.power }}  # {{ (scaling_factors.power / scaling_factors.co2_cost) | unit("tCO2/MWh") }} | https://www.ipcc-nggip.iges.or.jp/public/2006gl/pdf/2_Volume2/V2_2_Ch2_Stationary_Combustion.pdf

        locations:
            {% for region, capacity in capacity_per_region.iteritems() %}
            {% if capacity > 0 %}
            {{ region }}:
                techs.coal_power_plant.constraints.energy_cap_max: {{ capacity * scaling_factors.power }} # {{ (1 / scaling_factors.power) | unit("MW") }}
            {% endif %}
            {% endfor %}
"""


def generate_emissions_scenarios(
    path_to_fuel_costs, path_to_power_plants, path_to_units, scaling_factors,
    fuel_cost_source, fuel_cost_year, path_to_result

):
    """Generate a file that represents links in Calliope."""
    fuel_cost_df = pd.read_excel(path_to_fuel_costs, sheet_name="Appendix 1", header=0, index_col=[0, 1, 2, 3])
    fuel_cost_series = fuel_cost_df.xs((fuel_cost_source, fuel_cost_year, "2015 REF"), level=(1, 2, 3)).loc[:, "EU"]

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
        scaling_factors=scaling_factors,
        fuel_cost=fuel_cost_series
    )
    with open(path_to_result, "w") as result_file:
        result_file.write(scenarios)


if __name__ == "__main__":
    generate_emissions_scenarios(
        path_to_fuel_costs=snakemake.input.fuel_costs,
        path_to_power_plants=snakemake.input.power_plants,
        path_to_units=snakemake.input.units,
        scaling_factors=snakemake.params.scaling_factors,
        fuel_cost_source=snakemake.params["fuel_cost_source"],
        fuel_cost_year=snakemake.params["fuel_cost_year"],
        path_to_result=snakemake.output[0],
    )
