import jinja2
import pandas as pd

from eurocalliopelib import filters

import util

TEMPLATE = """
overrides:
    {% for scenario in ["2030_current", "2030_neutral", "2050_current", "2050_neutral"] %}
    {{ scenario }}:
        group_constraints:
            {% for row_id, row in emissions_targets.iterrows() %}
            systemwide_co2_max_{{ row_id }}:
                locs: {{ per_target_regions[row_id] }}
                cost_max.co2: {{ (row[starting_point] - less_coal[row_id]) * (1 - row[scenario]) * 1e6 * scaling_factors.co2_cost }}  # {{ (1 / scaling_factors.co2_cost) | unit("tCO2") }}
            {% endfor %}
    {% endfor %}
"""


def generate_emissions_scenarios(path_to_emissions_targets, path_to_regions, path_to_annual_demand, scaling_factors, year, projection_year, path_to_result):
    """Generate a file that represents links in Calliope."""
    emissions_targets = pd.read_csv(path_to_emissions_targets, header=0)
    regions = pd.read_csv(path_to_regions, header=0, index_col=0, squeeze=True)
    annual_demand = util.read_tdf(path_to_annual_demand).xs(year, level="year")
    per_target_regions = {}
    less_coal = {}
    if projection_year == "current":
        starting_point = "1990_energy_mtCO2eq"
    elif projection_year in ["2050", 2050]:
        starting_point = "1990_energy_steel_chemical_mtCO2eq"
    for _idx in emissions_targets.index:
        per_target_regions[_idx] = [
            i for i in regions.index
            if regions.loc[i] in emissions_targets.loc[_idx, "region"].split(",")
        ]
        try:
            less_coal[_idx] = (
                annual_demand.xs(
                    ("industry_demand", "industry", "coal"),
                    level=("dataset", "cat_name", "end_use")
                )
                .droplevel("unit")
                .loc[per_target_regions[_idx]]
                .mul(0.034)  # emissions factor MtCO2/0.1TWh
                .sum()
            )
        except KeyError:
            less_coal[_idx] = 0
    env = jinja2.Environment(lstrip_blocks=True, trim_blocks=True)
    env.filters["unit"] = filters.unit

    scenarios = env.from_string(TEMPLATE).render(
        emissions_targets=emissions_targets,
        scaling_factors=scaling_factors,
        per_target_regions=per_target_regions,
        starting_point=starting_point,
        less_coal=less_coal
    )
    with open(path_to_result, "w") as result_file:
        result_file.write(scenarios)


if __name__ == "__main__":
    generate_emissions_scenarios(
        path_to_emissions_targets=snakemake.input.emissions_targets,
        path_to_regions=snakemake.input.regions,
        path_to_annual_demand=snakemake.input.annual_demand,
        scaling_factors=snakemake.params.scaling_factors,
        year=snakemake.params.year,
        projection_year=snakemake.params.projection_year,
        path_to_result=snakemake.output[0],
    )
