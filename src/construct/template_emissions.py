import jinja2
import pandas as pd

from eurocalliopelib import filters

TEMPLATE = """
overrides:
    {% for scenario in ["2030_current", "2030_neutral", "2050_current", "2050_neutral"] %}
    {{ scenario }}:
        group_constraints:
            {% for row_id, row in emissions_targets.iterrows() %}
            systemwide_co2_max_{{ row_id }}:
                locs: {{ per_target_regions[row_id] }}
                cost_max.co2: {{ row["1990_mtCO2eq"] * (1 - row[scenario]) * 1e6 * scaling_factors.co2 }}  # {{ (1 / scaling_factors.co2) | unit("tCO2") }}
            {% endfor %}
    {% endfor %}
"""


def generate_emissions_scenarios(path_to_emissions_targets, scaling_factors, path_to_regions, path_to_result):
    """Generate a file that represents links in Calliope."""
    emissions_targets = pd.read_csv(path_to_emissions_targets, header=0)
    regions = pd.read_csv(path_to_regions, header=0, index_col=0, squeeze=True)
    per_target_regions = {}
    for _idx in emissions_targets.index:
        per_target_regions[_idx] = [
            i for i in regions.index
            if regions.loc[i] in emissions_targets.loc[_idx, "region"].split(",")
        ]
    env = jinja2.Environment(lstrip_blocks=True, trim_blocks=True)
    env.filters["unit"] = filters.unit

    scenarios = env.from_string(TEMPLATE).render(
        emissions_targets=emissions_targets,
        scaling_factors=scaling_factors,
        per_target_regions=per_target_regions
    )
    with open(path_to_result, "w") as result_file:
        result_file.write(scenarios)


if __name__ == "__main__":
    generate_emissions_scenarios(
        path_to_emissions_targets=snakemake.input.emissions_targets,
        path_to_regions=snakemake.input.regions,
        path_to_result=snakemake.output[0],
        scaling_factors=snakemake.params.scaling_factors
    )
