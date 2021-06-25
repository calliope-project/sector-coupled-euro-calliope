"""Applies config parameters to biofuel template."""
import numpy as np
import pandas as pd

import jinja2

from eurocalliopelib import filters

import util

TEMPLATE = """
overrides:
    biofuel-supply:
        techs:
            biofuel_supply:
                essentials:
                    name: Biofuel
                    parent: supply
                    carrier: biofuel
                constraints:
                    lifetime: 20  # arbritrarily chosen to avoid Calliope errors
                costs.monetary:
                    om_prod: {{ cost * scaling_factors.specific_costs }} # {{ (1 / scaling_factors.specific_costs) | unit("EUR/MWh") }}

        locations:
            {% for id in potentials.index %}
            {{ id }}.techs.biofuel_supply:
            {% endfor %}

    biofuel_maximum:
        group_constraints:
            {% for id, potential in potentials.iteritems() %}
            biofuel_{{ id }}:
                locs: [{{ id }}]
                techs: [biofuel_supply]
                carrier_prod_max:
                    biofuel: {{ potential * scaling_factors.power }} # {{ (1 / scaling_factors.power) | unit("MW") }}
            {% endfor %}

scenarios:
    add-biofuel: [biofuel-supply, biofuel_maximum]
"""


def parameterise_template(
    path_to_annual_demand, path_to_biofuel_potential_mwh, path_to_biofuel_costs,
    year, scaling_factors, path_to_result
):
    """Applies config parameters to template files."""

    annual_demand = util.read_tdf(path_to_annual_demand)
    potentials = (
        pd.read_csv(path_to_biofuel_potential_mwh, index_col=0, squeeze=True)
        .fillna(0)
    )
    cost = float(pd.read_csv(path_to_biofuel_costs).columns[0])
    try:
        annual_biofuel_demand = (
            annual_demand.xs(("biofuel", year), level=("end_use", "year"))
            .groupby(level="id").sum()
            .div(scaling_factors["power"])
        )
    except KeyError:
        annual_biofuel_demand = 0
    potentials = potentials.sub(annual_biofuel_demand).clip(lower=0)
    scaling_factors["specific_costs"] = scaling_factors["monetary"] / scaling_factors["power"]

    env = jinja2.Environment(lstrip_blocks=True, trim_blocks=True)
    env.filters["unit"] = filters.unit
    env.globals["mean"] = np.mean
    rendered = env.from_string(TEMPLATE).render(
        potentials=potentials,
        scaling_factors=scaling_factors,
        cost=cost
    )
    with open(path_to_result, "w") as result_file:
        result_file.write(rendered)


if __name__ == "__main__":
    parameterise_template(
        path_to_annual_demand=snakemake.input.annual_demand,
        path_to_biofuel_potential_mwh=snakemake.input.biofuel_potential,
        path_to_biofuel_costs=snakemake.input.biofuel_costs,
        scaling_factors=snakemake.params["scaling_factors"],
        year=snakemake.params["year"],
        path_to_result=snakemake.output[0]
    )
