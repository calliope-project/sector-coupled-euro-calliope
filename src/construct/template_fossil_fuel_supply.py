"""Applies config parameters to 2030 basis template files."""
import numpy as np
import pandas as pd

import jinja2

from eurocalliopelib import filters

TEMPLATE = """
overrides:
    fossil-fuel-supply:
        techs:  # emissions from [IPCC2006] https://www.ipcc-nggip.iges.or.jp/public/2006gl/pdf/2_Volume2/V2_2_Ch2_Stationary_Combustion.pdf
            methane_supply:
                essentials:
                    name: Natural gas
                    carrier: methane
                    parent: supply
                costs:
                    monetary:
                        om_prod: {{ fuel_cost.Gas * scaling_factors.specific_costs }}  # {{ (1 / scaling_factors.specific_costs) | unit("EUR2015/MWh") }}
                    co2:
                        om_prod: {{ 0.2 * scaling_factors.co2_cost / scaling_factors.power }}   # {{ (scaling_factors.power / scaling_factors.co2_cost) | unit("tCO2/MWh") }}
            diesel_supply:
                essentials:
                    name: Diesel
                    carrier: diesel
                    parent: supply
                costs:
                    monetary:
                        om_prod: {{ fuel_cost.Oil * scaling_factors.specific_costs }}  # {{ (1 / scaling_factors.specific_costs) | unit("EUR2015/MWh") }}
                    co2:
                        om_prod: {{ 0.26 * scaling_factors.co2_cost / scaling_factors.power }}   # {{ (scaling_factors.power / scaling_factors.co2_cost) | unit("tCO2/MWh") }}
            kerosene_supply:
                essentials:
                    name: Kerosene
                    carrier: kerosene
                    parent: supply
                costs:
                    monetary:
                        om_prod: {{ fuel_cost.Oil * scaling_factors.specific_costs }}  # {{ (1 / scaling_factors.specific_costs) | unit("EUR2015/MWh") }}
                    co2:
                        om_prod: {{ 0.27 * scaling_factors.co2_cost / scaling_factors.power }}   # {{ (scaling_factors.power / scaling_factors.co2_cost) | unit("tCO2/MWh") }}
            methanol_supply:
                essentials:
                    name: Methanol
                    carrier: methanol
                    parent: supply
                costs:
                    monetary:
                        om_prod: {{ fuel_cost.Oil * scaling_factors.specific_costs }}  # {{ (1 / scaling_factors.specific_costs) | unit("EUR2015/MWh") }}
                    co2:
                        om_prod: {{ 0.26 * scaling_factors.co2_cost / scaling_factors.power }}   # {{ (scaling_factors.power / scaling_factors.co2_cost) | unit("tCO2/MWh") }} | ASSUME: same as emissions for crude oil
            coal_supply:
                essentials:
                    name: Coal
                    carrier: coal
                    parent: supply
                costs:
                    monetary:
                        om_prod: {{ fuel_cost.Coal * scaling_factors.specific_costs }}  # {{ (1 / scaling_factors.specific_costs) | unit("EUR2015/MWh") }}
                    co2:
                        om_prod: {{ 0.34 * scaling_factors.co2_cost / scaling_factors.power }}   # {{ (scaling_factors.power / scaling_factors.co2_cost) | unit("tCO2/MWh") }}


        locations:
            {% for id in regions.index %}
            {{ id }}.techs:
                diesel_supply:
                kerosene_supply:
                methanol_supply:
                methane_supply:
                coal_supply:
            {% endfor %}
"""


def parameterise_template(path_to_fuel_costs, path_to_regions, scaling_factors, fuel_cost_source, fuel_cost_year, path_to_result):
    """Applies config parameters to template files."""

    fuel_cost_df = pd.read_excel(path_to_fuel_costs, sheet_name="Appendix 1", header=0, index_col=[0, 1, 2, 3])
    fuel_cost_series = fuel_cost_df.xs((fuel_cost_source, fuel_cost_year, "2015 REF"), level=(1, 2, 3)).loc[:, "EU"]
    assert isinstance(fuel_cost_series, pd.Series)
    scaling_factors["specific_costs"] = scaling_factors["monetary"] / scaling_factors["power"]
    regions = pd.read_csv(path_to_regions, index_col=0)

    env = jinja2.Environment(lstrip_blocks=True, trim_blocks=True)
    env.filters["unit"] = filters.unit
    env.globals["mean"] = np.mean
    rendered = env.from_string(TEMPLATE).render(
        regions=regions,
        scaling_factors=scaling_factors,
        fuel_cost=fuel_cost_series
    )
    with open(path_to_result, "w") as result_file:
        result_file.write(rendered)


if __name__ == "__main__":
    parameterise_template(
        path_to_fuel_costs=snakemake.input.fuel_costs,
        path_to_regions=snakemake.input.regions,
        scaling_factors=snakemake.params["scaling_factors"],
        fuel_cost_source=snakemake.params["fuel_cost_source"],
        fuel_cost_year=snakemake.params["fuel_cost_year"],
        path_to_result=snakemake.output[0]
    )
