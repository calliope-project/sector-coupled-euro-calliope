import jinja2
import pandas as pd

import util


TEMPLATE = """
overrides:
    heat_storage_max:
        locations:
        {% for idx in storage_requirement.index %}
            {{ idx }}.techs:
                heat_storage_small.constraints.storage_cap_max: {{ storage_requirement[idx] }}
                heat_storage_big.constraints.storage_cap_max: {{ storage_requirement[idx] }}
        {% endfor %}
    heat_techs:
        locations:
        {% for idx in storage_requirement.index %}
            {{ idx }}.techs:
                heat_storage_small:
                heat_storage_big:
                methane_boiler:
                biofuel_boiler:
                ashp:
                gshp:
                solar_thermal_collector:
                solar_thermal_energy:
                electric_heater:
                electric_hob:
                gas_hob:
                chp_biofuel:
                chp_biofuel_extraction:
                chp_methane_extraction:
                chp_methane_back_pressure_simple:
                chp_methane_back_pressure_combined:
                building_heat_to_storage:
                district_to_building_heat:
                demand_space_heat:
                demand_water_heat:
                demand_cooking:
        {% endfor %}

    heat_tech_grouping:
        group_constraints:
            cooking_demand_share:
                techs: [electric_hob, gas_hob]
                demand_share_per_timestep_decision.cooking: 1

            heating_demand_share:
                techs: [biofuel_boiler, methane_boiler, ashp, gshp, solar_thermal_collector, district_to_building_heat]  # TODO: add in building_heat_to_storage?
                demand_share_per_timestep_decision:
                    space_heat: 1
                    water_heat: 1

            district_heat_share:
                techs: [chp_biofuel, chp_biofuel_extraction, chp_methane_extraction, chp_methane_back_pressure_simple, chp_methane_back_pressure_combined]
                demand_share_per_timestep_decision:
                    district_heat: 1


scenarios:
    heat: [heat_storage_max, heat_techs, heat_tech_grouping]
"""


def fill_constraint_template(path_to_space_demand, path_to_water_demand, path_to_result, model_year, period):
    """Generate a file that represents links in Calliope."""
    space_heat = pd.read_csv(path_to_space_demand, parse_dates=True, index_col=0)
    water_heat = pd.read_csv(path_to_water_demand, parse_dates=True, index_col=0)
    storage_requirement = space_heat.add(water_heat).rolling(window=period).sum().abs().max()

    template = jinja2.Environment(lstrip_blocks=True, trim_blocks=True).from_string(TEMPLATE)
    heat_constraints = template.render(
        storage_requirement=storage_requirement
    )
    with open(path_to_result, "w") as result_file:
        result_file.write(heat_constraints)


if __name__ == "__main__":
    fill_constraint_template(
        path_to_space_demand=snakemake.input.space_heat_demand,
        path_to_water_demand=snakemake.input.water_heat_demand,
        path_to_result=snakemake.output[0],
        model_year=snakemake.params.model_year,
        period=snakemake.params.storage_period
    )
