
import xarray as xr
import pandas as pd

idx = pd.IndexSlice


def get_annual_heat_demand(
    annual_consumption, heat_tech_params, demand_out_path
):
    """
    demand = consumption of energy carrier * technology efficiency
    The only difference is for direct heat (already = demand)
    and heat pumps, because demand = ambient heat + electricity consumption
    """

    annual_consumption_df = pd.read_csv(annual_consumption)

    efficiencies = {
        'biogas': heat_tech_params['biogas_eff'],
        'biomass': heat_tech_params['biomass_eff'],
        'coal': heat_tech_params['coal_eff'],
        'direct_electric': heat_tech_params['direct_electric_eff'],
        'natural_gas': heat_tech_params['natural_gas_eff'],
        'oil': heat_tech_params['oil_eff'],
        'solar_thermal': heat_tech_params['solar_thermal_eff'],
        'heat': 1,
        # heat demand met by heat pumps = heat pump electricity + ambient heat
        'heat_pump': 1,
        'ambient_heat': 1
    }

    demand = pd.Series(index=annual_consumption.index)
    for k, v in efficiencies.items():
        demand += annual_consumption[k] * v

    return demand

    demand.to_csv(demand_out_path)


if __name__ == '__main__':
    get_annual_heat_demand(
        annual_consumption = snakemake.input.annual_consumption,
        heat_tech_params = snakemake.params.heat_tech_params,
        demand_out_path = snakemake.output.demand,
    )