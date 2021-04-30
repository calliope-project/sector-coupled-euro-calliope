# coding: utf-8
import calliope
from calliope.core.util.dataset import split_loc_techs
import pandas as pd
import numpy as np

m = calliope.read_netcdf('../outputs/resolution_test/link_cap_dynamic.nc')

prod = split_loc_techs(m._model_data.carrier_prod.sum("timesteps"), return_as="Series").where(lambda x: x > 1e-5).dropna().sum(level=["techs", "carriers"])
prod = prod.unstack('carriers')
prod.index = prod.index.str.split(":", expand=True)
prod = prod.groupby(level=0).sum(min_count=1)
prod = prod.stack().rename_axis(index=['techs', 'carrier_out'])

con = split_loc_techs(m._model_data.carrier_con.sum("timesteps"), return_as="Series").abs().where(lambda x: x > 1e-5).dropna().sum(level=["techs", "carriers"])
con = con.unstack('carriers')
con.index = con.index.str.split(":", expand=True)
con.groupby(level=0).sum()
con = con.groupby(level=0).sum(min_count=1)
con = con.stack().rename_axis(index=['techs', 'carrier_in'])

energy_flow = (
    pd.concat(prod.align(con), keys=['energy_out', 'energy_in'], axis=1)
    .div(10)
    .assign(unit_in='twh_per_year', unit_out='twh_per_year')
    .set_index(['unit_out', 'unit_in'], append=True)
    .reset_index()
)

energy_flow.loc[energy_flow.techs.str.find('transport_') > -1, 'unit_out'] = 'billion_km_per_year'
energy_flow.loc[energy_flow.techs.isin(['demand_heavy_transport', 'demand_light_transport']), 'unit_in'] = 'billion_km_per_year'
energy_flow = energy_flow.drop(energy_flow[energy_flow.techs.str.find('tech_heat') > -1].index)

for _flow in ["in", "out"]:
    energy_flow.loc[energy_flow[f"carrier_{_flow}"].isna(), f"unit_{_flow}"] = np.nan
    energy_flow.loc[energy_flow[f"carrier_{_flow}"].str.find('_heat') > -1, f'carrier_{_flow}'] = "heat"
    energy_flow.loc[energy_flow[f"carrier_{_flow}"] == 'co2', f"unit_{_flow}"] = '100kt_per_year'


energy_flow.to_csv('outputs/data_for_sentinel_partners/annual_energy_flows.csv')
