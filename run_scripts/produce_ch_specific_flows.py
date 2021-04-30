import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

import calliope
from calliope.core.util.dataset import split_loc_techs

GROUPS = {
     'ac_ohl_mountain_transmission': 'transmission',
     'ac_ohl_transmission': 'transmission',
     'battery': 'storage',
     'chp_biofuel_extraction': 'production',
     'chp_methane_extraction': 'production',
     'chp_wte_back_pressure': 'production',
     'dac': 'consumption',
     'dc_underground_transmission': 'transmission',
     'demand_elec': 'consumption',
     'electric_heater': 'consumption',
     'electric_hob': 'consumption',
     'electrolysis': 'consumption',
     'heavy_transport_ev': 'consumption',
     'hp': 'consumption',
     'hydro_reservoir': 'consumption',
     'hydro_run_of_river': 'consumption',
     'hydrogen_to_liquids': 'consumption',
     'hydrogen_to_methanol': 'consumption',
     'light_transport_ev': 'consumption',
     'open_field_pv': 'production',
     'pumped_hydro': 'storage',
     'roof_mounted_pv': 'production',
     'wind_onshore_competing': 'production',
     'wind_onshore_monopoly': 'production'
}


def get_flow(path_to_model):
    m = calliope.read_netcdf(path_to_model)  # I used "outputs/resolution_test_month/link_cap_dynamic.nc"
    flows = {}
    for flow in ["prod", "con"]:
        df = split_loc_techs(m._model_data[f"carrier_{flow}"]).loc[{"carriers": "electricity", "locs": ["CHE_1", "CHE_2"]}].sum("locs").to_series()
        df = df.unstack("techs")
        df.columns = df.columns.str.split(":", expand=True)
        df = df.groupby(axis=1, level=0).sum()
        df = df.where(lambda x: abs(x) > 1e-5).dropna(axis=1, how='all')
        flows[flow] = df

    aligned = pd.concat(
        [flow.fillna(0).rename_axis(columns='techs').stack() for flow in flows.values()],
        sort=True, axis=1, keys=flows.keys()
    )
    flow = aligned["prod"].fillna(0) + aligned["con"].fillna(0)
    flow = flow.unstack('techs').groupby(pd.Series(GROUPS), axis=1).sum().rename_axis(columns='flow')

    flow_scaled = flow.copy()
    flow_scaled.production *= flow_scaled.consumption.abs().sum() / flow_scaled.production.sum()
    flow_scaled.transmission = flow_scaled.production + flow_scaled.consumption + flow_scaled.storage

    flow["consumption"] = flow["consumption"].abs()
    flow_scaled["consumption"] = flow_scaled["consumption"].abs()


def plot_7d_average(df):
    fig, ax = plt.subplots(1, 1, figsize=(10, 10))
    df_to_plot = df.resample('7D').mean().stack().to_frame('Electricity flow (100 GW)').reset_index()
    sns.lineplot(data=df_to_plot, x='timesteps', y='Electricity flow (100 GW)', hue='flow', ax=ax)
    return fig, ax
