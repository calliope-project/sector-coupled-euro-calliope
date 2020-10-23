import sys
import os

import pandas as pd
import geopandas as gpd
import matplotlib.pyplot as plt
import matplotlib as mpl
from shapely.geometry import box

sys.path.append(os.getcwd())
from src.construct import util

idx = pd.IndexSlice

cmaps = {
    'Air/shipping fuel': 'GnBu',
    'CO$_2$ feedstock': 'GnBu',
    'Cooking heat': 'Oranges',
    'Electricity': 'Blues',
    'Energy feedstock': 'GnBu',
    'High enthalpy heat': 'YlOrRd',
    'Hot water': 'Oranges',
    'Other vehicles mileage': 'BuPu',
    'Passenger car mileage': 'BuPu',
    'Rail electricity': 'Blues',
    'Space Heat': 'Oranges',
}
timeseries_colors = {
    'Fuel': '#abd9e9',
    'Building heat': '#fdae61',
    'Electricity': '#2c7bb6',
}

GDF_SIMPLIFY = 0.2

plt.rcParams.update({
    "text.usetex": True,
    "font.family": "sans-serif",
    "font.sans-serif": ["Helvetica"]
})
mpl.rcParams['text.latex.preamble'] = [
    #r'\usepackage{siunitx}',   # i need upright \micro symbols, but you need...
    #r'\sisetup{detect-all}',   # ...this to force siunitx to actually use your fonts
    r'\usepackage{helvet}',    # set the normal font here
    r'\usepackage{sansmath}',  # load up the sansmath so that math -> helvet
    r'\sansmath'               # <- tricky! -- gotta actually tell tex to use!
]


def plot_figure_1(
    path_to_units, path_to_annual_demand, path_to_national_units,
    path_to_electricity_demand, path_to_space_heat_demand, path_to_water_heat_demand,
    path_to_cooking_demand, energy_scaling_factor, model_year, out_path
):
    units = gpd.read_file(path_to_units).set_index('id')
    annual_demand = util.read_tdf(path_to_annual_demand)
    electricity_demand = pd.read_csv(path_to_electricity_demand, index_col=0, parse_dates=True)
    space_heat_demand = pd.read_csv(path_to_space_heat_demand, index_col=0, parse_dates=True)
    water_heat_demand = pd.read_csv(path_to_water_heat_demand, index_col=0, parse_dates=True)
    cooking_demand = pd.read_csv(path_to_cooking_demand, index_col=0, parse_dates=True)

    mean_annual_demand = (
        annual_demand
        .drop(annual_demand.filter(regex=r"bau", axis=0).index)
        .mean(level=['dataset', 'cat_name', 'id', 'end_use'])
    )

    annual_electricity_demand = -1 * electricity_demand.sum()
    building_electricity = annual_electricity_demand - mean_annual_demand.xs('electricity', level='end_use').sum(level='id')

    figure_data_dict = {
        ('Building demand', 'Space Heat', 'TWh'): mean_annual_demand[idx[:, :, :, 'space_heat']].sum(level='id') / energy_scaling_factor,
        ('Building demand', 'Hot water', 'TWh'): mean_annual_demand[idx[:, :, :, 'water_heat']].sum(level='id') / energy_scaling_factor,
        ('Building demand', 'Cooking heat', 'TWh'): mean_annual_demand[idx[:, :, :, 'cooking']].sum(level='id') / energy_scaling_factor,
        ('Building demand', 'Electricity', 'TWh'): building_electricity / energy_scaling_factor,
        ('Industry demand', 'Energy feedstock', 'TWh LHV'): mean_annual_demand.loc[idx[:, :, :, ['hydrogen', 'methanol']]].sum(level='id') / energy_scaling_factor,
        ('Industry demand', 'CO$_2$ feedstock', 'kt'): mean_annual_demand[idx[:, :, :, 'co2']].sum(level='id'),
        ('Industry demand', 'High enthalpy heat', 'TWh methane'): mean_annual_demand[idx[:, :, :, 'methane']].sum(level='id') / energy_scaling_factor,
        ('Industry demand', 'Electricity', 'TWh'): mean_annual_demand.loc[idx['industry_demand', 'industry', :, 'electricity']].sum(level='id') / energy_scaling_factor,
        ('Transport demand', 'Passenger car mileage', 'x10$^{10}$ km'): mean_annual_demand[idx[:, :, :, 'passenger_car']].sum(level='id'),
        ('Transport demand', 'Other vehicles mileage', 'x10$^{10}$ km'): mean_annual_demand.loc[idx[:, :, :, ['motorcycle', 'bus', 'ldv', 'hdv']]].sum(level='id'),
        ('Transport demand', 'Air/shipping fuel', 'TWh LHV'): mean_annual_demand.loc[idx[:, ['air', 'shipping'], :, :]].sum(level='id') / energy_scaling_factor,
        ('Transport demand', 'Rail electricity', 'TWh'): mean_annual_demand[idx[:, 'rail', :, 'electricity']].sum(level='id') / energy_scaling_factor
    }
    figure_data_df = pd.concat(figure_data_dict.values(), axis=1, keys=figure_data_dict.keys())

    timeseries = pd.concat(
        [electricity_demand.sum(axis=1),
         space_heat_demand.sum(axis=1) + water_heat_demand.sum(axis=1) + cooking_demand.sum(axis=1)],
        keys=['Electricity', 'Building heat'],
        axis=1
    ).div(-1 * 10)
    timeseries['Fuel'] = (
        figure_data_df
        .loc[:, idx[:, ['Energy feedstock', 'High enthalpy heat', 'Air/shipping fuel'], :]]
        .sum().sum()
        / len(timeseries)
    )

    plot_gdf = (
        units[['country_code', 'geometry']]
        .merge(figure_data_df, right_index=True, left_index=True)
        .to_crs('epsg:3035')
    )

    national_gdf_geom = gpd.read_file(path_to_national_units).set_index('CNTR_CODE').to_crs(plot_gdf.crs)
    national_gdf = (
        gpd.GeoDataFrame(
            gpd.clip(national_gdf_geom.geometry.explode(), box(*plot_gdf.total_bounds))
            .reset_index('CNTR_CODE')
        )
        .dissolve('CNTR_CODE')
    )
    national_gdf.index = national_gdf.index.map(util.get_alpha3)
    national_gdf = national_gdf.merge(
        plot_gdf.groupby('country_code').sum(), right_index=True, left_index=True
    )

    ncols = 32
    fig = plt.figure(figsize=(11, 15))
    g = plt.GridSpec(
        nrows=10, ncols=ncols, figure=fig,
        height_ratios=[1, 1, 10, 2.5, 10, 2.5, 10, 3, 6, 4], hspace=0.12
    )
    axes = {}
    axes['annual_title'] = plt.subplot(g[0, :], frameon=False)
    axes['annual_title'].axis('off')
    axes['annual_title'].set_title('\\textbf{a. Annual demand}', fontweight='bold', loc='left', y=0.1)
    r = 1
    for row in figure_data_df.columns.levels[0]:
        c = 0
        axes[row] = plt.subplot(g[r, :], frameon=False)
        axes[row].axis('off')
        axes[row].set_title(f'\\textit{{{row}}}', loc='center', y=0.1)
        r += 1
        for col in figure_data_df.xs(row, axis=1).columns:
            axes[(row, ) + col] = plt.subplot(g[r, int(c * (ncols / 4)):int(c * (ncols / 4) + (ncols / 4))], frameon=False)
            axes[(row, ) + col].axis('off')
            axes[(row, ) + col].set_title(col[0])
            c += 1
        r += 1
    for col in figure_data_df.columns:
        if col[1] == 'Air/shipping fuel':
            national_gdf.plot(
                col, ax=axes[col], cmap=cmaps[col[1]], legend=True,
                legend_kwds={'label': f'\\small{{{col[2]}}}', 'orientation': 'horizontal', 'pad': 0.03}
            )
        else:
            plot_gdf.plot(
                col, ax=axes[col], cmap=cmaps[col[1]], legend=True,
                legend_kwds={'label': f'\\small{{{col[2]}}}', 'orientation': 'horizontal', 'pad': 0.03}
            )
        national_gdf.plot(fc='None', ec='black', lw=0.1, ax=axes[col], legend=False, linestyle='dashed')

    timeseries_axes = {
        'Title': plt.subplot(g[-3, :], frameon=False),
        'Winter': plt.subplot(g[-2, 2:int(ncols / 2) - 1], frameon=False),
        'Summer': plt.subplot(g[-2, int(ncols / 2) + 1:ncols - 2], frameon=False),
        'Legend': plt.subplot(g[-1, :], frameon=False),
    }
    timeseries_dates = {
        'Winter': slice(f'{model_year}-02-01', f'{model_year}-02-07'),
        'Summer': slice(f'{model_year}-08-01', f'{model_year}-08-07')
    }
    timeseries_axes['Title'].axis('off')
    timeseries_axes['Legend'].axis('off')
    timeseries_axes['Title'].set_title('\\textbf{b. Hourly demand}', fontweight='bold', loc='left', y=0)
    for i in ['Winter', 'Summer']:
        timeseries_axes[i].set_title(f'\\textit{{{i}}}', loc='center')
        timeseries_axes[i].set_ylim(0, max(timeseries.loc[v].sum(axis=1).max() for v in timeseries_dates.values()))
        timeseries.loc[timeseries_dates[i], timeseries_colors.keys()].plot.area(ax=timeseries_axes[i], legend=False, alpha=0.8, lw=0, color=timeseries_colors.values())
        timeseries_axes[i].yaxis.grid(True, zorder=-1)
        timeseries_axes[i].set_axisbelow(True)
        timeseries_axes[i].set_ylabel('Energy demand (TWh)')
    timeseries_axes['Summer'].yaxis.set_label_position("right")
    timeseries_axes['Summer'].yaxis.tick_right()

    handles, labels = timeseries_axes['Winter'].get_legend_handles_labels()
    legend = timeseries_axes['Legend'].legend(
        labels=labels, handles=handles, loc='upper center', bbox_to_anchor=[0.5, 0.5], ncol=5, frameon=False
    )

    fig.savefig(out_path, pad_inches=0, bbox_inches='tight')


if __name__ == '__main__':
    plot_figure_1(
        path_to_units=snakemake.input.units,
        path_to_annual_demand=snakemake.input.annual_demand,
        path_to_national_units=snakemake.input.national_units,
        path_to_electricity_demand=snakemake.input.electricity_demand,
        path_to_space_heat_demand=snakemake.input.space_heat_demand,
        path_to_water_heat_demand=snakemake.input.water_heat_demand,
        path_to_cooking_demand=snakemake.input.cooking_demand,
        energy_scaling_factor=snakemake.params.energy_scaling_factor,
        model_year=snakemake.params.model_year,
        out_path=snakemake.output.out_path
    )