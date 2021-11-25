import sys
import os

import pandas as pd
import geopandas as gpd
import matplotlib.pyplot as plt
import matplotlib as mpl
import seaborn as sns

sys.path.append(os.getcwd())
from src.construct import util

idx = pd.IndexSlice

CMAPS = {
    'Road vehicle mileage': 'BuPu',
    'Electricity': 'Blues',
    'Building heat': 'Oranges', #'OrRd'
    'Synthetic fuel': 'YlGnBu',
}
TIMESERIES_COLORS = {
    'Synthetic fuel': '#a1dab4',
    'Building heat': '#fdae61',
    'Electricity': '#2c7bb6',
    'Road transport electricity': '#894da3B3',
    'Road transport diesel': '#894da3B3',
}

GDF_SIMPLIFY = 0.01

plt.rcParams.update({
    "svg.fonttype": 'none'
})


def plot_figure_1(
    path_to_units, path_to_annual_demand, path_to_electricity_demand,
    path_to_current_electricity_demand, path_to_space_heat_demand,
    path_to_water_heat_demand, path_to_cooking_demand,
    transport_efficiency, energy_scaling_factor, model_year, path_to_output
):
    units = gpd.read_file(path_to_units).set_index('id')
    annual_demand = util.read_tdf(path_to_annual_demand)
    electricity_demand = pd.read_csv(path_to_electricity_demand, index_col=0, parse_dates=True)
    space_heat_demand = pd.read_csv(path_to_space_heat_demand, index_col=0, parse_dates=True)
    water_heat_demand = pd.read_csv(path_to_water_heat_demand, index_col=0, parse_dates=True)
    cooking_demand = pd.read_csv(path_to_cooking_demand, index_col=0, parse_dates=True)
    current_electricity_demand = pd.read_csv(path_to_current_electricity_demand, index_col=0, parse_dates=True)
    mean_annual_demand = (
        annual_demand
        .drop(annual_demand.filter(regex=r"bau", axis=0).index)
        .xs(model_year, level="year")
        .mean(level=['dataset', 'cat_name', 'id', 'end_use'])
    )

    annual_electricity_demand = -1 * electricity_demand.loc[str(model_year)].sum()

    figure_data_dict = {
        ('Building heat', 'TWh'): (
            mean_annual_demand
            .loc[idx[:, :, :, ['space_heat', 'water_heat', 'cooking']]]
            .sum(level='id')
            / energy_scaling_factor
        ),
        ('Electricity', 'TWh'): (
            annual_electricity_demand / energy_scaling_factor
        ),
        ('Synthetic fuel', 'TWh LHV'): (
            mean_annual_demand
            .loc[idx[:, :, :, ['hydrogen', 'methanol', 'methane', 'kerosene', 'diesel']]]
            .sum(level='id')
            / energy_scaling_factor
        ),
        ('Road vehicle mileage', 'billion vehicle km'): (
            mean_annual_demand
            .loc[idx["transport_demand", :, :, ['motorcycle', 'bus', 'ldv', 'hdv', 'passenger_car']]]
            .sum(level='id')
            / 10
        ),
    }
    figure_data_df = pd.concat(figure_data_dict.values(), axis=1, keys=figure_data_dict.keys())

    timeseries = (
        pd.concat(
            [electricity_demand.sum(axis=1),
             space_heat_demand.sum(axis=1) + water_heat_demand.sum(axis=1) + cooking_demand.sum(axis=1)],
            keys=['Electricity', 'Building heat'],
            axis=1
        )
        .div(-1 * 10)
        .loc[str(model_year)]
    )
    #timeseries['Fuel'] = figure_data_df['Synthetic fuel'].sum().sum() / len(timeseries)
    current_timeseries = (
        current_electricity_demand
        .sum(axis=1)
        .div(-1 * 10)
        .rename('Current electricity load')
        .rename_axis(None)
        .loc[str(model_year)]
    )
    plot_gdf = (
        units[['country_code', 'geometry']]
        .merge(figure_data_df, right_index=True, left_index=True)
        .to_crs('epsg:3035')
    )
    with sns.plotting_context("paper", font_scale=1.4):
        ncols = 32
        fig = plt.figure(figsize=(15, 10))
        g = plt.GridSpec(
            nrows=5, ncols=ncols, figure=fig,
            height_ratios=[1, 10, 1, 6, 2], hspace=0.1
        )
        spatial_axes = {}
        spatial_axes['spatial_title'] = plt.subplot(g[0, :], frameon=False)
        c = 0
        for col in figure_data_df.columns:
            spatial_axes[col] = plt.subplot(
                g[1, int(c * (ncols / 4)):int(c * (ncols / 4) + (ncols / 4))],
                frameon=False
            )
            c += 1
        plot_spatial(plot_gdf, spatial_axes, figure_data_df)

        for carrier in ["electricity", "diesel"]:
            figure_data_df[(f'Road transport {carrier}', 'TWh')] = sum(
                mean_annual_demand
                .loc[idx["transport_demand", :, :, vehicle]]
                .sum(level="id")
                * 100  # 100 mio_km to mio_km
                * efficiency[carrier]  # mio_km to MWh
                / 1e6  # MWh to TWh
                for vehicle, efficiency in transport_efficiency.items()
            )
        annual_total_axes = {
            'Title': plt.subplot(g[-3, :int(ncols / 4)], frameon=False),
            'Plot': plt.subplot(g[-2:, 2:int(ncols / 4)], frameon=False),
        }
        plot_annual_bar(annual_total_axes, figure_data_df.sum())

        timeseries_axes = {
            'Title': plt.subplot(g[-3, 1 + int(ncols / 4):], frameon=False),
            'Plot': plt.subplot(g[-2, 3 + int(ncols / 4):ncols], frameon=False),
            'Legend': plt.subplot(g[-1, 1 + int(ncols / 4):], frameon=False),
        }
        plot_timeseries(timeseries_axes, timeseries, current_timeseries)

        if path_to_output.endswith(".png"):
            kwargs = {"dpi": 300}
        else:
            kwargs = {}

        fig.savefig(path_to_output, bbox_inches='tight', **kwargs)


def plot_spatial(plot_gdf, axes, figure_data_df):
    national_gdf = plot_gdf.dissolve("country_code")
    axes['spatial_title'].axis('off')
    axes['spatial_title'].set_title('a. Regional annual demand', fontweight='bold', loc='left', y=0.5)
    for col in figure_data_df.columns:
        axes[col].axis('off')

        plot_gdf.plot(
            col, ax=axes[col], cmap=CMAPS[col[0]], legend=True,
            legend_kwds={'label': f'{col[1]}', 'orientation': 'horizontal', 'pad': 0.03}
        )
        national_gdf.plot(fc='None', ec='black', lw=0.1, ax=axes[col], legend=False, linestyle='dashed')
        axes[col].set_title(col[0])


def plot_annual_bar(axes, annual_data):
    axes['Title'].axis('off')
    #axes['Legend'].axis('off')
    axes['Title'].set_title('b. Total annual demand', fontweight='bold', loc='left', y=0)

    to_plot_stacked = (
        annual_data
        .drop([
            'Road vehicle mileage',
            'Road transport electricity',
            'Road transport diesel'
        ], level=0)
        .div(1000)
        .droplevel(1)  # units
        .to_frame()
        .T
    )
    to_plot_grouped = (
        annual_data
        .loc[idx[['Road transport electricity', 'Road transport diesel'], :]]
        .div(1000)
        .add(to_plot_stacked.sum().sum())
        .droplevel(1)  # units
        .to_frame()
        .T
    )

    to_plot_stacked.plot.bar(
        stacked=True, ax=axes["Plot"],
        color=[TIMESERIES_COLORS[i] for i in to_plot_stacked.columns]
    )
    _, labels = axes["Plot"].get_legend_handles_labels()
    for i in range(len(axes["Plot"].patches)):
        patch = axes["Plot"].patches[i]
        label = labels[i]
        bl = patch.get_xy()
        x = 0.5 * patch.get_width() + bl[0]
        y = 0.5 * patch.get_height() + bl[1]
        axes["Plot"].annotate(
            label, xy=(x, y), xycoords="data",
            verticalalignment='center', horizontalalignment='center',
            color="white"
        )

    to_plot_grouped.plot.bar(
        ax=axes["Plot"],
        color=[TIMESERIES_COLORS[i] for i in to_plot_grouped.columns],
        ec="black", lw=0.2,
        zorder=-1
    )
    _, labels = axes["Plot"].get_legend_handles_labels()

    for j in range(i + 1, len(axes["Plot"].patches)):
        patch = axes["Plot"].patches[j]
        label = "EV" if "electricity" in labels[j] else "ICE"
        bl = patch.get_xy()
        ymin = to_plot_stacked.sum().sum()
        x = 0.5 * patch.get_width() + bl[0]
        y = (ymin + (patch.get_height() - ymin) / 2) + bl[1]
        axes["Plot"].annotate(
            label, xy=(x, y), xycoords="data",
            verticalalignment='center', horizontalalignment='center',
            color="white", zorder=10
        )

    axes["Plot"].legend().set_visible(False)
    axes["Plot"].tick_params(axis="x", which="both", bottom=False, labelbottom=False)
    axes["Plot"].spines['top'].set_visible(False)
    axes["Plot"].spines['right'].set_visible(False)
    axes["Plot"].spines['bottom'].set_visible(False)
    axes["Plot"].set_ylabel('Energy demand (1000 TWh)')
    axes["Plot"].set_xlim(-0.3, 0.4)

    #move_legend_between_axes(axes)


def plot_timeseries(axes, model_timeseries_data, current_demand):
    axes['Title'].axis('off')
    axes['Legend'].axis('off')
    axes['Title'].set_title('c. Hourly demand', fontweight='bold', loc='left', y=0)

    current_demand.rolling(24 * 7, center=True).mean().plot(
        ax=axes["Plot"], lw=2, c='black', linestyle='--'
    )
    handles, labels = axes["Plot"].get_legend_handles_labels()

    for ts in model_timeseries_data.columns:
        model_timeseries_data[ts].rolling(24 * 7, center=True).mean().plot(
            ax=axes["Plot"], color=TIMESERIES_COLORS, label=ts
        )
        model_timeseries_data[ts].plot(
            ax=axes["Plot"], alpha=0.5, color=TIMESERIES_COLORS, legend=False
        )
        handles.append(
            mpl.lines.Line2D(
                xdata=[0, 1], ydata=[0, 1], color=TIMESERIES_COLORS[ts], linestyle='-'
            )
        )
        labels.append(ts)

    sns.despine(ax=axes["Plot"])
    axes["Plot"].set_ylabel('Energy demand (TWh)')
    move_legend_between_axes(axes, handles, labels)


def move_legend_between_axes(axes, handles=None, labels=None):
    if handles is None and labels is None:
        handles, labels = axes['Plot'].get_legend_handles_labels()
    axes['Plot'].legend().set_visible(False)
    axes['Legend'].legend(
        labels=labels, handles=handles,
        loc='upper center',
        bbox_to_anchor=[0.5, 0.5],
        ncol=5,
        frameon=False
    )


if __name__ == '__main__':
    plot_figure_1(
        path_to_units=snakemake.input.units,
        path_to_annual_demand=snakemake.input.annual_demand,
        path_to_electricity_demand=snakemake.input.electricity_demand,
        path_to_current_electricity_demand=snakemake.input.current_electricity_demand,
        path_to_space_heat_demand=snakemake.input.space_heat_demand,
        path_to_water_heat_demand=snakemake.input.water_heat_demand,
        path_to_cooking_demand=snakemake.input.cooking_demand,
        transport_efficiency=snakemake.params.transport_efficiency,
        energy_scaling_factor=snakemake.params.energy_scaling_factor,
        model_year=snakemake.params.model_year,
        path_to_output=snakemake.output[0]
    )