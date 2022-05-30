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
SUBTITLES = {
    'Electricity': "(appliances, cooling, rail, industry processes)",
    'Synthetic fuel': "(industry processes & feedstocks, marine, aviation)",
    'Building heat': "(space heating, hot water, cooking)",
    'Road vehicle mileage': "(passenger, commercial, freight)",
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
    "svg.fonttype": 'none',
    'font.family':'sans-serif',
    'font.sans-serif':'Arial',
    "font.size": 7
})

FIGWIDTH = 6.77165

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
        ('Road vehicle mileage', 'Billion vehicle km'): (
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
        .rename('Actual 2018 electricity load profile')
        .rename_axis(None)
        .loc[str(model_year)]
    )
    plot_gdf = (
        units[['country_code', 'geometry']]
        .merge(figure_data_df, right_index=True, left_index=True)
        .to_crs('epsg:3035')
    )


    ncols = 32
    fig = plt.figure(figsize=(FIGWIDTH, 10 * FIGWIDTH/15))
    g = plt.GridSpec(
        nrows=7, ncols=ncols, figure=fig,
        height_ratios=[2, 16, 1, 2, 2, 12, 4], hspace=0.1, wspace=0.25
    )
    fig.subplots_adjust(left=0, right=1, top=1, bottom=0)
    spatial_axes = {}
    spatial_axes['spatial_title'] = plt.subplot(g[0, :], frameon=False)
    c = 0
    for col in figure_data_df.columns:
        spatial_axes[col] = plt.subplot(
            g[1, int(c * (ncols / 4)):int(c * (ncols / 4) + (ncols / 4))],
            frameon=False
        )
        spatial_axes[col[0] + "_legend"] = plt.subplot(
            g[2, int(c * (ncols / 4)):int(c * (ncols / 4) + (ncols / 4))],
            frameon=False
        )
        c += 1
    plot_spatial(plot_gdf, spatial_axes, fig, figure_data_df)

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
        'Title': plt.subplot(g[-3, 2 + int(ncols / 4):], frameon=False),
        'Plot': plt.subplot(g[-2, 4 + int(ncols / 4):ncols], frameon=False),
        'Legend': plt.subplot(g[-1, 2 + int(ncols / 4):], frameon=False),
    }
    plot_timeseries(timeseries_axes, timeseries, current_timeseries)

    if path_to_output.endswith(".png") or path_to_output.endswith(".tif"):
        kwargs = {"dpi": 300}
    else:
        kwargs = {}

    fig.savefig(path_to_output, bbox_inches='tight', **kwargs)


def plot_spatial(plot_gdf, axes, fig, figure_data_df):
    national_gdf = plot_gdf.dissolve("country_code")
    axes['spatial_title'].axis('off')
    axes['spatial_title'].set_title('$\\bf{A}$ Regional annual service demands', loc='left', y=0.5)
    for col in figure_data_df.columns:
        axes[col].axis('off')

        plot_gdf.plot(
            col, ax=axes[col], cmap=CMAPS[col[0]], legend=False,
            legend_kwds={'label': f'{col[1]}', 'orientation': 'horizontal', 'pad': 0.03}
        )

        sm = plt.cm.ScalarMappable(
            cmap=CMAPS[col[0]],
            norm=plt.Normalize(vmin=plot_gdf[col].min(), vmax=plot_gdf[col].max()),
        )
        sm._A = []
        cbar = fig.colorbar(
            sm, cax=axes[col[0] + "_legend"],
            orientation="horizontal", label=f'{col[1]}', shrink=0.65, pad=0.05
        )
        cbar.outline.set_linewidth(0.1)
        cbar.ax.tick_params(axis='both', which='both', labelsize=6, pad=0, width=0.5)

        national_gdf.plot(fc='None', ec='black', lw=0.05, ax=axes[col], legend=False, linestyle='dashed')
        title_fontdict = {"fontsize": 5.5, "fontstyle": "italic"}
        axes[col].set_title(col[0])
        axes[col].text(0.5, 0.98,
            SUBTITLES[col[0]],
            verticalalignment="bottom", horizontalalignment="center",
            transform=axes[col].transAxes, **title_fontdict
        )

def plot_annual_bar(axes, annual_data):
    axes['Title'].axis('off')
    axes['Title'].set_title(
        '$\\bf{B}$ Total annual service demands', loc='left', y=0
    )
    axes['Title'].text(0.1, 0.25,
        "(road vehicle mileage translated to energy\nbased on model input vehicle efficiency)",
        verticalalignment="top", horizontalalignment="left",
        transform=axes['Title'].transAxes, fontsize=6
    )

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
    axes['Title'].set_title('$\\bf{C}$ Hourly load profiles for fixed timeseries demands', loc='left', y=0)

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
        if ts == "Electricity":
            labels.append("Modelled electricity demand")
        elif ts == "Building heat":
            labels.append("Modelled building heat demand")
        else:
            labels.append(ts)

    sns.despine(ax=axes["Plot"])
    axes["Plot"].set_ylabel('Energy demand (TWh)')
    axes["Plot"].tick_params(axis='x',which='minor',bottom=False)
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
        frameon=False,
        fontsize=7
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
