import os

import pandas as pd
import geopandas as gpd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import seaborn as sns

from friendly_data.converters import to_df
from frictionless.package import Package

COLORS = {
    "Oil": "#5d5d5d",
    "Natural gas": "#b9b9b9",
    "Other fossils": "#181818",
    "Nuclear heat": "#cc0000",
    "Biofuels": "#8fce00",
    "PV": "#ffd966",
    "pv": "#ffd966",
    "Onshore wind": "#674ea7",
    "onshore wind": "#674ea7",
    "Offshore wind": "#e062db",
    "offshore wind": "#e062db",
    "Hydro": "#2986cc",
    "Waste": "#ce7e00",
    "Direct heat": "#f6b26b",
    "Electricity": "#2986cc"
}
ENERGY_PRODUCERS = {
    "waste_supply": "Waste",
    "biofuel_supply": "Biofuels",
    "hydro_reservoir": "Hydro",
    "hydro_run_of_river": "Hydro",
    "nuclear": "Nuclear heat",
    "open_field_pv": "PV",
    "roof_mounted_pv": "PV",
    "wind_offshore": "Offshore wind",
    "wind_onshore": "Onshore wind",
}
NUCLEAR_HEAT_MUTIPLIER = 1 / 0.4  # our model uses an efficiency of 40% for nuclear

plt.rcParams.update({
    "svg.fonttype": 'none',
    'font.family':'sans-serif',
    'font.sans-serif':'Arial'
})


def make_plot(path_to_friendly_data, path_to_units, path_to_unit_groups, spore_num, path_to_output):
    units = gpd.read_file(path_to_units).set_index("id")
    units_m = units.to_crs("EPSG:3035")

    unit_groups = pd.read_csv(path_to_unit_groups, index_col=0, squeeze=True).reindex(units.index)

    units["Group"] = unit_groups.values
    units_m["Group"] = unit_groups.values

    units_grouped = units_m.dissolve("Group")
    dpkg_filepath = os.path.join(path_to_friendly_data, "datapackage.json")
    spores_datapackage = Package(dpkg_filepath)

    spores_data = {}

    for resource in spores_datapackage["resources"]:
        spores_data[resource["name"]] = pd.read_csv(
            os.path.join(path_to_friendly_data, resource["path"]),
            index_col=resource["schema"]["primaryKey"], squeeze=True
        )

    spore_name_dict = {spore_num: f'SPORE {spore_num}'}
    prod_regional, prod_regional_ymin, prod_regional_ymax = (
        get_pv_wind_grouped_region_prod(spores_data, units, spore_name_dict)
    )
    prod_sum, prod_sum_max = get_primary_energy_supply(spores_data, spore_name_dict)
    net_import, import_cbar_max = get_region_net_import(spores_data)
    synfuel_contribution, overprod_cbar_max = get_synfuel_contribution(spores_data)
    link_caps, link_cap_max, link_cap_mean = get_link_caps(spores_data)

    plot_comparison_maps(
        prod_regional,
        prod_regional_ymin,
        prod_regional_ymax,
        prod_sum,
        prod_sum_max,
        net_import,
        import_cbar_max,
        synfuel_contribution,
        overprod_cbar_max,
        link_caps,
        link_cap_max,
        link_cap_mean,
        units_m,
        units_grouped,
        spore_name_dict,
        path_to_output
    )

def get_output_energy(spores_data, spores=None):
    flow_out = spores_data["flow_out_sum"]
    flow_out_summed = (
        flow_out
        .unstack("spore")
        .sum(level=["techs", "carriers"])
        .groupby(ENERGY_PRODUCERS, level="techs").sum()
    )
    flow_out_summed.loc["Nuclear heat"] *= NUCLEAR_HEAT_MUTIPLIER
    if spores is None:
        return flow_out_summed
    else:
        return (
            flow_out_summed
            [spores]
            .loc[flow_out_summed[spores].mean(axis=1).sort_values(ascending=False).index]
        )


def get_pv_wind_grouped_region_prod(result_dict, units_grouped, spore_names):
    prod = (
        result_dict["flow_out_sum"]
        .unstack("spore")
        .groupby([ENERGY_PRODUCERS, units_grouped.Group.to_dict()], level=["techs", "locs"])
        .sum()
        .loc[["Onshore wind", "Offshore wind", "PV"]]
        .unstack("techs")
        .fillna(0)
        .stack("techs")
        .div(1000)
    )
    prod_ymin = prod.min().min()
    prod_ymax = prod.max().max()
    prod = prod[spore_names.keys()].where(prod.div(prod_ymax) > 0.01)

    return prod, prod_ymin, prod_ymax

def get_primary_energy_supply(result_dict, spore_names):
    prod_sum = get_output_energy(result_dict, spore_names.keys()).div(1000)
    prod_sum = prod_sum.where(prod_sum.div(prod_sum.sum()) > 0.01).dropna(how="all")
    prod_sum_max = get_output_energy(result_dict, None).div(1000).sum().max()

    return prod_sum, prod_sum_max


def get_region_net_import(result_dict):
    net_import = (
        result_dict
        ["net_import_sum"]
        .unstack("spore")
    )
    net_import = net_import[
        net_import.index.get_level_values("importing_region") !=
        net_import.index.get_level_values("exporting_region")
    ].sum(level="importing_region")
    import_cbar_max = net_import.abs().max().max()

    return net_import, import_cbar_max


def get_synfuel_contribution(result_dict):
    synfuel_contribution = (
        result_dict
        ["flow_out_sum"]
        .xs("electrolysis", level="techs")
        .unstack("spore")
        .sum(level="locs")
        .div(result_dict["flow_out_sum"].xs("electrolysis", level="techs").sum(level="spore"))
    )
    overprod_cbar_max = synfuel_contribution.max().max()

    return synfuel_contribution, overprod_cbar_max


def get_link_caps(result_dict):
    link_caps = (
        result_dict["grid_capacity_expansion"]
        .unstack("spore")
        .sum(level=["importing_region", "exporting_region"])
        .where(lambda x: x > 1e-5)
        .fillna(0)
    )
    link_cap_max = link_caps.max().max()
    link_cap_mean = link_caps.mean().mean()

    return link_caps, link_cap_max, link_cap_mean


def plot_comparison_maps(
    prod_regional,
    prod_regional_ymin,
    prod_regional_ymax,
    prod_sum,
    prod_sum_max,
    net_import,
    import_cbar_max,
    synfuel_contribution,
    overprod_cbar_max,
    link_caps,
    link_cap_max,
    link_cap_mean,
    units_m,
    units_grouped,
    spore_name_dict,
    path_to_output
):
    nspores = len(spore_name_dict.keys())
    nrows = nspores * 2 + 2
    centroids = units_m.centroid

    with sns.plotting_context("paper", font_scale=1.5):
        ax = {}
        fig = plt.figure(figsize=(20, nspores * 7 + 4))
        gs = mpl.gridspec.GridSpec(
            nrows=nrows,
            ncols=5,
            figure=fig,
            hspace=0.1,
            wspace=0.1,
            width_ratios=[1, 1, 23, 25, 25],
            height_ratios=[2, 18] * nspores + [2, 2]
        )
        _row = 0
        alpha_idx = 0

        for spore_num, spore_name in spore_name_dict.items():
            ax[spore_num] = {}
            ax[spore_num]["title"] = plt.subplot(gs[_row, :], frameon=False)
            plot_title(ax[spore_num]["title"], spore_name)

            ax[spore_num]["total_bar"] = plt.subplot(gs[_row + 1, 1], frameon=False)
            handles, labels = plot_total_energy_supply(
                ax[spore_num]["total_bar"],
                prod_sum, prod_sum_max,
                spore_num, spore_name,
                COLORS
            )


            ax[spore_num]["wind_pv_map"] = plt.subplot(gs[_row + 1, 2], frameon=False)
            plot_regional_energy_supply(
                ax[spore_num]["wind_pv_map"],
                prod_regional, prod_regional_ymax,
                units_m, units_grouped,
                spore_num, spore_name,
                COLORS
            )
            ax[spore_num]["wind_pv_map"].annotate(
                "Annual primary energy supply (bar)\n& annual regional PV & wind generation (map)",
                fontweight="bold",
                xy=(0.5, 1.1),
                xycoords='axes fraction',
                horizontalalignment="center",
                fontsize="small"
            )

            alpha_idx += 1
            ax[spore_num]["import_synfuel_map"] = plt.subplot(gs[_row + 1, 3], frameon=False)
            cmap = mpl.cm.coolwarm(np.linspace(0, 1, 200))
            import_cmap = mpl.colors.ListedColormap(cmap)
            synfuel_cmap = mpl.colors.LinearSegmentedColormap.from_list("", ["yellow", "grey", "black"])
            plot_imports_and_synfuel_overprod(
                ax[spore_num]["import_synfuel_map"],
                net_import, synfuel_contribution,
                overprod_cbar_max, import_cbar_max,
                spore_num, spore_name,
                units_m, centroids,
                import_cmap, synfuel_cmap
            )

            alpha_idx += 1
            ax[spore_num]["transmission_expansion_map"] = plt.subplot(gs[_row + 1, 4], frameon=False)
            plot_transmission_expansion(
                ax[spore_num]["transmission_expansion_map"],
                link_caps, link_cap_max,
                spore_num, spore_name,
                units_m, centroids,
            )
            _row += 2
            alpha_idx += 1

        ax[spore_num]["primary_energy_legend"] = plt.subplot(gs[-2:, 0:3], frameon=False)
        plot_regional_and_total_primary_energy_supply_legend(ax[spore_num]["primary_energy_legend"], handles, labels, ncol=3)

        ax[spore_num]["import_legend"] = plt.subplot(gs[-2, 3], frameon=False)
        ax[spore_num]["synfuel_legend"] = plt.subplot(gs[-1, 3], frameon=False)
        plot_imports_and_synfuel_overprod_legend(
            ax[spore_num]["import_legend"],
            ax[spore_num]["synfuel_legend"],
            import_cbar_max,
            overprod_cbar_max,
            import_cmap,
            synfuel_cmap
        )

        ax[spore_num]["transmission_legend"] = plt.subplot(gs[-2:, 4], frameon=False)
        plot_transmission_expansion_legend(ax[spore_num]["transmission_legend"], link_cap_mean, link_cap_max)

        fig.savefig(path_to_output, bbox_inches="tight", dpi=120)


def plot_regional_energy_supply(
    ax_map,
    prod_regional, prod_regional_ymax,
    units_m, grouped_units_m,
    spore_num, spore_name,
    tech_colors
):

    def _build_bar(mapx, mapy, ax, width, hideaxis=True):
        ax_h = ax.inset_axes(
            [mapx - width / 2, mapy - width / 2, width, width], transform=ax.transData,
        )
        if hideaxis:
            ax_h.set_xlabel("")
            ax_h.set_ylabel("")
        else:
            sns.despine(ax=ax_h)
        return ax_h

    ax_map.set_xlim((2.6e6, 6.23e6))

    group_data = prod_regional[spore_num].unstack()

    _selected_data = 0
    for group in group_data.index:
        if group_data.loc[group].sum() > prod_regional_ymax * 0.05:
            ax_h = _build_bar(
                grouped_units_m.centroid.loc[group].x,
                grouped_units_m.centroid.loc[group].y,
                ax_map,
                4e5,
            )
            group_data.loc[group].plot.bar(
                ax=ax_h,
                color=[tech_colors[i] for i in group_data.loc[group].index],
                width=1
            )
            ax_h.set_ylim(ymax=prod_regional_ymax)
            ax_h.tick_params(axis="x", which="both", bottom=False, labelbottom=False)
            ax_h.spines['right'].set_visible(False)
            ax_h.spines['top'].set_visible(False)
            ax_h.spines['left'].set_visible(False)
            ax_h.set_facecolor("None")
            ax_h.set_xlabel("")
            ax_h.set_ylabel("")
            grouped_units_m.loc[[group]].plot(ax=ax_map, fc="#E0ECF8", ec="white", lw=0.5)
            _selected_data += group_data.loc[group].sum()
            if group == 16:
                ax_h.spines['left'].set_visible(True)
                ax_h.set_ylabel("1000 TWh", fontweight="bold", fontsize="x-small")
            else:
                ax_h.tick_params(left=False)
                ax_h.set_yticklabels(['' for i in ax_h.get_yticklabels()])
        else:
            grouped_units_m.loc[[group]].plot(ax=ax_map, fc="lightgrey", ec="white", lw=0.5)
    units_m.plot(ax=ax_map, fc="None", linestyle="-", ec="white", lw=0.1)
    print(f"{100 * _selected_data / group_data.sum().sum():.2f}% Production plotted for SPORE {spore_num}")

    ax_map.axis("off")


def plot_total_energy_supply(
    ax_total_bar,
    prod_sum, prod_sum_max,
    spore_num, spore_name,
    tech_colors
):

    prod_sum[[spore_num]].T.plot.bar(
        stacked=True,
        ax=ax_total_bar,
        color=[tech_colors[i] for i in prod_sum.index],
        width=1
    )
    ax_total_bar.spines['right'].set_visible(False)
    ax_total_bar.spines['left'].set_visible(True)
    ax_total_bar.spines['top'].set_visible(False)
    ax_total_bar.spines['bottom'].set_visible(False)
    ax_total_bar.set_xticklabels([])
    ax_total_bar.set_xticks([])
    ax_total_bar.set_ylabel("1000 TWh")
    ax_total_bar.set_ylim(ymax=prod_sum_max)
    ax_total_bar.set_xlim((-0.51, 0.5))
    handles, labels = ax_total_bar.get_legend_handles_labels()
    ax_total_bar.legend().set_visible(False)
    ax_total_bar.set_xlabel("")

    return handles, labels


def plot_regional_and_total_primary_energy_supply_legend(ax_legend, handles, labels, ncol=4):
    ax_legend.axis("off")
    ax_legend.legend(handles=handles, labels=labels, ncol=ncol, frameon=False, loc="center")


def plot_imports_and_synfuel_overprod(
    ax_map,
    net_import, synfuel_contribution,
    overprod_cbar_max, import_cbar_max,
    spore_num, spore_name,
    units_m, centroids,
    import_cmap, synfuel_cmap
):
    units_m["net_import"] = net_import[spore_num].reindex(units_m.index)
    ax_map.set_xlim((2.6e6, 6.55e6))
    units_m.plot(
        "net_import",
        ax=ax_map,
        cmap=import_cmap,
        ec="white",
        lw=0.5,
        vmin=-import_cbar_max,
        vmax=import_cbar_max
    )
    unit_centroids = (
        centroids
        .buffer(5e4)
        .to_frame("geometry")
        .assign(overprod=synfuel_contribution[spore_num].reindex(units_m.index))
    )
    unit_centroids[unit_centroids.overprod > 0.05].plot(
        "overprod",
        cmap=synfuel_cmap,
        ax=ax_map,
        lw=0.2,
        ec="black",
        vmin=0.05, vmax=overprod_cbar_max
    )

    ax_map.annotate(
        "Regional electricity imports (choropleth)\n& synfuel production hubs (points)",
        fontweight="bold",
        xy=(0.5, 1.1),
        xycoords='axes fraction',
        horizontalalignment="center",
        fontsize="small"
    )
    ax_map.axis("off")


def plot_imports_and_synfuel_overprod_legend(
    ax_import_legend,
    ax_synfuel_legend,
    import_cbar_max,
    overprod_cbar_max,
    import_cmap,
    synfuel_cmap
):
    sm = plt.cm.ScalarMappable(cmap=import_cmap, norm=plt.Normalize(vmin=-import_cbar_max / 10000, vmax=import_cbar_max / 10000))
    cbar = plt.colorbar(sm, ax=ax_import_legend, orientation="horizontal", pad=0.1, anchor=(0.5, 1.5), aspect=40, fraction=1)
    cbar.ax.set_xlabel("Net electricity import (1000 TWh)")

    sm = plt.cm.ScalarMappable(cmap=synfuel_cmap, norm=plt.Normalize(vmin=0.05, vmax=overprod_cbar_max))
    cbar = plt.colorbar(sm, ax=ax_synfuel_legend, orientation="horizontal", pad=0.1, anchor=(0.5, 1.5), aspect=40, fraction=1)
    cbar.ax.set_xlabel("Fraction of total European hydrogen production")
    ax_import_legend.axis("off")
    ax_synfuel_legend.axis("off")


def plot_transmission_expansion(
    ax_map,
    link_caps, link_cap_max,
    spore_num, spore_name,
    units_m, centroids,
):
    def _add_line(ax, point_from, point_to, color, lw, dashes=(None, None), linestyle='-'):
        ax.add_line(mpl.lines.Line2D(
            xdata=(point_from[0], point_to[0]), ydata=(point_from[1], point_to[1]),
            color=color, lw=lw,
            linestyle=linestyle, dashes=dashes
        ))

    def _add_links(ax, link_caps, cap_max):
        links_completed = []
        for link in link_caps.index:
            if sorted(link) in links_completed:  # every link comes up twice
                continue
            cap = link_caps.loc[link]
            point_from = (centroids.loc[link[0]].x, centroids.loc[link[0]].y)
            point_to = (centroids.loc[link[1]].x, centroids.loc[link[1]].y)
            _add_line(
                ax, point_from, point_to,
                color='#a9a9a999',
                lw=0.5
            )
            _add_line(
                ax, point_from, point_to,
                color='blue',
                lw= 5 * cap / cap_max
            )

            links_completed.append(sorted(link))

    ax_map.set_xlim((2.6e6, 6.55e6))
    units_m.plot(ax=ax_map, fc='#add8e6', ec="white", lw=0.5)
    _add_links(ax_map, link_caps[spore_num], link_cap_max)

    ax_map.annotate(
        f"Transmission capacity expansion\n(Total: + {link_caps[spore_num].sum() / 10:.1f} TW)",
        fontweight="bold",
        xy=(0.5, 1.1),
        xycoords='axes fraction',
        horizontalalignment="center",
        fontsize="small"
    )
    ax_map.axis("off")


def plot_transmission_expansion_legend(ax_legend, link_cap_mean, link_cap_max, ncol=1):
    existing = mpl.lines.Line2D(
        xdata=[0, 1], ydata=[0, 1],
        color="#a9a9a999",
        lw=0.5,
        label="Existing link"
    )
    mean_line = mpl.lines.Line2D(
        xdata=[0, 1], ydata=[0, 1],
        color="blue",
        lw=5 * link_cap_mean / link_cap_max,
        label=f"+ {link_cap_mean * 1000:.0f} GW"
    )
    biggest_line = mpl.lines.Line2D(
        xdata=[0, 1], ydata=[0, 1],
        color="blue",
        lw=5 * link_cap_max / link_cap_max,
        label=f"+ {link_cap_max * 1000:.0f} GW"
    )
    ax_legend.axis("off")
    ax_legend.legend(
        handles=[existing, mean_line, biggest_line],
        frameon=False, loc='center', fontsize=12, ncol=ncol
    )


def plot_title(ax_title, spore_name):
    ax_title.axis("off")
    ax_title.text(
        0, 0.5,
        spore_name,
        va='bottom',
        ha='left',
        fontweight="bold"
    )


if __name__ == "__main__":
    make_plot(
        path_to_friendly_data=snakemake.input.friendly_data,
        path_to_units=snakemake.input.units,
        path_to_unit_groups=snakemake.input.unit_groups,
        spore_num=int(snakemake.wildcards.spore),
        path_to_output=snakemake.output[0]
    )
