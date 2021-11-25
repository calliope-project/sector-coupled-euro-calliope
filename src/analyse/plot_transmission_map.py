import os

import geopandas as gpd
import pandas as pd
import numpy as np

from friendly_data.converters import to_df
from frictionless.package import Package

import cartopy.crs as ccrs
from cartopy.io import shapereader

import matplotlib.pyplot as plt
from matplotlib import gridspec
import matplotlib.lines as mlines
import matplotlib.patches as mpatches
plt.rcParams['svg.fonttype'] = 'none'

import calliope

OUTER_COUNTRIES = ['RUS', 'UKR', 'TUR', 'MAR', 'DZA', 'TUN', 'LBY', 'EGY']

COLORS = {
    'outer': '#a9a9a9', 'ac': '#003366', 'dc': '#ffdf00',
    'eu_bg': '#add8e6', 'outer_bg': '#DCDCDC'
}


def plot_input_transmission(path_to_model, path_to_units, path_to_output, bounds):
    fig = plt.figure(figsize=(15, 16))
    gs = gridspec.GridSpec(2, 1, figure=fig, height_ratios=[20, 1], hspace=0)  # rows, columns

    extent = (bounds['x_min'], bounds['x_max'], bounds['y_min'], bounds['y_max'])

    ax_eu = plt.subplot(gs[0, 0])
    ax_eu.set_extent(extent)
    ax_eu.axis("off")  #outline_patch.set_linewidth(0)

    gdf = gpd.read_file(path_to_units).to_crs("EPSG:3035")
    plot_polygons(gdf, ax_eu, include_neighbours=True)

    link_caps = get_link_caps_from_input(path_to_model)
    cap_max = link_caps.max().max()
    add_links(link_caps, cap_max, ax_eu)

    ax_leg = plt.subplot(gs[1, 0])
    make_legend(ax_leg)

    if path_to_output.endswith(".png"):
        kwargs = {"dpi": 300}
    else:
        kwargs = {}

    fig.savefig(path_to_output, bbox_inches="tight", **kwargs)


def plot_spores_transmission(path_to_spores, path_to_units, spores, titles, path_to_output):
    fig = plt.figure(figsize=(20, 20))
    gs = gridspec.GridSpec(
        3, 2,  # rows, columns
        figure=fig, height_ratios=[20, 20, 1], hspace=0.1
    )
    link_caps = get_link_caps_from_output(path_to_spores)
    link_caps_min = link_caps.loc[spores].min()
    link_caps_max = link_caps.loc[spores].max()

    gdf = gpd.read_file(path_to_units).to_crs("EPSG:3035")
    ax = {}
    _row = 0
    _col = 0
    for spore in spores:
        ax[spore] = plt.subplot(gs[_row, _col])
        ax[spore].axis("off")  #outline_patch.set_linewidth(0)
        plot_polygons(gdf, ax[spore], include_neighbours=False)
        add_links(
            link_caps.loc[spore], link_caps_max, gdf, ax,
            distinguish_ac_dc=False, width_by_cap=False
        )
        add_links(
            link_caps.loc[spore], link_caps_max, gdf, ax,
            distinguish_ac_dc=False, width_by_cap=True
        )
        ax[spore].set_title(titles[spore])
        if _col == 1:
            _col = 0
            _row += 1
        else:
            _col += 1

    ax_leg = plt.subplot(gs[2, :])
    smallest_line = mlines.Line2D(
        xdata=(0, 1), ydata=(0, 0),
        color=COLORS["ac"],
        lw=1 + 1.5 * link_caps_min / link_caps_max,
        label=f"+ {link_caps_min / 10} TWh"
    )
    biggest_line = mlines.Line2D(
        xdata=(0, 1), ydata=(0, 0),
        color=COLORS["ac"],
        lw=1 + 1.5 * link_caps_max / link_caps_max,
        label=f"+ {link_caps_max / 10} TWh"
    )
    ax_leg.legend(
        handles=[smallest_line, biggest_line],
        frameon=False, loc='center', fontsize=12, ncol=3
    )

    if path_to_output.endswith(".png"):
        kwargs = {"dpi": 300}
    else:
        kwargs = {}

    fig.savefig(path_to_output, bbox_inches="tight", **kwargs)


def plot_polygons(gdf, ax, include_neighbours=True):
    gdf.plot(
        ax=ax, facecolor=COLORS['eu_bg'], edgecolor='gray',
        linestyle="--", linewidth=0.3
    )
    gdf.dissolve('country_code').plot(
        facecolor='None', edgecolor='white', linewidth=0.5
    )
    if include_neighbours:
        # get country borders
        shpfilename = shapereader.natural_earth(
            resolution='10m', category='cultural', name='admin_0_countries',
        )

        # read the shapefile using geopandas
        outer_gdf = gpd.read_file(shpfilename).to_crs("EPSG:3035")

        # read the outer country borders
        outer_gdf = outer_gdf.loc[~outer_gdf.ADM0_A3.isin(gdf.country_code.unique())]

        outer_gdf.plot(
            facecolor=COLORS['outer_bg'], edgecolor='white', linewidth=0.5
        )


def get_link_caps_from_input(path_to_model):
    m = calliope.read_netcdf(path_to_model)
    links = m.inputs.energy_cap_min.loc[
        {'loc_techs': m._model_data.loc_techs_transmission}
    ].to_pandas()
    links.index = (
        links.index.str.split(':', expand=True).droplevel(1).rename(('to', 'tech', 'from'))
    )
    return links.unstack(1)  # empty column from splitting '::'


def get_link_caps_from_output(path_to_spores_dpkg):
    spores_datapackage = Package(
        os.path.join(path_to_spores_dpkg, "spores", "datapackage.json")
    )
    spores_data_dicts = {
        resource["name"]: to_df(resource).squeeze().unstack("scenario")
        for resource in spores_datapackage["resources"]
    }
    return (
        spores_data_dicts["transmission_cap_expansion"]
        .rename_axis(index=('spore', 'from', 'to', 'tech'))
        .reorder_levels(('spore', 'to', 'tech', 'from'))
    )


def add_links(link_caps, cap_max, gdf, ax, distinguish_ac_dc=True, width_by_cap=True):

    def _add_line(ax, point_from, point_to, color, lw_num, dashes=(None, None), linestyle='-'):
        ax.add_line(mlines.Line2D(
            xdata=(point_from[0], point_to[0]), ydata=(point_from[1], point_to[1]),
            color=color, lw=1 + 1.5 * lw_num / cap_max,
            linestyle=linestyle, dashes=dashes
        ))

    def _get_line_color(ac_dc, link):
        if not distinguish_ac_dc and not width_by_cap:
            return COLORS['outer']
        elif not distinguish_ac_dc and width_by_cap:
            return COLORS['ac']

        if ac_dc == "ac":
            if any([f'{i}_1' in link for i in OUTER_COUNTRIES]):
                return COLORS['outer']
            else:
                return COLORS['ac']
        elif ac_dc == "dc":
            return COLORS["dc"]

    centroids = gdf.centroid

    links_completed = []
    for link, cap in link_caps.iterrows():
        if sorted(link) in links_completed:  # every link comes up twice
            continue

        point_from = (centroids.loc[link[0]].x, centroids.loc[link[0]].y)
        point_to = (centroids.loc[link[1]].x, centroids.loc[link[1]].y)

        cap_ac = cap.filter(regex='ac_').sum(min_count=1)
        cap_dc = cap.filter(regex='dc_').sum(min_count=1)
        if distinguish_ac_dc:
            if ~np.isnan(cap_dc) and ~np.isnan(cap_ac):
                _add_line(
                    ax, point_from, point_to,
                    color=_get_line_color("ac", link),
                    lw_num=lw, linestyle='--', dashes=(10, 1)
                )
                _add_line(
                    ax, point_from, point_to,
                    color=_get_line_color("dc", link),
                    lw_num=lw, linestyle='-', dashes=(5, 4)
                )
            elif np.isnan(cap_dc) and ~np.isnan(cap_ac):
                _add_line(
                    ax, point_from, point_to, color=_get_line_color("ac", link),
                    lw_num=cap_ac
                )
            elif ~np.isnan(cap_dc) and np.isnan(cap_ac):
                _add_line(
                    ax, point_from, point_to, color=_get_line_color("dc", link),
                    lw_num=cap_dc,
                )
        else:
            cap_all = cap.sum(min_count=1)
            if np.isnan(cap_all) or cap_all == 0:
                continue
            if width_by_cap:
                lw = cap_all
            else:
                lw = 0.5
            _add_line(
                ax, point_from, point_to,
                color=_get_line_color("ac", link),
                lw_num=lw
            )

        links_completed.append(sorted(link))


def make_legend(ax):
    ax.axis('off')
    _handles = [
        mlines.Line2D(xdata=[0, 1], ydata=[0, 1], color=COLORS['ac'], linestyle='-'),
        mlines.Line2D(xdata=[0, 1], ydata=[0, 1], color=COLORS['dc'], linestyle='-'),
        (mlines.Line2D(xdata=[0, 1], ydata=[0, 1],
                       color=COLORS['ac'], linestyle='-', dashes=(10, 1)),
         mlines.Line2D(xdata=[0, 1], ydata=[0, 1],
                       color=COLORS['dc'], linestyle='--', dashes=(5, 4))),
        mlines.Line2D(xdata=[0, 1], ydata=[0, 1], color=COLORS['outer'], linestyle='-'),
        mpatches.Patch(facecolor=COLORS['eu_bg']),
        mpatches.Patch(facecolor=COLORS['outer_bg']),
    ]

    _labels = [
        'AC transmission', 'DC transmission', 'AC and DC transmission',
        'Transmission to/from outside model scope', 'Model region',
        'Regions outside model scope'
    ]

    ax.legend(
        handles=_handles, labels=_labels, frameon=False, loc='center', fontsize=12, ncol=3
    )


if __name__ == '__main__':
    plot_system(
        path_to_model=snakemake.input.model,
        path_to_units=snakemake.input.units,
        path_to_output=snakemake.output[0],
        bounds=snakemake.params.bounds
)
