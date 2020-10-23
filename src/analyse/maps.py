import geopandas as gpd
import pandas as pd
import numpy as np

import cartopy
import cartopy.crs as ccrs
from cartopy.io import shapereader

import matplotlib.pyplot as plt
from matplotlib import gridspec
import matplotlib.lines as mlines
import matplotlib.patches as mpatches
plt.rcParams['svg.fonttype'] = 'none'

import calliope

OUTER_COUNTRIES = ['RUS', 'UKR', 'TUR', 'MAR', 'DZA', 'TUN', 'LBY', 'EGY']

PROJ = ccrs.AlbersEqualArea()
CRS = ccrs.PlateCarree()

COLORS = {
    'outer': '#a9a9a9', 'ac': '#003366', 'dc': '#ffdf00',
    'eu_bg': '#add8e6', 'outer_bg': '#DCDCDC'
}


def plot_system(path_to_model, path_to_units, path_to_output):
    gdf = gpd.read_file(path_to_units)

    fig = plt.figure(figsize=(15, 16))
    gs = gridspec.GridSpec(2, 1, figure=fig, height_ratios=[20, 1], hspace=0)  # rows, columns

    extent = (-11, 35, 30, 72) # (minx, maxx, miny, maxy)
    centroid = [np.mean((extent[0], extent[1])), np.mean((extent[2], extent[3]))]

    ax_eu = plt.subplot(gs[0, 0], projection=PROJ)
    ax_eu.set_extent(extent)
    ax_eu.outline_patch.set_linewidth(0)

    ax_eu.add_geometries(
        gdf.geometry.values,
        crs=CRS, facecolor=COLORS['eu_bg'], edgecolor='None', linewidth=0.3
    )
    ax_eu.add_geometries(
        gdf.dissolve('country_code').geometry.values,
        crs=CRS, facecolor='None', edgecolor='white', linewidth=0.5
    )

    # get country borders
    shpfilename = shapereader.natural_earth(
        resolution='10m', category='cultural', name='admin_0_countries',
    )

    # read the shapefile using geopandas
    df = gpd.read_file(shpfilename)

    # read the outer country borders
    outer_gdf = df.loc[~df.ADM0_A3.isin(gdf.country_code.unique())]

    ax_eu.add_geometries(
        outer_gdf.geometry.values,
        crs=CRS, facecolor=COLORS['outer_bg'], edgecolor='white', linewidth=0.5
    )

    add_links(path_to_model, ax_eu)

    ax_leg = plt.subplot(gs[1, 0])
    make_legend(ax_leg)

    fig.savefig(path_to_output, bbox_inches='tight', pad_inches=1)


def add_links(path_to_model, ax):

    def _xy_from_longlat(lon, lat, ax):
        """convert long-lat (degrees) to grid values"""
        return ax.projection.transform_point(lon, lat, CRS)

    def _add_line(ax, point_from, point_to, color, lw_num, dashes=(None, None), linestyle='-'):
        ax.add_line(mlines.Line2D(
            xdata=(point_from[0], point_to[0]), ydata=(point_from[1], point_to[1]),
            color=color, lw=1 + 1.5 * lw_num / cap_max,
            linestyle=linestyle, dashes=dashes
        ))

    m = calliope.read_netcdf(path_to_model)
    links = m.inputs.energy_cap_equals.loc[
        {'loc_techs': m._model_data.loc_techs_transmission}
    ].to_pandas()
    links.index = (
        links.index.str.split(':', expand=True).droplevel(1).rename(('to', 'tech', 'from'))
    )
    links = links.unstack(1)  # empty column from splitting '::'

    link_coords = m._model_data.loc_coordinates.to_pandas()
    link_lats = link_coords.loc['lat']
    link_lons = link_coords.loc['lon']
    cap_max = links.max().max()

    links_completed = []
    for link, cap in links.iterrows():
        if sorted(link) in links_completed:  # every link comes up twice
            continue

        point_from = _xy_from_longlat(link_lons[link[0]], link_lats[link[0]], ax)
        point_to = _xy_from_longlat(link_lons[link[1]], link_lats[link[1]], ax)

        cap_ac = cap['ac_transmission']
        cap_dc = cap['dc_transmission']

        if any([f'{i}_1' in link for i in OUTER_COUNTRIES]):
            ac_color = COLORS['outer']
        else:
            ac_color = COLORS['ac']

        if ~np.isnan(cap_dc) and ~np.isnan(cap_ac):
            lw = (cap_ac + cap_dc)
            _add_line(
                ax, point_from, point_to,
                color=ac_color, lw_num=lw, linestyle='--', dashes=(10, 1)
            )
            _add_line(
                ax, point_from, point_to,
                color=COLORS['dc'], lw_num=lw, linestyle='-', dashes=(5, 4)
            )
        elif np.isnan(cap_dc) and ~np.isnan(cap_ac):
            _add_line(
                ax, point_from, point_to, color=ac_color, lw_num=cap_ac
            )
        elif ~np.isnan(cap_dc) and np.isnan(cap_ac):
            _add_line(
                ax, point_from, point_to, color=COLORS['dc'], lw_num=cap_dc,
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
)