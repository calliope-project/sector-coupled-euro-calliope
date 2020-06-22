import os

import pandas as pd
import numpy as np

import matplotlib
import matplotlib.cm as cm
import matplotlib.pyplot as plt
import seaborn as sns
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

idx = pd.IndexSlice


def heatmap(correlation_df, palette=None, n_colors=256, colorscale_bounds=[-1, 1], fig=None, colorbar_width_scale=1, fontsize=14, xrotation=45, yrotation=0, xalignment='right', yalignment='center', reverse_y=True):
    """
    Hack of https://towardsdatascience.com/better-heatmaps-and-correlation-matrix-plots-in-python-41445d0f2bec

    Parameters
    ----------
    correlation_df : pandas.DataFrame
        The result of running 'df.corr()' on a pandas dataframe.
    palette : seaborn palette or None
        If None, will default to sns.diverging_palette(220, 20, n=n_colors)
    """
    _df = correlation_df.stack().reset_index()
    _df.columns = ['x', 'y', 'value']
    x = _df.x
    y = _df.y
    size = _df['value'].abs()
    color = _df['value']

    if palette is None:
        palette = sns.diverging_palette(220, 20, n=n_colors)  # Create the palette
    color_min, color_max = colorscale_bounds  # Range of values that will be mapped to the palette, i.e. min and max possible correlation

    def value_to_color(val):
        if np.isnan(val):
            return (0, 0, 0, 0)
        val_position = float((val - color_min)) / (color_max - color_min)  # position of value in the input range, relative to the length of the input range
        ind = int(val_position * (n_colors - 1))  # target index in the color palette
        return palette[ind]

    with sns.axes_style('whitegrid'):
        if fig is None:
            fig = plt.figure(figsize=(10, 10))
        gs = plt.GridSpec(1, 2, hspace=0.2, wspace=0.1, width_ratios=(20, 1 * colorbar_width_scale))  # Setup a 1x2 grid
        ax = plt.subplot(gs[:, 0])  # Use the leftmost column of the grid for the main plot

        # Mapping from column names to integer coordinates
        x_labels = x.drop_duplicates().reset_index().drop('index', axis=1).x
        if reverse_y is True:
            y_labels = y.drop_duplicates().sort_index(ascending=False).reset_index().drop('index', axis=1).y
        else:
            y_labels = y.drop_duplicates().reset_index().drop('index', axis=1).y
        x_to_num = {v: k for k, v in x_labels.to_dict().items()}
        y_to_num = {v: k for k, v in y_labels.to_dict().items()}

        size_scale = 100 / color_max
        ax.scatter(
            x=x.map(x_to_num),  # Use mapping for x
            y=y.map(y_to_num),  # Use mapping for y
            s=size * size_scale,  # Vector of square sizes, proportional to size parameter
            marker='s',  # Use square as scatterplot marker
            c=color.apply(value_to_color),
        )

        # Show column labels on the axes
        ax.set_xticks([x_to_num[v] for v in x_labels.values])
        ax.set_xticklabels(x_labels, rotation=xrotation, horizontalalignment=xalignment, fontsize=fontsize)
        ax.set_yticks([y_to_num[v] for v in y_labels.values])
        ax.set_yticklabels(y_labels, rotation=yrotation, fontsize=fontsize, verticalalignment=yalignment)
        ax.grid(False, 'major')
        ax.grid(True, 'minor')
        ax.set_xticks([t + 0.5 for t in ax.get_xticks()], minor=True)
        ax.set_yticks([t + 0.5 for t in ax.get_yticks()], minor=True)
        ax.set_xlim([-0.5, max([v for v in x_to_num.values()]) + 0.5])
        ax.set_ylim([-0.5, max([v for v in y_to_num.values()]) + 0.5])

    with sns.axes_style('darkgrid'):
        # Add color legend on the right side of the plot
        ax = plt.subplot(gs[:, 1])  # Use the rightmost column of the plot
        # Widths match the square sizes in the actual correlation plot (note that e size matching has been configured empirically)
        widths = np.concatenate((np.linspace(1, 0, int(len(palette) / 2)), np.linspace(0, 1, int(len(palette) / 2))))
        zero_x = [1.5 - i / 2 for i in widths]
        bar_y = np.linspace(color_min, color_max, n_colors)  # y coordinates for each of the n_colors bars

        bar_height = bar_y[1] - bar_y[0]

        ax.barh(
            y=bar_y,
            width=widths,
            left=zero_x,
            height=bar_height,
            color=palette,
            linewidth=0
        )
        ax.set_xlim(1, 2)  # Bars are going from 0 to 5, so lets crop the plot somewhere in the middle
        ax.grid(False)  # Hide grid
        ax.set_facecolor('white')  # Make background white
        ax.set_xticks([])  # Remove horizontal ticks
        ax.set_yticks(np.linspace(min(bar_y), max(bar_y), 3))  # Show vertical ticks for min, middle and max
        ax.tick_params(axis='y', labelsize=fontsize)
        ax.yaxis.tick_right()  # Show vertical ticks on the right

    return fig
