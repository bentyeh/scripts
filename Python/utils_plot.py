import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt

# region --- matplotlib

def hist_colorbar(values, names, bins=10, colors='viridis', nticks=None,
                  implementation='bar', kwargs=None, colorbar_args=None, ax=None):
    '''
    Plot histogram with each data point colored based on its name.

    Args
    - values: float
        Values to plot on the histogram
    - names: str
        Names of each value defining the group (and color) to which each value belongs
    - bins: int, or list of float. default=10
        int: Number of bins
        list: Bin edges (see matplotlib.pyplot.hist() documentation)
        - if using 'imshow' `implementation`, bin widths must be identical
    - colors: str or numpy.ndarray of shape (N, 3) or (N, 4). default='viridis'
        Colormap name, or RGB(A) array
    - nticks: int. default=None
        Number of ticks on the colorbar
    - implementation: str. default='bar'
        Implement using matplotlib.axes.Axes.bar ('bar'),
        matplotlib.axes.Axes.imshow ('imshow'), or pandas.DataFrame.plot.scatter ('scatter')
    - kwargs: dict. default=None
        `implementation` == 'bar': arguments to pass to matplotlib.axes.Axes.bar
          Common arguments: align ('center', 'edge'), log (False, True)
        `implementation` == 'imshow': arguments to pass to matplotlib.axes.Axes.imshow
          Common arguments: align ('center', 'edge'), log (False, True)
        `implementation` == 'scatter': arguments to pass to pandas.DataFrame.plot.scatter
    - colorbar_args: dict. default=None
        Keyword arguments to pass to matplotlib.figure.Figure.colorbar
        Common arguments: label (str), orientation ('vertical', 'horizontal')
    - ax: matplotlib.axes.Axes. default=None
        Axes on which to make histogram

    Returns: matplotlib.axes.Axes, matplotlib.colorbar.Colorbar
      Axes on which histogram is plotted, and colorbar object (to get colorbar axes: cbar.ax).
      The colorbar object is returned as None if `implementation` is 'scatter'.
    '''
    assert len(values) == len(names)
    assert np.issubdtype(type(bins), np.integer) or isinstance(bins, list)
    assert implementation in ('bar', 'imshow', 'scatter')

    # create dataframe with x and y bar values
    df = pd.DataFrame({'values': values, 'names': names, 'x': None, 'y': None})
    bin_map = {}
    if np.issubdtype(type(bins), np.integer):
        nbins = bins
        min_value = df['values'].min()
        max_value = df['values'].max()
        width = (max_value - min_value) / bins
        df['bin'] = ((df['values'] - min_value) / (max_value - min_value) // (1 / bins)).astype(int)
        df['x'] = df['bin'] * width + min_value
        bin_map = df[['bin', 'x']].drop_duplicates().set_index('x').squeeze().to_dict()
    else:
        width = np.diff(bins)
        nbins = len(width)
        for i in range(nbins):
            bin_center = (bins[i] + bins[i+1])/2
            bin_map[bin_center] = i
            if i == nbins:
                df.loc[(df['values'] >= bins[i]) & (df['values'] < bins[i+1]), 'x'] = bin_center
            else:
                df.loc[(df['values'] >= bins[i]) & (df['values'] <= bins[i+1]), 'x'] = bin_center
        if df['x'].isna().sum() > 0:
            print('Some values were outside the provided bin edges and will not be plotted.')
            df.dropna(subset='x', inplace=True)
    df['y'] = df.groupby('x').transform(lambda group: np.arange(group.shape[0]))

    if implementation == 'imshow':
        assert len(np.unique(width)) == 1

    unames = df['names'].unique() # pandas.Series.unique() returns unique names in order
    n_unames = len(unames)
    if ax is None and implementation in ('bar', 'imshow'):
        fig, ax = plt.subplots(constrained_layout=True)
    if kwargs is None:
        kwargs = {}
    if colorbar_args is None:
        colorbar_args = {}
    if isinstance(colors, str):
        colors = plt.get_cmap(colors, n_unames).colors
        colors_map = {}
        for i in range(n_unames):
            colors_map[unames[i]] = colors[i,:]
        colors = np.vstack(df['names'].map(colors_map))

    if implementation == 'bar':
        patches = ax.bar(x=df['x'], height=1, width=width, bottom=df['y'], color=colors, **kwargs)
        p = matplotlib.collections.PatchCollection(patches)
        p.set_array(colors)
    elif implementation == 'imshow':
        mat = np.ones([df['y'].max() + 1, nbins]) * np.nan
        mat[df['y'], df['x'].map(bin_map)] = np.arange(n_unames)
        p = ax.imshow(mat, origin='lower', aspect='auto', **kwargs)
        transformer = matplotlib.transforms.Affine2D().scale(np.unique(width)[0], 1).translate(min(bin_map.keys()), 0.5)
        p.set_transform(transformer + ax.transData)
        ax.set_ylim(0, df['y'].max()+1.5)
        ax.set_xlim(transformer.transform_point((-1, 0))[0], transformer.transform_point((nbins, 0))[0])
    else:
        df['y'] += 1
        ax = df.plot.scatter('x', 'y', c='names', ax=ax, colormap=matplotlib.colors.ListedColormap(colors), **kwargs)
        ax.set(ylabel=None)
        ax.tick_params(reset=True, top=False, right=False)
        # df.plot.scatter appends the colorbar axes to the end of the figure axes
        ax.figure.axes[-1].set_ylabel(None)
        return ax, None
    cbar = fig.colorbar(p, ax=ax, **colorbar_args)
    if nticks is None:
        nticks = len(cbar.get_ticks())
    ticklabels = df['names'].iloc[::int(np.ceil(len(df['names']) / nticks))]
    if implementation == 'bar':
        cbar.set_ticks(np.arange(0, 1, 1/nticks))
    if cbar.orientation == 'vertical':
        cbar.ax.set_yticklabels(ticklabels)
    else:
        cbar.ax.set_xticklabels(ticklabels)
    return ax, cbar

# endregion --- matplotlib

# region ------ seaborn

def fix_colorbar_height(axs):
    '''
    Set colorbars to same height as their corresponding heatmaps.

    Arg: list of matplotlib.axes.Axes
    - First n axes are heatmaps, last n axes are colorbars
    - Example: axes in a figure in which a heatmap was drawn with seaborn.heatmap()

    Returns: None
      Axes objects are modified inplace.
    '''
    n_heatmaps = int(len(axs) / 2)
    for k in range(n_heatmaps):
        x0 = axs[n_heatmaps + k].get_position().bounds[0]
        y0 = axs[k].get_position().bounds[1]
        width = axs[n_heatmaps + k].get_position().bounds[2]
        height = axs[k].get_position().bounds[3]
        axs[n_heatmaps + k].set_position([x0, y0, width, height])

# endregion --- seaborn

# region ------ pandas

def background_gradient(self, cvals=None, cmin=None, cmax=None, cmap='viridis', **css):
    '''
    Apply a heatmap color gradient elementwise to a dataframe.
    Usage: df.style.apply(background_gradient, axis=None, cmap='Oranges', cmin=0, cmax=1)

    Args
    - self: pandas.DataFrame
        The calling DataFrame. This argument is automatically passed in by the
        `pandas.DataFrame.style.apply` method
    - cvals: pandas.DataFrame. default=None
        If specified this DataFrame is used to determine the color gradient
    - cmin: float. default=None
        If specified, any values below this will be clipped to the bottom of the cmap
    - cmax: float. default=None
        If specified, any values above this will be clipped to the top of the cmap
    - cmap: colormap or str. default='viridis'
        The colormap to use
    - css: dict
        Extra inline css key/value pairs to pass to the styler

    Returns: pd.DataFrame

    Example: 

    Source: https://github.com/pandas-dev/pandas/issues/15204#issuecomment-274730740
    '''
    if cvals is None:
        cvals = self.values.ravel().copy()
    else:
        assert cvals.shape == self.shape
        cvals = cvals.values.ravel().copy()
    cvals -= cmin or cvals.min()
    cvals /= cmax or cvals.max()
    cvals = cvals.clip(0, 1)
    styles = []
    for rgb in plt.get_cmap(cmap)(cvals):
        style = [
            "{}: {}".format(key.replace('_', '-'), value)
            for key, value in css.items()
        ]
        style.append("background-color: {}".format(matplotlib.colors.rgb2hex(rgb)))
        styles.append('; '.join(style))
    styles = np.asarray(styles).reshape(self.shape)
    return pd.DataFrame(styles, index=self.index, columns=self.columns)

# endregion --- pandas
