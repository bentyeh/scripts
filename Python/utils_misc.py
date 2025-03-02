import multiprocessing, os, sys
sys.path.append(os.path.abspath(os.path.dirname(__file__)))
import numpy as np
import pandas as pd
import utils

def intervals_merge(intervals):
    '''
    Merge intervals.

    Args
    - intervals: list of 2-tuples of real numbers
        List of (start, stop) intervals

    Returns: list of 2-tupes of real numbers
    '''
    n_intervals = len(intervals)
    if n_intervals < 2:
        return intervals

    start_dtype, end_dtype = int, int
    if any([np.issubdtype(type(interval[0]), np.floating) for interval in intervals]):
        start_dtype = float
    if any([np.issubdtype(type(interval[1]), np.floating) for interval in intervals]):
        end_dtype = float
    values = np.array(intervals, dtype=[('start', start_dtype), ('end', end_dtype)])
    values.sort(axis=0, order=['start', 'end'])
    merged = []
    start, end = values[0]
    for i in range(1, n_intervals):
        if values[i][0] > end:
            merged.append((start, end))
            start = values[i][0]
        end = values[i][1]
        if i == n_intervals - 1:
            merged.append((start, end))
    return merged

def intervals_weightedOverlap(intervals):
    '''
    Return non-overlapping set of regions whose weights are the sum of the weights of
    overlapping input intervals spanning that region

    Args
    - intervals: list of 2- or 3-tuples of real numbers
        (start, end[, weight]) where start < end
        weight assumed to be 1 if not provided

    Returns: 3-tuple of (list of 3-tuples of real numbers, list of list of real numbers, list of list of integers)
    - List of non-overlapping intervals (start, end, weight) where start < end
      - These non-overlapping intervals span the entire range of the input intervals.
        Regions in this range not overlapped by at least 1 input interval are included
        with weight 0.
    - List of list of weights of original intervals contributing to new interval weights
    - List of list of indices of original intervals contributing to new intervals
    '''
    n = len(intervals)

    values = np.array(intervals)
    s = values[:,0]
    e = values[:,1]
    if values.shape[1] != 3:
        w = np.ones(s.shape)
    else:
        w = values[:,2]

    s_argsort = np.argsort(s)
    s_sorted = s[s_argsort]
    e_argsort = np.argsort(e)
    e_sorted = e[e_argsort]

    start = s_sorted[0]        # current start
    weight = [w[s_argsort[0]]] # current weight
    s_idx = 1
    e_idx = 0

    newIntervals = []     # list of non-overlapping intervals
    weights = []          # list of list of weights of original intervals contributing to new interval weights
    overlapIntervals = [] # list of list of indices of original intervals contributing to new intervals
    curOverlap = [s_argsort[0]]

    while e_idx < n:
        if (s_idx < n) and (s_sorted[s_idx] < e_sorted[e_idx]):
            end = s_sorted[s_idx]
            newIntervals.append((start, end))
            weights.append(list(weight))
            overlapIntervals.append(list(curOverlap))
            weight.append(w[s_argsort[s_idx]])
            curOverlap.append(s_argsort[s_idx])
            s_idx += 1
        else:
            end = e_sorted[e_idx]
            overlapIntervals.append(list(curOverlap))
            weights.append(list(weight))
            newIntervals.append((start, end))
            weight.remove(w[e_argsort[e_idx]])
            curOverlap.remove(e_argsort[e_idx])
            e_idx += 1
        start = end

    # post-processing: remove length-0 intervals where start == end
    validIntervals = [i for i in range(len(newIntervals)) if newIntervals[i][0] < newIntervals[i][1]]
    weights = [weights[i] for i in validIntervals]
    newIntervals = [newIntervals[i] for i in validIntervals]
    overlapIntervals = [overlapIntervals[i] for i in validIntervals]

    return (newIntervals, weights, overlapIntervals)

def intervals_value(intervals, pos, check_intervals=True):
    '''
    Get value at a specific position.

    Args
    - intervals: list of 3-tuples of float/int, or pandas.DataFrame
        Non-overlapping intervals (start, end, value) where start < end
    - pos: numpy.ndarray or list of float/int
        positions at which to look-up value
    - check_intervals: bool. default=True
        Perform basic checks that intervals are non-overlapping and start < end

    Returns: list of float/int or None
      Values at specified positions. If a position is not in the range of intervals,
      its corresponding value is None.

    Note: If intervals are overlapping, run intervals_weightedOverlap(intervals) first.
    '''
    # wrangle input arguments to desired form
    if isinstance(intervals, pd.DataFrame):
        df = intervals
    else:
        df = pd.DataFrame(intervals)
    df = df.sort_values(by=[df.columns[0], df.columns[1]], axis=0)
    start = df.iloc[:, 0].values
    end = df.iloc[:, 1].values
    value = df.iloc[:, 2].values

    pos = np.asarray(pos).reshape(-1,1) # reshape to column vector

    # verify that intervals are non-overlapping
    if check_intervals:
        assert np.all(start < end)
        assert np.all(end[:-1] <= start[1:])

    array_idx = (start <= pos) & (end > pos)
    in_intervals = np.sum(array_idx, axis=1, dtype=np.bool)

    pos_values = np.full(len(pos), None)
    pos_values[in_intervals] = value[np.where(array_idx)[1]]
    return pos_values

def df_split_row(df, col, sep, keep=False):
    '''
    Split the values of a column and expand so the new DataFrame has one split
    value per row. Drops rows where the values in the column is NA.

    Args
    - df: pandas.DataFrame
        Dataframe with the column to split and expand
    - col: str
        Name of the column to split and expand.
    - sep: str
        The string used to split the column's values
    - keep: bool
        Retain the presplit value as its own row

    Returns: pandas.DataFrame

    Source: https://stackoverflow.com/a/39946744

    Note: There is no universal inverse operation. For example, if the original DataFrame
    was already partially "split", it is impossible to know how to recombine that
    partially split structure.

          col1 col2                                    col1 col2            col1  col2
        0    a  b,c   df_split_row(df, 'col2', ',')  0    a    b  invert  0    a b,c,d
        1    a    d   ---------------------------->  1    a    c  ----->
                                                     2    a    d
    A general combine can be performed as follows:
      df.groupby(list(set(df.columns) - set([col]))) \
        .agg(lambda x: sep.join(x)) \
        .reset_index()
    '''
    indexes = list()
    new_values = list()
    df = df.dropna(subset=[col])
    for i, presplit in enumerate(df[col].astype(str)):
        values = presplit.split(sep)
        if keep and len(values) > 1:
            indexes.append(i)
            new_values.append(presplit)
        for value in values:
            indexes.append(i)
            new_values.append(value)
    new_df = df.iloc[indexes, :].copy()
    new_df[col] = new_values
    return new_df

def pandas_apply_parallel(grouped, func, nproc=None, concat_axis=None, **kwargs):
    '''
    Parallel pandas.DataFrame.apply

    Args
    - grouped: pandas.Series or pandas.DataFrame
    - func: function
        Accepts a pandas.Series or pandas.DataFrame
    - nproc: int. default=None
        Number of processes to use. If None, defaults to the number of usable CPUs.
    - concat_axis: 0 or 1. default=None
        None: Return a dict where names are the group names and values are the result of
          applying `func` to a group.
        0 or 1: Return a pandas DataFrame consisting of the results concatenated along the axis
          given by `concat_axis'
    - **kwargs
        Additional keyword arguments (e.g., kwds, callback, error_callback) to pass to
        multiprocessing.Pool.apply_async().

    Returns: pandas.DataFrame or dict
      See the `concat_axis` argument.

    Based on https://stackoverflow.com/a/29281494
    '''
    if nproc is None:
        nproc = utils.get_available_cpus()

    with multiprocessing.Pool(nproc) as pool:
        result_dict = {name: pool.apply_async(func, (group,), **kwargs) for name, group in grouped}
        pool.close()
        pool.join()

    results = {name: result.get() for name, result in result_dict.items()}
    if concat_axis is None:
        return results
    else:
        df = pd.concat(results, axis=concat_axis)
        if concat_axis:
            df = df.transpose()
        return df
