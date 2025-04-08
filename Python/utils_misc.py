import itertools
import collections
import collections.abc
import multiprocessing
import os
import sys
import typing

import numpy as np
import pandas as pd

sys.path.append(os.path.abspath(os.path.dirname(__file__)))
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
    in_intervals = np.sum(array_idx, axis=1, dtype=bool)

    pos_values = np.full(len(pos), None)
    pos_values[in_intervals] = value[np.where(array_idx)[1]]
    return pos_values


def _advance(iterator):
    """Get the next item from an iterator, returning None if exhausted."""
    try:
        return next(iterator)
    except StopIteration:
        return None


def overlap_fraction(start1: int, end1: int, start2: int, end2: int) -> tuple[float, float]:
    """
    Return the fraction of overlap between two intervals, with respect to the length of the first and second intervals,
    respectively. Arguments are assumed to represent 0-based, start-closed, end-open intervals.
    """
    overlap = max(0, min(end1, end2) - max(start1, start2))
    len1 = end1 - start1
    len2 = end2 - start2
    return (overlap / len1), (overlap / len2)


def overlap_length(start1: int, end1: int, start2: int, end2: int, allow_negative=False) -> int:
    """
    Return the length of the overlap between two intervals.
    Arguments are assumed to represent 0-based, start-closed, end-open intervals.
    """
    if allow_negative:
        return min(end1, end2) - max(start1, start2)
    else:
        return max(0, min(end1, end2) - max(start1, start2))


def simple_overlap(start1: int, end1: int, start2: int, end2: int) -> bool:
    """
    Check if two intervals overlap. Arguments are assumed to represent 0-based, start-closed, end-open intervals.
    """
    return (start1 < end2) and (start2 < end1)


def intervals_overlap(
    intervals1: collections.abc.Iterable[tuple[typing.Any, int, int]],
    intervals2: collections.abc.Iterable[tuple[typing.Any, int, int]],
    values1: collections.abc.Iterable | None = None,
    values2: collections.abc.Iterable | None = None,
    overlap_func: collections.abc.Callable | None = None,
    how: str = "inner",
    max_distance: int = 0,
) -> collections.abc.Iterable[tuple[tuple[typing.Any, int, int] | None, tuple[typing.Any, int, int] | None, typing.Any, typing.Any]]:
    """
    Find overlapping intervals between two sets of intervals.
    Implemented lazily as a generator to be as memory efficient as possible.

    Args
    - intervals1: The first set of intervals, where each interval is represented as a tuple (group, start, end).
        Assumptions: intervals are sorted by start and end; group is ignored.
        start and end are 0-based, start-closed, end-open.
        Intervals can be overlapping; each interval is treated independently.
    - intervals2: The second set of intervals, where each interval is represented as a tuple (group, start, end).
        Assumptions: intervals are sorted by start and end; group is ignored.
        start and end are 0-based, start-closed, end-open.
        Intervals can be overlapping; each interval is treated independently.
    - values1: The data associated with the first set of intervals. If None, treated as a repeating iterable of None.
    - values2: The data associated with the second set of intervals. If None, treated as a repeating iterable of None.
    - overlap_func: A function to determine if two intervals overlap. If None, defaults to simple_overlap().
        Must take 4 arguments (start1, end1, start2, end2) and return a boolean.
        Examples:
          - Require minimum overlap of 4 positions:
              overlap_func=lambda *args: overlap_length(*args) >= 4
          - Require minimum overlap of 50%, mutually:
              overlap_func=lambda *args: min(overlap_fraction(*args)) >= 0.5
          - Consider intervals overlapping even if they are separated by up to distance 5 positions:
              overlap_func=lambda *args: overlap_length(*args, allow_negative=True) > -5, max_distance=5
                This will consider the following intervals as overlapping:
                - interval1: (0, 0, 10):  ----------
                                                    |||| 
                - interval2: (0, 14, 20):               ------
    - how: Which intervals to report
        - 'left': report all intervals in intervals1, including those without overlapping intervals in intervals2
            - The same interval in interval1 may be reported multiple times if it overlaps with multiple intervals in
            intervals2.
            - Intervals without overlaps will have None for the other interval and value: (interval1, value1, None, None)
        - 'right': report all intervals in intervals2, including those without overlapping intervals in intervals1
            - The same interval in interval2 may be reported multiple times if it overlaps with multiple intervals in
            intervals1.
            - Intervals without overlaps will have None for the other interval and value: (None, None, interval2, value2)
        - 'inner': only report overlapping intervals
        - 'outer': report all intervals, including those without overlapping intervals in the other set
    - max_distance: maximum distance between intervals such that they could be considered overlapping
        - Only relevant if a custom overlap_func is provided that may consider two overlapping intervals as overlapping
          even if separated by up to max_distance positions.
        - Two intervals that are adjacent (end1 == start2) but non-overlapping are considered to be separated by
          distance 1.

    Returns: iterable of (interval1, interval2, value1, value2)
    """
    if how not in ('left', 'right', 'inner', 'outer'):
        raise ValueError(f"how='{how}' not supported; must be one of 'left', 'right', 'inner', 'outer'")
    if values1 is None:
        values1 = itertools.repeat(None)
    if values2 is None:
        values2 = itertools.repeat(None)
    if overlap_func is None:
        overlap_func = simple_overlap

    iter1 = enumerate(zip(intervals1, values1)) # yields id1, ((group1, start1, end1), value1)
    iter2 = enumerate(zip(intervals2, values2)) # yields id2, ((group2, start2, end2), value2)

    id1, (interval1, value1) = _advance(iter1) or (None, (None, None))

    d = collections.deque() # deque of id2, ((group2, start2, end2), value2)
    id2s_seen = set() # track which interval2s in the deque have previously been yielded

    while id1 is not None:
        group1, start1, end1 = interval1
        id1_seen = False

        # pop all lesser non-overlapping intervals (end2 <= start1) off of the deque
        while d and d[0][1][0][2] <= (start1 - max_distance):
            id2, (interval2, value2) = d.popleft()
            if id2 in id2s_seen:
                id2s_seen.remove(id2) # id2 will never be seen again; can remove it from the set to free up memory
            elif how in ('right', 'outer'):
                yield None, interval2, None, value2

        # Check for overlaps with remaining intervals in the deque
        last_deque_interval2_exceeds_end1 = False
        for id2, (interval2, value2) in d:
            group2, start2, end2 = interval2
            if overlap_func(start1, end1, start2, end2):
                yield interval1, interval2, value1, value2
                id2s_seen.add(id2)
                id1_seen = True
            elif start2 >= end1 + max_distance:
                last_deque_interval2_exceeds_end1 = True
                break

        # Add new interval2s to the deque until one exceeds the end of interval1
        # - If a new interval2 overlap with interval1, yield them
        # - Note that it is possible for one interval2 to overlap with interval1, a subsequent interval2 not to overlap
        #   with interval1, and then a third interval2 to overlap with interval1 again.
        #   - Example
        #       id1, interval1 = 0, (0, 10, 20)           ----------
        #       id2, interval2 = 0, (0, 5, 15)       ----------
        #       id2, interval2 = 1, (0, 6, 7)         -
        #       id2, interval2 = 2, (0, 19, 25)                    ------
        if not last_deque_interval2_exceeds_end1:
            id2, (interval2, value2) = _advance(iter2) or (None, (None, None))
            while id2 is not None:
                group2, start2, end2 = interval2
                d.append((id2, (interval2, value2)))
                if start2 >= end1 + max_distance:
                    break
                elif overlap_func(start1, end1, start2, end2):
                    yield interval1, interval2, value1, value2
                    id2s_seen.add(id2)
                    id1_seen = True
                id2, (interval2, value2) = _advance(iter2) or (None, (None, None))

        # Advance the iterator for intervals1
        if not id1_seen and how in ('left', 'outer'):
            yield interval1, None, value1, None
        id1, (interval1, value1) = _advance(iter1) or (None, (None, None))

    if how in ('right', 'outer'):
        # process remaining intervals2 in the deque
        while d:
            id2, (interval2, value2) = d.popleft()
            if id2 not in id2s_seen:
                yield None, interval2, None, value2

        # process remaining intervals2 in iter2
        for id2, (interval2, value2) in iter2:
            yield None, interval2, None, value2


def intervals_filter_by_overlap(
    intervals1: collections.abc.Iterable[tuple[typing.Any, int, int]],
    intervals2: collections.abc.Iterable[tuple[typing.Any, int, int]],
    values1: collections.abc.Iterable | None = None,
    how: str = "semi",
    **kwargs
) -> collections.abc.Iterable[tuple[tuple[typing.Any, int, int], typing.Any]]:
    """
    Filtering joins: filter intervals in intervals1 based on the presence or absence of overlaps in intervals2

    This implementation is not optimally efficient, as it is a simple wrapper around intervals_overlap_grouped(), which
    tries to consider all pairwise overlaps, whereas semi-join and anti-join effectively only care about the first
    overlap.

    Args
    - intervals1: see intervals_overlap_grouped()
    - intervals2: see intervals_overlap_grouped()
    - values1: see intervals_overlap_grouped()
    - how: which intervals to report
        - 'semi': only intervals in intervals1 that overlap with intervals2
        - 'anti': only intervals in intervals1 that do not overlap with intervals2
    - **kwargs: passed to intervals_overlap_grouped()

    Yields: interval1, value1
    """
    if how == 'semi':
        if values1 is None:
            values1 = itertools.repeat(None)
        # add id to each value1 to track which intervals have been yielded
        last_id1_seen = -1
        for interval1, _, value1, _ in intervals_overlap_grouped(intervals1, intervals2, values1=enumerate(values1), how='inner', **kwargs):
            id1, value = value1
            if id1 > last_id1_seen:
                last_id1_seen = id1
                yield interval1, value
    elif how == 'anti':
        for interval1, interval2, value1, _ in intervals_overlap_grouped(intervals1, intervals2, values1=values1, how='left', **kwargs):
            if interval2 is None:
                yield interval1, value1
    else:
        raise ValueError(f"how='{how}' not supported; must be 'semi' or 'anti'")


def intervals_overlap_grouped(
    intervals1: collections.abc.Iterable[tuple[typing.Any, int, int]],
    intervals2: collections.abc.Iterable[tuple[typing.Any, int, int]],
    values1: collections.abc.Iterable | None = None,
    values2: collections.abc.Iterable | None = None,
    how: str = "inner",
    **kwargs
) -> collections.abc.Iterable[tuple[tuple[typing.Any, int, int] | None, tuple[typing.Any, int, int] | None, typing.Any, typing.Any]]:
    """
    Find overlapping intervals between two sets of intervals. Only intervals in the same group are compared.

    Args: see intervals_overlap()
    - Only additional assumptions are that intervals1 and intervals2 are sorted by group, start and end.

    Returns: iterable of (interval1, value1, interval2, value2)

    Example usage
    1. Overlap between 2 sets of intervals in BED6 format, using in-memory pandas DataFrames. Only intervals on the
       same chromosome and strand are compared.

        CHROMS = [f'chr{i}' for i in range(1, 22)] + ['chrX', 'chrY', 'chrM']
        DTYPE_CHROMS = pd.CategoricalDtype(categories=CHROMS, ordered=True)
        DTYPE_STRANDS = pd.CategoricalDtype(categories=['+', '-'], ordered=True)
        COLS_BED6 = ['chrom', 'start', 'end', 'name', 'score', 'strand']
        dtypes = dict(chrom=pd.CategoricalDtype(), start=int, end=int, name=str, score=int, strand=DTYPE_STRANDS)
        df1 = pd.read_csv('input1.bed', sep='\t', header=None, names=COLS_BED6, dtype=dtypes).sort_values(['chrom', 'strand', 'start', 'end'])
        df2 = pd.read_csv('input2.bed', sep='\t', header=None, names=COLS_BED6, dtype=dtypes).sort_values(['chrom', 'strand', 'start', 'end'])
        intervals_overlap_grouped(
            intervals1=zip(zip(df1['chrom'], df1['strand']), df1['start'], df1['end']),
            intervals2=zip(zip(df2['chrom'], df2['strand']), df2['start'], df2['end']),
            values1=zip(df1['name'], df1['score']),
            values2=zip(df2['name'], df2['score']),
            overlap_func=simple_overlap,
            how='inner'
        )

    2. Intersect 2 sets of intervals in BED3 format, streaming from disk. Write intersecting regions to standard output.
       Analogous to bedtools intersect -a input1.bed -b input2.bed.

        def parse_bed3(path):
            with open(path) as f:
                for line in f:
                    line = line.strip()
                    if line == '':
                        continue
                    chrom, start, end = line.split('\t')
                    yield (chrom, int(start), int(end))
        intervals1 = parse_bed3('input1.bed')
        intervals2 = parse_bed3('input2.bed')
        for (interval1, interval2, _, _) in intervals_overlap_grouped(intervals1, intervals2, how='inner'):
            chrom1, start1, end1 = interval1
            chrom2, start2, end2 = interval2
            print(chrom1, max(start1, start2), min(end1, end2), sep='\t')
    """
    if how not in ('left', 'right', 'inner', 'outer'):
        raise ValueError(f"how='{how}' not supported; how must be one of 'left', 'right', 'inner', 'outer'")
    if values1 is None:
        values1 = itertools.repeat(None)
    if values2 is None:
        values2 = itertools.repeat(None)

    iter1 = zip(intervals1, values1) # yields (group1, start1, end1), value1
    iter2 = zip(intervals2, values2)

    grouped1 = itertools.groupby(iter1, key=lambda x: x[0][0]) # yields group1, iterable((group1, start1, end1), value1)
    grouped2 = itertools.groupby(iter2, key=lambda x: x[0][0])

    group1, sub_iter1 = _advance(grouped1) or (None, None)
    group2, sub_iter2 = _advance(grouped2) or (None, None)

    while group1 is not None or group2 is not None:
        if group2 is None or (group1 is not None and group1 < group2):
            # Group only exists in intervals1
            if how in ('left', 'outer'):
                for interval1, value1 in sub_iter1:
                    yield interval1, None, value1, None
            group1, sub_iter1 = _advance(grouped1) or (None, None)
        elif group1 is None or (group2 is not None and group2 < group1):
            # Group only exists in intervals2
            if how in ('right', 'outer'):
                for interval2, value2 in sub_iter2:
                    yield None, interval2, None, value2
            group2, sub_iter2 = _advance(grouped2) or (None, None)
        else: # group1 == group2, groups match
            # Process the matched groups using the helper function
            intervals1, values1 = zip(*sub_iter1)
            intervals2, values2 = zip(*sub_iter2)
            yield from intervals_overlap(intervals1, intervals2, values1=values1, values2=values2, how=how, **kwargs)

            # Advance both group iterators
            group1, sub_iter1 = _advance(grouped1) or (None, None)
            group2, sub_iter2 = _advance(grouped2) or (None, None)


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
