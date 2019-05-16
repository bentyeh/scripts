import itertools
import numpy as np
import pandas as pd

def isint(s):
    '''
    Check if a string represents an int.
    Source: https://stackoverflow.com/a/1265696
    '''
    if s[0] in ('-', '+'):
        return s[1:].isdigit()
    return s.isdigit()

def isfloat(s):
    '''
    Check if a string represents a float.

    Notes
    - If the string represents an int, this function will still return True.
    - To determine whether a string `s` represents an int or a float, consider the following options:
      1. Use `isint(s)` first, then `isfloat(s)` if the former returns False.
      2. Use `isfloat()` first, then `s.is_integer()` if the former returns True.

    Source: https://stackoverflow.com/a/15357477
    '''
    try:
        float(s)
    except ValueError:
        return False
    return True

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
