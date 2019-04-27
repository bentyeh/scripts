import itertools
import numpy as np

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

def weightedOverlapIntervals(intervals):
    '''
    Return non-overlapping set of regions whose weights are the sum of the weights of
    overlapping input intervals spanning that region

    Args
    - intervals: list of 2- or 3-tuples of real numbers
        (start, end[, weight]) where start < end
        weight assumed to be 1 if not provided

    Returns: 2-tuple of (list of 3-tuples of real numbers, list of sets of integers)
    - Element 0: list of non-overlapping intervals (start, end, weight) where start < end
      - These non-overlapping intervals span the entire range of the input intervals.
        Regions in this range not overlapped by at least 1 input interval are included
        with weight 0.
    - Element 1: list of sets of indices of original intervals contributing to new intervals
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

    start = s_sorted[0]      # current start
    weight = w[s_argsort[0]] # current weight
    s_idx = 1
    e_idx = 0

    newIntervals = []     # list of non-overlapping intervals
    overlapIntervals = [] # list of sets of indices of original intervals contributing to new intervals
    curOverlap = set([s_argsort[0]])

    while e_idx < n:
        if (s_idx < n) and (s_sorted[s_idx] < e_sorted[e_idx]):
            end = s_sorted[s_idx]
            newIntervals.append((start, end, weight))
            overlapIntervals.append(set(curOverlap))
            weight += w[s_argsort[s_idx]]
            curOverlap.add(s_argsort[s_idx])
            s_idx += 1
        else:
            end = e_sorted[e_idx]
            overlapIntervals.append(set(curOverlap))
            newIntervals.append((start, end, weight))
            weight -= w[e_argsort[e_idx]]
            curOverlap.remove(e_argsort[e_idx])
            e_idx += 1
        start = end
    
    # post-processing: remove length-0 intervals where start == end
    validIntervals = [i for i in range(len(newIntervals)) if newIntervals[i][0] < newIntervals[i][1]]
    newIntervals = [newIntervals[i] for i in validIntervals]
    overlapIntervals = [overlapIntervals[i] for i in validIntervals]
    
    return (newIntervals, overlapIntervals)