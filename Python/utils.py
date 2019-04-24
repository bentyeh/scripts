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

def findOverlapIntervals(intervals, n):
    '''
    Find regions overlapped by n or more intervals.

    Args
    - intervals: list of 2-tuples of real numbers (float or int)
        Pythonic intervals [start, end), where start < end.
    - n: int
        The output regions will be overlapped by at least `n` of the intervals given in `Intervals`.

    Returns: list of 2-tuples of real numbers (float or int)

    Algorithm based on checking for balanced parentheses. Credits to Zachary Taylor for the idea.
    '''

    values = np.array(list(itertools.chain.from_iterable(intervals)))
    starts = np.array([1, -1]*len(intervals)) # 1: start value; -1: end value
    idx = np.argsort(values)
    values = values[idx]
    starts = starts[idx]
    counts = np.cumsum(starts)
    highCountsIdx = np.flatnonzero(counts >= n)

    regions = []
    start = None
    for i in range(len(highCountsIdx)):
        if start is None:
            start = highCountsIdx[i]

        if (i == len(highCountsIdx) - 1) or (highCountsIdx[i+1] - highCountsIdx[i] > 1):
            regions.append((values[start], values[highCountsIdx[i]+1]))
            start = None
    return regions