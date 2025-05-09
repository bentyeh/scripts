from pathlib import Path
import sys

sys.path.append(str(Path(__file__).resolve(strict=True).parent))
import utils_misc

INTERVALS1 = [
    ('chr1', 100, 200),
    ('chr1', 300, 400),
    ('chr2', 100, 200),
    ('chr2', 300, 400)
]
INTERVALS2 = [
    ('chr1', 150, 250),
    ('chr1', 350, 450),
    ('chr2', 150, 250),
    ('chr2', 450, 550)
]

def test1():
    """Interval intersection"""
    expected_results = [
        ('chr1', 150, 200),
        ('chr1', 350, 400),
        ('chr2', 150, 200)
    ]
    results = []
    for (interval1, interval2, _, _) in utils_misc.intervals_overlap_grouped(INTERVALS1, INTERVALS2, how='inner'):
        chrom1, start1, end1 = interval1
        chrom2, start2, end2 = interval2
        results.append((chrom1, max(start1, start2), min(end1, end2)))
    assert results == expected_results, f"Expected {expected_results}, but got {results}"
    print("Test 1 passed!")

def test2():
    """Semi-join intervals"""
    expected_results = [
        ('chr1', 100, 200),
        ('chr1', 300, 400),
        ('chr2', 100, 200),
    ]
    results = []
    for (interval1, _) in utils_misc.intervals_filter_by_overlap(INTERVALS1, INTERVALS2, how='semi'):
        results.append(interval1)
    assert results == expected_results, f"Expected {expected_results}, but got {results}"
    print("Test 2 passed!")

def test3():
    """Anti-join intervals"""
    expected_results = [
        ('chr2', 300, 400)
    ]
    results = []
    for (interval1, _) in utils_misc.intervals_filter_by_overlap(INTERVALS1, INTERVALS2, how='anti'):
        results.append(interval1)
    assert results == expected_results, f"Expected {expected_results}, but got {results}"
    print("Test 3 passed!")

def test4():
    """Overlapping at distance"""
    expected_results = [
        (('chr1', 100, 200), ('chr1', 150, 250)),
        (('chr1', 300, 400), ('chr1', 350, 450)),
        (('chr2', 100, 200), ('chr2', 150, 250)),
        (('chr2', 300, 400), None)
    ]
    results = []
    for interval1, interval2, _, _ in utils_misc.intervals_overlap_grouped(
        INTERVALS1,
        INTERVALS2,
        how='left',
        overlap_func=lambda *args: utils_misc.overlap_length(*args, allow_negative=True) > -49,
        max_distance=49
    ):
        results.append((interval1, interval2))
    assert results == expected_results, f"Expected {expected_results}, but got {results}"

    expected_results = [
        (('chr1', 100, 200), ('chr1', 150, 250)),
        (('chr1', 300, 400), ('chr1', 350, 450)),
        (('chr2', 100, 200), ('chr2', 150, 250)),
        (('chr2', 300, 400), ('chr2', 450, 550))
    ]
    results = []
    for interval1, interval2, _, _ in utils_misc.intervals_overlap_grouped(
        INTERVALS1,
        INTERVALS2,
        how='left',
        overlap_func=lambda *args: utils_misc.overlap_length(*args, allow_negative=True) > -50,
        max_distance=50
    ):
        results.append((interval1, interval2))
    assert results == expected_results, f"Expected {expected_results}, but got {results}"
    print("Test 4 passed!")

def main():
    test1()
    test2()
    test3()
    print("All tests passed!")

if __name__ == '__main__':
    main()
