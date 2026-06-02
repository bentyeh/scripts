#!/usr/bin/env python3
"""
test_intersect.py

Unit tests for intersect.py, verifying alignment grouping, coordinate calculations,
gzip streaming detection, sweep-line sorting validation, and intersection modes.

Written by Gemini 3.5 Flash.
"""

import os
import io
import gzip
import sys
import tempfile
import unittest
from typing import Any
from unittest.mock import patch, MagicMock
import pysam

# Import components from the intersect module
sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)), '..'))
from intersect import (
    ConcatStream,
    EntryA,
    EntryB,
    open_bed_stream,
    get_a_iterator,
    process_intersections,
    open_output
)


class MockAlignedSegment:
    """
    A lightweight mock representing a pysam.AlignedSegment.
    Enables isolated unit testing of EntryA and intersection logic without full BAM files.
    """
    def __init__(self, chrom: str, start: int, end: int, name: str, is_reverse: bool = False):
        """
        Initializes the MockAlignedSegment.

        Args:
            chrom: The chromosome name.
            start: The 0-based start coordinate.
            end: The 0-based end coordinate.
            name: The query/read name.
            is_reverse: True if aligned to reverse strand, False otherwise.
        """
        self.reference_name = chrom
        self.reference_start = start
        self.reference_end = end
        self.query_name = name
        self.is_reverse = is_reverse
        self.tags: dict[str, Any] = {}

    def set_tag(self, tag: str, val: Any) -> None:
        """
        Simulates setting a SAM attribute tag.

        Args:
            tag: Name of the SAM tag.
            val: Value of the tag.
        """
        self.tags[tag] = val


class TestConcatStream(unittest.TestCase):
    """
    Tests for the ConcatStream class, ensuring transparent prepend buffering.
    """
    def test_concat_stream_reading(self) -> None:
        """
        Tests that magic bytes and normal payload are correctly stitched together.
        """
        prefix = b"\x1f\x8b"
        underlying = io.BytesIO(b"compressed_data_here")
        stream = ConcatStream(prefix, underlying)
        
        self.assertTrue(stream.readable())
        # Partial read containing prefix and part of payload
        self.assertEqual(stream.read(4), b"\x1f\x8bco")
        # Remaining read of payload
        self.assertEqual(stream.read(8), b"mpressed")
        # Read to EOF
        self.assertEqual(stream.read(), b"_data_here")


class TestEntryClasses(unittest.TestCase):
    """
    Tests for EntryA and EntryB structures, focusing on coordinate calculations and outputs.
    """
    def test_entry_a_single_end(self) -> None:
        """
        Validates coordinate and strand generation for single-end alignments.
        """
        read = MockAlignedSegment("chr1", 100, 150, "read_se", is_reverse=True)
        entry = EntryA([read])
        self.assertEqual(entry.chrom, "chr1")
        self.assertEqual(entry.start, 100)
        self.assertEqual(entry.end, 150)
        self.assertEqual(entry.strand, "-")

    def test_entry_a_paired_end(self) -> None:
        """
        Validates span coordinate and strand generation for paired-end alignments.
        """
        read1 = MockAlignedSegment("chr1", 100, 150, "read_pe")
        read2 = MockAlignedSegment("chr1", 300, 350, "read_pe", is_reverse=True)
        entry = EntryA([read1, read2])
        self.assertEqual(entry.chrom, "chr1")
        # Span should represent the entire template fragment (100 to 350)
        self.assertEqual(entry.start, 100)
        self.assertEqual(entry.end, 350)
        self.assertEqual(entry.strand, ".")

    def test_entry_b_parsing(self) -> None:
        """
        Verifies correct BED interval conversion of standard and extended BED records.
        """
        entry = EntryB(["chr1", "100", "200", "featX", "1000", "+"])
        self.assertEqual(entry.chrom, "chr1")
        self.assertEqual(entry.start, 100)
        self.assertEqual(entry.end, 200)


class TestBedStreaming(unittest.TestCase):
    """
    Tests the file reading and automated compression type detection mechanics of open_bed_stream.
    """
    def setUp(self) -> None:
        self.temp_files: list[str] = []

    def tearDown(self) -> None:
        for path in self.temp_files:
            if os.path.exists(path):
                os.remove(path)

    def test_open_bed_stream_plain(self) -> None:
        """
        Tests opening and processing plain-text BED files.
        """
        with tempfile.NamedTemporaryFile(mode='w', delete=False) as f:
            f.write("chr1\t100\t200\n")
            f.write("# comment block\n")
            f.write("chr1\t300\t400\tfeatA\n")
            self.temp_files.append(f.name)
            path = f.name

        entries = list(open_bed_stream(path))
        self.assertEqual(len(entries), 2)
        self.assertEqual(entries[0].start, 100)
        self.assertEqual(entries[1].fields[3], "featA")

    def test_open_bed_stream_gzipped(self) -> None:
        """
        Tests auto-detecting and reading gzip-compressed BED files.
        """
        content = b"chr1\t500\t600\nchr2\t10\t50\n"
        compressed = gzip.compress(content)
        
        with tempfile.NamedTemporaryFile(mode='wb', delete=False) as f:
            f.write(compressed)
            self.temp_files.append(f.name)
            path = f.name

        entries = list(open_bed_stream(path))
        self.assertEqual(len(entries), 2)
        self.assertEqual(entries[0].chrom, "chr1")
        self.assertEqual(entries[1].chrom, "chr2")


class TestGetAIterator(unittest.TestCase):
    """
    Tests iterator wrapping of SAM/BAM alignments, mapping coordinates and sorting.
    """
    def setUp(self) -> None:
        self.temp_files: list[str] = []

    def tearDown(self) -> None:
        for path in self.temp_files:
            if os.path.exists(path):
                os.remove(path)

    def create_sam_file(self, sam_lines: list[str]) -> str:
        """
        Helper to write a custom mock SAM file with coordinate headers.
        """
        header = (
            "@HD\tVN:1.6\tSO:coordinate\n"
            "@SQ\tSN:chr1\tLN:10000\n"
            "@SQ\tSN:chr2\tLN:10000\n"
        )
        with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix=".sam") as f:
            f.write(header)
            for line in sam_lines:
                f.write(line + "\n")
            self.temp_files.append(f.name)
            return f.name

    def test_single_end_alignments(self) -> None:
        """
        Ensures coordinate parsing handles 1-based SAM to 0-based Python conversions.
        """
        lines = [
            "read1\t0\tchr1\t101\t60\t50M\t*\t0\t0\t*\t*",  # 0-based: 100-150
            "read2\t16\tchr1\t201\t60\t50M\t*\t0\t0\t*\t*"  # 0-based: 200-250
        ]
        sam_path = self.create_sam_file(lines)
        with pysam.AlignmentFile(sam_path, "r") as f:
            entries = list(get_a_iterator(f))
            self.assertEqual(len(entries), 2)
            self.assertEqual(entries[0].start, 100)
            self.assertEqual(entries[0].strand, "+")
            self.assertEqual(entries[1].start, 200)
            self.assertEqual(entries[1].strand, "-")

    def test_paired_end_alignments(self) -> None:
        """
        Verifies read paired grouping and template span resolution.
        """
        lines = [
            "pair1\t99\tchr1\t101\t60\t50M\t=\t301\t250\t*\t*",  # Read 1 Forward, start=100
            "pair1\t147\tchr1\t301\t60\t50M\t=\t101\t-250\t*\t*" # Read 2 Reverse, start=300, end=350
        ]
        sam_path = self.create_sam_file(lines)
        with pysam.AlignmentFile(sam_path, "r") as f:
            entries = list(get_a_iterator(f))
            self.assertEqual(len(entries), 1)
            self.assertEqual(entries[0].start, 100)
            self.assertEqual(entries[0].end, 350)
            self.assertEqual(entries[0].strand, ".")

    def test_sorting_assertion(self) -> None:
        """
        Verifies that out-of-order coordinate records cause an immediate ValueError.
        """
        lines = [
            "read1\t0\tchr1\t501\t60\t50M\t*\t0\t0\t*\t*",  # 500-550
            "read2\t0\tchr1\t201\t60\t50M\t*\t0\t0\t*\t*"   # 200-250 (out of order!)
        ]
        sam_path = self.create_sam_file(lines)
        with pysam.AlignmentFile(sam_path, "r") as f:
            with self.assertRaises(ValueError):
                list(get_a_iterator(f))


class TestProcessIntersections(unittest.TestCase):
    """
    Comprehensive suite evaluating the core sweep-line intersection configurations.
    """
    def setUp(self) -> None:
        self.ref_dict = {"chr1": 0, "chr2": 1}

    def test_default_intersection(self) -> None:
        """
        Tests the default action, printing overlapping BED3 intervals.
        """
        a = EntryA([MockAlignedSegment("chr1", 100, 200, "rA")])
        b = EntryB(["chr1", "150", "250"])
        out = io.StringIO()
        
        process_intersections(iter([a]), iter([b]), self.ref_dict, None, out, None)
        self.assertEqual(out.getvalue(), "chr1\t150\t200\n")

    def test_fractional_overlaps_f_and_F(self) -> None:
        """
        Validates overlap fraction rules for both input A and input B.
        """
        # A length = 100. Overlaps B by 50bp (fraction of A = 0.5)
        # B length = 100. Overlaps A by 50bp (fraction of B = 0.5)
        a = EntryA([MockAlignedSegment("chr1", 100, 200, "rA")])
        b = EntryB(["chr1", "150", "250"])

        # 1. 50% overlap should clear default thresholds (1e-9)
        out = io.StringIO()
        process_intersections(iter([a]), iter([b]), self.ref_dict, None, out, None)
        self.assertEqual(out.getvalue(), "chr1\t150\t200\n")

        # 2. Minimum fraction of A = 60%. Overlap (50%) should fail.
        out = io.StringIO()
        process_intersections(iter([a]), iter([b]), self.ref_dict, None, out, None, f=0.6)
        self.assertEqual(out.getvalue(), "")

        # 3. Minimum fraction of B = 60%. Overlap (50%) should fail.
        out = io.StringIO()
        process_intersections(iter([a]), iter([b]), self.ref_dict, None, out, None, F=0.6)
        self.assertEqual(out.getvalue(), "")

    def test_mode_wa(self) -> None:
        """
        Validates -wa option (write original A coordinate).
        """
        a = EntryA([MockAlignedSegment("chr1", 100, 200, "rA")])
        b = EntryB(["chr1", "150", "250"])
        out = io.StringIO()
        
        process_intersections(iter([a]), iter([b]), self.ref_dict, "wa", out, "BED3")
        self.assertEqual(out.getvalue(), "chr1\t100\t200\n")

    def test_mode_wb_exact_once(self) -> None:
        """
        Ensures -wb dumps matched B features, preventing duplicates even if hit by multiple As.
        """
        a1 = EntryA([MockAlignedSegment("chr1", 100, 150, "rA1")])
        a2 = EntryA([MockAlignedSegment("chr1", 170, 220, "rA2")])
        b = EntryB(["chr1", "120", "200", "featB"])
        out = io.StringIO()
        
        process_intersections(iter([a1, a2]), iter([b]), self.ref_dict, "wb", out, None)
        # Should only report featB once, despite overlapping with both a1 and a2
        self.assertEqual(out.getvalue(), "chr1\t120\t200\tfeatB\n")

    def test_mode_wab(self) -> None:
        """
        Validates -wab option, yielding standard 6-column BEDPE output.
        """
        a = EntryA([MockAlignedSegment("chr1", 100, 200, "rA")])
        b = EntryB(["chr1", "150", "250"])
        out = io.StringIO()
        
        process_intersections(iter([a]), iter([b]), self.ref_dict, "wab", out, None)
        self.assertEqual(out.getvalue(), "chr1\t100\t200\tchr1\t150\t250\n")

    def test_mode_v(self) -> None:
        """
        Validates -v option, only displaying segments in A lacking any overlap in B.
        """
        a1 = EntryA([MockAlignedSegment("chr1", 100, 200, "rA1")])  # overlaps
        a2 = EntryA([MockAlignedSegment("chr1", 300, 400, "rA2")])  # does not overlap
        b = EntryB(["chr1", "150", "250"])
        out = io.StringIO()
        
        process_intersections(iter([a1, a2]), iter([b]), self.ref_dict, "v", out, "BED3")
        self.assertEqual(out.getvalue(), "chr1\t300\t400\n")

    def test_mode_c(self) -> None:
        """
        Validates -c option, generating columns showing hit tallies.
        """
        a = EntryA([MockAlignedSegment("chr1", 100, 300, "rA")])
        b1 = EntryB(["chr1", "120", "150"])
        b2 = EntryB(["chr1", "200", "250"])
        out = io.StringIO()
        
        process_intersections(iter([a]), iter([b1, b2]), self.ref_dict, "c", out, "BED3")
        self.assertEqual(out.getvalue(), "chr1\t100\t300\t2\n")

    def test_chrom_transitions(self) -> None:
        """
        Confirms streaming functions correctly across multiple chromosome segments.
        """
        a1 = EntryA([MockAlignedSegment("chr1", 100, 200, "rA1")])
        a2 = EntryA([MockAlignedSegment("chr2", 100, 200, "rA2")])
        
        b1 = EntryB(["chr1", "150", "250"])
        b2 = EntryB(["chr2", "150", "250"])
        out = io.StringIO()
        
        process_intersections(iter([a1, a2]), iter([b1, b2]), self.ref_dict, None, out, None)
        self.assertEqual(out.getvalue(), "chr1\t150\t200\nchr2\t150\t200\n")


class TestOpenOutput(unittest.TestCase):
    """
    Tests for open_output context manager using mocks.
    """
    @patch('pysam.AlignmentFile')
    def test_open_output_sam(self, mock_alignment_file: MagicMock) -> None:
        """
        Ensures open_output correctly binds system stdout for SAM format targets.
        """
        template = MagicMock()
        with open_output('SAM', template) as _:
            pass
        mock_alignment_file.assert_called_once_with(sys.stdout, 'w', template=template)

    def test_open_output_bed(self) -> None:
        """
        Ensures BED3/BED6 formats default straightforwardly to sys.stdout.
        """
        with open_output('BED3', None) as f:
            self.assertIs(f, sys.stdout)


if __name__ == '__main__':
    unittest.main()