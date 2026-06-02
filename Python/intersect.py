#!/usr/bin/env python3
"""
intersect.py

Determines whether 2 sets of genomic features intersect.
It accepts a template-coordinate-sorted SAM/BAM file (-a) and a position-sorted 
BED file (-b), returning intervals or original entries that overlap based on constraints.

It processes interval comparison in a streaming manner using a sweep-line algorithm.
This module can be executed via CLI or imported into other Python projects.

Generated with help from Gemini 3.1 Pro.
"""

import argparse
import contextlib
import gzip
import sys
import io
from typing import Any, Literal, BinaryIO
from collections.abc import Iterator
import pysam


class ConcatStream(io.BufferedIOBase):
    """
    A readable stream that concatenates a byte prefix with an underlying stream.
    Used to prepend peeked magic bytes back onto a stream without consuming them.
    """
    def __init__(self, prefix: bytes, stream: BinaryIO):
        """
        Initializes the ConcatStream.

        Args:
            prefix: The byte prefix to prepend to the stream.
            stream: The underlying buffered IO stream.
        """
        self.prefix = prefix
        self.stream = stream
        self.pos = 0

    def read(self, size: int = -1) -> bytes:
        """
        Reads from the concatenated stream.

        Args:
            size: The number of bytes to read. Defaults to -1 (read all).

        Returns:
            The bytes read from the stream.
        """
        if self.pos < len(self.prefix):
            if size == -1 or size >= len(self.prefix) - self.pos:
                res = self.prefix[self.pos:]
                self.pos = len(self.prefix)
                
                if size != -1:
                    extra = self.stream.read(size - len(res))
                    if extra: res += extra
                else:
                    extra = self.stream.read()
                    if extra: res += extra
                return res
            else:
                res = self.prefix[self.pos : self.pos + size]
                self.pos += size
                return res
        return self.stream.read(size)

    def read1(self, size: int = -1) -> bytes:
        """
        Reads and returns up to size bytes, delegating to read().
        This is required by io.TextIOWrapper under certain conditions.

        Args:
            size: The maximum number of bytes to read. Defaults to -1.

        Returns:
            The bytes read.
        """
        return self.read(size)

    def readable(self) -> bool:
        """
        Returns whether the stream is readable.

        Returns:
            True, as the stream is always readable.
        """
        return True


class EntryA:
    """
    Represents an entry from File A (either a single read or a paired fragment).
    """
    def __init__(self, reads: list[pysam.AlignedSegment]):
        """
        Initializes an EntryA object.

        Args:
            reads: A list of pysam.AlignedSegment objects representing the entry.
        """
        self.reads = reads
        self.chrom = reads[0].reference_name
        # Template start/end for the full span (0-based, half-open)
        self.start = min(r.reference_start for r in reads)
        self.end = max(r.reference_end for r in reads)
        self.name = reads[0].query_name
        
        if len(reads) == 1:
            self.strand = '-' if reads[0].is_reverse else '+'
        else:
            self.strand = '.'

    def write(self, out_file: Any, output_type: str | None, hits: int | None = None, tag: str = 'XC') -> None:
        """
        Writes the original entry in A to the specified output.

        Args:
            out_file: The output file-like object to write to.
            output_type: The format to output ('SAM', 'BAM', 'UBAM', 'BED3', 'BED6', or None).
            hits: The number of intersection hits, optionally reported.
            tag: The SAM tag used for reporting hits (default: 'XC').

        Returns:
            None
        """
        if output_type in ('SAM', 'BAM', 'UBAM'):
            for r in self.reads:
                if hits is not None:
                    r.set_tag(tag, hits)
                out_file.write(r)
        elif output_type in ('BED3', 'BED6'):
            row = [self.chrom, str(self.start), str(self.end)]
            if output_type == 'BED6' or hits is not None:
                # Need to extend array if BED6 requested
                if output_type == 'BED6':
                    row.extend([self.name, '0', self.strand])
            if hits is not None:
                row.append(str(hits))
            out_file.write("\t".join(row) + "\n")


class EntryB:
    """
    Represents an entry from File B (a BED row).
    """
    def __init__(self, fields: list[str]):
        """
        Initializes an EntryB object.

        Args:
            fields: A list of strings representing the tab-separated fields of a BED row.
        """
        self.fields = fields
        self.chrom = fields[0]
        self.start = int(fields[1])
        self.end = int(fields[2])
        self.written = False

    def write(self, out_file: Any, output_type: str | None = None) -> None:
        """
        Writes the original entry in B to the specified output.

        Args:
            out_file: The output file-like object to write to.
            output_type: The format to output ('BED3', 'BED6', or None).

        Returns:
            None
        """
        if output_type in ('BED3', 'BED6'):
            row = [self.chrom, str(self.start), str(self.end)]
            if output_type == 'BED6':
                name = self.fields[3] if len(self.fields) > 3 else '.'
                score = self.fields[4] if len(self.fields) > 4 else '0'
                strand = self.fields[5] if len(self.fields) > 5 else '.'
                row.extend([name, score, strand])
            out_file.write("\t".join(row) + "\n")
        else:
            out_file.write("\t".join(self.fields) + "\n")


@contextlib.contextmanager
def open_default(file: str | Any, opener=open, default=None, **kwargs) -> Iterator[Any]:
    """
    Context manager to open a file if given as a string path, or yield the file object if already opened.
    """
    if isinstance(file, str):
        with opener(file, **kwargs) as f:
            yield f
    else:
        yield default


@contextlib.contextmanager
def open_gzip_auto(concat_stream, magic) -> Iterator[Any]:
    """
    Context manager to open a gzip stream if magic bytes indicate gzip compression.

    Args:
        concat_stream: The ConcatStream object to read from.
        magic: The first two bytes read from the stream to check for gzip magic numbers.
    
    Yields:
        A file-like object for reading the BED data.
    """
    if magic == b'\x1f\x8b':
        with gzip.open(concat_stream, 'rt') as f:
            yield f
    else:
        with io.TextIOWrapper(concat_stream) as f:
            yield f


def open_bed_stream(filepath: str) -> Iterator[EntryB]:
    """
    Opens a BED file, detecting gzip compression via magic bytes.
    
    Args:
        filepath: Path to the BED file, or '-' to read from standard input.
        
    Yields:
        Parsed EntryB objects.
    """
    with open_default(filepath if filepath != '-' else None, mode='rb', default=sys.stdin.buffer) as stream:
        magic = stream.read(2)
        concat_stream = ConcatStream(magic, stream)

        with open_gzip_auto(concat_stream, magic) as f:    
            for line in f:
                line = line.strip()
                if line and not line.startswith('#'):
                    yield EntryB(line.split('\t'))


def get_a_iterator(alignment_file: pysam.AlignmentFile) -> Iterator[EntryA]:
    """
    Iterates over a coordinate-sorted SAM/BAM file and yields EntryA objects 
    grouped by fragments. Invalid pairs are skipped.
    
    Validates that the yielded entries are strictly template-coordinate sorted.

    Args:
        alignment_file: The opened SAM/BAM alignment file object.

    Yields:
        Parsed EntryA objects representing single reads or fragments.
    """
    unmapped_count = 0
    ignored_pairs_count = 0
    waiting_mates: dict[str, pysam.AlignedSegment] = {}
    
    last_chrom = None
    last_start = -1

    for read in alignment_file.fetch(until_eof=True):
        if read.is_unmapped:
            unmapped_count += 1
            continue

        if not read.is_paired:
            entry = EntryA([read])
        else:
            # Pair validations: unmapped mate, diff chromosomes, non-FR orientation
            if read.mate_is_unmapped or read.reference_id != read.next_reference_id or read.is_reverse == read.mate_is_reverse:
                if read.is_read1 or not read.is_read2: # roughly avoid double-counting pairs
                    ignored_pairs_count += 1
                continue

            qname = read.query_name
            if qname in waiting_mates:
                mate = waiting_mates.pop(qname)
                entry = EntryA([mate, read])
            else:
                waiting_mates[qname] = read
                continue

        # Validate sorting dynamically
        if entry.chrom != last_chrom:
            last_chrom = entry.chrom
            last_start = -1
            
        if entry.start < last_start:
            raise ValueError(
                f"File A is not sorted by template start! "
                f"Entry starting at {entry.start} found after {last_start}."
            )
        last_start = entry.start
        yield entry

    print(f"Intersect.py Warning: Ignored {ignored_pairs_count} incorrectly oriented/mapped read pairs "
          f"and {unmapped_count} unmapped single reads.", file=sys.stderr)


def process_intersections(
    iter_A: Iterator[EntryA],
    iter_B: Iterator[EntryB],
    ref_dict: dict[str, int],
    mode: Literal["wa", "wb", "wab", "c", "v"] | None,
    out_file: Any,
    output_type: Literal['SAM', 'BAM', 'UBAM', 'BED3', 'BED6'] | None,
    f: float = 1e-9,
    F: float = 1e-9,
    tag: str = 'XC'
) -> None:
    """
    Performs streaming interval intersections using a sweep-line algorithm.
    Only holds active, overlapping B intervals in memory.

    Time Complexity: O(N_A + N_B) 
    Memory Complexity: O(K), where K is the max number of overlapping B intervals at a time.

    Args:
        iter_A: Iterator of EntryA objects.
        iter_B: Iterator of EntryB objects.
        ref_dict: Dictionary mapping reference sequence names to their integer indices.
        mode: The intersection reporting mode ('wa', 'wb', 'wab', 'c', 'v', or None).
        out_file: The output file-like object to write to.
        output_type: The format to output.
        f: Minimum overlap required as a fraction of A (default: 1e-9).
        F: Minimum overlap required as a fraction of B (default: 1e-9).
        tag: The SAM tag used for reporting hits (default: 'XC').

    Returns:
        None
    """
    try:
        current_b: EntryB | None = next(iter_B)
    except StopIteration:
        current_b = None

    active_B: list[EntryB] = []

    for a in iter_A:
        a_chrom_idx = ref_dict.get(a.chrom, -1)
        if a_chrom_idx == -1:
            continue

        # 1. Prune obsolete B entries out of the active buffer
        new_active = []
        for b in active_B:
            b_chrom_idx = ref_dict.get(b.chrom, -1)
            # Retain if it's on the same chromosome and ends after A starts
            if b_chrom_idx == a_chrom_idx and b.end > a.start:
                new_active.append(b)
        active_B = new_active

        # 2. Advance the B stream to catch up with A
        while current_b is not None:
            b_chrom_idx = ref_dict.get(current_b.chrom, float('inf'))
            if b_chrom_idx == float('inf'):
                current_b = next(iter_B, None)
                continue

            if b_chrom_idx < a_chrom_idx:
                current_b = next(iter_B, None)
            elif b_chrom_idx == a_chrom_idx:
                if current_b.end <= a.start:
                    current_b = next(iter_B, None)
                elif current_b.start < a.end:
                    active_B.append(current_b)
                    current_b = next(iter_B, None)
                else:
                    break  # B is ahead of A, hold it
            else:
                break  # B's chromosome is ahead of A's, hold it

        # 3. Find specific overlaps within active_B
        hits = 0
        for b in active_B:
            overlap_start = max(a.start, b.start)
            overlap_end = min(a.end, b.end)
            
            if overlap_start < overlap_end:
                overlap_len = overlap_end - overlap_start
                frac_a = overlap_len / (a.end - a.start) if a.end > a.start else 0
                frac_b = overlap_len / (b.end - b.start) if b.end > b.start else 0
                
                if frac_a >= f and frac_b >= F:
                    hits += 1
                    
                    if mode == 'wab':
                        out_file.write(f"{a.chrom}\t{a.start}\t{a.end}\t{b.chrom}\t{b.start}\t{b.end}\n")
                    elif mode == 'wb' and not b.written:
                        b.write(out_file, output_type)
                        b.written = True
                    elif mode is None:
                        # Default intersection behavior (BED3)
                        out_file.write(f"{a.chrom}\t{overlap_start}\t{overlap_end}\n")

        # 4. Handle modifiers attached to entry A
        if mode == 'wa' and hits > 0:
            a.write(out_file, output_type)
        elif mode == 'v' and hits == 0:
            a.write(out_file, output_type)
        elif mode == 'c':
            a.write(out_file, output_type, hits=hits, tag=tag)


def parse_args() -> argparse.Namespace:
    """
    Parses command-line arguments.

    Returns:
        Parsed arguments.
    """
    parser = argparse.ArgumentParser(description="Streaming intersection of genomic features.")
    parser.add_argument("-a", required=True, help="Path to template-coordinate-sorted File A (SAM/BAM)")
    parser.add_argument("-b", required=True, help="Path to position-sorted File B (BED/GZ-BED)")
    
    group = parser.add_mutually_exclusive_group()
    group.add_argument("-wa", action="store_true", help="Write original A entry on overlap")
    group.add_argument("-wb", action="store_true", help="Write original B entry on overlap")
    group.add_argument("-wab", action="store_true", help="Write both A and B entries (BEDPE)")
    group.add_argument("-v", action="store_true", help="Only report entries in A with no overlaps in B")
    group.add_argument("-c", action="store_true", help="Count hits in B for each A")
    
    parser.add_argument("--output_type", choices=["SAM", "BAM", "UBAM", "BED3", "BED6"], 
                        help="Format of the output.")
    parser.add_argument("-f", type=float, default=1e-9, help="Minimum overlap fraction of A (default 1e-9)")
    parser.add_argument("-F", type=float, default=1e-9, help="Minimum overlap fraction of B (default 1e-9)")
    parser.add_argument("--tag", default="XC", help="SAM tag for counts (-c). Defaults to XC.")
    
    args = parser.parse_args()
    
    if args.a == '-' and args.b == '-':
        parser.error("Either -a or -b (but not both) can be passed the value '-' to read from stdin")
    
    if args.wab and args.output_type is not None:
        parser.error("-wab is mutually exclusive with --output_type")
        
    if args.wb and args.output_type not in (None, "BED3", "BED6"):
        parser.error("-wb is only compatible with BED3 or BED6 --output_type")
        
    # Defaults
    if args.output_type is None:
        if args.wa or args.v or args.c:
            args.output_type = "SAM"

    return args


@contextlib.contextmanager
def open_output(output_type: str | None, template: pysam.AlignmentFile | None) -> Iterator[Any]:
    """
    Context manager to open a BAM file if given as a string path, or yield the file object if already opened.

    Args:
        output_type: The desired output format, or None.
        template: The template alignment file for headers if outputting BAM/SAM.

    Yields: 
        Opened output file object (BAM file or sys.stdout).
    """
    if output_type in ('SAM', 'BAM', 'UBAM'):
        mode_map: dict[str, Literal['w', 'wb', 'wbu']] = {'SAM': 'w', 'BAM': 'wb', 'UBAM': 'wbu'}
        with pysam.AlignmentFile(
            sys.stdout.buffer if output_type in ('BAM', 'UBAM') else sys.stdout,
            mode_map[output_type],
            template=template
        ) as f:
            yield f
    else:
        yield sys.stdout


def main() -> None:
    """
    Main execution entry point for the script.
    """
    args = parse_args()

    if args.wa:
        mode = 'wa'
    elif args.wb:
        mode = 'wb'
    elif args.wab:
        mode = 'wab'
    elif args.c:
        mode = 'c'
    elif args.v:
        mode = 'v'
    else:
        mode = None

    # Determine input A configuration
    a_path = sys.stdin if args.a == '-' else args.a
    with pysam.AlignmentFile(a_path) as alignment_file:
        with open_output(args.output_type, template=alignment_file) as out_file:    
            # Track Canonical Order using A's Sequence Header Dictionary
            ref_dict = {ref: i for i, ref in enumerate(alignment_file.references)}

            iter_A = get_a_iterator(alignment_file)
            iter_B = open_bed_stream(args.b)
            process_intersections(
                iter_A,
                iter_B,
                ref_dict,
                mode,
                out_file,
                args.output_type,
                f=args.f,
                F=args.F,
                tag=args.tag
            )


if __name__ == "__main__":
    main()