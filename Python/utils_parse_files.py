# add path of this file (e.g., a scripts directory) to sys.path
import collections, multiprocessing, os, re, sys, time
from pathlib import Path

import pandas as pd

sys.path.append(str(Path(__file__).resolve(strict=True).parent))
import utils, utils_files

def _to_csv(x, col, folder_out, ext='', **kwargs):
    '''
    Save DataFrame to file. This function is necessary for shard_df_pandas(),
    since local lambda functions cannot be pickled and used in a multiprocessing
    context.
    
    Args
    - x: pandas.DataFrame or pandas.Series
        Sharded table
    - col: str
        Name of column containing a single key on which `x` was sharded
    - folder_out
        Directory where to save DataFrame. Filename is {key}{ext} where {key} is the
        unique key on which `x` was sharded, and {ext} is given by the `ext` argument
    - ext: str. default=''
        Extension of output files, including the period. Compression is inferred.
    - **kwargs (sep, header, index, mode, ...):
        See pandas.DataFrame.to_csv().
    
    Returns: None
    '''
    x.to_csv(Path(folder_out, x[col].unique()[0] + ext), **kwargs)

def shard_df_pandas(file_in, folder_out, col, sep, header=None, chunksize=None,
                    nproc=None, ext=None, verbose=True, **kwargs):
    '''
    Shard a table by keys in a given column.

    Args
    - file_in: str
        Path to file to shard
    - folder_out: str
        Directory to save sharded file.
        Filename of output is {key}{ext} - see the `ext` argument.
    - col: str
        Name of column containing keys on which to shard
    - sep: str. default=None
        Delimiter. If None, splits on consecutive whitespace.
    - header: int. default=None.
        Row number to use as the column names.
    - chunksize: int. default=None
        Number of rows to read for each chunk. If None, reads the entire `file_in` into memory.
        Memory usage will be primarily determined by this argument.
    - nproc: int. default=None
        Number of processes to use. If None, defaults to the number of usable CPUs.
    - ext: str. default=None
        Extension of output files, including the period. Compression is inferred.
        If None, automatically detects extension based on `file_in`.
    - verbose: bool. default=True
        Print out progress
    - **kwargs:
        Additional arguments to pass to pandas.read_csv() when reading `file_in`

    Returns: None

    Notes
    - Unlike shard_df(), this function writes out a single file per key, rather than multiple files
      per key that have to be concatenated later.
    - The groupby() operation and subsequent iteration over groups is performed by a single process
      and can therefore be a performance bottleneck. Therefore, a moderate chunksize (e.g., 1e6 or 1e7)
      is recommended.
    - Parallelization is achieved in saving grouped dataframes to disk. Each forked process
      occupies the same amount of memory as the main process. Therefore, memory usage can be substantial
      if nproc is large and/or chunksize is large.
    '''
    if ext is None:
        m = re.search(r'\.[^.\s]+(\.(gz|xz|zip|bz2))?$', file_in)
        ext = '' if m is None else m.group()
    
    if chunksize is not None and chunksize > os.path.getsize(file_in):
        chunksize = None

    reader = pd.read_csv(file_in, sep=sep, header=header, chunksize=chunksize, **kwargs)

    # special case if the entire file is read in
    if isinstance(reader, (pd.Series, pd.DataFrame)):
        chunksize = reader.shape[0]
        reader = [reader]

    if verbose:
        start_time = time.time()
        nChunks = 0

    for chunk in reader:
        grouped = chunk.groupby(col)
        utils.pandas_apply_parallel(
            grouped, _to_csv, nproc=nproc, concat_axis=None,
            col=col, folder_out=folder_out, ext=ext,
            sep=sep, header=(nChunks == 0) and (header is not None),
            index=False, mode='w' if nChunks == 0 else 'a')
        if verbose:
            nChunks += 1
            nLines = int(nChunks * chunksize)
            elapsed = time.time() - start_time
            print(f'Lines processed: {nLines}',
                  f'Time elapsed: {elapsed:.2g} s',
                  f'Average speed: {nLines/elapsed:.2g} lines/s',
                  sep=', ')

def shard_df(file_in, folder_out, col, start=0, end=None, sep=None, flush=None, ext=None, close=True, verbose=True):
    '''
    Shard a table by keys in a given column.

    Args
    - file_in: str
        Path to file to shard
    - folder_out: str
        Directory to save sharded file.
        Filename of output is {key}_{start}{ext}, where the extension matches that of file_in
    - col: int
        Column containing keys on which to shard
    - sep: str. default=None
        Delimiter. If None, splits on consecutive whitespace.
    - start: int. default=0
        Offset in bytes from the start of the file to start sharding
    - end: int. default=None
        Offset in bytes from the start of the file to end sharding
        If None, defaults to the end of the file
    - flush: int. default=None
        Number of lines processed after which to write flush and print out progress if verbose
    - ext: str. default=None
        Extension of output files, including the period. Compression is inferred.
        If None, automatically detects extension based on `file_in`.
    - close: bool. default=True
        Close all file objects pointing to sharded file outputs, even if there is an error.
    - verbose: bool. default=True
        Print out progress

    Note
    - This function shards from the first full line after the start offset (e.g., after a newline)
      through the line "intersected" by the end offset.
    - This function can be thought of as the "map" operation in a MapReduce framework.

    Returns: None
    '''
    # automatically detect extension
    if ext is None:
        m = re.search(r'\.[^.\s]+(\.(gz|xz|zip|bz2))?$', file_in)
        ext = '' if m is None else m.group()
    
    file_size = os.path.getsize(file_in)
    if end is None or end > file_size:
        end = file_size

    # validate arguments
    assert start >= 0 and start <= end
    assert flush is None or flush > 0
    assert close in (True, False)
    assert Path(folder_out).exists() and Path(folder_out).is_dir()

    # dictionary of keys -> file objects opened for writing
    files_out = {}

    with open(file_in) as f:
        # seek to next full line
        f.seek(start)
        if start != 0:
            f.readline()

        if flush is not None:
            nLines = 0
        if verbose:
            start_time = time.time()

        try:
            while f.tell() <= end:
                line = f.readline()
                key = line.split(sep)[col]
                files_out \
                    .setdefault(
                        key,
                        utils_files.createFileObject(Path(folder_out, f'{key}_{start}{ext}'), 'wt')) \
                    .write(line)

                if flush is not None:
                    nLines += 1
                    if nLines % flush == 0:
                        for g in files_out.values():
                            g.flush()
                        if verbose:
                            nBytes = f.tell() - start
                            print(
                                f'Offset: {start}',
                                f'Bytes processed: {nBytes}',
                                f'Bytes remaining: {end - start - nBytes}',
                                f'Lines processed: {nLines}',
                                f'Lines per second: {nLines/(time.time() - start_time):.3g}',
                                sep=', ')
        finally:
            if close:
                for g in files_out.values():
                    g.close()

def shard_df_mp(file_in, folder_out, col, nproc=None, shard_size=None, **kwargs):
    '''
    Shard a table by keys in a given column. Multiprocessing wrapper around shard_df().

    Args
    - file_in, folder_out, col, sep, flush, ext, verbose: see shard_df()
    - nproc: int. default=None
        Number of processes to use. If None, defaults to the number of usable CPUs.
    - shard_size: int. default=None
        Approximate size of each shard in bytes. If None, defaults to size of the file / nproc.
    - **kwargs (sep, flush, ext, close, verbose):
        Arguments to shard_df

    Returns: None
    '''
    if nproc is None:
        nproc = len(os.sched_getaffinity(0))
    file_size = os.path.getsize(file_in)
    if shard_size is None:
        shard_size = int(file_size / nproc)
    with open(file_in) as f:
        with multiprocessing.Pool(nproc) as pool:
            for shard_start in range(0, file_size, shard_size):
                shard_end = shard_start + shard_size
                pool.apply_async(shard_df, (file_in, folder_out, col, shard_start, shard_end), kwargs)
            pool.close()
            pool.join()

def etree_to_dict(t):
    '''
    Convert xml.etree.ElementTree.Element to dictionary.
    Source: https://stackoverflow.com/a/10077069
    '''
    d = {t.tag: {} if t.attrib else None}
    children = list(t)
    if children:
        dd = collections.defaultdict(list)
        for dc in map(etree_to_dict, children):
            for k, v in dc.items():
                dd[k].append(v)
        d = {t.tag: {k: v[0] if len(v) == 1 else v for k, v in dd.items()}}
    if t.attrib:
        d[t.tag].update(('@' + k, v) for k, v in t.attrib.items())
    if t.text:
        text = t.text.strip()
        if children or t.attrib:
            if text:
                d[t.tag]['#text'] = text
        else:
            d[t.tag] = text
    return d
