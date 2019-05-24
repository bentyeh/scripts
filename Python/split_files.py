# add path of this file (e.g., a scripts directory) to sys.path
import sys, os
sys.path.append(os.path.dirname(__file__))

import multiprocessing, re, time
import pandas as pd

import utils, utils_files

def _to_csv(x, col, sep, header, folder_out, ext='', mode='a', **kwargs):
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
    - sep, header, mode, **kwargs:
        See pandas.to_csv()
    
    Returns: None
    '''
    x.to_csv(os.path.join(folder_out, x[col].unique()[0] + ext),
             sep=sep, header=header, index=False, mode=mode, **kwargs)

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
      and can therefore by a performance bottleneck. Therefore, a moderate chunksize (e.g., 1e6 or 1e7)
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
            grouped, _to_csv, nproc=nproc, concat=False,
            col=col, sep=sep, header=(nChunks == 0) and (header is not None), folder_out=folder_out,
            ext=ext, mode='w' if nChunks == 0 else 'a'
        )
        if verbose:
            nChunks += 1
            nLines = int(nChunks * chunksize)
            elapsed = time.time() - start_time
            print(f'Lines processed: {nLines}, ' \
                  f'Time elapsed: {round(elapsed, 2)} s, '\
                  f'Average speed: {round(nLines/elapsed, 2)} lines/s')

def shard_df(file_in, folder_out, col, sep=None, start=0, end=None, flush=None, ext=None, verbose=True):
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
    assert os.path.exists(folder_out)

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

        while f.tell() <= end:
            line = f.readline()
            key = line.split(sep)[col]
            files_out.setdefault(key,
                                 utils_files.createFileObject(os.path.join(folder_out, f'{key}_{start}{ext}'), 'wt')) \
                     .write(line)

            if flush is not None:
                nLines += 1
                if nLines % flush == 0:
                    for g in files_out.values():
                        g.flush()
                    if verbose:
                        nBytes = f.tell() - start
                        print(f'Offset: {start}, Bytes processed: {nBytes}, Bytes remaining: {end - start - nBytes}, ' \
                              f'Lines processed: {nLines}, Lines per second: {nLines/(time.time() - start_time)}')

def shard_df_mp(file_in, folder_out, col, nproc=None, shard_size=None, sep=None, flush=None, ext=None, verbose=True):
    '''
    Shard a table by keys in a given column. Multiprocessing wrapper around shard_df().

    Args
    - file_in, folder_out, col, sep, flush, ext, verbose: see shard_df()
    - nproc: int. default=None
        Number of processes to use. If None, defaults to the number of usable CPUs.
    - shard_size: int. default=None
        Approximate size of each shard in bytes. If None, defaults to size of the file / nproc.

    Returns: None
    '''
    if nproc is None:
        nproc = len(os.sched_getaffinity(0))
    file_size = os.path.getsize(file_in)
    if shard_size is None:
        shard_size = int(file_size / nproc)
    with open(file_in) as f:
        with multiprocessing.Pool(nproc) as pool:
            for shard_start in range(start, file_size, shard_size):
                shard_end = shard_start + shard_size
                pool.apply_async(shard_df, (file_in, folder_out, col, sep, shard_start, shard_end, flush, verbose))
            pool.close()
            pool.join()
