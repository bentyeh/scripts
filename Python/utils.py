import multiprocessing, os, queue, signal, subprocess, threading, time
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

def wrap_signal_handler(fun, sig, handler=signal.SIG_DFL):
    '''
    Wrap a function with a signal handler.

    Args
    - fun: function
        Function to to wrap
    - sig: signal
        Signal to handle with handler
    - handler: function or signal.SIG_IGN or signal.SIG_DFL. default=signal.SIG_DFL
        Signal handler

    Returns: function

    Examples
    - Wrap a function such that it ignores SIGINT:
        wrap_signal_handler(fun, signal.SIGINT, handler=signal.SIG_IGN)
    '''
    def newfun(*args, **kwargs):
        signal.signal(sig, handler)
        return fun(*args, **kwargs)
    return newfun

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
      df.groupby(list(set(df.columns) - set([col])))
        .agg(lambda x: sep.join(x))
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

def pandas_apply_parallel(grouped, func, nproc=None, concat=True, **kwargs):
    '''
    Parallel pandas.DataFrame.apply

    Args
    - grouped: pandas.Series or pandas.DataFrame
    - func: function
        Accepts a pandas.Series or pandas.DataFrame
    - nproc: int. default=None
        Number of processes to use. If None, defaults to the number of usable CPUs.
    - concat: bool or None. default=True
        Specifies the return type
          True: Concatenate results (outputs of `func`) with pandas.concat().
          False: List where each element is the result of applying `func` to a group.
          None: Returns None
    - **kwargs:
        Additional arguments to pass to func

    Returns: pandas.DataFrame, list, or None
      See the `concat` argument.

    Based on https://stackoverflow.com/a/29281494
    '''
    if nproc is None:
        nproc = len(os.sched_getaffinity(0))

    with multiprocessing.Pool(nproc) as pool:
        # ret_list = pool.map(func, [group for name, group in grouped])
        result_list = []
        for name, group in grouped:
            result_list.append(pool.apply_async(func, (group,), kwargs))
        pool.close()
        pool.join()

    if concat is None:
        return None

    results = [result.get() for result in result_list]
    if concat:
        return pd.concat(results)
    else:
        return results

class ThreadPool:
    '''
    Thread pool to drive subprocesses.

    Usage
      pool = ThreadPool(nproc)
      pool.schedule(name, cmd)
      pool.join()

    Example
      pool = ThreadPool(4)
      pool.schedule('example1', {'args': ['ls', '-l']})
      pool.join()

    Notes
    - This implementation relies on subprocess.Popen() to spawn subprocesses. Compare with the following:
      - subprocess.run()

    Compare with the standard multiprocessing library
    - SIGINT signals (e.g., from Ctrl-C) will terminate process directly created
      by multiprocessing.Process or multiprocessing.[dummy.]Pool, but not descendant processes
    - Default multiprocessing.Pool on Unix systems uses a 'fork' start method, which will
      duplicate the entire parent's memory.
    '''

    class __worker_data__:
        def __init__(self):
            self.available = True
            self.name = None
            self.cmd = None
            self.execute = threading.Semaphore(0)

    def __init__(self, nproc):
        assert nproc > 0
        self.nproc = nproc
        self.names = set()
        self.done = False
        self.terminated = False
        self.count = 0
        self.thunk_scheduler = threading.Semaphore(0)
        self.worker_available = threading.Semaphore(0)
        self.worker_data = []
        self.workers = []
        self.queue = []
        self.queue_lock = threading.Lock()
        self.cv = threading.Condition(lock=self.queue_lock)
        self.dt = threading.Thread(target=self.dispatcher, name='dispatcher')
        self.running = {}
        self.finished = {}

        # setup worker threads and worker data
        for i in range(nproc):
            self.workers.append(threading.Thread(target=self.worker, name=f'worker_{i}', args=(i,)))
            self.worker_data.append(self.__worker_data__())

        # start threads
        for w in self.workers:
            w.start()
        self.dt.start()

    def schedule(self, name, cmd, ignore_sigint=True):
        '''
        Args
        - name: hashable object (e.g., str, int, tuple, etc.)
        - cmd: dict
            Keyword arguments to pass to subprocess.Popen()
            Keys (str) are argument names. Some of the args listed below:
            - args: iterable of str, or str
                Sequence of program arguments, where the program to execute is the first item in args
                On POSIX, if args is a string, the string is interpreted as the name or path of the program to execute.
            - preexec_fn: callable
                Called in the child process just before the child is executed
        - ignore_sigint: bool. default=True
            Ignore SIGINT (e.g., Ctrl-C) signals sent to worker processes. Useful to prevent
            keyboard interrupts from interruping ThreadPool in interactive sessions (e.g., Jupyter Notebook).
            Call ThreadPool.terminate() to terminate subprocesses.
        '''
        assert name not in self.names
        assert not self.terminated
        self.names.add(name)
        if ignore_sigint:
            cmd.setdefault('preexec_fn', lambda: signal.signal(signal.SIGINT, signal.SIG_IGN))
        with self.queue_lock:
            if self.terminated:
                return
            self.queue.append((name, cmd))
            self.count += 1
            self.thunk_scheduler.release()

    def get_available_worker(self):
        for i in range(self.nproc):
            if self.worker_data[i].available:
                return i
        return None

    def dispatcher(self):
        while True:
            # block until the queue of outstanding functions is nonempty
            self.thunk_scheduler.acquire()

            if self.done:
                return

            # wait for a worker thread to become available and select it
            self.worker_available.acquire()
            i = self.get_available_worker()
            assert i is not None
            self.worker_data[i].available = False

            # dequeue the least recently scheduled function
            # put a copy of that function in a place where the selected worker (and only that worker) can find it

            with self.queue_lock:
                if not self.terminated:
                    name, cmd = self.queue.pop(0)
                else:
                    self.worker_data[i].available = True
                    continue
            self.worker_data[i].name = name
            self.worker_data[i].cmd = cmd

            # signal the worker thread to execute the thunk
            self.worker_data[i].execute.release()

    def worker(self, i):
        self.worker_available.release()
        while True:
            # block until the dispatcher thread signals worker to execute an assigned function
            self.worker_data[i].execute.acquire()

            # if done but not available, that means the worker was assigned a thunk but not yet executed
            # before the ThreadPool destructor was called
            if self.done and self.worker_data[i].available:
                return

            # invoke the function, wait for it to execute
            name = self.worker_data[i].name
            cmd = self.worker_data[i].cmd
            p = subprocess.Popen(**cmd)
            self.running[name] = p

            # if terminate() was called before adding p to self.running,
            # manually terminate the subprocess
            if self.terminated and p.poll() is None:
                p.terminate()
            p.wait()

            with self.cv:
                self.finished[name] = self.running.pop(name)
                self.count -= 1
                if self.count == 0:
                    self.cv.notify_all()

            self.worker_data[i].available = True
            self.worker_data[i].name = None
            self.worker_data[i].cmd = None
            self.worker_available.release()

    def status(self, _print=True, _return=False):
        '''
        Args
        - _print: bool. default=True
            Print progress to stdout
        - _return: bool. default=False
            Return progress as (never_run, queued, running, error, success)
            where each element of the tuple is an iterable of process names

        Returns: 5-tuple or None
        '''
        success = []
        error = []
        with self.queue_lock:
            running = list(self.running.keys())
            for name, p in self.finished.items():
                if p.poll() == 0:
                    success.append(name)
                else:
                    error.append(name)
            queued = [name for name, _ in self.queue]
        never_run = self.names - set(running + success + error + queued)
        if _print:
            print('Never run:', never_run, end='\n\n')
            print('Queued:', queued, end='\n\n')
            print('Running:', running, end='\n\n')
            print('Finished with error:', error, end='\n\n')
            print('Finished successfully:', success)
        if _return:
            return (never_run, queued, running, error, success)

    def terminate(self):
        with self.queue_lock:
            # clear queue
            self.queue.clear()

            # terminate processes
            for p in self.running.values():
                p.terminate()

            # clear thunk_scheduler
            while self.thunk_scheduler.acquire(blocking=False):
                pass

            self.terminated = True

    def join(self):
        # wait until all scheduled functions have executed to completion
        with self.cv:
            if not self.terminated:
                self.cv.wait_for(lambda: self.count == 0)
        self.done = True

        # inform dispatcher and worker threads to exit
        self.thunk_scheduler.release()
        for i in range(len(self.worker_data)):
            self.worker_data[i].execute.release()

        # wait for dispatcher and worker threads to exit
        self.dt.join()
        for w in self.workers:
            w.join()
