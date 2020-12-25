import os, signal, subprocess, sys, threading

def isint(s):
    '''
    Check if a string represents an int.
    Note: Numbers in scientific notation (e.g., 1e3) are floats in Python,
          and `isint('1e3')` will return False.
    Source: https://stackoverflow.com/a/1265696
    '''
    if len(s) < 1:
        return False
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
      2. Use `isfloat()` first, then `float(s).is_integer()` if the former returns True.

    Source: https://stackoverflow.com/a/15357477
    '''
    try:
        float(s)
    except ValueError:
        return False
    return True

def wrap_structure(x, structure):
    '''
    Enclose object in a data structure if it is not already of that type.

    Args
    - x: object
    - structure: str or data structure
        'list', list, 'tuple', tuple, 'set', or set
        Raises ValueError if unsupported data structure given.

    Returns: list, tuple, or set
    '''
    if structure in ('list', list):
        return x if type(x) == list else [x]
    elif structure in ('tuple', tuple):
        return x if type(x) == tuple else (x,)
    elif structure in ('set', set):
        return x if type(x) == set else {x}
    else:
        raise ValueError(f'{str(structure)} not supported. Returning original object.')

def get_available_cpus(job_scheduler=None):
    '''
    Get the number of CPUs the current process can use.

    Args
    - job_scheduler: str. default=None
        Job scheduling environment. Currently supports 'SLURM'

    Returns: int or None
      May return None if the number of CPUs in the system is undetermined.

    See https://docs.python.org/3/library/os.html#os.cpu_count.
    '''
    try:
        return len(os.sched_getaffinity(0))
    except:
        try:
            if job_scheduler and job_scheduler.upper() == 'SLURM':
                return os.environ['SLURM_CPUS_PER_TASK']
        except:
            pass
    print('Unable to detect number of available CPUs. Returning number of total CPUs.', file=sys.stderr)
    return os.cpu_count()

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

class empty_context:
    def __enter__(self):
        pass
    def __exit__(self, exc_type, exc_value, traceback):
        pass

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

    def __del__(self):
        # timeout is necessary to prevent deadlock - see https://docs.python.org/3/reference/datamodel.html#object.__del__
        self.terminate(timeout=10)

    def schedule(self, name, cmd, ignore_sigint=True):
        '''
        Args
        - name: hashable object (e.g., str, int, tuple, etc.)
        - cmd: dict
            Keyword arguments to pass to subprocess.Popen()
            Keys (str) are argument names. Some of the args listed below:
            - args: iterable of str, or str
                Sequence of program arguments, where the program to execute is the first item in args.
                On POSIX, if args is a string, the string is interpreted as the name or path of the program to execute.
            - preexec_fn: callable
                Called in the child process just before the child is executed
        - ignore_sigint: bool. default=True
            Ignore SIGINT (e.g., Ctrl-C) signals sent to worker processes. Useful to prevent keyboard interrupts from
            interruping ThreadPool in interactive sessions (e.g., Jupyter Notebook).
            Call ThreadPool.terminate() to terminate subprocesses.
        '''
        with self.queue_lock:
            if self.terminated:
                print('Pool has been terminated and can no longer schedule commands.', file=sys.stderr)
                return
            if name in self.names:
                print(f'{name} has already been used as the name of a process. Not scheduling ...', file=sys.stderr)
                return
            self.names.add(name)
            if ignore_sigint:
                cmd.setdefault('preexec_fn', lambda: signal.signal(signal.SIGINT, signal.SIG_IGN))
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
                    self.worker_data[i].name = name
                    self.worker_data[i].cmd = cmd
                else:
                    self.worker_data[i].available = True
                    self.worker_available.release()
                    continue

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
            # acquiring the queue_lock is necessary to prevent a race condition with terminate()
            with self.queue_lock:
                if self.terminated:
                    self.worker_data[i].available = True
                    self.worker_data[i].name = None
                    self.worker_data[i].cmd = None
                    self.worker_available.release()
                    continue
                name = self.worker_data[i].name
                cmd = self.worker_data[i].cmd
                p = subprocess.Popen(**cmd)
                self.running[name] = p

            # Calling p.wait() without any lock assumes that it is thread-safe, which appears to be the case based on
            # https://bugs.python.org/issue21291. The subprocess module also only notes a lack of thread-safety once
            # when using the preexec_fn parameter with subprocess.Popen(). Specifically, we assume that there is no race
            # condition involved between calling p.wait() and p.terminate().
            # 
            # Note: subprocess.Popen.wait() is currently implemented by CPython using busy waiting :(
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
                print('Finished successfully:', success, flush=True)
            if _return:
                return (never_run, queued, running, error, success)

    def terminate(self, timeout=-1):
        lock_acquired = self.queue_lock.acquire(timeout=timeout)
        if lock_acquired:
            try:
                # clear queue
                self.queue.clear()

                # terminate processes
                for p in self.running.values():
                    p.terminate()

                # clear thunk_scheduler
                while self.thunk_scheduler.acquire(blocking=False):
                    pass

                self.terminated = True
            finally:
                self.queue_lock.release()
        else:
            print('terminate() called unsuccessfully. Could not acquire queue lock.', file=sys.stderr)

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
