import gzip

def createFileObject(file, mode=None):
    '''
    Create file object.
    
    Args
    - file: str
        Path to file. If extension ends with '.gz', gzip compression is assumed.
    - mode: str. default=None
        Mode to open file with.
    
    Returns: file object (io.TextIOBase)
    '''
    
    assert(type(file) is str)
    if file.endswith('.gz'):
        mode = mode if mode is not None else 'rt'
        f = gzip.open(file, mode=mode)
    else:
        mode = mode if mode is not None else 'r'
        f = open(file, mode=mode)
    return f

def readFile(file, strip=None):
    '''
    Reads entire file and returns lines as a list.

    Args
    - file: str
        Path to file. If extension ends with '.gz', gzip compression is assumed.
    - strip: str. default=None.
        Parameter to str.strip() applied to each read line.
        - None: strip whitespace
        - '': do not strip
    
    Returns: list of str
      Each element in the list is a line from the file.
    '''
    
    with createFileObject(file) as f:
        contents = f.read().strip(strip).splitlines()
    return contents
