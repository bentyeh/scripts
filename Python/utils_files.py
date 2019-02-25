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

def readFile(file):
    '''
    Reads entire file and returns lines as a list.
    '''
    
    with createFileObject(file) as f:
        contents = f.read().splitlines()
    return contents
