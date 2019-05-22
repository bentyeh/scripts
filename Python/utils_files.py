import bz2, gzip, lzma, os, re

def createValidPath(s, sub_spaces=None):
    '''
    Create valid path name.
    - Remove characters that are not alphanumeric, space, underscore, dash, or period
    - [Optional] Remove whitespace
    - Normalize pathname

    Args
    - s: str
        String to convert to valid path
    - sub_spaces: str. default=None
        String with which to replace whitespace, e.g., '_'

    Return: str

    Note: Inspired by https://stackoverflow.com/a/295466 and django's get_valid_filename()
      https://github.com/django/django/blob/97d3321e89c8d4434927bdbc308db1ccffa99d3b/django/utils/text.py#L221
    '''
    s = re.sub(r'(?u)[^-\w. ]', '', s.strip()).strip()
    if sub_spaces is not None:
        s = re.sub(r' ', sub_spaces, s).strip()
    return os.path.normpath(s)

def createFileObject(file, mode=None):
    '''
    Create file object.

    Args
    - file: str
        Path to file. If ends with '.gz', '.bz2', or '.xz', appropriate compression is assumed.
    - mode: str. default=None
        Mode to open file with.

    Returns: file object (io.TextIOBase)
    '''
    assert(type(file) is str)
    if file.endswith('.gz'):
        mode = mode if mode is not None else 'rt'
        f = gzip.open(file, mode=mode)
    elif file.endswith('.bz2'):
        mode = mode if mode is not None else 'rt'
        f = bz2.open(file, mode=mode)
    elif file.endswith('.xz'):
        mode = mode if mode is not None else 'rt'
        f = lzma.open(file, mode=mode)
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
