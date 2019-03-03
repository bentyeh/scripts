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