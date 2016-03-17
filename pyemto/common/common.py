def lat_to_ibz(lat):
    """Returns the Bravais lattice ibz code based on the input string code.

    e.g. lat='bcc' or lat='fco'

    :param lat: lattice string code
    :type lat: str
    :returns: ibz code corresponding to the Bravais lattice
            given by the 'lat' key
    :rtype: int
    """

    # One should add some error hangling here.
    ltoi = {'sc': 1, 'fcc': 2, 'bcc': 3, 'hcp': 4, 'st': 5, 'bct': 6, 'trig': 7, 'so': 8,
            'baco': 9, 'bco': 10, 'fco': 11, 'sm': 12, 'bacm': 13, 'bcm': 13, 'stric': 14,
            'B2':1,'L12':1}

    return ltoi[lat]

def ibz_to_lat(lat):
    """Returns the Bravais lattice ibz code based on the input string code.

    e.g. lat='bcc' or lat='fco'

    :param lat: lattice string code
    :type lat: str
    :returns: ibz code corresponding to the Bravais lattice
            given by the 'lat' key
    :rtype: int
    """

    # One should add some error hangling here.
    ltoi = {1:'sc', 2:'fcc', 3:'bcc', 4:'hcp', 5:'st', 6:'bct', 7:'trig', 8:'so',
            9:'baco', 10:'bco', 11:'fco', 12:'sm', 13:'bacm', 13:'bcm', 14:'stric'}

    return ltoi[lat]


def check_folders(*args):
    """Checks whether or not given folders exist.

    :param *args: The name of the folder or a list of names of folders
    :type *args: str or list(str)
    :returns: None
    :rtype: None
    """

    import os

    for arg in args:
        if not os.path.exists(arg):
            os.makedirs(arg)
    return

def cleanup_path(string):
    """Cleans up directory path strings from double forward slashes.

    :param string: Cleaned up path string
    :type string: str
    :returns: Cleaned up directory path string
    :rtype: str
    """

    import re

    string = re.sub('\/+', '/', string)
    #string = string.lstrip('/')

    return string
