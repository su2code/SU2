# SU2/__init__.py

from . import run
from . import io
from . import mesh
from . import eval
from . import opt
from . import util

try:
    import readline
    import rlcompleter
    if 'libedit' in readline.__doc__:
        readline.parse_and_bind("bind ^I rl_complete")
    else:
        readline.parse_and_bind("tab: complete")
except ImportError:
    pass
