# SU2/__init__.py

class EvaluationFailure(RuntimeError):
    pass
class DivergenceFailure(EvaluationFailure):
    pass


# Please do not remove next imports
# this is in place to save the need for additional import lines in user scripts
# It's important for the whole python package to be accessible with one import SU2
# See issue #246

from SU2 import run
from SU2 import io
from SU2 import eval
from SU2 import opt
from SU2 import util

try:
    import readline
    import rlcompleter
    if readline.__doc__ and 'libedit' in readline.__doc__:
        readline.parse_and_bind("bind ^I rl_complete")
    else:
        readline.parse_and_bind("tab: complete")
except:
    pass




