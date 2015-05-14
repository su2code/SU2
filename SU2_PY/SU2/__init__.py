# SU2/__init__.py

class EvaluationFailure(RuntimeError):
    pass
class DivergenceFailure(EvaluationFailure):
    pass


import run
import io
import mesh
import eval
import opt
import util

try:
    import readline
    import rlcompleter
    if readline.__doc__ and 'libedit' in readline.__doc__:
        readline.parse_and_bind("bind ^I rl_complete")
    else:
        readline.parse_and_bind("tab: complete")
except:
    pass




