# SU2/amginria/__init__.py

from importlib import import_module
from distutils.sysconfig import get_config_var
_amgio = import_module("_amgio"+get_config_var('EXT_SUFFIX'))
sys.modules["_amgio"+get_config_var('EXT_SUFFIX')] = _amgio
# import _amgio

from tools     import *
from interface import *