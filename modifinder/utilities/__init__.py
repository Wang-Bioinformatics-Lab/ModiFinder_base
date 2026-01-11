from importlib import import_module

from . import general_utils as gu
from .general_utils import read_mgf, write_mgf

__all__ = [
    "vis",
    "gu",
    "read_mgf",
    "write_mgf",
]


def __getattr__(name):
    if name == "vis":
        return import_module("modifinder.utilities.visualizer")
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
