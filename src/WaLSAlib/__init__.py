"""Auxiliary routines for spectroscopy and solar feature tracking.

LineFit and linefit refer to the same unchanged fitting function. Tracking
namespaces load on demand; importing WaLSAlib never imports the ML stack.
"""
from importlib import import_module
from importlib.metadata import PackageNotFoundError, version

try:
    __version__ = version("WaLSAlib")
except PackageNotFoundError:
    __version__ = "0+unknown"

__all__ = ["LineFit", "linefit", "SoFT", "SoFTg"]


def __getattr__(name):
    if name in {"LineFit", "linefit"}:
        function = import_module(".linefit", __name__).linefit
        globals().update(LineFit=function, linefit=function)
        return function
    if name in {"SoFT", "SoFTg"}:
        module = import_module("." + name, __name__)
        globals()[name] = module
        return module
    if name == "WaLSAtools":
        module = import_module("WaLSAtools")
        globals()[name] = module
        return module
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
