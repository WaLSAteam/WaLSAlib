"""
LineFit module (adaptive Voigt spectral line fitting).

The implementation lives in WaLSA_LineFit.py.
We provide a simple alias:
    linefit.linefit  -> WaLSA_LineFit
"""

from .WaLSA_LineFit import WaLSA_LineFit as linefit

__all__ = ["linefit"]