"""Solar Feature Tracking; algorithms maintained by mib-unitn/SoFT.

Example: SoFT.detection(image, l_thr=20, sign="pos").
The full upstream core is available as SoFT.core; optional helpers as SoFT.extras.
"""
from importlib import import_module

__all__ = ["track_all", "detection", "identification", "associate", "core", "extras"]


def __getattr__(name):
    if name.startswith("_"):
        raise AttributeError(name)
    if name == "extras":
        try:
            return import_module("._vendor.soft.extras", __package__)
        except ModuleNotFoundError as exc:
            raise ModuleNotFoundError(
                "SoFT extras require: pip install 'WaLSAlib[soft-extras]'"
            ) from exc
    core = import_module("._vendor.soft.soft", __package__)
    if name == "core":
        return core
    return getattr(core, name)
