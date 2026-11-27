"""Granule segmentation and optional UNet++ workflows from mib-unitn/SoFTg.

Classical segmentation works with the base installation. Model, prediction,
and training workflows require pip install 'WaLSAlib[softg]'. Trained weights
are supplied by the user; no data or checkpoints are downloaded automatically.
"""
from importlib import import_module

__all__ = ["segment_solar_granules", "granule_detection_lib", "model", "predict", "train", "dataset", "losses", "utils"]
_MODULES = set(__all__) - {"segment_solar_granules"}
_MODULES.add("label_tool")


def __getattr__(name):
    if name == "segment_solar_granules":
        return import_module("._vendor.softg.granule_detection_lib", __package__).segment_solar_granules
    if name in _MODULES:
        try:
            return import_module("._vendor.softg." + name, __package__)
        except ModuleNotFoundError as exc:
            extra = "softg-gui" if name == "label_tool" else "softg"
            raise ModuleNotFoundError(
                f"SoFTg {name} requires: pip install 'WaLSAlib[{extra}]'"
            ) from exc
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
