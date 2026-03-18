<div style="display: block; margin-bottom: 20px;">
    <img align="left" src="https://walsa.team/WaLSAlib_logo.svg" alt="WaLSAlib" width="400" height="auto">
</div>

<br><br><br>

# WaLSAlib &ndash; WaLSA Library of Auxiliary Routines
<p align="left">
    <a href="#"><img src="https://img.shields.io/badge/WaLSAlib-v1.0.0-0066cc"></a> 
    <a href="https://walsa.team" target="_blank"><img src="https://img.shields.io/badge/powered%20by-WaLSA%20Team-000d1a"></a>
    <a href="#"><img src="https://img.shields.io/badge/license-Apache%202.0-green"></a>
</p>

**WaLSAlib** is a companion repository to [**WaLSAtools**](https://github.com/WaLSAteam/WaLSAtools), collecting practical workflow routines developed within the WaLSA Team to prepare *analysis-ready* products and time series for wave/oscillation studies (and related diagnostics). The focus here is on upstream and auxiliary steps that often sit between calibrated observations and the actual wave analysis.

## Scope and philosophy
Many routines in WaLSAlib are motivated by (and initially developed for) **Solar Physics** use cases — e.g., spectroscopy/spectro-polarimetry and time-series imaging — although they may be useful more broadly.

## Current status
Early development. **The first module included is `LineFit`**, an adaptive multi-line fitting routine for dense, evolving spectra. It provides robust line-centre tracking and extracts line-core intensity and LOS-velocity measurements from either a single spectrum or a sequence, combining stable coarse seeding near expected centres, conservative per-line windowing with safety bounds, and bounded Voigt-family fitting (with optional asymmetry handling). It is particularly aimed at dense spectral windows where profiles evolve, become asymmetric, blend (including close pairs), or intermittently show split/self-reversal morphologies.

More routines will be added in the near future (e.g., motion/registration utilities, additional spectroscopy helpers, and pipeline-level glue code).

## Relationship to WaLSAtools
- **WaLSAlib**: preprocessing / preparation workflows (e.g., fitting, motion magnification, feature detection)
- **WaLSAtools**: wave and oscillation analysis methods (FFT/Welch, wavelets, cross-spectra, significance, etc.)

## Installation (development)
For now, install from a local clone:

```bash
git clone https://github.com/WaLSAteam/WaLSAlib.git
cd WaLSAlib
python -m pip install -e .
```

## Quick start

```python
from WaLSAlib import linefit
```

(See `examples/` for current demos)

## Contributing

Contributions are welcome. If you plan to add a new routine/module, please open an issue first so we can agree on scope, dependencies, and how it fits the overall structure.

## License

Apache License 2.0 (see LICENSE).
