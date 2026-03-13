# WaLSA-Pipelines

WaLSA-Pipelines is a companion repository to [**WaLSAtools**](https://github.com/WaLSAteam/WaLSAtools), collecting practical workflow routines developed within the WaLSA Team to prepare *analysis-ready* products and time series for wave/oscillation studies (and related diagnostics). The focus here is on upstream and auxiliary steps that often sit between calibrated observations and the actual wave analysis.

## Scope and philosophy
Many routines in WaLSA-Pipelines are motivated by (and initially developed for) **Solar Physics** use cases — e.g., spectroscopy/spectro-polarimetry and time-series imaging — although they may be useful more broadly.

## Current status
Early development. **The first module included is `LineFit`**, an adaptive **Voigt** spectral line fitting routine designed for robust extraction of line parameters (e.g., intensity and LOS velocity proxies) from high-resolution spectra.

More routines will be added in the near future (e.g., motion/registration utilities, additional spectroscopy helpers, and pipeline-level glue code).

## Relationship to WaLSAtools
- **WaLSA-Pipelines**: preprocessing / preparation workflows (e.g., fitting, motion magnification, feature detection)
- **WaLSAtools**: wave and oscillation analysis methods (FFT/Welch, wavelets, cross-spectra, significance, etc.)

## Installation (development)
For now, install from a local clone:

```bash
git clone https://github.com/WaLSAteam/WaLSA-Pipelines.git
cd WaLSA-Pipelines
python -m pip install -e .
```

## Quick start

```python
from WaLSA_pipelines import linefit
```

(See examples/ for current demos)

## Contributing

Contributions are welcome. If you plan to add a new routine/module, please open an issue first so we can agree on scope, dependencies, and how it fits the overall structure.

## License

Apache License 2.0 (see LICENSE).
