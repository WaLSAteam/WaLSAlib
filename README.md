<div style="display: block; margin-bottom: 20px;">
    <img align="left" src="https://raw.githubusercontent.com/WaLSAteam/WaLSAlib/main/docs/images/WaLSAlib_logo.svg" alt="WaLSAlib" width="400" height="auto">
</div>

<br><br><br>

# WaLSAlib &ndash; WaLSA Library of Auxiliary Routines
<p align="left">
    <a href="#"><img src="https://img.shields.io/pypi/v/WaLSAlib"></a>
    <a href="https://walsa.team" target="_blank"><img src="https://img.shields.io/badge/powered%20by-WaLSA%20Team-000d1a"></a>
    <a href="#"><img src="https://img.shields.io/badge/license-Apache%202.0%20%2B%20GPLv3-green"></a>
    <a href="https://doi.org/10.5281/zenodo.23043078"><img src="https://img.shields.io/badge/DOI-10.5281%2Fzenodo.23043078-blue.svg" alt="DOI">
</a>
</p>

**WaLSAlib** is a companion repository to [**WaLSAtools**](https://github.com/WaLSAteam/WaLSAtools), collecting practical workflow routines developed within the WaLSA Team to prepare *analysis-ready* products and time series for wave/oscillation studies (and related diagnostics). The focus here is on upstream and auxiliary steps that often sit between calibrated observations and the actual wave analysis.

## Scope and philosophy
Many routines in WaLSAlib are motivated by (and initially developed for) **Solar Physics** use cases — e.g., spectroscopy/spectro-polarimetry and time-series imaging — although they may be useful more broadly.

## Included modules

- **LineFit**: adaptive multi-line spectral fitting, line-core intensities, and LOS velocities. Its published implementation and original proof-of-principle notebooks are unchanged.
- **SoFT**: solar feature detection, association, and tabulation, maintained independently at [mib-unitn/SoFT](https://github.com/mib-unitn/SoFT).
- **SoFTg**: classical granule segmentation and optional UNet++ training/inference, maintained independently at [mib-unitn/SoFTg](https://github.com/mib-unitn/SoFTg).

Upstream algorithms are pinned to reviewed commits, bundled with the package, and adapted for Python namespace imports and bounded temporary-directory cleanup. Installation never fetches Git repositories. See [upstream integration](https://github.com/WaLSAteam/WaLSAlib/blob/main/docs/upstream.md).

## Relationship to WaLSAtools
- **WaLSAlib**: preprocessing / preparation workflows (e.g., fitting, motion magnification, feature detection)
- **WaLSAtools**: wave and oscillation analysis methods (FFT/Welch, wavelets, cross-spectra, significance, etc.)

## Installation

Python 3.10 or newer is required. After the first PyPI publication of this release:

```bash
python -m pip install WaLSAlib
```

Until publication, install this checkout:

```bash
python -m pip install .
# For development:
python -m pip install -e '.[test,build]'
```

The base installation includes LineFit, SoFT tracking, and classical SoFTg segmentation. Optional workflows:

```bash
python -m pip install 'WaLSAlib[softg]'        # UNet++ training and prediction
python -m pip install 'WaLSAlib[softg-gui]'    # interactive labeling (desktop/display required)
python -m pip install 'WaLSAlib[soft-extras]' # SoFT SunPy / EMD helpers
python -m pip install 'WaLSAlib[walsatools]'   # companion wave analysis library
```

## Quick start

```python
from WaLSAlib import LineFit, SoFT, SoFTg
# Existing code also works: from WaLSAlib import linefit

# LineFit is the existing linefit callable, with the same arguments and results.
# result = LineFit(spectra, wavelengths=wavelengths, refwavelength=centres, ...)

# image is a 2D NumPy array; parameters depend on instrument and data units.
# magnetic_labels = SoFT.detection(image, l_thr=20, sign="pos")
# labels, seeds, skeleton = SoFTg.segment_solar_granules(image, plate_scale=50)
```

Use `SoFT.track_all(...)` for the upstream FITS tracking pipeline and `SoFT.extras` for optional analysis helpers. Give tracking a dedicated working directory containing `00-data/` and a trailing path separator: the pipeline resets its `01-mask/`, `02-id/`, and `03-assoc/` outputs on each run. ML tools live at `SoFTg.model`, `SoFTg.predict`, and `SoFTg.train`. For command-line workflows:

```bash
python -m WaLSAlib._vendor.softg.train --help
python -m WaLSAlib._vendor.softg.predict --help
```

Provide your own training data and trained checkpoint; model weights are not bundled. Importing WaLSAlib does not load PyTorch, download weights, or choose a GUI backend.

Original LineFit publication examples remain in [examples/LineFit_Frontiers](https://github.com/WaLSAteam/WaLSAlib/tree/main/examples/LineFit_Frontiers), including the original PoP CNN notebook and outputs. They are repository examples and are excluded from PyPI artifacts.

## Contributing and releases

See [CONTRIBUTING.md](https://github.com/WaLSAteam/WaLSAlib/blob/main/CONTRIBUTING.md), the [repository review](https://github.com/WaLSAteam/WaLSAlib/blob/main/docs/repository-review.md), and [release instructions](https://github.com/WaLSAteam/WaLSAlib/blob/main/docs/release.md). Upstream refreshes propose pull requests; scientific review and acceptance remain with maintainers.

## Citation and licenses

Please cite the original methods as well as WaLSAlib ([Zenodo](https://doi.org/10.5281/zenodo.23043078)). SoFT's reference is [Berretti et al. (2025)](https://doi.org/10.1051/0004-6361/202452665). SoFTg's upstream README requests citation of Berretti et al. (2026); consult the upstream project for its final bibliographic details.

Native WaLSAlib code remains Apache-2.0 (see LICENSE). SoFT and SoFTg retain GPLv3. **The combined distribution must be redistributed in compliance with GPLv3**, while retaining the Apache notices for the native code. The SPDX expression `Apache-2.0 AND GPL-3.0-only` records both sets of obligations; it is not a choice of licenses. See [THIRD_PARTY_NOTICES.md](https://github.com/WaLSAteam/WaLSAlib/blob/main/THIRD_PARTY_NOTICES.md) for provenance and adaptation details.
