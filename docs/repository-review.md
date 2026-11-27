# General repository review — 2026-11-27

Scope: public baseline dbfed37, release infrastructure, and SoFT/SoFTg integration. The LineFit algorithm was not revised. Local convex-hull, hybrid acceleration, research, and hybrid_runs work is excluded. Existing published Frontiers examples, including PoP_CNN_emulator, are retained in Git.

| Finding | Revision |
| --- | --- |
| pyproject version 0.1.0 disagreed with README/tag v1.0.0 | Prepare next release as 1.1.0; tag/version check before publication |
| Python >=3.9 claim conflicted with evaluated `list[...] | None` annotations | Declare >=3.10; CI matrix covers 3.10/3.12/3.13 |
| Wildcard package discovery could ship local experiments | Explicit packages plus build-time module allowlist; sdist pruning and artifact audit |
| Requested LineFit name was missing | Add capitalised alias while preserving lowercase linefit; implementation unchanged |
| Eager imports loaded optional WaLSAtools and masked all exceptions | Lazy public imports; WaLSAtools becomes an explicit optional extra |
| Repository badge/logo/links were unsuitable for PyPI README | Remote logo, dynamic version badge, absolute documentation links |
| Little release metadata or build validation | SPDX license metadata/files, authors, classifiers, URLs; build/twine checks; wheel and sdist installation tests |
| SoFT shell cleanup interpolated paths and deleted temp* indiscriminately | Generated adapter uses Python cleanup limited to numbered temporary directories; tested paths with spaces and unrelated-folder preservation |
| Independent upstream code was not Python-package ready | Locked original sources, deterministic relative-import adapters, public namespaces, tests |
| GPL upstream code cannot be described as Apache-only | Preserve Apache-native terms and GPL source/notices; combined distribution documents GPLv3 obligations |
| Large ML/solar-download/GUI stacks would burden all users | Base scientific tracking dependencies; separate ML, SoFT helpers, desktop, WaLSAtools extras |
| No controlled update/release process | Weekly upstream PR proposals with validation; manual tag-only Trusted Publishing workflow |
| No contributor instructions | Add contribution/release/upstream guides and changelog |

## Remaining maintainer actions

Configure PyPI/TestPyPI ownership and Trusted Publishers, protected publishing environments, Actions PR creation, and main-branch review requirements. These account settings are not part of source code. Review the combined-distribution GPL obligations before publishing. No PyPI upload has been performed.

Scientific acceptance of tracking/segmentation changes requires reference observations and domain review beyond packaging tests. SoFTg supplies no trained checkpoint in this snapshot; no trained-model accuracy claim is made. Desktop GUI and GPU workflows need platform-specific validation. Dependency lower bounds are declared compatibility floors, not a fully tested minimum-version matrix.

Local installed-artifact verification on macOS arm64 / Python 3.12 passed all nine API/scientific contract tests, including the ML forward/backward, packaged CLI, and four-frame FITS/multiprocessing tracking checks. scikit-image 0.26 warns that the upstream erosion/dilation calls will be removed in 0.28; the dependency is therefore capped below 0.28 until upstream adapts them. The algorithms were left unchanged.

Additional checks passed: clean sdist rebuild and base-only install (eight passed, one ML skip); SoFT optional analysis imports; dependency consistency; actionlint validation; and exclusion canaries covering known/unknown local modules and stale build outputs. PyPI returned HTTP 404 for the project name during preparation; this does not guarantee name availability or account ownership.
