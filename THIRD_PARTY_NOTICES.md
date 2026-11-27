# Third-party source and distribution terms

WaLSAlib's original code, including the unchanged LineFit implementation, remains under the Apache License 2.0 in LICENSE. New WaLSAlib integration/build/test code is also Apache-2.0. Vendored upstream code is not relicensed.

| Component | Author / source | License | Source snapshot |
| --- | --- | --- | --- |
| SoFT | Michele Berretti and upstream contributors, https://github.com/mib-unitn/SoFT | GNU GPL v3 | vendor/SoFT |
| SoFTg | mib-unitn/SoFTg contributors, https://github.com/mib-unitn/SoFTg | GNU GPL v3 | vendor/SoFTg |

Exact commit IDs, selected source files, and SHA-256 hashes are in vendor/upstream-lock.json. Original license texts, READMEs, and dependency declarations accompany the snapshots. The `GPL-3.0-only` identifier is used conservatively because the repository-level licenses do not explicitly grant an “or later” option.

Bundling these GPL components means the combined distribution must satisfy GPLv3, including provision of corresponding source. It cannot be advertised as Apache-only. Native Apache grants and notices remain applicable to their files. See the Apache Software Foundation's explanation: https://www.apache.org/licenses/GPL-compatibility.html.

The source distribution contains the original snapshots, generated runtime sources, adapter generator, dependency declarations, and build instructions. Release both the wheel and source distribution together. No upstream NOTICE files exist in the pinned snapshots. Preserve any new notices when accepting future updates.

## Modifications to upstream code

For the 1.1.0 integration (2026-11-27), tools/sync_upstreams.py generates runtime copies under src/WaLSAlib/_vendor. Absolute sibling imports are changed to relative imports; upstream SoFT's initializer is replaced with a minimal private initializer. SoFT's shell-based `rm -rf {datapath}temp*` cleanup is replaced with Python removal limited to non-symlink directories temp0 through temp9. This supports spaces in paths and preserves unrelated temp-prefixed folders. No scientific calculations, parameters, or numerical expressions are changed. Original copies remain byte-for-byte in vendor/. Wrapper APIs outside _vendor are WaLSAlib code.

Dependency packages and user-supplied model weights retain their own licenses. No datasets, trained weights, upstream bytecode, or upstream example products are redistributed in the Python package.
