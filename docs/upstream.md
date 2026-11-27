# Independent upstream integration

## Decision

Use **pinned, filtered source vendoring**, with reviewed update PRs. This is an explicit alternative to a literal `git subtree` import, not a Git submodule or subtree history merge. Like subtrees, files are present in ordinary clones, wheels, and sdists. It avoids importing upstream tracked __pycache__, examples/data, and unrelated history; it also makes SoFTg's script-to-package import and bounded filesystem cleanup adaptations reproducible. Neither independent repository needs to change, and users do not need Git installed.

Upstream originals live in vendor/SoFT and vendor/SoFTg. A lock records repository, branch, immutable SHA, file selection, and hashes. Generated runtime copies live in WaLSAlib._vendor; do not edit them. Public adapters are WaLSAlib.SoFT and WaLSAlib.SoFTg. We do not install a global `soft`, `utils`, `dataset`, or `model` package, and do not manipulate sys.path or alias global modules.

## Refresh

```bash
python tools/sync_upstreams.py --check  # offline verification
python tools/sync_upstreams.py --update # fetch upstream main, refresh pins, generate
python tools/sync_upstreams.py         # regenerate only, without network
```

An update fetches into disposable Git repositories and reads selected blobs without executing upstream code. All existing selected paths must remain available. License changes and unrecognised import layouts fail for manual review. New upstream files/dependencies require a maintainer to extend the lock/adapters/dependency extras and build allowlist; selection does not silently grow. A changed upstream README or requirements file appears in the PR for review.

The weekly `upstream-sync.yml` workflow (also manually dispatchable) generates a candidate, checks hashes/adapters, builds artifacts, runs base and ML contract tests in a job with read-only credentials, and uploads a patch. A separate job with write permissions commits only that patch to the reserved `upstream-sync` branch and opens/updates a PR. It never merges. Enable “Allow GitHub Actions to create and approve pull requests” in repository Actions settings. The workflow must be on the default branch for scheduled/manual runs to be available.

Ordinary PR CI also runs for upstream PRs; GitHub may require a maintainer to approve runs initiated by GITHUB_TOKEN. Candidate validation already runs before the proposal job. Protect main and require CI plus human review. Do not use the reserved automation branch for hand edits.

## Validation and scientific review

Tests cover imports, simple synthetic feature separation and size filtering, a four-frame FITS/multiprocessing tracking pipeline with area/flux/velocity checks, SoFTg segmentation output contracts, upstream/adapted numerical parity, a small UNet++ forward/loss step, and artifact contents. They do not establish instrument-independent scientific validity or trained-model accuracy. Before accepting algorithm changes, review the upstream diff between old/new commits and compare established reference observations and scientific metrics. No approved observational reference set or trained SoFTg checkpoint is supplied by the upstream snapshots.

## Future modules

Native modules need entries in build_support.py and pyproject.toml, public exports, dependencies, tests, citation, and documentation. Another imported component needs an independent lock entry, preserved notices, a deterministic adapter if necessary, dependency/license review, and validation. Explicit lists keep unreleased work out of distributions.
