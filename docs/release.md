# Publishing WaLSAlib

The prepared version is **1.1.0**, following the repository's existing v1.0.0 tag (the old pyproject incorrectly said 0.1.0). This work does not upload to PyPI or create a release tag.

## One-time account setup

1. Confirm ownership/availability of `WaLSAlib` on PyPI and TestPyPI. They use separate accounts. Configure a pending publisher if the project does not exist.
2. Register a GitHub Trusted Publisher on each index: owner `WaLSAteam`, repository `WaLSAlib`, workflow filename `publish.yml`, environment `pypi` or `testpypi` respectively. No long-lived API token is needed.
3. Create the matching GitHub environments. Add required reviewers for `pypi`, and restrict deployment refs to reviewed release tags. Protect main and require the release CI checks.
4. Enable Actions PR creation for upstream sync. Merge the reviewed preparation branch before enabling the schedule.

References: [PyPA publishing guide](https://packaging.python.org/en/latest/guides/publishing-package-distribution-releases-using-github-actions-ci-cd-workflows/) and [PyPI Trusted Publishers](https://docs.pypi.org/trusted-publishers/adding-a-publisher/).

## Validate locally

```bash
python -m venv .venv
. .venv/bin/activate
python -m pip install '.[build,test,softg]'
python tools/sync_upstreams.py --check
python -m build
python -m twine check --strict dist/*
python tools/check_dist.py dist
python -m pip install --force-reinstall --no-deps dist/*.whl
python -m pytest tests/release
```

The CI workflow additionally rebuilds the wheel from the sdist, installs it, and tests that installation. CI tests base installs on Python 3.10/3.12/3.13 and ML integration on 3.12. Linux is the automated target; local verification records any additional platform tested. GPU acceleration and desktop labeling require separate platform checks.

## Release

1. Review licensing/provenance and scientific validation. Ensure the README, version, and changelog describe the release.
2. Merge to main. Create and push `v1.1.0` only when ready; the workflow checks that the tag equals pyproject.toml's version.
3. Manually run `publish.yml` on that tag with target `testpypi`. Inspect artifacts and validate an install in a fresh environment (TestPyPI dependencies may need PyPI as a separate dependency source).
4. Run the same workflow on the same tag with target `pypi`. The build/test job must succeed, then the `pypi` environment reviewer approves upload. Both sdist and wheel publish together.
5. Confirm `python -m pip install WaLSAlib==1.1.0` and all three public imports from an empty directory. Add GitHub release notes and confirm archive/citation metadata.

Publishing is manual (`workflow_dispatch`) and accepts only version tags. A branch push or upstream PR cannot publish. PyPI versions are immutable: correct a failed/revised release with a new version/tag rather than replacing published files.
