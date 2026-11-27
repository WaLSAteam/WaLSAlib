# Contributing

Discuss new methods/modules in a GitHub issue before submitting. Include documentation, dependency and license information, citation, and meaningful numerical/interface tests. Native LineFit changes are outside this packaging/integration revision.

Use Python 3.10+, create a virtual environment, and install `python -m pip install -e '.[test,build]'`. Add `[softg]` for ML validation. Run `python tools/sync_upstreams.py --check`, `python -m pytest tests/release`, `python -m build`, `python -m twine check --strict dist/*`, and `python tools/check_dist.py dist`.

Never edit generated _vendor files directly. See docs/upstream.md for updates and docs/release.md for publishing. New public packages/modules require changes to both the package and module release allowlists. Keep private research, observational data, checkpoints, convex-hull work, and new hybrid runs out of release commits. Stage explicit paths, inspect `git diff --cached --stat` and `git diff --cached`, and preserve the original Frontiers proof-of-principle examples.

Keep API and scientific behavior changes documented. Follow .github/CODE_OF_CONDUCT.md and report security concerns as described in .github/SECURITY.md.
