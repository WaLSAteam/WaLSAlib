#!/usr/bin/env python3
"""Audit wheel and sdist contents, licenses, and PyPI metadata before upload."""
from email.parser import BytesParser
from pathlib import Path, PurePosixPath
import sys
import tarfile
import zipfile

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from build_support import ReleaseBuildPy

EXPECTED_MODULES = {
    package.replace(".", "/") + "/" + module + ".py"
    for package, modules in ReleaseBuildPy._native.items() for module in modules
}
FORBIDDEN = ("convex", "hybrid", "research", "__pycache__", "examples/", ".ds_store")


def check_names(names):
    for name in names:
        if any(word in name.lower() for word in FORBIDDEN) or name.endswith((".pyc", ".pt", ".fits", ".npy", ".npz")):
            raise ValueError(f"Unexpected release content: {name}")
        path = PurePosixPath(name)
        if path.is_absolute() or ".." in path.parts:
            raise ValueError(f"Unsafe archive path: {name}")


def check_metadata(content):
    metadata = BytesParser().parsebytes(content)
    assert metadata["Name"] == "WaLSAlib"
    assert metadata["Requires-Python"] == ">=3.10"
    assert metadata["License-Expression"] == "Apache-2.0 AND GPL-3.0-only"
    assert len(metadata.get_all("License-File", [])) == 4
    assert {"softg", "soft-extras", "softg-gui", "walsatools"} <= set(metadata.get_all("Provides-Extra", []))
    return metadata


def audit(directory):
    wheels = list(Path(directory).glob("*.whl"))
    sdists = list(Path(directory).glob("*.tar.gz"))
    if len(wheels) != 1 or len(sdists) != 1:
        raise ValueError("Expected exactly one wheel and one source distribution")
    with zipfile.ZipFile(wheels[0]) as archive:
        names = archive.namelist()
        check_names(names)
        assert {name for name in names if name.endswith(".py")} == EXPECTED_MODULES
        metadata_name = next(n for n in names if n.endswith(".dist-info/METADATA"))
        metadata = check_metadata(archive.read(metadata_name))
        prefix = metadata_name.removesuffix("METADATA")
        for license_file in metadata.get_all("License-File"):
            assert prefix + "licenses/" + license_file in names
    with tarfile.open(sdists[0]) as archive:
        members = archive.getmembers()
        assert all(m.isfile() or m.isdir() for m in members), "No archive links permitted"
        names = [m.name.partition("/")[2] for m in members if m.isfile()]
        check_names(names)
        assert {n.removeprefix("src/") for n in names if n.startswith("src/WaLSAlib/") and n.endswith(".py")} == EXPECTED_MODULES
        assert {"build_support.py", "tools/sync_upstreams.py", "vendor/upstream-lock.json", "THIRD_PARTY_NOTICES.md", "vendor/SoFT/LICENSE", "vendor/SoFTg/LICENSE"} <= set(names)
        pkg_info = next(m for m in members if m.name.count("/") == 1 and m.name.endswith("/PKG-INFO"))
        assert check_metadata(archive.extractfile(pkg_info).read())["Version"] == metadata["Version"]
    print(f"Audited {wheels[0].name} and {sdists[0].name}: release modules, licenses, metadata, exclusions OK")


if __name__ == "__main__":
    audit(sys.argv[1] if len(sys.argv) > 1 else ROOT / "dist")
