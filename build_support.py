"""Release-only module selection, even when building in a research checkout."""
from pathlib import Path
import shutil
from setuptools.command.build_py import build_py


class ReleaseBuildPy(build_py):
    # Adding a native public module is an explicit release decision.
    _native = {
        "WaLSAlib": {"__init__", "SoFT", "SoFTg"},
        "WaLSAlib.linefit": {"__init__", "WaLSA_LineFit"},
        "WaLSAlib._vendor": {"__init__"},
        "WaLSAlib._vendor.soft": {"__init__", "soft", "extras"},
        "WaLSAlib._vendor.softg": {
            "__init__", "dataset", "granule_detection_lib", "label_tool",
            "losses", "model", "predict", "train", "utils",
        },
    }

    def find_package_modules(self, package, package_dir):
        return [item for item in super().find_package_modules(package, package_dir)
                if item[1] in self._native.get(package, set())]

    def run(self):
        # A previous build may contain an experimental module that is no longer
        # selected. Start the package staging directory afresh on every build.
        staged = Path(self.build_lib) / "WaLSAlib"
        if staged.resolve() == Path("src/WaLSAlib").resolve():
            raise ValueError("Build output must not replace the source package")
        if staged.exists():
            shutil.rmtree(staged)
        super().run()
