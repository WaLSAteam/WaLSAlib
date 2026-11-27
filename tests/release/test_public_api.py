"""Installed public API and small, deterministic scientific contract checks."""
import importlib.util
import os
from pathlib import Path
import subprocess
import sys
import types

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[2]


def test_public_imports_and_legacy_alias():
    from WaLSAlib import LineFit, linefit, SoFT, SoFTg
    from WaLSAlib.linefit.WaLSA_LineFit import WaLSA_LineFit
    assert LineFit is linefit is WaLSA_LineFit
    assert callable(SoFT.detection)
    assert callable(SoFTg.segment_solar_granules)


def test_imports_do_not_load_optional_stacks():
    code = """
import sys
from WaLSAlib import LineFit, SoFT, SoFTg
assert not {'torch', 'lightning', 'cv2', 'sunpy', 'WaLSAtools'} & sys.modules.keys()
assert not {'soft', 'utils', 'model', 'dataset'} & sys.modules.keys()
"""
    subprocess.run([sys.executable, "-c", code], cwd=ROOT.parent, check=True)


def magnetic_image():
    image = np.zeros((64, 64), dtype=float)
    image[8:18, 10:20] = 50
    image[38:50, 40:52] = -70
    image[2, 2] = 50  # rejected by min_size
    return image


def test_soft_detection_and_size_filter():
    from WaLSAlib import SoFT
    image = magnetic_image()
    detected = SoFT.detection(image, l_thr=20, sign="both")
    assert detected.shape == image.shape
    assert detected[12, 14] > 0 and detected[43, 44] < 0
    filtered = SoFT.identification(detected, min_size=4)
    assert filtered[2, 2] == 0
    assert filtered[12, 14] > 0 and filtered[43, 44] < 0
    assert np.all(filtered[image == 0] == 0)


@pytest.mark.parametrize("separation", [False, True])
def test_soft_matches_original_snapshot(separation):
    from WaLSAlib import SoFT
    path = ROOT / "vendor/SoFT/soft/soft.py"
    spec = importlib.util.spec_from_file_location("original_soft_core", path)
    upstream = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(upstream)
    image = magnetic_image()
    actual = SoFT.detection(image, 20, sign="both", separation=separation)
    expected = upstream.detection(image, 20, sign="both", separation=separation)
    np.testing.assert_array_equal(actual, expected)
    np.testing.assert_array_equal(SoFT.identification(actual, 4), upstream.identification(expected, 4))


def granulation_image():
    y, x = np.mgrid[:64, :64]
    return 1 + .3 * np.cos(x * np.pi / 12) * np.cos(y * np.pi / 12)


def test_softg_segmentation_and_original_parity():
    from WaLSAlib import SoFT, SoFTg
    source = (ROOT / "vendor/SoFTg/granule_detection_lib.py").read_text()
    # Bind the original module's only global-package import to the reviewed SoFT
    # core. No numerical code is changed in this independent reference load.
    source = source.replace("import soft.soft as st\n", "")
    upstream = types.ModuleType("original_softg_segmentation")
    upstream.st = SoFT.core
    exec(compile(source, "original_softg_segmentation", "exec"), upstream.__dict__)
    image = granulation_image()
    actual = SoFTg.segment_solar_granules(image, plate_scale=50)
    expected = upstream.segment_solar_granules(image, plate_scale=50)
    assert len(actual) == 3
    for output, reference in zip(actual, expected):
        assert output.shape == image.shape
        assert np.issubdtype(output.dtype, np.integer)
        np.testing.assert_array_equal(output, reference)
    labels, seeds, skeleton = actual
    assert labels.max() > 0 and seeds.max() > 0
    assert set(np.unique(skeleton)) <= {0, 1}
    assert not np.any(labels < 0)


def test_ml_forward_loss_and_packaged_script_imports():
    if os.environ.get("WALSALIB_REQUIRE_ML"):
        import torch
        import torchvision
        import lightning
        import cv2
    else:
        torch = pytest.importorskip("torch")
        pytest.importorskip("torchvision")
        pytest.importorskip("lightning")
        pytest.importorskip("cv2")
    from WaLSAlib import SoFTg
    torch.set_num_threads(1)
    net = SoFTg.model.UNetPlusPlus(n_channels=1, n_classes=2, scale=16)
    inputs = torch.linspace(0, 1, 32 * 32).reshape(1, 1, 32, 32)
    logits = net(inputs)
    assert logits.shape == (1, 2, 32, 32)
    loss = SoFTg.losses.SoftIoULoss()(logits, torch.zeros((1, 32, 32), dtype=torch.long))
    assert torch.isfinite(loss)
    loss.backward()
    assert net.final.weight.grad is not None
    assert callable(SoFTg.predict.load_trained_model)
    assert callable(SoFTg.train.main)
    for module in ("train", "predict"):
        subprocess.run([sys.executable, "-m", f"WaLSAlib._vendor.softg.{module}", "--help"], cwd=ROOT.parent, check=True, capture_output=True)


def test_missing_ml_dependency_has_install_hint():
    code = """
import importlib.abc
import sys
class NoTorch(importlib.abc.MetaPathFinder):
    def find_spec(self, fullname, path=None, target=None):
        if fullname == 'torch':
            raise ModuleNotFoundError('blocked torch', name='torch')
sys.meta_path.insert(0, NoTorch())
from WaLSAlib import SoFTg
try:
    SoFTg.model
except ModuleNotFoundError as exc:
    assert 'WaLSAlib[softg]' in str(exc)
else:
    raise AssertionError('Expected optional dependency error')
"""
    subprocess.run([sys.executable, "-c", code], cwd=ROOT.parent, check=True)


def test_soft_fits_tracking_pipeline_in_worker_processes(tmp_path):
    """Exercise installed module names through real multiprocessing and FITS I/O."""
    import json
    from astropy.io import fits
    tmp_path = tmp_path / "work with spaces"
    tmp_path.mkdir()
    keep = tmp_path / "temp_keep"
    keep.mkdir()
    (keep / "keep.txt").write_text("preserve unrelated files")
    data = tmp_path / "00-data"
    data.mkdir()
    for frame in range(4):
        image = np.zeros((32, 32), dtype=np.float32)
        image[10:16, 8 + frame:14 + frame] = 50
        fits.writeto(data / f"frame{frame:03d}.fits", image)
    runner = tmp_path / "run_tracking.py"
    runner.write_text(
        "from WaLSAlib import SoFT\n"
        "if __name__ == '__main__':\n"
        f"    SoFT.track_all({str(tmp_path) + os.sep!r}, cores=1, l_thr=20, "
        "min_size=4, dx=50, dt=1, sign='pos', separation=False)\n"
    )
    completed = subprocess.run([sys.executable, str(runner)], cwd=tmp_path, capture_output=True, text=True, timeout=90)
    assert completed.returncode == 0, completed.stdout + completed.stderr
    table = json.loads((tmp_path / "dataframe.json").read_text())
    assert list(table["Lifetime"].values()) == [4]
    assert len(list((tmp_path / "03-assoc").glob("*.fits"))) == 4
    assert list(table["Area"].values()) == [[36, 36, 36, 36]]
    np.testing.assert_allclose(list(table["Vx"].values()), [[50, 50, 50, 50]])
    np.testing.assert_allclose(list(table["Vy"].values()), [[0, 0, 0, 0]])
    np.testing.assert_allclose(list(table["Flux"].values()), [[50, 50, 50, 50]])
    assert (keep / "keep.txt").read_text() == "preserve unrelated files"
    assert not (tmp_path / "temp0").exists()
