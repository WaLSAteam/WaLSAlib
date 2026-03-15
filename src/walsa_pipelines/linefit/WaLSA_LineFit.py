"""
LineFit: Adaptive multi-line centre extraction for dense solar spectra
======================================================================

This module provides an adaptive multi-line fitting routine (WaLSA_LineFit) for
robust extraction of:

  - line-centre wavelength  λ̂_i
  - LOS Doppler velocity    v̂_i  (computed relative to rest wavelength λ0_i)
  - line-core intensity     I(λ̂_i)  (sampled from the *original* spectrum)

The implementation is tailored to dense spectral windows (tens to hundreds of lines)
where profiles can evolve rapidly and become asymmetric, blended, or intermittently
multi-lobed (split/self-reversal). The core design goals are:

  1) Stable coarse seeding near expected rest wavelengths (optionally close-pair safe)
  2) Conservative per-line adaptive windowing (small-by-default, expand only if justified)
  3) Parametric centre refinement by bounded non-linear least squares
  4) Special handling for declared reversal candidates using an inverted-space
     split-core detector, with a “weak-but-real” acceptance path for shallow reversals.

IMPORTANT: 
-----------------------------------------
The numerical behaviour of the fitter depends on:
  - per-line window sizes and min/max bounds
  - reversal detector thresholds
  - smoothing window and normalisation
  - choice of model family (voigt vs asymmetric_voigt)

Debug prints (e.g. REVDBG/USED_HALFWIN) exist for development and
are not required for scientific use; they can be silenced with `silent=True`.

Data conventions
----------------
`spectra` may be:
  - 1D array (nw,) : a single spectrum
  - 2D array (ns, nw): a stack of spectra (e.g. multiple spatial pixels, or time frames)

Dependencies
------------
numpy, scipy (curve_fit, find_peaks, uniform_filter1d, wofz), astropy (optional FITS writing), csv.

"""

import numpy as np  # type: ignore
from scipy.optimize import curve_fit  # type: ignore
from scipy.signal import find_peaks  # type: ignore
from scipy.special import wofz  # type: ignore
from scipy.ndimage import uniform_filter1d  # type: ignore
from astropy.io import fits  # type: ignore
import csv  # type: ignore

# 🔹 Constants
c_speed = 299792.458  # Speed of light in km/s


# 🔹 Voigt function
def voigt(x, A0, A, x0, sigma, gamma):
    """
    Symmetric Voigt profile used in *inverted spectrum space*.

    Parameters
    ----------
    x : array
        Wavelength axis (nm).
    A0 : float
        Baseline offset (in inverted space).
    A : float
        Amplitude (in inverted space; positive values correspond to absorption cores).
    x0 : float
        Line centre (nm).
    sigma : float
        Gaussian width parameter (nm-like scale in the current setup).
    gamma : float
        Lorentzian width parameter.

    Notes
    -----
    - sigma and gamma are lower-bounded to keep the fit numerically stable.
    - The function returns a *peak* at x0 in inverted-space (i.e. absorption core).
    """
    sigma = max(float(sigma), 0.001)  # Gaussian width
    gamma = max(float(gamma), 0.001)  # Lorentzian width
    z = ((x - x0) + 1j * gamma) / (sigma * np.sqrt(2.0))
    return A0 + A * np.real(wofz(z)) / sigma / np.sqrt(2.0 * np.pi)


# 🔹 Asymmetric Voigt (wing asymmetry)
def asymmetric_voigt(x, A0, A, x0, sigma_blue, gamma_blue, alpha, red_asymmetry=True):
    """
    Asymmetric Voigt profile with controllable wing broadening.

    The asymmetry is implemented by using different (sigma, gamma) values on either
    side of the centre x0. If `red_asymmetry=True`, the red wing is broadened by
    a factor (1+alpha); otherwise the blue wing is broadened.

    Parameters
    ----------
    x : array
        Wavelength axis (nm).
    A0, A, x0 : floats
        As in `voigt`.
    sigma_blue, gamma_blue : float
        Base (blue-side) widths.
    alpha : float
        Wing asymmetry factor (allowed to be negative or positive; bounded by caller).
    red_asymmetry : bool
        If True, apply asymmetry to the red wing; else to the blue wing.

    Notes
    -----
    - Width parameters are clipped to [0.002, 0.05] for stability.
    """
    sigma_blue = float(sigma_blue)
    gamma_blue = float(gamma_blue)
    alpha = float(alpha)

    if red_asymmetry:
        sigma_red = sigma_blue * (1.0 + alpha)
        gamma_red = gamma_blue * (1.0 + alpha)
    else:
        sigma_red = sigma_blue
        gamma_red = gamma_blue
        sigma_blue = sigma_blue * (1.0 + alpha)
        gamma_blue = gamma_blue * (1.0 + alpha)

    sigma = np.where(x < x0, sigma_blue, sigma_red)
    gamma = np.where(x < x0, gamma_blue, gamma_red)

    sigma = np.clip(sigma, 0.002, 0.05)
    gamma = np.clip(gamma, 0.002, 0.05)

    z = ((x - x0) + 1j * gamma) / (sigma * np.sqrt(2.0))
    return A0 + A * np.real(wofz(z)) / sigma / np.sqrt(2.0 * np.pi)


# =============================================================================
# Close-pair-safe coarse localisation (ownership split by valley)
# =============================================================================
def coarse_centers_closepair_safe(
    wl, I, lambda0_nm,
    search_half_win=120,
    smooth=3,
    close_pairs=None,
    peak_distance=2,
    height_rel=0.05,
    valley_pad_bins=0,
):
    """
    Coarse localisation of absorption-core peaks in inverted space.

    Starting from expected rest wavelengths `lambda0_nm`, this routine searches for
    prominent absorption-core peaks (peaks in inv = 1 - smoothed_norm_spectrum).
    For declared close pairs, it enforces "ownership" by splitting the shared search
    interval at the inter-line valley (minimum of inv between the pair), preventing
    swapping between the two members.

    Returns
    -------
    centers_nm : (n_lines,) float
        Coarse centre estimates in nm.
    idx : (n_lines,) int
        Corresponding indices into `wl`.
    """
    wl = np.asarray(wl, float)
    I  = np.asarray(I, float)
    lam0 = np.asarray(lambda0_nm, float)
    nlines = lam0.size
    n = wl.size

    if close_pairs is None:
        close_pairs = []

    # sanitize + normalize
    y = np.asarray(I, float)
    if np.any(~np.isfinite(y)):
        m = np.nanmedian(y)
        y = np.where(np.isfinite(y), y, m if np.isfinite(m) else 0.0)

    y01 = (y - np.min(y)) / (np.max(y) - np.min(y) + 1e-12)
    sm = uniform_filter1d(y01, size=int(smooth))
    inv = 1.0 - sm  # peaks correspond to absorption cores (in inverted space)

    def nearest_index(w, x):
        return int(np.argmin(np.abs(np.asarray(w, float) - float(x))))

    def pick_peak_in_interval(i, lo, hi):
        lo = int(max(0, lo))
        hi = int(min(n - 1, hi))
        if hi <= lo + 4:
            return int(np.clip(nearest_index(wl, lam0[i]), lo, hi))

        seg = inv[lo:hi+1]
        if seg.size < 10 or (not np.isfinite(seg).all()):
            return int(lo + np.nanargmax(seg))

        thr = float(height_rel * np.nanmax(seg))
        pks, props = find_peaks(seg, distance=int(peak_distance), height=thr)
        if pks.size == 0:
            return int(lo + np.nanargmax(seg))

        pks_abs = lo + pks
        dist = np.abs(wl[pks_abs] - lam0[i])
        height = inv[pks_abs]
        score = dist + 0.10 * (1.0 - height / (np.nanmax(height) + 1e-12))
        return int(pks_abs[np.argmin(score)])

    # independent search
    idx = np.full(nlines, -1, dtype=int)
    for i in range(nlines):
        idx0 = nearest_index(wl, lam0[i])
        lo = max(0, idx0 - int(search_half_win))
        hi = min(n - 1, idx0 + int(search_half_win))
        idx[i] = pick_peak_in_interval(i, lo, hi)

    # enforce ownership for declared close pairs
    for (i, j) in close_pairs:
        ii = nearest_index(wl, lam0[i])
        jj = nearest_index(wl, lam0[j])
        a = min(ii, jj)
        b = max(ii, jj)

        if b <= a + 2:
            valley = (a + b) // 2
        else:
            seg = inv[a:b+1]
            valley = int(a + np.argmin(seg))

        valley_L = max(0, valley - int(valley_pad_bins))
        valley_R = min(n - 1, valley + int(valley_pad_bins))

        lo_pair = max(0, a - int(search_half_win))
        hi_pair = min(n - 1, b + int(search_half_win))

        idx[i] = pick_peak_in_interval(i, lo_pair, valley_L)
        idx[j] = pick_peak_in_interval(j, valley_R, hi_pair)

    centers_nm = wl[idx].astype(float)
    return centers_nm, idx


# =============================================================================
# Peak candidate helpers
# =============================================================================
def _find_candidate_peaks(inv, lo, hi, peak_distance=2, height_rel=0.05):
    """
    Find candidate absorption-core peaks in an interval of inverted spectrum `inv`.
    Returns absolute indices (into inv/wavelength array) and the raw `find_peaks` props.
    """
    seg = inv[lo:hi+1]
    if seg.size < 10 or (not np.isfinite(seg).all()):
        return np.array([], dtype=int), {}

    thr = float(height_rel * np.nanmax(seg))
    pks, props = find_peaks(seg, distance=int(peak_distance), height=thr)
    if pks.size == 0:
        return np.array([], dtype=int), {}

    pks_abs = lo + pks
    return pks_abs.astype(int), props


def _score_candidates(wavelengths, inv, cand_idx, x_prior, w_dist_nm=1.0, w_height=0.10):
    """
    Score candidate peaks by (distance to prior centre) + a small penalty for low peak height.
    Lower is better.
    """
    x = wavelengths[cand_idx]
    dist = np.abs(x - float(x_prior))

    h = inv[cand_idx]
    hmax = np.nanmax(h) if np.isfinite(np.nanmax(h)) and np.nanmax(h) > 0 else 1.0
    height_term = 1.0 - (h / hmax)

    score = w_dist_nm * dist + w_height * height_term
    return score


def local_continuum_estimate(y):
    """
    Robust local continuum estimate from the top 20% of samples in a window.

    Notes
    -----
    In deep-line regions or crowded windows, this provides a stable local Ic without
    requiring the spectrum to return fully to the true continuum.
    """
    y = np.asarray(y, float)
    if y.size == 0 or np.all(~np.isfinite(y)):
        return np.nan
    y = y[np.isfinite(y)]
    if y.size == 0:
        return np.nan
    top = np.sort(y)[int(0.8 * len(y)):]
    return np.median(top) if top.size else np.max(y)


# =============================================================================
# Conservative adaptive window helpers (small-by-default, expand only if justified)
# =============================================================================
def _local_depression(wavelengths, original_spectrum, idx_center, half_search=120):
    """
    Compute a local absorption depression profile around idx_center:
        D = 1 - I/Ic   (clipped to >= 0)
    using a robust local continuum Ic from `local_continuum_estimate`.

    Returns
    -------
    D : array or None
    Ic : float or None
    lo, hi : int
        window bounds (inclusive indices in original_spectrum)
    """
    n = len(wavelengths)
    lo = max(0, idx_center - int(half_search))
    hi = min(n - 1, idx_center + int(half_search))
    y = original_spectrum[lo:hi+1]

    if y.size < 10 or (not np.isfinite(y).all()):
        return None, None, lo, hi

    Ic = local_continuum_estimate(y)
    if (not np.isfinite(Ic)) or Ic <= 0:
        return None, None, lo, hi

    D = 1.0 - (y / Ic)
    D = np.clip(D, 0.0, None)
    return D, Ic, lo, hi


def _halfwidth_from_depression(D, k0, frac=0.35, margin_bins=4):
    """
    Estimate a half-width (in index bins) from a depression profile D:

    Define a threshold thr = frac * max(D), then walk left/right from k0 until D <= thr.
    The returned half-width is the max distance to either side plus `margin_bins`.

    Returns None if D is invalid or has no measurable depth.
    """
    Dmax = float(np.nanmax(D)) if D is not None else np.nan
    if (not np.isfinite(Dmax)) or Dmax <= 0:
        return None

    thr = float(frac) * Dmax

    L = int(k0)
    while L > 0 and D[L] > thr:
        L -= 1
    R = int(k0)
    while R < D.size - 1 and D[R] > thr:
        R += 1

    half = int(max(k0 - L, R - k0) + int(margin_bins))
    return int(max(2, half))


def _secondary_peak_guard(inv, lo, hi, primary_idx, peak_distance=2, height_rel=0.05, secondary_peak_rel=0.35):
    """
    Guard against unsafe window expansion: if a second peak is strong relative to the primary
    and is not close to the primary index, expanding the window risks swallowing a blend or
    a competing structure.
    """
    # If there is a second strong peak away from the core → expanding is dangerous (blend / extra lobe / wing issue).
    seg = inv[lo:hi+1]
    if seg.size < 10 or (not np.isfinite(seg).all()):
        return False

    thr = float(height_rel * np.nanmax(seg))
    pks, props = find_peaks(seg, distance=int(peak_distance), height=thr)
    if pks.size < 2:
        return False

    pabs = lo + pks
    h = inv[pabs]
    if h.size < 2:
        return False

    order = np.argsort(h)[::-1]
    p1 = int(pabs[order[0]])
    p2 = int(pabs[order[1]])
    h1 = float(inv[p1])
    h2 = float(inv[p2])

    if h1 > 0 and (h2 / h1) >= float(secondary_peak_rel):
        if abs(p2 - int(primary_idx)) > 3:
            return True
    return False


def _is_reversal_present(
    inv,
    idx_anchor,
    search_half=160,
    peak_distance=3,
    height_rel=0.03,
    prominence_rel=0.02,
    min_sep_bins=10,
    valley_depth_rel=0.15,
    peak_min_rel=0.18,
    peak_balance_rel=0.25,
    require_straddle=False,
    anchor_pad_bins=10,
    max_lobe_dist_bins=60,
    debug=False,
    metrics_out=None,

    # ------------------------------------------------------------------
    # NEW: weak-but-real reversal acceptance path
    # ------------------------------------------------------------------
    weak_valley_depth_rel=0.05,
    weak_peak_min_rel=0.14,
    weak_balance_rel=0.20,
    weak_sep_min_bins=12,
):
    """
    Detect split-core / self-reversal morphology in inverted absorption-core space.

    We work in inv-space (peaks correspond to absorption lobes) and search for a
    *pair of nearby lobes* around idx_anchor separated by a valley (central brightening
    in intensity space). The key metric is:
        depth_metric = (0.5*(hL+hR) - v) / (0.5*(hL+hR))
    where hL/hR are the inv-peak heights and v is the inv value at the valley.

    Acceptance modes
    ----------------
    Strong acceptance:
        depth_metric >= valley_depth_rel
        and lobes pass minimum strength and balance constraints.

    Weak acceptance (to catch shallow-but-clear reversals):
        depth_metric >= weak_valley_depth_rel
        and sep >= weak_sep_min_bins
        and lobes pass weaker strength/balance constraints.

    Locality constraints
    --------------------
    - Both lobes must be within `max_lobe_dist_bins` of the anchor to avoid selecting
      distant wings or unrelated blends.
    - Optional anchor rules: require_straddle / anchor_pad_bins.

    Parameters
    ----------
    metrics_out : dict or None
        If provided, the best candidate metrics are written into this dict. This is
        used downstream for the "huge reversal override" logic and for debugging prints.

    Returns
    -------
    bool
        True if a split-core/reversal is detected, else False.
    """
    inv = np.asarray(inv, float)
    n = inv.size

    lo = max(0, int(idx_anchor) - int(search_half))
    hi = min(n - 1, int(idx_anchor) + int(search_half))
    seg = inv[lo:hi + 1]

    if seg.size < 30 or (not np.isfinite(seg).all()):
        return False

    smax = float(np.nanmax(seg))
    if (not np.isfinite(smax)) or smax <= 0:
        return False

    thr_h = float(height_rel) * smax
    thr_p = float(prominence_rel) * smax

    pks, props = find_peaks(
        seg,
        distance=int(peak_distance),
        height=thr_h,
        prominence=thr_p,
    )
    if pks.size < 2:
        return False

    pabs = (lo + pks).astype(int)
    heights = inv[pabs]
    prom = props.get("prominences", np.zeros_like(pks, float))

    score = 0.7 * prom + 0.3 * heights
    order = np.argsort(score)[::-1]
    pabs = pabs[order]

    # Track best candidate even if it fails the threshold
    best = {
        "depth": -np.inf,
        "L": -1, "R": -1, "valley": -1,
        "hL": np.nan, "hR": np.nan, "v": np.nan,
        "smax": smax,
        "sep": -1,
        "bal": np.nan,
        "accepted_mode": "none",
    }

    for a in range(min(6, pabs.size)):
        for b in range(a + 1, min(6, pabs.size)):
            p1 = int(pabs[a])
            p2 = int(pabs[b])

            if abs(p2 - p1) < int(min_sep_bins):
                continue

            left, right = (p1, p2) if p1 < p2 else (p2, p1)

            # Anchor condition
            if require_straddle:
                ok_anchor = (left < int(idx_anchor) < right)
            else:
                ok_anchor = (
                    (left <= int(idx_anchor) <= right)
                    or (abs(int(idx_anchor) - left) <= int(anchor_pad_bins))
                    or (abs(int(idx_anchor) - right) <= int(anchor_pad_bins))
                )
            if not ok_anchor:
                continue

            # Both lobes must stay local to the anchor
            if max(abs(left - int(idx_anchor)), abs(right - int(idx_anchor))) > int(max_lobe_dist_bins):
                if debug:
                    print(f"[REVDBG2] reject: far lobes  L={left} R={right} idx={int(idx_anchor)}")
                continue

            hL = float(inv[left])
            hR = float(inv[right])

            # Strength gate
            if (min(hL, hR) < float(peak_min_rel) * smax):
                if debug:
                    print(f"[REVDBG2] reject: weak lobe  hL={hL:.5g} hR={hR:.5g}  smax={smax:.5g}")
                continue

            # Balance gate
            bal = min(hL, hR) / max(hL, hR)
            if bal < float(peak_balance_rel):
                if debug:
                    print(f"[REVDBG2] reject: imbalance bal={bal:.3f}  hL={hL:.5g} hR={hR:.5g}")
                continue

            seg_lr = inv[left:right + 1]
            if seg_lr.size < 5 or (not np.isfinite(seg_lr).all()):
                continue

            valley_rel = int(np.nanargmin(seg_lr))
            valley_idx = left + valley_rel
            valley_val = float(inv[valley_idx])

            hmean = 0.5 * (hL + hR)
            if hmean <= 0:
                continue

            depth_metric = (hmean - valley_val) / hmean
            sep = int(right - left)

            # Track best candidate even if it doesn't pass
            if depth_metric > best["depth"]:
                best.update({
                    "depth": depth_metric,
                    "L": left, "R": right, "valley": valley_idx,
                    "hL": hL, "hR": hR, "v": valley_val,
                    "sep": sep,
                    "bal": bal,
                    "accepted_mode": "candidate",
                })

            if debug:
                print(
                    f"[REVDBG2] idx={int(idx_anchor)} L={left} R={right} valley={valley_idx} "
                    f"sep={sep} depth={depth_metric:.3f} (thr={float(valley_depth_rel):.3f}) "
                    f"hL={hL:.5g} hR={hR:.5g} v={valley_val:.5g} bal={bal:.3f}"
                )

            # Strong acceptance
            if depth_metric >= float(valley_depth_rel):
                best.update({"accepted_mode": "strong"})
                if metrics_out is not None and isinstance(metrics_out, dict):
                    metrics_out.clear()
                    metrics_out.update(best)
                if debug:
                    print(f"[REVDBG2] ACCEPT_STRONG  idx={int(idx_anchor)} L={left} R={right}")
                return True

            # Weak acceptance (shallow-but-clear reversals)
            weak_reversal = (
                (depth_metric >= float(weak_valley_depth_rel)) and
                (sep >= int(weak_sep_min_bins)) and
                (min(hL, hR) >= float(weak_peak_min_rel) * smax) and
                (bal >= float(weak_balance_rel))
            )

            if weak_reversal:
                best.update({"accepted_mode": "weak"})
                if metrics_out is not None and isinstance(metrics_out, dict):
                    metrics_out.clear()
                    metrics_out.update(best)
                if debug:
                    print(
                        f"[REVDBG2] ACCEPT_WEAK  idx={int(idx_anchor)} L={left} R={right} "
                        f"depth={depth_metric:.3f} sep={sep} bal={bal:.3f}"
                    )
                return True

    if metrics_out is not None and isinstance(metrics_out, dict):
        metrics_out.clear()
        metrics_out.update(best)

    return False


def _reversal_window_halfwidth_bins(
    wavelengths,
    inv,
    idx_anchor,
    base_half=20,
    search_half=160,
    peak_distance=3,
    height_rel=0.03,
    prominence_rel=0.02,
    margin_bins=10,
    max_half_cap=80,
):
    """
    Estimate a suitable half-window for potential reversal/multi-lobe lines.

    Uses the spread between the extreme detected peaks within a local search window,
    plus a margin. This provides a conservative bound for lines whose core may split.
    """
    n = len(wavelengths)
    lo = max(0, int(idx_anchor) - int(search_half))
    hi = min(n - 1, int(idx_anchor) + int(search_half))

    seg = np.asarray(inv[lo:hi+1], float)
    if seg.size < 12 or (not np.isfinite(seg).all()):
        return int(base_half)

    seg_max = float(np.nanmax(seg))
    if not np.isfinite(seg_max) or seg_max <= 0:
        return int(base_half)

    thr_h = float(height_rel * seg_max)
    thr_p = float(prominence_rel * seg_max)

    pks, props = find_peaks(
        seg,
        distance=int(peak_distance),
        height=thr_h,
        prominence=thr_p,
    )
    if pks.size < 2:
        return int(base_half)

    pks_abs = lo + pks
    left = int(np.min(pks_abs))
    right = int(np.max(pks_abs))

    half = int(max(idx_anchor - left, right - idx_anchor) + int(margin_bins))
    half = max(int(base_half), half)
    half = min(int(max_half_cap), half)
    return int(half)


def _local_emission(wavelengths, original_spectrum, idx_center, half_search=160):
    """
    Compute a local *emission* proxy profile around idx_center:
        E = I/Ic - 1   (clipped to >= 0)
    based on a robust local continuum estimate Ic.

    This is used to decide whether window expansion is warranted due to local
    brightening structures (for non-reversal lines) in crowded windows.
    """
    n = len(wavelengths)
    lo = max(0, idx_center - int(half_search))
    hi = min(n - 1, idx_center + int(half_search))
    y = original_spectrum[lo:hi+1]

    if y.size < 10 or (not np.isfinite(y).all()):
        return None, None, lo, hi

    Ic = local_continuum_estimate(y)
    if (not np.isfinite(Ic)) or Ic <= 0:
        return None, None, lo, hi

    E = (y / Ic) - 1.0
    E = np.clip(E, 0.0, None)
    return E, Ic, lo, hi


def _halfwidth_from_emission(E, k0, frac=0.20, margin_bins=6, min_width_bins=6):
    """
    Estimate an emission-based half-width from E (>=0) using a fractional threshold
    of the local peak emission.
    """
    if E is None:
        return None

    Emax = float(np.nanmax(E))
    if (not np.isfinite(Emax)) or Emax <= 0:
        return None

    thr = float(frac) * Emax

    L = int(k0)
    while L > 0 and E[L] > thr:
        L -= 1
    R = int(k0)
    while R < E.size - 1 and E[R] > thr:
        R += 1

    width = int(R - L)
    if width < int(min_width_bins):
        return None

    half = int(max(k0 - L, R - k0) + int(margin_bins))
    return int(max(2, half))


def _is_emission_present(
    wavelengths,
    original_spectrum,
    idx_anchor,
    half_search=160,
    emis_rel=0.010,
    emis_abs=0.003,
):
    """
    Boolean check for a detectable local emission bump around idx_anchor.

    The bump threshold is max(emis_abs, emis_rel * Ic) where Ic is a robust local
    continuum estimate from within the search window.
    """
    n = len(wavelengths)
    lo = max(0, int(idx_anchor) - int(half_search))
    hi = min(n - 1, int(idx_anchor) + int(half_search))
    y = np.asarray(original_spectrum[lo:hi+1], float)

    if y.size < 20 or (not np.isfinite(y).all()):
        return False

    Ic = local_continuum_estimate(y)
    if (not np.isfinite(Ic)) or Ic <= 0:
        return False

    bump = float(np.nanmax(y) - Ic)
    thr = max(float(emis_abs), float(emis_rel) * float(Ic))
    return (np.isfinite(bump) and bump > thr)


# =============================================================================
# Main WaLSA fitter
# =============================================================================
def WaLSA_LineFit(
    spectra,
    fit_func=voigt,
    num_iterations=10,
    save_fitted_spectra=False,
    date_identifier=None,
    refwavelength=None,
    smoothing_window=3,
    wavelengths=None,
    element=None,
    pixel=None,
    window_size=4,
    plot_result=False,
    alpha_init=0.3,
    red_asymmetry=True,
    return_vlos_intensity=True,
    write_csv=False,
    adaptive_window=True,
    per_line_window_size=None,
    init_centres_nm=None,
    use_coarse_centers: bool = False,
    close_pairs: list[tuple[int, int]] | None = None,
    reversal_lines: list[int] | None = None,
    coarse_search_half_win: int = 120,
    coarse_smooth: int = 3,
    coarse_peak_distance: int = 2,
    coarse_height_rel: float = 0.05,
    coarse_valley_pad_bins: int = 0,
    max_center_shift_nm: float = 0.015,
    reversal_window_floor: int = 20,
    silent: bool = False,
    DEBUG_PRINT = False,

    reversal_valley_depth_rel: float = 0.10,
    reversal_peak_min_rel: float = 0.16,
    reversal_peak_balance_rel: float = 0.20,
    reversal_weak_valley_depth_rel: float = 0.05,
    reversal_weak_peak_min_rel: float = 0.12,
    reversal_weak_balance_rel: float = 0.15,
    reversal_weak_sep_min_bins: int = 10,

    # -------------------------------------------------------------------------
    # Conservative adaptive window controls (small-by-default; expand only if justified)
    # -------------------------------------------------------------------------
    per_line_min_halfwin=None,
    per_line_max_halfwin=None,
    allow_shrink_below_base: bool = False,
    max_expand_bins: int = 6,
    max_expand_factor: float = 1.30,
    depth_min_for_expand: float = 0.03,
    frac_for_width: float = 0.35,
    reject_expand_if_secondary_peak: bool = True,
    secondary_peak_rel: float = 0.35,
    emission_rel_min: float = 0.03,
    emission_frac_for_width: float = 0.20,
    emission_margin_bins: int = 8,
    emission_min_width_bins: int = 8,
    emission_expand_bins: int = 6,
    emission_expand_factor: float = 1.60,

    # Huge-reversal rare-case trigger (NOT continuum-based):
    # It uses the inv-space depth metric: depth = (hmean - valley)/hmean
    huge_emission_rel_thr: float = 0.35,
    huge_emission_force_max: bool = True,
):
    """
    Adaptive multi-line line-centre fitting and time-series extraction.

    This routine fits parametric profiles to multiple spectral lines in each input
    spectrum and returns LOS velocity and line-core intensity per line.

    Inputs
    ------
    spectra : array
        Either a single spectrum (nw,) or a stack (ns, nw), where `ns` is the number
        of spectra to process (e.g. spatial pixels, time frames, or any other grouping).
    fit_func : callable
        Either `voigt` or `asymmetric_voigt` (or any compatible callable).
    refwavelength : (n_lines,) array
        Rest wavelengths λ0_i for each line, used for seeding and Doppler conversion.
    wavelengths : (nw,) array
        Wavelength axis (nm), same length as the spectral dimension.
    element : (n_lines,) array-like
        Labels for each line (strings), written to CSV if enabled.
    pixel : int | list[int] | range | None
        Select which rows in `spectra` to process. If None, process all rows.

    Core algorithm (per spectrum)
    -----------------------------
    1) Normalise to [0,1], smooth with a boxcar (`smoothing_window`), then invert:
           inv = 1 - smooth(norm(I))
       so absorption cores appear as peaks in `inv`.

    2) Initial seeding:
       - If `use_coarse_centers=True`, call `coarse_centers_closepair_safe` to obtain
         coarse indices and prevent swapping in declared close pairs.
       - Otherwise seed around `refwavelength` and/or local peaks.

       Additionally, for each line in iteration==0, a local search interval is used to
       select a candidate peak in inv-space, with optional close-pair “ownership split”
       and a maximum allowed centre shift (`max_center_shift_nm`) relative to λ0_i.

    3) Adaptive per-line windowing (iteration==0):
       - Each line starts from a per-line baseline half-window (from `per_line_window_size`
         or `window_size`), bounded by per-line minimum/maximum half-windows.
       - For non-reversal lines, the window can expand based on:
           (a) absorption depression width (D = 1 - I/Ic)
           (b) emission proxy width (E = I/Ic - 1) if a local brightening is present
         Expansion is capped by `max_expand_bins` and `max_expand_factor`, and can be
         rejected if a strong secondary peak suggests blends (`_secondary_peak_guard`).

       - For declared reversal candidates (`reversal_lines`), a dedicated detector
         `_is_reversal_present` is evaluated *per frame* in inv-space. If a reversal is
         detected, the half-window is pushed to the per-line maximum; otherwise the
         window may be conservatively shrunk to protect the core fit.

       - A rare-case override exists: if the best reversal candidate has a very large
         inv-space depth metric (>= `huge_emission_rel_thr`), reversal mode can be forced.

    4) Parametric fit and centre refinement:
       - Fit the chosen model within the final window using bounded `curve_fit`.
       - Update the fitted centre from the fitted x0 parameter.
       - Repeat for `num_iterations` (or stop early if both total residual and centre
         updates converge).

       Asymmetry handling:
       - If using `asymmetric_voigt`, the fit includes an asymmetry parameter `alpha`
         and the red/blue wing choice is controlled by `red_asymmetry`. This provides
         robustness when wings evolve asymmetrically.

    5) Outputs:
       - LOS velocity: v̂_i = c * (λ̂_i - λ0_i) / λ0_i
       - Line-core intensity: I(λ̂_i), sampled from the *original* (non-normalised) spectrum.

    Optional outputs / side effects
    -------------------------------
    write_csv : bool
        Writes extracted per-line parameters to:
            extracted_parameters_{date_identifier}_{pixel_identifier}.csv
    save_fitted_spectra : bool
        Writes fitted model segments (in inv-space) to:
            fitted_spectra_{date_identifier}_{pixel_identifier}.fits

    Returns
    -------
    vlos_out : (ns, n_lines) float32
        LOS velocities (km/s) for each processed spectrum and line.
    int_out : (ns, n_lines) float32
        Line-core intensity (original units) evaluated at the fitted centre.

    Notes on scientific interpretation
    ---------------------------------
    - “Reversal” here refers to a split/multi-lobed core structure (two inv peaks with a
      valley between them), not necessarily an emission feature above the continuum level.
    - Continuum estimation is local and robust; this is important for deep lines and
      crowded windows where the spectrum may not return to a true continuum within
      the fitting interval.

    """
    # ----------------------------- (NO CODE CHANGES BELOW THIS LINE) -----------------------------

    spectra = np.asarray(spectra)
    if spectra.ndim == 1:
        image_data = spectra[np.newaxis, :]
    elif spectra.ndim == 2:
        image_data = spectra
    else:
        raise ValueError("Input spectra must be 1D or 2D array.")

    if refwavelength is None or wavelengths is None or element is None:
        raise ValueError("refwavelength, wavelengths, and element must be provided.")

    wavelengths = np.asarray(wavelengths, dtype=float)
    refwavelength = np.asarray(refwavelength, dtype=float)
    element = np.asarray(element, dtype=object)

    n_samples = image_data.shape[0]
    n_lines = len(refwavelength)

    # --- DEBUG controls (edit) ---
    DEBUG_LINES = {5}
    DEBUG_ITERS = {0}
    DEBUG_PIXELS = None
    DEBUG_ONLY_IF_ALPHA_NEAR_BOUND = True
    ALPHA_BOUND_EPS = 0.02

    specific_spatial_pixel = pixel

    if n_samples == 1:
        spatial_pixel_range = [0]
        specific_spatial_pixel = 0
    else:
        if specific_spatial_pixel is None:
            spatial_pixel_range = range(n_samples)
        elif isinstance(specific_spatial_pixel, int):
            spatial_pixel_range = [specific_spatial_pixel]
        else:
            spatial_pixel_range = specific_spatial_pixel

    if per_line_window_size is not None:
        per_line_window_size = np.asarray(per_line_window_size, dtype=int)
        if per_line_window_size.shape[0] != n_lines:
            raise ValueError("per_line_window_size must have same length as refwavelength")
        if np.any(per_line_window_size < 1):
            raise ValueError("per_line_window_size entries must be >= 1")

    if per_line_min_halfwin is not None:
        per_line_min_halfwin = np.asarray(per_line_min_halfwin, dtype=int)
        if per_line_min_halfwin.shape[0] != n_lines:
            raise ValueError("per_line_min_halfwin must have same length as refwavelength")
    else:
        per_line_min_halfwin = np.ones(n_lines, dtype=int) * 4

    if per_line_max_halfwin is not None:
        per_line_max_halfwin = np.asarray(per_line_max_halfwin, dtype=int)
        if per_line_max_halfwin.shape[0] != n_lines:
            raise ValueError("per_line_max_halfwin must have same length as refwavelength")
    else:
        if per_line_window_size is not None:
            per_line_max_halfwin = per_line_window_size + 8
        else:
            per_line_max_halfwin = np.ones(n_lines, dtype=int) * (int(window_size) + 8)

    per_line_max_halfwin = np.maximum(per_line_max_halfwin, per_line_min_halfwin)

    vlos_out = np.full((n_samples, n_lines), np.nan, dtype=np.float32)
    int_out  = np.full((n_samples, n_lines), np.nan, dtype=np.float32)

    all_results = []
    all_fitted_spectra = []

    for spatial_pixel in spatial_pixel_range:
        original_spectrum = np.asarray(image_data[int(spatial_pixel), :], dtype=np.float64)

        bad = ~np.isfinite(original_spectrum)
        if bad.any():
            good_idx = np.flatnonzero(~bad)
            if good_idx.size >= 2:
                original_spectrum[bad] = np.interp(np.flatnonzero(bad), good_idx, original_spectrum[good_idx])
            elif good_idx.size == 1:
                original_spectrum[bad] = original_spectrum[good_idx[0]]
            else:
                original_spectrum[:] = 0.0

        smin = np.min(original_spectrum)
        smax = np.max(original_spectrum)
        den  = smax - smin
        if (not np.isfinite(den)) or (den <= 0.0):
            spectrum = np.zeros_like(original_spectrum, dtype=np.float64)
            spectrum[0] = 1e-6
        else:
            spectrum = (original_spectrum - smin) / den

        smoothed_spectrum = uniform_filter1d(spectrum, size=int(smoothing_window))
        inverted_spectrum = 1.0 - smoothed_spectrum

        pk, properties = find_peaks(inverted_spectrum, height=0.02, width=2)
        if pk.size == 0:
            pk = np.array([np.argmin(np.abs(wavelengths - rw)) for rw in refwavelength], dtype=int)
            fwhm = np.full(pk.shape, 3.0, dtype=float)
        else:
            fwhm = properties.get("widths", np.full(pk.shape, 6.0, dtype=float))

        if init_centres_nm is not None:
            init_centres_nm = np.asarray(init_centres_nm, dtype=float)
            if init_centres_nm.shape[0] != n_lines:
                raise ValueError("init_centres_nm must have same length as refwavelength")
            fitted_centres = init_centres_nm.copy()
        else:
            fitted_centres = np.copy(refwavelength)

        coarse_idx = None
        if use_coarse_centers:
            cpairs = close_pairs if close_pairs is not None else []
            coarse_centres_nm, coarse_idx = coarse_centers_closepair_safe(
                wl=wavelengths,
                I=original_spectrum,
                lambda0_nm=refwavelength,
                search_half_win=int(coarse_search_half_win),
                smooth=int(coarse_smooth),
                close_pairs=cpairs,
                peak_distance=int(coarse_peak_distance),
                height_rel=float(coarse_height_rel),
                valley_pad_bins=int(coarse_valley_pad_bins),
            )
            if init_centres_nm is None:
                fitted_centres = coarse_centres_nm.copy()

        prev_total_residual = np.inf
        residual_tol = 1e-6
        param_tol = 1e-4

        used_halfwin = np.full(n_lines, -1, dtype=int)
        rev_detected = np.zeros(n_lines, dtype=bool)  # gating fallback refit

        for iteration in range(int(num_iterations)):
            matched_lines = []
            final_fitted_spectra = np.full_like(inverted_spectrum, np.nan)
            total_residual = 0.0

            for i, ref_wl in enumerate(fitted_centres):
                is_rev = False
                rev_half = None

                if iteration == 0:
                    cpairs = close_pairs or []
                    reversal_set = set(reversal_lines or [])

                    x_prior = float(refwavelength[i])
                    idx_prior = int(np.argmin(np.abs(wavelengths - x_prior)))

                    if use_coarse_centers and (coarse_idx is not None):
                        idx_prior = int(coarse_idx[i])

                    search_half = 60
                    Lloc = max(0, idx_prior - int(search_half))
                    Rloc = min(len(wavelengths) - 1, idx_prior + int(search_half))

                    if cpairs:
                        for (a, b) in cpairs:
                            if i == a or i == b:
                                ia = int(np.argmin(np.abs(wavelengths - refwavelength[a])))
                                ib = int(np.argmin(np.abs(wavelengths - refwavelength[b])))
                                lo_pair = min(ia, ib)
                                hi_pair = max(ia, ib)
                                seg = inverted_spectrum[lo_pair:hi_pair+1]
                                valley = (lo_pair + hi_pair) // 2 if seg.size < 3 else (lo_pair + int(np.argmin(seg)))
                                if i == min(a, b):
                                    Rloc = min(Rloc, valley)
                                else:
                                    Lloc = max(Lloc, valley)

                    if i in reversal_set:
                        rev_half = _reversal_window_halfwidth_bins(
                            wavelengths=wavelengths,
                            inv=inverted_spectrum,
                            idx_anchor=idx_prior,
                            base_half=int(reversal_window_floor),
                            search_half=120,
                            peak_distance=2,
                            height_rel=0.08,
                            margin_bins=8,
                        )

                    cand, _ = _find_candidate_peaks(inverted_spectrum, Lloc, Rloc, peak_distance=2, height_rel=0.05)

                    if cand.size == 0:
                        peak_idx = idx_prior
                        peak_wavelength = float(wavelengths[peak_idx])
                        fwhm_peak = 6.0
                    else:
                        guard_nm = float(max_center_shift_nm)
                        cand = cand[np.abs(wavelengths[cand] - x_prior) <= guard_nm]
                        if cand.size == 0:
                            peak_idx = idx_prior
                            peak_wavelength = float(wavelengths[peak_idx])
                            fwhm_peak = 6.0
                        else:
                            score = _score_candidates(
                                wavelengths=wavelengths,
                                inv=inverted_spectrum,
                                cand_idx=cand,
                                x_prior=x_prior,
                                w_dist_nm=1.0,
                                w_height=0.08,
                            )
                            peak_idx = int(cand[np.argmin(score)])
                            peak_wavelength = float(wavelengths[peak_idx])
                            fwhm_peak = 6.0

                    fitted_centres[i] = float(peak_wavelength)

                else:
                    peak_idx = int(np.argmin(np.abs(wavelengths - float(ref_wl))))
                    peak_wavelength = float(ref_wl)
                    fwhm_peak = 6.0

                local_half_base = int(window_size)
                if per_line_window_size is not None:
                    local_half_base = int(per_line_window_size[i])

                local_half = int(local_half_base)

                min_half = int(per_line_min_halfwin[i])
                max_half = int(per_line_max_halfwin[i])

                if allow_shrink_below_base:
                    local_half = int(max(min_half, min(local_half, local_half_base)))

                is_rev = False

                if adaptive_window:
                    idx_c = int(np.argmin(np.abs(wavelengths - float(fitted_centres[i]))))

                    # reversal lines: detect split in inv-space
                    if silent:
                        debug = False
                    else:
                        debug = True
                        
                    if (reversal_lines is not None) and (i in set(reversal_lines)):
                        m = {}
                        is_rev = _is_reversal_present(
                            inv=inverted_spectrum,
                            idx_anchor=idx_c,
                            search_half=160,
                            peak_distance=3,
                            height_rel=0.03,
                            prominence_rel=0.02,
                            min_sep_bins=10,
                            valley_depth_rel=float(reversal_valley_depth_rel),
                            peak_min_rel=float(reversal_peak_min_rel),
                            peak_balance_rel=float(reversal_peak_balance_rel),
                            weak_valley_depth_rel=float(reversal_weak_valley_depth_rel),
                            weak_peak_min_rel=float(reversal_weak_peak_min_rel),
                            weak_balance_rel=float(reversal_weak_balance_rel),
                            weak_sep_min_bins=int(reversal_weak_sep_min_bins),
                            max_lobe_dist_bins=60,
                            require_straddle=False,
                            anchor_pad_bins=10,
                            debug=debug,
                            metrics_out=m,
                        )

                        # HUGE reversal override based on inv depth_metric, not continuum
                        # depth_metric = (hmean - valley)/hmean
                        if (not is_rev) and (m.get("depth", -np.inf) >= float(huge_emission_rel_thr)) and bool(huge_emission_force_max):
                            is_rev = True
                            if not silent:
                                print(
                                    f"[REVDBG_HUGE] frame={int(spatial_pixel)} line={i} "
                                    f"forcing is_rev=True due to huge reversal depth={m.get('depth', np.nan):.3f} "
                                    f"(thr={float(huge_emission_rel_thr):.3f})"
                                )

                        if (not silent) and (iteration == 0):
                            print(f"[REVDBG] frame={int(spatial_pixel)} line={i} is_rev={is_rev} idx_c={idx_c}")

                        if iteration == 0:
                            rev_detected[i] = bool(is_rev)

                        # If reversal -> push to max window
                        if is_rev:
                            local_half = max(local_half, int(max_half))
                        else:
                            # no reversal -> allow shrink (core-protect)
                            if iteration == 0:
                                D, Ic, loD, hiD = _local_depression(
                                    wavelengths=wavelengths,
                                    original_spectrum=original_spectrum,
                                    idx_center=idx_c,
                                    half_search=120,
                                )
                                if D is not None:
                                    k0 = int(idx_c - loD)
                                    hw_core = _halfwidth_from_depression(D, k0=k0, frac=0.55, margin_bins=4)
                                    if hw_core is not None:
                                        local_half = int(max(min_half, min(local_half, hw_core)))

                    else:
                        # Non-reversal lines: keep emission widening + absorption widening logic
                        if iteration == 0:
                            E, IcE, loE, hiE = _local_emission(
                                wavelengths=wavelengths,
                                original_spectrum=original_spectrum,
                                idx_center=idx_c,
                                half_search=160,
                            )
                            if E is not None:
                                Emax = float(np.nanmax(E))
                                if np.isfinite(Emax) and (Emax >= float(emission_rel_min)):
                                    k0e = int(idx_c - loE)
                                    hwE = _halfwidth_from_emission(
                                        E, k0=k0e,
                                        frac=float(emission_frac_for_width),
                                        margin_bins=int(emission_margin_bins),
                                        min_width_bins=int(emission_min_width_bins),
                                    )
                                    if hwE is not None:
                                        targetE = int(hwE)
                                        targetE = min(targetE, int(np.floor(local_half_base * float(emission_expand_factor))))
                                        targetE = min(targetE, int(local_half_base + int(emission_expand_bins)))
                                        local_half = max(local_half, int(targetE))

                        if iteration == 0:
                            D, Ic, loD, hiD = _local_depression(
                                wavelengths=wavelengths,
                                original_spectrum=original_spectrum,
                                idx_center=idx_c,
                                half_search=120,
                            )

                            if D is not None:
                                Dmax = float(np.nanmax(D))
                                if np.isfinite(Dmax) and (Dmax >= float(depth_min_for_expand)):
                                    k0 = int(idx_c - loD)
                                    hw = _halfwidth_from_depression(D, k0=k0, frac=float(frac_for_width), margin_bins=4)

                                    if hw is not None:
                                        target = int(hw)
                                        target = min(target, int(np.floor(local_half_base * float(max_expand_factor))))
                                        target = min(target, int(local_half_base + int(max_expand_bins)))

                                        if _is_emission_present(
                                            wavelengths=wavelengths,
                                            original_spectrum=original_spectrum,
                                            idx_anchor=idx_c,
                                            half_search=160,
                                            emis_rel=0.010,
                                            emis_abs=0.003,
                                        ):
                                            target = max(target, int(min(max_half, local_half_base + 12)))

                                        local_half = max(local_half, int(target))

                                        if reject_expand_if_secondary_peak:
                                            if _secondary_peak_guard(
                                                inv=inverted_spectrum,
                                                lo=max(0, idx_c - 120),
                                                hi=min(len(wavelengths) - 1, idx_c + 120),
                                                primary_idx=idx_c,
                                                peak_distance=2,
                                                height_rel=0.05,
                                                secondary_peak_rel=float(secondary_peak_rel),
                                            ):
                                                target = local_half_base

                                        local_half = max(local_half, int(target))

                local_half = int(np.clip(local_half, min_half, max_half))

                if iteration == 0:
                    used_halfwin[i] = int(local_half)

                idx_c = int(np.argmin(np.abs(wavelengths - float(fitted_centres[i]))))

                fit_range = slice(
                    max(0, idx_c - local_half),
                    min(len(wavelengths), idx_c + local_half + 1),
                )

                x_data = wavelengths[fit_range]
                y_data = inverted_spectrum[fit_range]

                peak_wavelength = float(wavelengths[idx_c])

                if (y_data.size < 5) or (not np.isfinite(y_data).all()):
                    continue

                if fit_func == voigt:
                    p0 = [np.min(y_data), np.max(y_data) - np.min(y_data), peak_wavelength, 0.008, 0.005]
                    bounds = (
                        [-np.inf, 0, peak_wavelength - 0.005, 0.002, 0.002],
                        [ np.inf, np.inf, peak_wavelength + 0.005, 0.02,  0.01],
                    )
                    try:
                        popt, _ = curve_fit(fit_func, x_data, y_data, p0=p0, bounds=bounds, maxfev=200000)
                        fitted_segment = fit_func(x_data, *popt)
                    except RuntimeError:
                        if not silent:
                            print(f"Fit failed for {element[i]} at {refwavelength[i]:.4f} nm in sample {spatial_pixel}")
                        continue
                else:
                    dw = float(np.median(np.diff(wavelengths)))
                    fwhm_nm = float(fwhm_peak) * dw
                    sigma_init = float(np.clip(fwhm_nm / 2.355, 0.002, 0.05))
                    gamma_init = float(np.clip(fwhm_nm / 2.0,   0.002, 0.05))
                    alpha0 = float(alpha_init)

                    def fit_asym(x, A0, A, x0, sigma_blue, gamma_blue, alpha):
                        return asymmetric_voigt(
                            x, A0, A, x0, sigma_blue, gamma_blue, alpha,
                            red_asymmetry=red_asymmetry
                        )

                    p0 = [np.min(y_data),
                          np.max(y_data) - np.min(y_data),
                          peak_wavelength,
                          sigma_init, gamma_init, alpha0]

                    bounds = (
                        [-np.inf, 0, peak_wavelength - 0.005, 0.002, 0.002, -0.5],
                        [ np.inf, np.inf, peak_wavelength + 0.005, 0.05,  0.05,   0.5],
                    )

                    try:
                        popt, _ = curve_fit(fit_asym, x_data, y_data, p0=p0, bounds=bounds, maxfev=200000)
                        fitted_segment = fit_asym(x_data, *popt)
                    except RuntimeError:
                        if not silent:
                            print(f"Fit failed for {element[i]} at {refwavelength[i]:.4f} nm in sample {spatial_pixel}")
                        continue

                if DEBUG_PRINT:
                    do_pixel = (DEBUG_PIXELS is None) or (int(spatial_pixel) in DEBUG_PIXELS)
                    do_line  = (int(i) in DEBUG_LINES)
                    do_iter  = (int(iteration) in DEBUG_ITERS)

                    if do_pixel and do_line and do_iter:
                        x0_fit = float(popt[2])
                        xL = float(x_data[0])
                        xR = float(x_data[-1])

                        alpha_fit = float(popt[5]) if (len(popt) >= 6) else np.nan

                        if (not DEBUG_ONLY_IF_ALPHA_NEAR_BOUND) or (
                            np.isfinite(alpha_fit) and (
                                abs(alpha_fit - 0.5) < ALPHA_BOUND_EPS or abs(alpha_fit + 0.5) < ALPHA_BOUND_EPS
                            )
                        ):
                            print(
                                f"[DEBUG] pix={int(spatial_pixel)} it={iteration} line={i} "
                                f"x0={x0_fit:.6f}  win=[{xL:.6f},{xR:.6f}]  "
                                f"A0={popt[0]:.3e} A={popt[1]:.3e} "
                                f"sigmaB={popt[3]:.4f} gammaB={popt[4]:.4f} alpha={alpha_fit:.4f}"
                            )

                def _needs_refit_asym(popt, bounds, x_data, peak_wavelength):
                    x0 = float(popt[2])
                    sigmaB = float(popt[3])
                    gammaB = float(popt[4])
                    alpha = float(popt[5])

                    if (x0 - x_data[0]) < 0.20 * (x_data[-1] - x_data[0]):
                        return True
                    if (x_data[-1] - x0) < 0.20 * (x_data[-1] - x_data[0]):
                        return True

                    alo, ahi = bounds[0][5], bounds[1][5]
                    if abs(alpha - alo) < 0.03 or abs(alpha - ahi) < 0.03:
                        return True

                    slo, shi = bounds[0][3], bounds[1][3]
                    glo, ghi = bounds[0][4], bounds[1][4]
                    if abs(sigmaB - slo) < 0.002 or abs(sigmaB - shi) < 0.002:
                        return True
                    if abs(gammaB - glo) < 0.002 or abs(gammaB - ghi) < 0.002:
                        return True

                    return False

                # prevent fallback refit when reversal was detected
                if (fit_func != voigt) and (not rev_detected[i]) and _needs_refit_asym(popt, bounds, x_data, peak_wavelength):
                    local_half_refit = max(min_half, int(np.floor(0.75 * local_half)))
                    idx_c_refit = int(np.argmin(np.abs(wavelengths - float(fitted_centres[i]))))
                    fit_range2 = slice(
                        max(0, idx_c_refit - local_half_refit),
                        min(len(wavelengths), idx_c_refit + local_half_refit + 1),
                    )
                    x2 = wavelengths[fit_range2]
                    y2 = inverted_spectrum[fit_range2]

                    bounds2 = (
                        [bounds[0][0], bounds[0][1], bounds[0][2], bounds[0][3], bounds[0][4], -0.15],
                        [bounds[1][0], bounds[1][1], bounds[1][2], bounds[1][3], bounds[1][4],  0.15],
                    )

                    p0_2 = list(p0)
                    p0_2[5] = float(np.clip(p0_2[5], -0.10, 0.10))

                    try:
                        popt2, _ = curve_fit(fit_asym, x2, y2, p0=p0_2, bounds=bounds2, maxfev=200000)
                        fitted_segment2 = fit_asym(x2, *popt2)

                        sse1 = float(np.sum((y2 - fit_asym(x2, *popt))**2))
                        sse2 = float(np.sum((y2 - fitted_segment2)**2))
                        if np.isfinite(sse2) and (sse2 < sse1):
                            popt = popt2
                            fit_range = fit_range2
                            x_data = x2
                            y_data = y2
                            fitted_segment = fitted_segment2
                    except Exception:
                        pass

                fitted_centres[i] = float(popt[2])

                velocity = c_speed * (fitted_centres[i] - refwavelength[i]) / refwavelength[i]

                peak_idx_int = int(np.argmin(np.abs(wavelengths - fitted_centres[i])))
                intensity_at_centre = float(original_spectrum[peak_idx_int])

                vlos_out[int(spatial_pixel), i] = float(velocity)
                int_out[int(spatial_pixel), i]  = float(intensity_at_centre)

                matched_lines.append([element[i], refwavelength[i], fitted_centres[i], velocity, intensity_at_centre])

                if y_data.shape != fitted_segment.shape:
                    y_data = inverted_spectrum[fit_range]

                final_fitted_spectra[fit_range] = fitted_segment

                residual = float(np.sum((y_data - fitted_segment) ** 2))
                total_residual += residual

            if (not silent) and (iteration == 0):
                print(
                    f"[USED_HALFWIN] frame={int(spatial_pixel)}  "
                    + " ".join([f"{j}:{used_halfwin[j]}" for j in range(n_lines)])
                )

            if iteration > 0:
                if (not silent) and (specific_spatial_pixel is not None):
                    print(f"\nIteration {iteration + 1} Summary:")
                    print(f"  Total residual: {total_residual:.6e}")
                    for i2, ref_wl2 in enumerate(refwavelength):
                        delta_centre = fitted_centres[i2] - prev_fitted_centres[i2]
                        print(f"  {element[i2]} | Δ centre: {delta_centre:.6e} nm | New centre: {fitted_centres[i2]:.6f} nm")
                else:
                    if not silent:
                        print(f"Pixel {int(spatial_pixel):3d} | Iteration {iteration + 1:2d} | Total residual: {total_residual:.3e}", end="\r")

                if (np.abs(total_residual - prev_total_residual) < residual_tol and
                    np.all(np.abs(fitted_centres - prev_fitted_centres) < param_tol)):
                    if specific_spatial_pixel is not None and (not silent):
                        print(f"Converged after {iteration+1} iterations (both residual and centre).")
                    break

            prev_fitted_centres = np.copy(fitted_centres)
            prev_total_residual = float(total_residual)

        all_results.extend([[int(spatial_pixel)] + line for line in matched_lines])
        if save_fitted_spectra:
            all_fitted_spectra.append(final_fitted_spectra)

    if plot_result and specific_spatial_pixel is not None:
        import matplotlib.pyplot as plt  # type: ignore
        from matplotlib.ticker import AutoMinorLocator  # type: ignore

        plt.figure(figsize=(30, 5))
        plt.plot(wavelengths, smoothed_spectrum, label="Smoothed Spectrum", color="blue", lw=2.5)

        valid_indices = np.where(~np.isnan(final_fitted_spectra))[0]
        if len(valid_indices) > 0:
            split_ranges = np.split(valid_indices, np.where(np.diff(valid_indices) > 1)[0] + 1)
            for plot_fit_range in split_ranges:
                x_fit = wavelengths[plot_fit_range]
                y_fit = 1.0 - final_fitted_spectra[plot_fit_range]
                if len(x_fit) > 0:
                    plt.plot(x_fit, y_fit, linestyle="--", color="red", lw=1.5)

        for line in matched_lines:
            plt.axvline(float(line[2]), linestyle="--", color="green", alpha=0.7)

        plt.xlim(wavelengths[0], wavelengths[-1])
        plt.tick_params(axis='both', which='major', length=10, width=1.4, pad=5)
        plt.tick_params(axis='x', which='major', pad=8)
        plt.tick_params(axis='both', which='minor', length=6, width=1.4)
        ax = plt.gca()
        ax.xaxis.set_minor_locator(AutoMinorLocator(5))
        ax.yaxis.set_minor_locator(AutoMinorLocator(4))
        plt.xlabel("Wavelength (nm)", fontsize=16)
        plt.ylabel("Normalised Intensity", fontsize=16)
        plt.title(f"Detected and Fitted Spectral Lines for Spatial Pixel {specific_spatial_pixel}", fontsize=18)
        plt.show()

    if write_csv:
        if specific_spatial_pixel is None:
            pixel_identifier = "all_pixels"
        elif isinstance(specific_spatial_pixel, int):
            pixel_identifier = f"pixel_{specific_spatial_pixel}"
        elif isinstance(specific_spatial_pixel, range):
            pixel_identifier = f"pixels_{specific_spatial_pixel.start}-{specific_spatial_pixel.stop-1}"
        elif isinstance(specific_spatial_pixel, list):
            pixel_identifier = f"pixels_{min(specific_spatial_pixel)}-{max(specific_spatial_pixel)}"
        else:
            pixel_identifier = "pixels_custom"

        if date_identifier is None:
            date_identifier = "data"

        output_file = f"extracted_parameters_{date_identifier}_{pixel_identifier}.csv"
        with open(output_file, mode="w", newline="") as file:
            writer = csv.writer(file)
            writer.writerow(["Spatial Pixel", "Element", "Ref Wavelength (nm)", "Fitted centre (nm)", "LOS Velocity (km/s)", "Line Core Intensity"])
            writer.writerows(all_results)

        if save_fitted_spectra and len(all_fitted_spectra) > 0:
            all_fitted_spectra = np.array(all_fitted_spectra)
            fits.writeto(f"fitted_spectra_{date_identifier}_{pixel_identifier}.fits", all_fitted_spectra, overwrite=True)

        if not silent:
            print("\nFinal Extracted Spectral Line Parameters:")
            for line in all_results:
                print(
                    f"Pixel {line[0]} - {line[1]}: Ref Wavelength = {line[2]:.4f} nm, "
                    f"Fitted centre = {line[3]:.6f} nm, LOS Velocity = {line[4]:.2f} km/s, Intensity = {line[5]:.2f}"
                )

    if return_vlos_intensity:
        return vlos_out, int_out