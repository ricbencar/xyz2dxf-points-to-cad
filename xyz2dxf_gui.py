#!/usr/bin/env python3
r"""XYZ to DXF Converter — point-cloud processing and surface interpolation.

GENERAL DESCRIPTION
    Read irregular Cartesian X Y Z observations, optionally screen vertical
    outliers, apply deterministic minimum-distance thinning, and construct a
    regular surface grid. Export the retained observations, interpolated grid,
    confidence diagnostics, a DXF drawing and a processing report. DXF contains
    source/grid POINT entities, elevation TEXT labels and an optional boundary.

    The four interpolation choices are local RBF (default), Clough–Tocher
    cubic, linear TIN, and voxel reduction + cubic + Gaussian smoothing.
    Local overshoot protection is enabled by default, with bounded reference
    fits for sparse or nearly collinear support. NumPy/SciPy provide the
    established numerical routines; Numba compiles numerical processing loops.
    Grid evaluation and file exports are batched to limit working memory.

    Defaults: minimum distance 1 m, grid spacing 10 m, output precision 2
    decimals, PDMODE 3, robust outlier removal OFF, outer boundary ON, and
    overshoot protection ON. Distances are metres when input XY uses metres.
    All coordinates must use consistent Cartesian units. No coordinate-system
    transformation is performed; internal arithmetic uses float64.

INSTALLATION — WINDOWS COMMAND PROMPT / ANACONDA PROMPT
    Use 64-bit Python 3.10 or newer (Python 3.12 recommended), with Tcl/Tk
    enabled. Open the terminal in the folder containing xyz2dxf_gui.py.

    Create and activate an isolated environment:

        py -m venv .venv
        .venv\Scripts\activate.bat
        python -m pip install --upgrade pip setuptools wheel

    If the Windows py launcher is unavailable, use python -m venv .venv.

    Install the supplied requirements:

        python -m pip install --upgrade -r requirements.txt

    If only this Python file was downloaded, install the dependencies directly:

        python -m pip install "numpy>=1.26,<3" "scipy>=1.11,<2" "numba>=0.62,<1" "ezdxf>=1.3,<2" "threadpoolctl>=3.2,<4"

    Tkinter is included in the standard python.org Windows installation; do not
    install a package named tkinter with pip. Numba installs its LLVM runtime
    through llvmlite; no GPU or separate C/C++ compiler is required. The first
    numerical run can take longer while JIT kernels compile. Normal Python
    runs cache compiled kernels for later use.

USAGE
    Open the graphical interface:

        python xyz2dxf_gui.py

    Select the XYZ file, review the settings and click Run Conversion. Outlier
    removal is optional and starts unchecked. The worker runs in the background;
    Cancel stops between processing batches. Successful runs end the log with
    "Finished in xxx minutes." using the measured time to two decimal places.

    Optional command-line examples:

        python xyz2dxf_gui.py --input "C:\Data\survey.xyz"
        python xyz2dxf_gui.py --input survey.xyz --method linear --outliers
        python xyz2dxf_gui.py --input survey.xyz --method rbf --smoothing auto
        python xyz2dxf_gui.py --help

    For survey.xyz, the five products beside the source file are:
        survey.xyz.filtered.xyz
        survey.xyz.grid.xyz
        survey.xyz.confidence.xyz
        survey.xyz.dxf
        survey.xyz.rpt.txt

    The input is preserved. Completed runs replace products with the same
    names; failed/cancelled runs discard incomplete new data products.

BUILD A SINGLE WINDOWS EXECUTABLE — PYINSTALLER
    Run these commands in the script folder. Keep the environment activated
    while installing dependencies and building the executable:

        py -m venv .venv
        .venv\Scripts\activate.bat
        python -m pip install --upgrade pip setuptools wheel
        python -m pip install --upgrade -r requirements.txt
        python -m pip install --upgrade pyinstaller
        python -m PyInstaller --noconfirm --clean --onefile --windowed --collect-all numba --collect-all llvmlite --name xyz2dxf_gui xyz2dxf_gui.py

    If requirements.txt is unavailable, use the direct dependency-install
    command above before installing PyInstaller. The executable is created at
    dist\xyz2dxf_gui.exe. PyInstaller supplies the .exe extension on Windows.
    The --collect-all options include Numba/llvmlite modules and runtime files.
    Build and test the Windows executable on Windows. In a frozen executable,
    persistent JIT caching is disabled; required kernels compile in each new
    process. The supplied deliverable is Python source, not a prebuilt .exe.

VERIFICATION
        python -m unittest -v test_xyz2dxf

    See README.md for numerical methods, diagnostics and file-format details.
"""

from __future__ import annotations

import argparse
from array import array
from collections import OrderedDict
from dataclasses import dataclass, replace
import math
import os
from pathlib import Path
import queue
import sys
import tempfile
import threading
import time
from typing import Callable, Iterator
import warnings

try:
    import numpy as np
    import scipy
    from scipy import ndimage
    from scipy.interpolate import (
        CloughTocher2DInterpolator, LinearNDInterpolator, RBFInterpolator,
    )
    from scipy.spatial import ConvexHull, Delaunay, QhullError, cKDTree
    import ezdxf
    from ezdxf.lldxf.tagwriter import TagWriter
    from threadpoolctl import threadpool_limits
    import numba
    from numba import njit, prange, types
    from numba.typed import Dict
except ImportError as exc:
    _dependency_error = (
        f"Missing dependency: {exc}.\n\nInstall the required libraries with:\n"
        "python -m pip install -r requirements.txt\n\n"
        "Tkinter is included with the standard python.org Windows installer."
    )
    if __name__ == "__main__":
        print(_dependency_error, file=sys.stderr)
        if len(sys.argv) == 1:
            try:
                import tkinter as _tk
                from tkinter import messagebox as _mb
                _root = _tk.Tk()
                _root.withdraw()
                _mb.showerror("XYZ to DXF — dependencies", _dependency_error)
                _root.destroy()
            except Exception:
                pass
        raise SystemExit(1)
    raise


APP_TITLE = "XYZ to DXF Converter"
METHODS = {
    "rbf": "Local RBF (default)",
    "cubic": "Clough–Tocher cubic",
    "linear": "Linear TIN",
    "voxel": "Voxel + cubic + Gaussian",
}
CONFIDENCE_HEADER = (
    "# x y z confidence residual lambda boundaryPenalty fallback "
    "outsideHull clampedLocal neighborCount patchCount\n"
)
SUFFIXES = (".filtered.xyz", ".grid.xyz", ".confidence.xyz", ".dxf", ".rpt.txt")
Progress = Callable[[str, float | None], None]
JIT_CACHE = not getattr(sys, "frozen", False)
CPU_THREADS = max(1, min(8, os.cpu_count() or 1, numba.config.NUMBA_NUM_THREADS))
CELL_KEY_TYPE = types.UniTuple(types.int64, 2)


def configure_acceleration() -> None:
    """Set a bounded native-thread count in the calling worker thread."""
    if not numba.config.DISABLE_JIT:
        numba.set_num_threads(CPU_THREADS)


@njit(cache=JIT_CACHE, nogil=True, parallel=True, fastmath=False)
def _outlier_keep_batch(points, distances, indices, start, radius):
    """Independent local plane/MAD tests; no shared mutable state across rows."""
    rows, k = indices.shape
    keep = np.ones(rows, dtype=np.bool_)
    for row in prange(rows):
        own = start + row
        count = 0
        mean_r2 = 0.0
        zcentre = 0.0
        for j in range(k):
            index = indices[row, j]
            distance = distances[row, j]
            if index != own and distance <= radius:
                count += 1
                mean_r2 += distance * distance
                zcentre += points[index, 2]
        if count < 6:
            continue
        scale = math.sqrt(max(mean_r2 / count, 1e-12))
        zcentre /= count
        gram = np.zeros((3, 3))
        rhs = np.zeros(3)
        for j in range(k):
            index = indices[row, j]
            distance = distances[row, j]
            if index == own or distance > radius:
                continue
            u = (points[index, 0] - points[own, 0]) / scale
            v = (points[index, 1] - points[own, 1]) / scale
            weight = 1.0 / (1.0 + (distance / scale) ** 2)
            value = points[index, 2] - zcentre
            gram[0, 0] += weight
            gram[0, 1] += weight * u
            gram[0, 2] += weight * v
            gram[1, 1] += u * weight * u
            gram[1, 2] += u * weight * v
            gram[2, 2] += v * weight * v
            rhs[0] += weight * value
            rhs[1] += u * weight * value
            rhs[2] += v * weight * value
        gram[1, 0] = gram[0, 1]
        gram[2, 0] = gram[0, 2]
        gram[2, 1] = gram[1, 2]
        eigenvalues = np.linalg.eigvalsh(gram)
        if eigenvalues[0] <= 1e-10 * max(eigenvalues[2], 1.0):
            continue
        for j in range(3):
            gram[j, j] += 1e-12
        coefficients = np.linalg.solve(gram, rhs)
        residuals = np.empty(count)
        offset = 0
        for j in range(k):
            index = indices[row, j]
            if index == own or distances[row, j] > radius:
                continue
            u = (points[index, 0] - points[own, 0]) / scale
            v = (points[index, 1] - points[own, 1]) / scale
            fitted = coefficients[0] + u * coefficients[1] + v * coefficients[2]
            residuals[offset] = points[index, 2] - zcentre - fitted
            offset += 1
        median = np.median(residuals)
        for j in range(count):
            residuals[j] = abs(residuals[j] - median)
        sigma = max(1.4826 * np.median(residuals), 1e-8)
        error = abs(points[own, 2] - zcentre - coefficients[0] - median)
        keep[row] = error <= 3.5 * sigma
    return keep


@njit(cache=JIT_CACHE, nogil=True, fastmath=False)
def _thin_batch(points, cells, order, distance_squared, heads, links, accepted, count):
    """Serial greedy ordering is intentional: parallel acceptance would change it."""
    for sequence in range(len(order)):
        index = order[sequence]
        x, y = points[index, 0], points[index, 1]
        cx, cy = cells[index, 0], cells[index, 1]
        too_close = False
        for dx in range(-1, 2):
            for dy in range(-1, 2):
                neighbour_key = (cx + dx, cy + dy)
                other = heads[neighbour_key] if neighbour_key in heads else -1
                while other >= 0:
                    ddx = x - points[other, 0]
                    ddy = y - points[other, 1]
                    if ddx * ddx + ddy * ddy < distance_squared:
                        too_close = True
                        break
                    other = links[other]
                if too_close:
                    break
            if too_close:
                break
        if not too_close:
            accepted[count] = index
            count += 1
            key = (cx, cy)
            links[index] = heads[key] if key in heads else -1
            heads[key] = index
    return count


@njit(cache=JIT_CACHE, nogil=True, fastmath=False)
def _classify_hull(xy, vertices, equations, rank):
    outside = np.zeros(len(xy), dtype=np.bool_)
    distance = np.zeros(len(xy))
    for row in range(len(xy)):
        x, y = xy[row, 0], xy[row, 1]
        is_outside = rank < 2
        if rank == 2:
            for j in range(len(equations)):
                if x * equations[j, 0] + y * equations[j, 1] + equations[j, 2] > 1e-10:
                    is_outside = True
                    break
        if is_outside:
            nearest_squared = np.inf
            for j in range(len(vertices)):
                a, b = vertices[j], vertices[(j + 1) % len(vertices)]
                ex, ey = b[0] - a[0], b[1] - a[1]
                norm_squared = ex * ex + ey * ey
                t = ((x - a[0]) * ex + (y - a[1]) * ey) / norm_squared if norm_squared else 0.0
                t = min(1.0, max(0.0, t))
                dx, dy = x - (a[0] + t * ex), y - (a[1] + t * ey)
                nearest_squared = min(nearest_squared, dx * dx + dy * dy)
            distance[row] = math.sqrt(nearest_squared)
        outside[row] = is_outside if rank == 2 else distance[row] > 1e-10
    return outside, distance


@njit(cache=JIT_CACHE, nogil=True, fastmath=False)
def _guard_support(xy, values, spacing, distances, indices):
    """Local range, positive-weight estimate and geometry/density quality."""
    result = np.empty((len(indices), 4))
    for row in range(len(indices)):
        count = indices.shape[1]
        first = indices[row, 0]
        lower, upper = values[first], values[first]
        mean_x = mean_y = total = norm = 0.0
        local_spacing = np.empty(count)
        reference_distance = max(distances[row, 0], 1e-12)
        for j in range(count):
            index = indices[row, j]
            value = values[index]
            lower, upper = min(lower, value), max(upper, value)
            mean_x += xy[index, 0]
            mean_y += xy[index, 1]
            local_spacing[j] = spacing[index]
            weight = (reference_distance / max(distances[row, j], 1e-12)) ** 2
            total += weight * value
            norm += weight
        mean_x /= count
        mean_y /= count
        xx = yy = cross = 0.0
        for j in range(count):
            index = indices[row, j]
            dx, dy = xy[index, 0] - mean_x, xy[index, 1] - mean_y
            xx += dx * dx
            yy += dy * dy
            cross += dx * dy
        trace = xx + yy
        root = math.sqrt(max((xx - yy) ** 2 + 4 * cross * cross, 0.0))
        aspect = max(trace - root, 0.0) / max(trace + root, 1e-24)
        distance_ratio = distances[row, 0] / max(np.median(local_spacing), 1e-12)
        sparse = min(1.0, max(0.0, (distance_ratio - 1.5) / 2.5))
        elongated = min(1.0, max(0.0, (0.03 - aspect) / 0.03))
        risk = max(sparse, elongated)
        risk = risk * risk * (3.0 - 2.0 * risk)
        if distances[row, 0] <= 1e-10:
            estimate = values[first]
            risk = 0.0
        else:
            estimate = total / norm
        result[row, 0], result[row, 1] = lower, upper
        result[row, 2], result[row, 3] = estimate, risk
    return result


@njit(cache=JIT_CACHE, nogil=True, fastmath=False)
def _canonical_simplex(queries, simplex, transforms, neighbours):
    """Choose the same containing triangle at shared edges, in any batch order."""
    result = simplex.copy()
    for row in range(len(queries)):
        if simplex[row] < 0:
            continue
        pending = [simplex[row]]
        offset = 0
        chosen = simplex[row]
        while offset < len(pending):
            face = pending[offset]
            offset += 1
            transform = transforms[face]
            dx = queries[row, 0] - transform[2, 0]
            dy = queries[row, 1] - transform[2, 1]
            a = transform[0, 0] * dx + transform[0, 1] * dy
            b = transform[1, 0] * dx + transform[1, 1] * dy
            c = 1.0 - a - b
            if min(a, b, c) < -1e-9:
                continue
            chosen = min(chosen, face)
            for edge in range(3):
                weight = a if edge == 0 else b if edge == 1 else c
                neighbour = neighbours[face, edge]
                if abs(weight) <= 1e-9 and neighbour >= 0 and neighbour not in pending:
                    pending.append(neighbour)
        result[row] = chosen
    return result


def containing_simplex(triangulation, queries):
    initial = triangulation.find_simplex(queries)
    return _canonical_simplex(queries, initial, triangulation.transform, triangulation.neighbors)


class CancelledError(RuntimeError):
    """Cooperative cancellation between processing batches."""


def check_cancel(cancel: threading.Event | None) -> None:
    if cancel is not None and cancel.is_set():
        raise CancelledError("Conversion cancelled. Existing data outputs were preserved.")


@dataclass(frozen=True)
class Settings:
    input_file: Path
    min_dist: float = 1.0
    precision: int = 2
    pdmode: int = 3
    grid_spacing: float = 10.0
    use_outliers: bool = False
    use_boundary: bool = True
    prevent_overshoot: bool = True
    method: str = "rbf"
    neighbors: int = 50
    smoothing: float | str = 0.1
    kernel: str = "thin_plate_spline"
    voxel_size: float | str = "auto"
    gaussian_sigma: float = 1.0
    batch_size: int = 4096
    validation_samples: int = 256
    max_boundary_cells: int = 2_000_000

    def validate(self) -> None:
        if not Path(self.input_file).is_file():
            raise ValueError("Please select an existing input XYZ file.")
        if not math.isfinite(self.min_dist) or self.min_dist < 0:
            raise ValueError("Min Dist must be a finite number >= 0.")
        if not isinstance(self.precision, int) or not 0 <= self.precision <= 15:
            raise ValueError("Precision must be an integer between 0 and 15.")
        if not isinstance(self.pdmode, int) or not -32768 <= self.pdmode <= 32767:
            raise ValueError("PDMODE must be a DXF 16-bit integer (-32768 to 32767).")
        if not math.isfinite(self.grid_spacing) or self.grid_spacing <= 0:
            raise ValueError("Grid Spacing must be a positive finite number.")
        if self.method not in METHODS:
            raise ValueError(f"Unknown interpolation method: {self.method}")
        if self.method == "rbf":
            if not isinstance(self.neighbors, int) or not 3 <= self.neighbors <= 512:
                raise ValueError("RBF Neighbours must be an integer between 3 and 512.")
            if self.kernel not in ("thin_plate_spline", "cubic"):
                raise ValueError("RBF kernel must be thin_plate_spline or cubic.")
            if self.smoothing != "auto":
                try:
                    value = float(self.smoothing)
                except (TypeError, ValueError) as exc:
                    raise ValueError("RBF Smoothing must be >= 0, or auto.") from exc
                if not math.isfinite(value) or value < 0:
                    raise ValueError("RBF Smoothing must be finite and >= 0, or auto.")
        if self.method == "voxel":
            if self.voxel_size != "auto":
                try:
                    value = float(self.voxel_size)
                except (TypeError, ValueError) as exc:
                    raise ValueError("Voxel Cell Size must be positive, or auto.") from exc
                if not math.isfinite(value) or value <= 0:
                    raise ValueError("Voxel Cell Size must be positive, or auto.")
            if not math.isfinite(self.gaussian_sigma) or not 0 <= self.gaussian_sigma <= 4:
                raise ValueError("Gaussian Sigma must be between 0 and 4 grid cells.")
        if not isinstance(self.batch_size, int) or self.batch_size < 64:
            raise ValueError("Batch size must be an integer >= 64.")
        if not isinstance(self.validation_samples, int) or self.validation_samples < 0:
            raise ValueError("Validation samples must be an integer >= 0.")
        if self.max_boundary_cells < 10000:
            raise ValueError("Boundary raster limit must be >= 10000 cells.")


@dataclass
class ReadStats:
    lines: int = 0
    accepted: int = 0
    ignored: int = 0
    nonfinite: int = 0


@dataclass
class Validation:
    xy: np.ndarray
    errors: np.ndarray
    fallback_count: int = 0

    @classmethod
    def empty(cls) -> Validation:
        return cls(np.empty((0, 2)), np.empty(0))

    @property
    def count(self) -> int:
        return len(self.errors)

    @property
    def rmse(self) -> float:
        return float(np.sqrt(np.mean(self.errors ** 2))) if self.count else math.nan

    @property
    def mae(self) -> float:
        return float(np.mean(np.abs(self.errors))) if self.count else math.nan

    @property
    def p95(self) -> float:
        return float(np.percentile(np.abs(self.errors), 95)) if self.count else math.nan

    def description(self) -> str:
        if not self.count:
            return "Predictive validation unavailable (too few points, disabled, or failed)."
        return (f"Independent holdout: RMSE={self.rmse:.8g}, MAE={self.mae:.8g}, "
                f"P95={self.p95:.8g}, n={self.count}, fallback predictions={self.fallback_count}.")


@dataclass(frozen=True)
class Grid:
    xmin: float
    ymin: float
    spacing: float
    nx: int
    ny: int

    @classmethod
    def from_points(cls, points: np.ndarray, spacing: float) -> Grid:
        lower = points[:, :2].min(axis=0) - 1.5 * spacing
        upper = points[:, :2].max(axis=0) + 1.5 * spacing
        widths = np.maximum(upper - lower, 1.0)
        if not np.all(np.isfinite(widths / spacing)):
            raise ValueError("Grid extent is too large for the requested spacing.")
        nx, ny = (math.ceil(float(w) / spacing) + 1 for w in widths)
        if nx * ny > np.iinfo(np.int64).max:
            raise ValueError("Grid dimensions exceed the 64-bit node count limit.")
        if lower[0] + spacing == lower[0] or lower[1] + spacing == lower[1]:
            raise ValueError("Grid spacing is below floating-point resolution at these coordinates.")
        return cls(float(lower[0]), float(lower[1]), spacing, nx, ny)

    @property
    def total(self) -> int:
        return self.nx * self.ny

    def coordinates(self, i0: int, i1: int, j0: int, j1: int) -> np.ndarray:
        x = self.xmin + np.arange(i0, i1, dtype=float) * self.spacing
        y = self.ymin + np.arange(j0, j1, dtype=float) * self.spacing
        return np.column_stack((np.repeat(x, len(y)), np.tile(y, len(x))))

    def contains(self, xy: np.ndarray) -> np.ndarray:
        tolerance = 1e-8 * self.spacing
        return ((xy[:, 0] >= self.xmin - tolerance) &
                (xy[:, 0] <= self.xmin + (self.nx - 1) * self.spacing + tolerance) &
                (xy[:, 1] >= self.ymin - tolerance) &
                (xy[:, 1] <= self.ymin + (self.ny - 1) * self.spacing + tolerance))


class Reporter:
    def __init__(self, path: Path, progress: Progress | None):
        self.stream = path.open("w", encoding="utf-8", newline="\n")
        self.progress = progress

    def emit(self, message: str, fraction: float | None = None) -> None:
        self.stream.write(message + "\n")
        self.stream.flush()
        if self.progress is not None:
            self.progress(message, fraction)

    def close(self) -> None:
        self.stream.close()


def read_xyz(path: Path, cancel: threading.Event | None = None,
             log: Callable[[str], None] | None = None) -> tuple[np.ndarray, ReadStats]:
    """Parse the first three numeric fields; ignore headers and malformed rows."""
    values = array("d")
    stats = ReadStats()
    with Path(path).open("r", encoding="utf-8-sig", errors="replace") as stream:
        for line in stream:
            stats.lines += 1
            if stats.lines % 10000 == 0:
                check_cancel(cancel)
            text = line.strip()
            if not text or text.startswith("#"):
                stats.ignored += 1
                continue
            fields = text.replace(",", " ").split()
            try:
                if len(fields) < 3:
                    raise ValueError
                xyz = (float(fields[0]), float(fields[1]), float(fields[2]))
            except ValueError:
                stats.ignored += 1
                continue
            if not all(math.isfinite(v) for v in xyz):
                stats.nonfinite += 1
                continue
            values.extend(xyz)
            stats.accepted += 1
            if log and stats.accepted % 250000 == 0:
                log(f"Read {stats.accepted:,} valid XYZ points ...")
    check_cancel(cancel)
    if not stats.accepted:
        raise ValueError("No valid finite XYZ points were found in the input file.")
    return np.frombuffer(values, dtype=np.float64).reshape(-1, 3), stats


def planar_spacing(points: np.ndarray) -> float:
    if len(points) < 2:
        return 1.0
    span = np.ptp(points[:, :2], axis=0)
    return math.sqrt(max(float(span[0] * span[1]), 1e-12) / len(points))


def remove_outliers(points: np.ndarray, min_dist: float,
                    cancel: threading.Event | None = None,
                    log: Callable[[str], None] | None = None) -> np.ndarray:
    """Batched local plane/MAD screening, before thinning, as in the C++ pipeline."""
    if len(points) < 7:
        return points
    # Stable ordering also makes equal-distance KD-tree ties reproducible.
    points = points[np.lexsort((points[:, 2], points[:, 1], points[:, 0]))]
    radius = float(max(5.0 * min_dist, 3.0 * planar_spacing(points), 0.01))
    if log:
        log(f"Outlier neighbourhood radius={radius:.8g}; MAD threshold factor=3.5.")
    tree = cKDTree(points[:, :2])
    keep = np.ones(len(points), dtype=bool)
    k = min(len(points), 96)
    last = time.monotonic()
    configure_acceleration()
    workers = CPU_THREADS if len(points) >= 65536 else 1
    for start in range(0, len(points), 8192):
        check_cancel(cancel)
        stop = min(start + 8192, len(points))
        distance, indices = tree.query(points[start:stop, :2], k=k, workers=workers)
        keep[start:stop] = _outlier_keep_batch(points, distance, indices, start, radius)
        if log and (time.monotonic() - last > 1.0 or stop == len(points)):
            log(f"Outlier scan: {stop:,}/{len(points):,}; removed {np.count_nonzero(~keep[:stop]):,}.")
            last = time.monotonic()
    result = points[keep]
    if not len(result):
        raise ValueError("All points were removed by robust outlier screening.")
    return result


def thin_points(points: np.ndarray, distance: float,
                cancel: threading.Event | None = None,
                log: Callable[[str], None] | None = None) -> np.ndarray:
    """Cell-centre ordering and strict XY distance test preserve C++ semantics."""
    check_cancel(cancel)
    if distance <= 0:
        return np.unique(points, axis=0)
    scaled = np.floor(points[:, :2] / distance)
    if not np.all(np.isfinite(scaled)) or np.max(np.abs(scaled)) >= 9e18:
        raise ValueError("Min Dist is too small for the coordinate magnitudes.")
    cells = scaled.astype(np.int64)
    centre_delta = points[:, :2] - (cells.astype(float) + 0.5) * distance
    score = np.sum(centre_delta ** 2, axis=1)
    order = np.lexsort((points[:, 2], points[:, 1], points[:, 0], score,
                       cells[:, 1], cells[:, 0]))
    heads = Dict.empty(key_type=CELL_KEY_TYPE, value_type=types.int64)
    links = np.full(len(points), -1, dtype=np.int64)
    accepted = np.empty(len(points), dtype=np.int64)
    count = 0
    last = time.monotonic()
    for start in range(0, len(order), 32768):
        check_cancel(cancel)
        stop = min(start + 32768, len(order))
        count = _thin_batch(points, cells, order[start:stop], float(distance * distance),
                            heads, links, accepted, count)
        if log and (time.monotonic() - last > 1.0 or stop == len(order)):
            log(f"Minimum-distance thinning: {stop:,}/{len(points):,}; retained={count:,}.")
            last = time.monotonic()
    return points[accepted[:count]]


def unique_xy_controls(points: np.ndarray) -> tuple[np.ndarray, int]:
    """Keep the median Z at coincident XY for a well-defined single-valued fit.

    This does not modify the filtered XYZ product (which preserves all distinct
    XYZ observations when Min Dist is zero).
    """
    order = np.lexsort((points[:, 2], points[:, 1], points[:, 0]))
    ordered = points[order]
    starts = np.r_[0, np.flatnonzero(np.any(np.diff(ordered[:, :2], axis=0), axis=1)) + 1]
    counts = np.diff(np.r_[starts, len(ordered)])
    controls = ordered[starts].copy()
    for i in np.flatnonzero(counts > 1):
        controls[i, 2] = np.median(ordered[starts[i]:starts[i] + counts[i], 2])
    return controls, int(len(points) - len(controls))


def voxel_reduce(points: np.ndarray, size: float) -> np.ndarray:
    """Sparse 2-D voxel means: O(N) arrays, never a dense bounding-box matrix."""
    scaled = np.floor(points[:, :2] / size)
    if not np.all(np.isfinite(scaled)) or np.max(np.abs(scaled)) >= 9e18:
        raise ValueError("Voxel cell size is too small for the coordinate magnitudes.")
    cells = scaled.astype(np.int64)
    order = np.lexsort((points[:, 2], points[:, 1], points[:, 0],
                       cells[:, 1], cells[:, 0]))
    cells = cells[order]
    ordered = points[order]
    starts = np.r_[0, np.flatnonzero(np.any(np.diff(cells, axis=0), axis=1)) + 1]
    counts = np.diff(np.r_[starts, len(ordered)])
    # Average actual XY locations, rather than displacing observations to cell centres.
    origin = ordered[0]
    return np.add.reduceat(ordered - origin, starts, axis=0) / counts[:, None] + origin


def characteristic_spacing(xy: np.ndarray) -> float:
    if len(xy) < 2:
        return 1.0
    sample = xy[np.linspace(0, len(xy) - 1, min(4096, len(xy)), dtype=int)]
    distance = cKDTree(xy).query(sample, k=2, workers=1)[0][:, 1]
    positive = distance[distance > 0]
    return max(float(np.median(positive)) if len(positive) else 1.0, 1e-9)


@dataclass
class Hull:
    """Convex support diagnostics without a full Delaunay triangulation for RBF."""
    vertices: np.ndarray
    equations: np.ndarray | None
    rank: int
    direction: np.ndarray | None = None

    @classmethod
    def build(cls, xy: np.ndarray) -> Hull:
        centre = xy.mean(axis=0)
        _, singular, vectors = np.linalg.svd(xy - centre, full_matrices=False)
        rank = int(np.count_nonzero(singular > max(float(singular[0]), 1e-12) * 1e-12))
        if len(xy) >= 3 and rank == 2:
            try:
                hull = ConvexHull(xy)
                return cls(xy[hull.vertices], hull.equations, 2)
            except QhullError:
                rank = 1
        if rank == 0:
            return cls(xy[:1], None, 0)
        direction = vectors[0]
        projection = xy @ direction
        return cls(xy[[np.argmin(projection), np.argmax(projection)]], None, 1, direction)

    def classify(self, xy: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        equations = self.equations if self.equations is not None else np.empty((0, 3))
        return _classify_hull(np.ascontiguousarray(xy), self.vertices, equations, self.rank)


def _outer_loop(mask: np.ndarray) -> np.ndarray:
    """Trace exposed cell edges, choosing the largest closed outer loop."""
    adjacency: dict[tuple[int, int], list[tuple[int, int]]] = {}
    below = np.zeros_like(mask)
    below[:, 1:] = mask[:, :-1]
    above = np.zeros_like(mask)
    above[:, :-1] = mask[:, 1:]
    left = np.zeros_like(mask)
    left[1:, :] = mask[:-1, :]
    right = np.zeros_like(mask)
    right[:-1, :] = mask[1:, :]
    for exposed, offset_a, offset_b in (
        (mask & ~below, (0, 0), (1, 0)),
        (mask & ~right, (1, 0), (1, 1)),
        (mask & ~above, (1, 1), (0, 1)),
        (mask & ~left, (0, 1), (0, 0)),
    ):
        ix, iy = np.nonzero(exposed)
        for x, y in zip(ix.tolist(), iy.tolist()):
            a = (x + offset_a[0], y + offset_a[1])
            b = (x + offset_b[0], y + offset_b[1])
            adjacency.setdefault(a, []).append(b)
    best = np.empty((0, 2))
    best_area = 0.0
    while adjacency:
        start = next(iter(adjacency))
        current = start
        loop = []
        while current in adjacency:
            loop.append(current)
            outgoing = adjacency[current]
            following = outgoing.pop()
            if not outgoing:
                del adjacency[current]
            current = following
            if current == start:
                break
        if current != start or len(loop) < 4:
            continue
        vertices = np.asarray(loop, dtype=float)
        edge_before = vertices - np.roll(vertices, 1, axis=0)
        edge_after = np.roll(vertices, -1, axis=0) - vertices
        cross = edge_before[:, 0] * edge_after[:, 1] - edge_before[:, 1] * edge_after[:, 0]
        vertices = vertices[cross != 0]
        area = abs(float(np.sum(vertices[:, 0] * np.roll(vertices[:, 1], -1)
                                - vertices[:, 1] * np.roll(vertices[:, 0], -1))))
        if area > best_area:
            best, best_area = vertices, area
    return best


@dataclass
class Boundary:
    origin: np.ndarray
    cell_size: float
    mask: np.ndarray | None
    loop: np.ndarray
    convex: Hull | None = None

    def contains(self, xy: np.ndarray) -> np.ndarray:
        if self.convex is not None:
            return ~self.convex.classify(xy - self.origin)[0]
        position = (xy - self.origin) / self.cell_size
        # Include the closed cell edges, including a point on an upper raster edge.
        near_edge = np.abs(position - np.rint(position)) < 1e-8
        position = np.where(near_edge, np.rint(position), position)
        cells = np.floor(position).astype(np.int64)
        inside = np.zeros(len(xy), dtype=bool)
        for dx, dy in ((0, 0), (-1, 0), (0, -1), (-1, -1)):
            indices = cells + np.array([dx, dy])
            valid = ((indices[:, 0] >= 0) & (indices[:, 0] < self.mask.shape[0]) &
                     (indices[:, 1] >= 0) & (indices[:, 1] < self.mask.shape[1]))
            if dx:
                valid &= near_edge[:, 0]
            if dy:
                valid &= near_edge[:, 1]
            ids = np.flatnonzero(valid)
            inside[ids] |= self.mask[indices[ids, 0], indices[ids, 1]]
        return inside


def build_boundary(points: np.ndarray, settings: Settings,
                   cancel: threading.Event | None,
                   log: Callable[[str], None]) -> Boundary:
    xy = points[:, :2]
    lower, upper = xy.min(axis=0), xy.max(axis=0)
    cell = max(settings.min_dist if settings.min_dist > 0 else planar_spacing(points), 1e-6)
    # The full closing halo prevents erosion of outermost observed cells.
    margin = 35
    for attempt in range(6):
        check_cancel(cancel)
        while True:
            dims_float = np.ceil((upper - lower) / cell) + 2 * margin + 2
            total = float(np.prod(dims_float))
            if np.isfinite(total) and total <= settings.max_boundary_cells:
                break
            factor = max(1.1, math.sqrt(total / settings.max_boundary_cells)) if math.isfinite(total) else 10.0
            cell *= factor
            if not math.isfinite(cell):
                raise ValueError("Coordinates span an unsupported boundary extent.")
        nx, ny = (int(v) for v in dims_float)
        origin = lower - margin * cell
        indices = np.floor((xy - origin) / cell).astype(int)
        raw = np.zeros((nx, ny), dtype=bool)
        raw[indices[:, 0], indices[:, 1]] = True
        mask = ndimage.binary_fill_holes(raw)
        components = ndimage.label(mask)[1]
        used_radius = 0
        for radius in range(1, 33):
            if components <= 1:
                break
            check_cancel(cancel)
            dilated = ndimage.maximum_filter(raw, size=2 * radius + 1, mode="constant", cval=0)
            closed = ndimage.minimum_filter(dilated, size=2 * radius + 1, mode="constant", cval=0)
            mask = ndimage.binary_fill_holes(closed | raw)
            components = ndimage.label(mask)[1]
            used_radius = radius
        log(f"Boundary iteration {attempt + 1}: cell={cell:.8g}, raster={nx} x {ny}, "
            f"components={components}, closing radius={used_radius}.")
        if components == 1:
            loop = _outer_loop(mask)
            if len(loop) >= 3:
                boundary = Boundary(origin, cell, mask, origin + loop * cell)
                if np.all(boundary.contains(xy)):
                    return boundary
        cell *= 1.35
    # A connected concave footprint cannot always describe disjoint surveys.
    # Keep every observation and document the convex-envelope fallback.
    centre = lower + 0.5 * (upper - lower)
    hull = Hull.build(xy - centre)
    if hull.rank != 2:
        raise ValueError("Unable to construct a connected outer boundary. Disable Use outer boundary.")
    log("Boundary fallback: convex envelope used to contain all disconnected survey groups.")
    return Boundary(centre, cell, None, hull.vertices + centre, hull)


@dataclass
class Prediction:
    z: np.ndarray
    fallback: np.ndarray
    neighbors: np.ndarray
    patches: np.ndarray
    clamped: np.ndarray | None = None

    def __post_init__(self):
        if self.clamped is None:
            self.clamped = np.zeros(len(self.z), dtype=bool)


class Surface:
    """Established SciPy interpolators, with explicit and traceable fallbacks."""

    def __init__(self, points: np.ndarray, settings: Settings,
                 log: Callable[[str], None] | None = None):
        if settings.method not in METHODS:
            raise ValueError(f"Unknown interpolation method: {settings.method}")
        self.settings = settings
        self.log = log or (lambda message: None)
        self.voxel_size: float | None = None
        if settings.method == "voxel":
            self.voxel_size = (max(settings.grid_spacing / 2.0, 1e-9)
                               if settings.voxel_size == "auto" else float(settings.voxel_size))
            original_count = len(points)
            points = voxel_reduce(points, self.voxel_size)
            self.log(f"Voxel aggregation: {original_count:,} XY controls -> {len(points):,} means; "
                     f"cell size={self.voxel_size:.8g}.")
        self.points = points
        self.origin = points[:, :2].min(axis=0) + 0.5 * np.ptp(points[:, :2], axis=0)
        self.scale = characteristic_spacing(points[:, :2])
        self.xy = np.ascontiguousarray((points[:, :2] - self.origin) / self.scale)
        self.zoffset = float(np.median(points[:, 2]))
        self.values = points[:, 2] - self.zoffset
        self.tree = cKDTree(self.xy)
        self.hull = Hull.build(self.xy)
        self.interpolator = None
        self.linear = None
        self.guard_cache = OrderedDict()
        self.control_spacing = np.ones(len(points))
        if settings.prevent_overshoot and settings.method in ("rbf", "cubic", "voxel") and len(points) > 1:
            for start in range(0, len(points), 32768):
                stop = min(start + 32768, len(points))
                self.control_spacing[start:stop] = self.tree.query(
                    self.xy[start:stop], k=2, workers=1)[0][:, 1]
        self.rbf_models: dict[int, RBFInterpolator] = {}
        self.effective_neighbors = min(settings.neighbors, len(points))
        self.global_fallback = self.hull.rank < 2
        self.smoothing = float(settings.smoothing) if settings.method == "rbf" else 0.0
        if self.global_fallback:
            self.log("Surface controls are collinear or insufficient for 2-D fitting; "
                     "nearest-neighbour fallback is flagged for every prediction.")
        elif settings.method == "rbf":
            self.interpolator = self._rbf(self.effective_neighbors)
        elif settings.method in ("cubic", "linear", "voxel"):
            # Delaunay storage is required only for the TIN-based methods.
            try:
                triangulation = Delaunay(self.xy)
                self.linear = LinearNDInterpolator(triangulation, self.values, fill_value=np.nan)
                if settings.method == "linear":
                    self.interpolator = self.linear
                else:
                    with warnings.catch_warnings(record=True) as caught:
                        warnings.simplefilter("always")
                        self.interpolator = CloughTocher2DInterpolator(
                            triangulation, self.values, fill_value=np.nan,
                            tol=1e-6, maxiter=400)
                    if caught:
                        self.log("Cubic gradient estimation warning: using linear TIN fallback.")
                        self.interpolator = self.linear
                        self.global_fallback = True
            except QhullError:
                self.log("Delaunay triangulation failed; nearest-neighbour fallback is flagged.")
                self.global_fallback = True
                self.interpolator = None
        self.log(f"Interpolation support: {len(points):,} controls, spatial rank={self.hull.rank}, "
                 f"isotropic XY normalization scale={self.scale:.8g}.")

    def _rbf(self, count: int) -> RBFInterpolator:
        if count not in self.rbf_models:
            self.rbf_models[count] = RBFInterpolator(
                self.xy, self.values, neighbors=count,
                kernel=self.settings.kernel, degree=1, smoothing=self.smoothing)
        return self.rbf_models[count]

    def _nearest_indices(self, q: np.ndarray) -> np.ndarray:
        if len(self.xy) == 1:
            return np.zeros(len(q), dtype=int)
        distance, indices = self.tree.query(q, k=2, workers=1)
        result = indices[:, 0].copy()
        tolerance = 1e-12 * np.maximum(distance[:, 0], 1.0)
        ties = np.flatnonzero(distance[:, 1] - distance[:, 0] <= tolerance)
        if len(ties):
            groups = self.tree.query_ball_point(q[ties], distance[ties, 0] + tolerance[ties],
                                                return_sorted=True, workers=1)
            for index, group in zip(ties, groups):
                result[index] = min(group)
        return result

    def predict(self, xy: np.ndarray,
                cancel: threading.Event | None = None) -> Prediction:
        xy = np.asarray(xy, dtype=np.float64).reshape(-1, 2)
        if not len(xy):
            return Prediction(np.empty(0), np.empty(0, bool), np.empty(0, int), np.empty(0, int))
        check_cancel(cancel)
        q = (xy - self.origin) / self.scale
        nearest_indices = self.tree.query(q, k=1, workers=1)[1]
        nearest_values = self.values[nearest_indices]
        z = nearest_values.copy()
        fallback = np.full(len(q), self.global_fallback, dtype=bool)
        neighbors = np.ones(len(q), dtype=int)
        patches = np.zeros(len(q), dtype=int)
        method = self.settings.method
        if method == "rbf" and not self.global_fallback:
            # A single singular local neighbourhood must not fail a whole batch.
            # Split failed batches, expand local support, then use explicit nearest.
            def evaluate(indices: np.ndarray, count: int) -> None:
                check_cancel(cancel)
                try:
                    result = np.asarray(self._rbf(count)(q[indices])).ravel()
                    good = np.isfinite(result)
                    ids = indices[good]
                    z[ids] = result[good]
                    neighbors[ids] = count
                    patches[ids] = 1
                    if count != self.effective_neighbors:
                        fallback[ids] = True
                    if np.all(good):
                        return
                    indices = indices[~good]
                except (np.linalg.LinAlgError, ValueError):
                    pass
                if len(indices) > 1:
                    mid = len(indices) // 2
                    evaluate(indices[:mid], count)
                    evaluate(indices[mid:], count)
                elif count < min(512, len(self.xy)):
                    evaluate(indices, min(max(count + 1, 2 * count), 512, len(self.xy)))
                else:
                    fallback[indices] = True
            evaluate(np.arange(len(q)), self.effective_neighbors)
        elif self.interpolator is not None:
            values = np.asarray(self.interpolator(q)).ravel()
            good = np.isfinite(values)
            z[good] = values[good]
            neighbors[good] = 3
            patches[good] = 1
            failed = ~good
            # Cubic numerical failures inside a triangle can still use the TIN.
            if np.any(failed) and self.linear is not None and self.interpolator is not self.linear:
                ids = np.flatnonzero(failed)
                linear = np.asarray(self.linear(q[ids])).ravel()
                usable = np.isfinite(linear)
                z[ids[usable]] = linear[usable]
                neighbors[ids[usable]] = 3
                patches[ids[usable]] = 1
            fallback[failed] = True
        else:
            fallback[:] = True
        nearest_used = (patches == 0) & fallback
        if np.any(nearest_used):
            # Stable equal-distance ties also prevent Gaussian tile seams caused
            # by equivalent coordinates differing in their final floating bit.
            ids = self._nearest_indices(q[nearest_used])
            z[nearest_used] = self.values[ids]
        z += self.zoffset
        if not np.all(np.isfinite(z)):
            raise ValueError("Interpolation produced a non-finite Z value.")
        return self.protect(xy, Prediction(z, fallback, neighbors, patches), cancel)

    def protect(self, xy: np.ndarray, prediction: Prediction,
                cancel: threading.Event | None = None) -> Prediction:
        """Bound smooth fits and favour convex interpolation in sparse support.

        The RBF path uses small cached local TINs only where needed, preserving
        its bounded-memory design. Linear TIN and nearest are already bounded.
        """
        if (not self.settings.prevent_overshoot or not len(xy) or self.hull.rank < 2
                or self.settings.method == "linear"):
            return prediction
        check_cancel(cancel)
        q = np.ascontiguousarray((xy - self.origin) / self.scale)
        count = min(len(self.xy), max(16, min(self.settings.neighbors, 64)))
        distance, indices = self.tree.query(q, k=count, workers=1)
        distance = distance.reshape(len(q), count)
        indices = indices.reshape(len(q), count)
        support = _guard_support(self.xy, self.values, self.control_spacing, distance, indices)
        lower, upper = support[:, 0].copy(), support[:, 1].copy()
        safe, weight = support[:, 2].copy(), support[:, 3].copy()
        raw = prediction.z - self.zoffset
        safe_neighbors = np.full(len(q), count, dtype=int)
        safe_patches = np.zeros(len(q), dtype=int)
        tolerance = 1e-10 * np.maximum(1.0, np.maximum(np.abs(lower), np.abs(upper)))
        if self.linear is not None:
            simplex = containing_simplex(self.linear.tri, q)
            inside = simplex >= 0
            ids = np.flatnonzero(inside)
            if len(ids):
                vertices = self.linear.tri.simplices[simplex[ids]]
                vertex_values = self.values[vertices]
                lower[ids], upper[ids] = vertex_values.min(axis=1), vertex_values.max(axis=1)
                estimate = np.asarray(self.linear(q[ids])).ravel()
                safe[ids] = np.where(np.isfinite(estimate), estimate, safe[ids])
                safe_neighbors[ids], safe_patches[ids] = 3, 1
            # TIN methods already use nearest outside their hull.
            weight[~inside] = 0.0
            safe[~inside] = raw[~inside]
        else:
            risk = (weight > 1e-8) | (raw < lower - tolerance) | (raw > upper + tolerance)
            candidates = np.flatnonzero(risk)
            if len(candidates):
                outside = self.hull.classify(q[candidates])[0]
                exterior = candidates[outside]
                safe[exterior] = self.values[self._nearest_indices(q[exterior])]
                weight[exterior], safe_neighbors[exterior] = 1.0, 1
                interior = candidates[~outside]
                if len(interior):
                    extended_count = min(512, len(self.xy))
                    extended = self.tree.query(q[interior], k=extended_count, workers=1)[1]
                    extended = extended.reshape(len(interior), extended_count)
                    for row, selected in enumerate(interior):
                        if row % 64 == 0:
                            check_cancel(cancel)
                        key = tuple(sorted(extended[row].tolist()))
                        if key not in self.guard_cache:
                            chosen = np.asarray(key, dtype=int)
                            try:
                                model = LinearNDInterpolator(self.xy[chosen], self.values[chosen],
                                                             fill_value=np.nan)
                            except QhullError:
                                model = None
                            self.guard_cache[key] = model
                            if len(self.guard_cache) > 32:
                                self.guard_cache.popitem(last=False)
                        else:
                            self.guard_cache.move_to_end(key)
                            model = self.guard_cache[key]
                        if model is not None:
                            simplex = int(containing_simplex(model.tri, q[selected:selected + 1])[0])
                            if simplex >= 0:
                                result = float(np.asarray(model(q[selected:selected + 1])).ravel()[0])
                                if math.isfinite(result):
                                    vertex_values = model.values[model.tri.simplices[simplex], 0]
                                    lower[selected], upper[selected] = vertex_values.min(), vertex_values.max()
                                    safe[selected] = result
                                    safe_neighbors[selected], safe_patches[selected] = 3, 1
                        # If the expanded controls cannot enclose the node, the
                        # positive-weight local estimate stays inside its range.
        outside_range = (raw < lower - tolerance) | (raw > upper + tolerance)
        weight[outside_range] = 1.0
        corrected = np.clip((1.0 - weight) * raw + weight * safe, lower, upper)
        changed = np.abs(corrected - raw) > tolerance
        prediction.z = corrected + self.zoffset
        prediction.clamped |= changed
        prediction.fallback |= changed
        replaced = changed & (weight >= 1.0 - 1e-8)
        prediction.neighbors[replaced] = safe_neighbors[replaced]
        prediction.patches[replaced] = safe_patches[replaced]
        return prediction


def _domain_mask(xy: np.ndarray, boundary: Boundary | None, grid: Grid) -> np.ndarray:
    valid = grid.contains(xy)
    if boundary is not None:
        valid &= boundary.contains(xy)
    return valid


def predict_smoothed(surface: Surface, xy: np.ndarray, boundary: Boundary | None,
                     grid: Grid, cancel: threading.Event | None = None) -> Prediction:
    """The voxel-mode Gaussian operator, also used on independent holdout points.

    Uses the same finite Gaussian stencil and normalized domain mask as the
    tiled grid exporter. Other methods and sigma=0 return the base prediction.
    """
    base = surface.predict(xy, cancel)
    sigma = surface.settings.gaussian_sigma
    if surface.settings.method != "voxel" or sigma <= 0:
        return base
    radius = int(4.0 * sigma + 0.5)
    offsets = np.arange(-radius, radius + 1)
    weight = np.exp(-0.5 * (offsets / sigma) ** 2)
    total = np.zeros(len(xy))
    norm = np.zeros(len(xy))
    fallback = base.fallback.copy()
    clamped = base.clamped.copy()
    for i, dx in enumerate(offsets):
        for j, dy in enumerate(offsets):
            q = xy + np.array([dx, dy]) * grid.spacing
            # Equivalent grid stencils use exactly the same coordinates, even
            # when floating-point additions would otherwise differ by one bit.
            origin = np.array([grid.xmin, grid.ymin])
            node = (q - origin) / grid.spacing
            rounded = np.rint(node)
            on_grid = np.all(np.abs(node - rounded) < 1e-8, axis=1)
            q[on_grid] = origin + rounded[on_grid] * grid.spacing
            valid = _domain_mask(q, boundary, grid)
            if np.any(valid):
                value = surface.predict(q[valid], cancel)
                w = weight[i] * weight[j]
                total[valid] += w * value.z
                norm[valid] += w
                fallback[valid] |= value.fallback
                clamped[valid] |= value.clamped
    usable = norm > 0
    base.z[usable] = total[usable] / norm[usable]
    base.fallback = fallback
    base.clamped = clamped
    return surface.protect(xy, base, cancel)


def tune_and_validate(points: np.ndarray, settings: Settings,
                      boundary: Boundary | None, grid: Grid,
                      cancel: threading.Event | None,
                      log: Callable[[str], None]) -> tuple[Settings, Validation]:
    n = len(points)
    count = min(settings.validation_samples, n // 10)
    # Always resolve auto before constructing a numerical interpolator.
    resolved = replace(settings, smoothing=0.1) if settings.smoothing == "auto" else settings
    if count < 2 or n < 12:
        log("Too few independent XY controls for predictive holdout; auto smoothing uses 0.1.")
        return resolved, Validation.empty()
    permutation = np.random.default_rng(1729).permutation(n)
    test = permutation[:count]
    train_mask = np.ones(n, dtype=bool)
    train_mask[test] = False
    if settings.method == "rbf" and settings.smoothing == "auto":
        tuning = permutation[count:2 * count]
        tuning_mask = train_mask.copy()
        tuning_mask[tuning] = False
        best = math.inf
        for smoothing in (0.0, 0.01, 0.1, 1.0):
            check_cancel(cancel)
            candidate = replace(settings, smoothing=smoothing)
            model = Surface(points[tuning_mask], candidate)
            prediction = model.predict(points[tuning, :2], cancel)
            errors = prediction.z - points[tuning, 2]
            metrics = Validation(points[tuning, :2], errors)
            score = metrics.rmse + 0.35 * metrics.p95 + 0.15 * metrics.mae
            log(f"RBF tuning holdout: smoothing={smoothing:g}, RMSE={metrics.rmse:.8g}, "
                f"score={score:.8g}, n={len(tuning)}.")
            # Prefer modest regularization for differences below numeric significance.
            if score < best - 1e-10 or (abs(score - best) <= 1e-10 and smoothing == 0.1):
                best, resolved = score, candidate
        log(f"Selected RBF smoothing={resolved.smoothing}; final validation uses a separate holdout.")
    check_cancel(cancel)
    training = Surface(points[train_mask], resolved)
    prediction = predict_smoothed(training, points[test, :2], boundary, grid, cancel)
    validation = Validation(points[test, :2], prediction.z - points[test, 2],
                            int(np.count_nonzero(prediction.fallback)))
    log(validation.description())
    return resolved, validation


def grid_blocks(grid: Grid, batch_size: int) -> Iterator[tuple[int, int, int, int]]:
    """X-major order, including grids whose Y dimension exceeds one batch."""
    columns = min(grid.ny, batch_size)
    rows = max(1, batch_size // columns)
    for i0 in range(0, grid.nx, rows):
        for j0 in range(0, grid.ny, columns):
            yield i0, min(i0 + rows, grid.nx), j0, min(j0 + columns, grid.ny)


def evaluate_grid_block(surface: Surface, grid: Grid, boundary: Boundary | None,
                        block: tuple[int, int, int, int],
                        cancel: threading.Event | None) -> tuple[np.ndarray, Prediction]:
    i0, i1, j0, j1 = block
    xy = grid.coordinates(i0, i1, j0, j1)
    valid = np.ones(len(xy), dtype=bool) if boundary is None else boundary.contains(xy)
    if not np.any(valid):
        return xy[valid], surface.predict(xy[valid], cancel)
    if surface.settings.method != "voxel" or surface.settings.gaussian_sigma <= 0:
        return xy[valid], surface.predict(xy[valid], cancel)
    sigma = surface.settings.gaussian_sigma
    radius = int(4.0 * sigma + 0.5)
    hi0, hi1 = max(0, i0 - radius), min(grid.nx, i1 + radius)
    hj0, hj1 = max(0, j0 - radius), min(grid.ny, j1 + radius)
    halo_xy = grid.coordinates(hi0, hi1, hj0, hj1)
    halo_valid = np.ones(len(halo_xy), dtype=bool) if boundary is None else boundary.contains(halo_xy)
    shape = (hi1 - hi0, hj1 - hj0)
    z = np.zeros(len(halo_xy))
    fallbacks = np.zeros(len(halo_xy), dtype=bool)
    clamps = np.zeros(len(halo_xy), dtype=bool)
    neighbors = np.zeros(len(halo_xy), dtype=int)
    patches = np.zeros(len(halo_xy), dtype=int)
    ids = np.flatnonzero(halo_valid)
    for start in range(0, len(ids), surface.settings.batch_size):
        part = ids[start:start + surface.settings.batch_size]
        prediction = surface.predict(halo_xy[part], cancel)
        z[part], fallbacks[part] = prediction.z, prediction.fallback
        clamps[part] = prediction.clamped
        neighbors[part], patches[part] = prediction.neighbors, prediction.patches
    numerator = ndimage.gaussian_filter(z.reshape(shape), sigma, mode="constant", cval=0,
                                        radius=radius)
    denominator = ndimage.gaussian_filter(halo_valid.reshape(shape).astype(float), sigma,
                                          mode="constant", cval=0, radius=radius)
    smooth = np.divide(numerator, denominator, out=np.zeros_like(numerator), where=denominator > 0)
    propagated = ndimage.maximum_filter(fallbacks.reshape(shape), size=2 * radius + 1,
                                         mode="constant", cval=0)
    propagated_clamps = ndimage.maximum_filter(clamps.reshape(shape), size=2 * radius + 1,
                                               mode="constant", cval=0)
    cut = (slice(i0 - hi0, i1 - hi0), slice(j0 - hj0, j1 - hj0))
    prediction = Prediction(
        smooth[cut].ravel()[valid], propagated[cut].ravel()[valid],
        neighbors.reshape(shape)[cut].ravel()[valid], patches.reshape(shape)[cut].ravel()[valid],
        propagated_clamps[cut].ravel()[valid])
    return xy[valid], surface.protect(xy[valid], prediction, cancel)


class Diagnostics:
    def __init__(self, surface: Surface, validation: Validation, grid: Grid):
        self.surface, self.validation, self.grid = surface, validation, grid
        self.tree = cKDTree((validation.xy - surface.origin) / surface.scale) if validation.count else None
        self.zscale = max(float(np.std(surface.values)), 1e-8)

    def evaluate(self, xy: np.ndarray, prediction: Prediction) -> np.ndarray:
        q = (xy - self.surface.origin) / self.surface.scale
        nearest = self.surface.tree.query(q, k=1, workers=1)[0]
        outside, distance = self.surface.hull.classify(q)
        penalty = np.exp(-distance * self.surface.scale /
                         max(self.grid.spacing, self.surface.scale))
        if self.tree is not None:
            d, ids = self.tree.query(q, k=min(8, self.validation.count), workers=1)
            if d.ndim == 1:
                d, ids = d[:, None], ids[:, None]
            weights = 1.0 / np.maximum(d, 0.05) ** 2
            residual = np.sqrt(np.sum(weights * self.validation.errors[ids] ** 2, axis=1)
                               / weights.sum(axis=1))
            error_scale = max(self.validation.rmse, 0.01 * self.zscale, 1e-8)
            error_quality = 1.0 / (1.0 + residual / error_scale)
        else:
            # Unknown predictive error is deliberately not represented as zero.
            residual = np.full(len(xy), np.nan)
            error_quality = 0.25
        confidence = np.exp(-nearest / 2.0) * penalty * error_quality
        confidence *= np.where(prediction.fallback, 0.35, 1.0)
        if self.surface.hull.rank < 2:
            confidence *= 0.25
        confidence = np.clip(confidence, 0.0, 1.0)
        return np.column_stack((
            xy, prediction.z, confidence, residual,
            np.full(len(xy), self.surface.smoothing), penalty,
            prediction.fallback.astype(int), outside.astype(int), prediction.clamped.astype(int),
            prediction.neighbors, prediction.patches))


class DXFStream:
    """Minimal R12 point/label drawing using ezdxf's streaming tag serializer.

    Explicit fixed-decimal values preserve the requested 0–15 places instead
    of the six-place rounding used by ezdxf.addons.r12writer for some entities.
    The drawing is never accumulated in memory.
    """
    LAYERS = (("xyz_points", 7), ("xyz_labels", 3), ("boundary", 1),
              ("grid_points", 5), ("grid_labels", 4))

    def __init__(self, stream, settings: Settings, grid: Grid, points: np.ndarray,
                 boundary: Boundary | None):
        self.writer = TagWriter(stream, dxfversion="AC1009", write_handles=False)
        self.precision = settings.precision
        xmin = min(float(points[:, 0].min()), grid.xmin)
        xmax = max(float(points[:, 0].max()), grid.xmin + (grid.nx - 1) * grid.spacing)
        ymin = min(float(points[:, 1].min()), grid.ymin)
        ymax = max(float(points[:, 1].max()), grid.ymin + (grid.ny - 1) * grid.spacing)
        if boundary is not None:
            xmin, ymin = np.minimum([xmin, ymin], boundary.loop.min(axis=0))
            xmax, ymax = np.maximum([xmax, ymax], boundary.loop.max(axis=0))
        centre = (0.5 * (xmin + xmax), 0.5 * (ymin + ymax))
        view_size = max(1.0, max(xmax - xmin, ymax - ymin) * 1.1)
        self.tags((0, "SECTION"), (2, "HEADER"), (9, "$ACADVER"), (1, "AC1009"),
                  (9, "$PDMODE"), (70, settings.pdmode), (9, "$PDSIZE"), (40, "0.5"),
                  (0, "ENDSEC"), (0, "SECTION"), (2, "TABLES"),
                  (0, "TABLE"), (2, "LTYPE"), (70, 1),
                  (0, "LTYPE"), (2, "CONTINUOUS"), (70, 0), (3, "Solid line"),
                  (72, 65), (73, 0), (40, "0.0"), (0, "ENDTAB"),
                  (0, "TABLE"), (2, "LAYER"), (70, len(self.LAYERS) + 1))
        for layer, color in (("0", 7), *self.LAYERS):
            self.tags((0, "LAYER"), (2, layer), (70, 0), (62, color), (6, "CONTINUOUS"))
        self.tags((0, "ENDTAB"), (0, "TABLE"), (2, "STYLE"), (70, 1),
                  (0, "STYLE"), (2, "STANDARD"), (70, 0), (40, "0.0"), (41, "1.0"),
                  (50, "0.0"), (71, 0), (42, "1.0"), (3, "txt"), (4, ""), (0, "ENDTAB"),
                  (0, "TABLE"), (2, "VPORT"), (70, 1), (0, "VPORT"), (2, "*ACTIVE"),
                  (70, 0), (10, "0.0"), (20, "0.0"), (11, "1.0"), (21, "1.0"),
                  (12, self.fmt(centre[0])), (22, self.fmt(centre[1])), (40, self.fmt(view_size)),
                  (41, "1.0"), (0, "ENDTAB"), (0, "ENDSEC"),
                  (0, "SECTION"), (2, "ENTITIES"))
        if boundary is not None:
            self.tags((0, "POLYLINE"), (8, "boundary"), (66, 1), (70, 1),
                      (10, "0.0"), (20, "0.0"), (30, "0.0"))
            for x, y in boundary.loop:
                self.tags((0, "VERTEX"), (8, "boundary"), (10, self.fmt(x)),
                          (20, self.fmt(y)), (30, "0.0"), (70, 0))
            self.tags((0, "SEQEND"), (8, "boundary"))

    def tags(self, *tags: tuple[int, object]) -> None:
        for code, value in tags:
            self.writer.write_tag2(code, value)

    def fmt(self, value: float) -> str:
        return f"{float(value):.{self.precision}f}"

    def points(self, xyz: np.ndarray, grid: bool = False) -> None:
        point_layer, label_layer = ("grid_points", "grid_labels") if grid else ("xyz_points", "xyz_labels")
        precision = self.precision
        point_prefix = f"  0\nPOINT\n  8\n{point_layer}\n"
        label_prefix = f"  0\nTEXT\n  8\n{label_layer}\n"
        records = []
        for x, y, z in xyz:
            x, y, z = float(x), float(y), float(z)
            sx, sy, sz = f"{x:.{precision}f}", f"{y:.{precision}f}", f"{z:.{precision}f}"
            lx, ly = f"{x + 0.2:.{precision}f}", f"{y + 0.2:.{precision}f}"
            records.append(f"{point_prefix} 10\n{sx}\n 20\n{sy}\n 30\n{sz}\n"
                           f"{label_prefix} 10\n{lx}\n 20\n{ly}\n 30\n0.0\n 40\n1.0\n  1\n{sz}\n")
        self.writer.write_str("".join(records))

    def close(self) -> None:
        self.tags((0, "ENDSEC"), (0, "EOF"))


def commit_outputs(staged: list[Path], destinations: list[Path], staging: Path) -> None:
    """Replace products only after success; roll back a failed multi-file commit."""
    backups: dict[Path, Path] = {}
    installed: list[Path] = []
    try:
        for index, (source, destination) in enumerate(zip(staged, destinations)):
            if destination.exists():
                backup = staging / f"previous-{index}"
                os.replace(destination, backup)
                backups[destination] = backup
            os.replace(source, destination)
            installed.append(destination)
    except OSError:
        for destination in reversed(destinations):
            if destination in installed:
                destination.unlink(missing_ok=True)
            if destination in backups:
                os.replace(backups[destination], destination)
        raise


@dataclass
class RunResult:
    outputs: tuple[Path, ...]
    input_points: int
    filtered_points: int
    grid_points: int
    validation: Validation
    elapsed: float


def process_xyz_to_dxf(settings: Settings, progress: Progress | None = None,
                        cancel: threading.Event | None = None) -> RunResult:
    """Run the complete workflow; GUI and CLI share this exact processing path."""
    settings.validate()
    source = Path(settings.input_file).expanduser().resolve()
    settings = replace(settings, input_file=source)
    destinations = [Path(str(source) + suffix) for suffix in SUFFIXES]
    started = time.monotonic()
    # Same filesystem permits atomic renaming. No output data file is published
    # until all computations and writes have succeeded.
    with tempfile.TemporaryDirectory(prefix=".xyz2dxf-", dir=source.parent) as directory:
        staging = Path(directory)
        staged = [staging / f"result{suffix}" for suffix in SUFFIXES]
        reporter = Reporter(staged[-1], progress)
        report_closed = False
        try:
            with threadpool_limits(limits=1, user_api="blas"):
                log = reporter.emit
                log(APP_TITLE + " — Tkinter / SciPy")
                log(f"Python {sys.version.split()[0]}; NumPy {np.__version__}; "
                    f"SciPy {scipy.__version__}; ezdxf {ezdxf.__version__}; Numba {numba.__version__}.")
                configure_acceleration()
                log(f"Numba CPU acceleration: up to {CPU_THREADS} threads; strict float64 arithmetic; "
                    "first-use compilation is cached for subsequent runs.")
                log(f"Input: {source}")
                log(f"Settings: Min Dist={settings.min_dist:g}, Precision={settings.precision}, "
                    f"PDMODE={settings.pdmode}, Grid Spacing={settings.grid_spacing:g}, "
                    f"outliers={settings.use_outliers}, boundary={settings.use_boundary}, "
                    f"method={METHODS[settings.method]}, overshoot protection={settings.prevent_overshoot}.")
                if settings.method == "rbf":
                    log(f"RBF: kernel={settings.kernel}, neighbours={settings.neighbors}, "
                        f"smoothing={settings.smoothing}, polynomial degree=1.")
                if settings.method == "voxel":
                    log(f"Voxel: cell size={settings.voxel_size}, Gaussian sigma={settings.gaussian_sigma:g} grid cells.")
                log("Phase 1/7: reading original XYZ input ...", 0.01)
                points, read_stats = read_xyz(source, cancel, log)
                log(f"Input rows: {read_stats.lines:,}; valid={read_stats.accepted:,}; "
                    f"ignored/header={read_stats.ignored:,}; non-finite={read_stats.nonfinite:,}.")
                check_cancel(cancel)
                log("Phase 2/7: robust local residual screening on original cloud ...", 0.08)
                if settings.use_outliers:
                    before = len(points)
                    points = remove_outliers(points, settings.min_dist, cancel, log)
                    log(f"After outlier removal: {len(points):,}; removed={before - len(points):,}.")
                else:
                    log("Robust outlier removal disabled.")
                log("Phase 3/7: deterministic minimum-distance thinning ...", 0.22)
                points = thin_points(points, settings.min_dist, cancel, log)
                log(f"After minimum-distance filter: {len(points):,} points.")
                grid = Grid.from_points(points, settings.grid_spacing)
                log(f"Grid planned: {grid.nx:,} x {grid.ny:,} = {grid.total:,} nodes; "
                    "1.5-spacing margin; X-major / Y-minor order.")
                check_cancel(cancel)
                log("Phase 4/7: constructing outer boundary ...", 0.30)
                boundary = build_boundary(points, settings, cancel, log) if settings.use_boundary else None
                if boundary is None:
                    log("Outer boundary disabled: all rectangular-grid nodes are kept.")
                else:
                    log(f"Outer boundary: one closed loop, {len(boundary.loop):,} vertices; all filtered points contained.")
                log("Phase 5/7: writing filtered.xyz ...", 0.38)
                fmt = f"%.{settings.precision}f"
                with staged[0].open("w", encoding="ascii", newline="\n", buffering=1024 * 1024) as stream:
                    for start in range(0, len(points), settings.batch_size):
                        check_cancel(cancel)
                        np.savetxt(stream, points[start:start + settings.batch_size], fmt=fmt)
                log("Phase 6/7: fitting interpolator and predictive validation ...", 0.42)
                controls, merged = unique_xy_controls(points)
                if merged:
                    log(f"Interpolation-only duplicate XY merge: {merged:,} observations combined using median Z; "
                        "filtered.xyz remains the retained observation cloud.")
                tuned, validation = tune_and_validate(controls, settings, boundary, grid, cancel, log)
                if not validation.count:
                    log(validation.description())
                check_cancel(cancel)
                surface = Surface(controls, tuned, log)
                diagnostics = Diagnostics(surface, validation, grid)
                log("Confidence is a support/error diagnostic, not a probability or confidence interval.")
                log("Residual is an inverse-distance weighted RMS of up to 8 independent holdout errors; "
                    "nan means validation unavailable. lambda is normalized-coordinate RBF smoothing (0 otherwise).")
                log("outsideHull uses the interpolator-control convex support; boundaryPenalty=exp(-outsideDistance/"
                    "max(gridSpacing, normalizationScale)). Fallbacks reduce confidence by a factor of 0.35.")
                log("clampedLocal=1 marks a range-limited or sparse-support correction. "
                    "Safeguards use bounded linear triangles, or positive weights if no triangle encloses the node. "
                    "Neighbour/patch counts describe the base fit; Gaussian flags include the stencil.")
                log("Phase 6/7: streaming grid, confidence and DXF entities ...", 0.52)
                count = fallback_count = outside_count = clamped_count = 0
                confidence_sum = residual_sum = 0.0
                residual_count = 0
                last_update = time.monotonic()
                with (staged[1].open("w", encoding="ascii", newline="\n", buffering=1024 * 1024) as grid_stream,
                      staged[2].open("w", encoding="ascii", newline="\n", buffering=1024 * 1024) as conf_stream,
                      staged[3].open("w", encoding="ascii", newline="\n", buffering=1024 * 1024) as dxf_stream):
                    conf_stream.write(CONFIDENCE_HEADER)
                    dxf = DXFStream(dxf_stream, tuned, grid, points, boundary)
                    # The C++ DXF source points are read back from filtered.xyz.
                    # Preserve that rounding behaviour without holding another cloud.
                    with staged[0].open("r", encoding="ascii") as filtered:
                        block_values: list[list[float]] = []
                        for line in filtered:
                            block_values.append([float(value) for value in line.split()])
                            if len(block_values) >= settings.batch_size:
                                check_cancel(cancel)
                                dxf.points(np.asarray(block_values))
                                block_values.clear()
                        if block_values:
                            dxf.points(np.asarray(block_values))
                    for block in grid_blocks(grid, settings.batch_size):
                        check_cancel(cancel)
                        xy, prediction = evaluate_grid_block(surface, grid, boundary, block, cancel)
                        if len(xy):
                            xyz = np.column_stack((xy, prediction.z))
                            np.savetxt(grid_stream, xyz, fmt=fmt)
                            table = diagnostics.evaluate(xy, prediction)
                            np.savetxt(conf_stream, table, fmt=[fmt] * 7 + ["%d"] * 5)
                            dxf.points(xyz, grid=True)
                            count += len(xy)
                            confidence_sum += float(table[:, 3].sum())
                            good_residual = np.isfinite(table[:, 4])
                            residual_sum += float(table[good_residual, 4].sum())
                            residual_count += int(good_residual.sum())
                            fallback_count += int(np.count_nonzero(prediction.fallback))
                            outside_count += int(np.count_nonzero(table[:, 8]))
                            clamped_count += int(np.count_nonzero(prediction.clamped))
                        i0, i1, j0, j1 = block
                        completed = (i1 - 1) * grid.ny + j1
                        if time.monotonic() - last_update > 0.7 or completed == grid.total:
                            log(f"Grid processed: {completed:,}/{grid.total:,}; exported nodes={count:,}.",
                                0.52 + 0.43 * completed / grid.total)
                            last_update = time.monotonic()
                    if count == 0:
                        raise ValueError("No regular-grid node falls inside the outer boundary. "
                                         "Reduce Grid Spacing or disable Use outer boundary.")
                    dxf.close()
                check_cancel(cancel)
                log("Phase 7/7: finalizing DXF and publishing completed products ...", 0.97)
                mean_residual = residual_sum / residual_count if residual_count else math.nan
                log(f"Confidence summary: nodes={count:,}, meanConfidence={confidence_sum / count:.8g}, "
                    f"meanResidual={mean_residual:.8g}, fallbacks={fallback_count:,}, "
                    f"outsideHull={outside_count:,}, clampedLocal={clamped_count:,}.")
                log(f"DXF: {len(points):,} source POINT/TEXT pairs; {count:,} grid POINT/TEXT pairs; "
                    f"boundary polylines={int(boundary is not None)}; PDSIZE=0.5.")
                log("Predictive validation is conditional on the cleaned cloud; it is not external survey "
                    "verification or a spatial-block extrapolation test.")
                elapsed = time.monotonic() - started
                log(f"Done. Total time: {elapsed:.2f} sec.")
                reporter.stream.write(f"Finished in {elapsed / 60:.2f} minutes.\n")
                reporter.close()
                report_closed = True
                commit_outputs(staged, destinations, staging)
                if progress:
                    progress(f"Done. {count:,} grid nodes exported in {elapsed:.2f} sec.", 1.0)
                return RunResult(tuple(destinations), read_stats.accepted, len(points), count, validation, elapsed)
        except Exception as exc:
            message = ("Insufficient memory. Increase Min Dist or select local RBF / voxel reduction."
                       if isinstance(exc, MemoryError) else str(exc))
            if not report_closed:
                reporter.emit(f"{'Cancelled' if isinstance(exc, CancelledError) else 'Error'}: {message}")
                reporter.emit("Incomplete new data products were discarded; previous data outputs were preserved.")
                reporter.close()
                report_closed = True
                os.replace(staged[-1], destinations[-1])
            raise
        finally:
            if not report_closed:
                reporter.close()


class ConverterApp:
    """Fixed-size Tkinter interface; only the main thread accesses Tk widgets."""

    def __init__(self, root):
        import tkinter as tk
        from tkinter import ttk
        from tkinter.scrolledtext import ScrolledText
        import tkinter.font as tkfont

        self.root, self.tk, self.ttk = root, tk, ttk
        self.events: queue.SimpleQueue = queue.SimpleQueue()
        self.cancel_event = threading.Event()
        self.worker: threading.Thread | None = None
        self.busy = False
        self.closing = False
        self.destroyed = False
        self.input_widgets = []
        root.title(APP_TITLE)
        root.resizable(False, False)
        root.protocol("WM_DELETE_WINDOW", self.close)
        background, panel = "#eaf0f7", "#ffffff"
        text, secondary, border = "#101827", "#27364b", "#b5c4d6"
        accent, selected = "#174bc4", "#123c9d"
        family = "Segoe UI" if sys.platform == "win32" else "DejaVu Sans"
        for name in ("TkDefaultFont", "TkTextFont", "TkMenuFont"):
            tkfont.nametofont(name).configure(family=family, size=11)
        root.configure(background=background)
        root.columnconfigure(0, weight=1)
        root.rowconfigure(0, weight=1)
        style = ttk.Style(root)
        if "clam" in style.theme_names():
            style.theme_use("clam")
        style.configure(".", font=(family, 11), foreground=text, background=panel)
        style.configure("App.TFrame", background=background)
        style.configure("Card.TLabelframe", background=panel, bordercolor=border,
                        lightcolor=panel, darkcolor=panel, borderwidth=1, relief="solid")
        style.configure("Card.TLabelframe.Label", font=(family, 12, "bold"),
                        foreground="#0c223f", background=panel)
        style.configure("TLabel", background=panel, foreground=text)
        style.configure("Hint.TLabel", foreground=secondary, background=panel, font=(family, 11))
        style.configure("Bold.TLabel", foreground="#0c223f", background=panel, font=(family, 11, "bold"))
        style.configure("TEntry", foreground=text, fieldbackground=panel, bordercolor=border,
                        lightcolor=panel, darkcolor=panel, padding=(6, 5))
        style.map("TEntry", foreground=[("disabled", "#475569")],
                  fieldbackground=[("disabled", "#f1f5f9")], bordercolor=[("focus", accent)])
        style.configure("TCombobox", foreground=text, fieldbackground=panel, background="#e3ebf5",
                        arrowcolor="#0c223f", bordercolor=border, padding=(5, 4))
        style.map("TCombobox", foreground=[("disabled", "#475569"), ("readonly", text)],
                  fieldbackground=[("disabled", "#f1f5f9"), ("readonly", panel)],
                  selectforeground=[("readonly", text)], selectbackground=[("readonly", panel)])
        for name in ("TCheckbutton", "TRadiobutton"):
            style.configure(name, foreground=text, background=panel, padding=(0, 2),
                            indicatorbackground=panel, indicatorforeground=accent)
            style.map(name, foreground=[("disabled", "#475569"), ("active", text)],
                      background=[("active", panel)], indicatorbackground=[("selected", accent)],
                      indicatorforeground=[("selected", panel), ("!selected", text)])
        style.configure("TButton", foreground="#10243e", background="#e2eaf5", bordercolor=border,
                        padding=(16, 7), font=(family, 11, "bold"))
        style.map("TButton", foreground=[("disabled", "#475569")],
                  background=[("disabled", "#e8edf3"), ("active", "#cbd9eb")])
        style.configure("Accent.TButton", foreground=panel, background=accent, bordercolor=accent)
        style.map("Accent.TButton", foreground=[("disabled", "#475569"), ("!disabled", panel)],
                  background=[("disabled", "#d8e1ed"), ("pressed", selected), ("active", selected)])
        style.configure("Horizontal.TProgressbar", background=accent, troughcolor="#e3ebf5",
                        bordercolor=border, lightcolor=accent, darkcolor=accent, thickness=12)
        self.path = tk.StringVar()
        self.min_dist = tk.StringVar(value="1")
        self.precision = tk.StringVar(value="2")
        self.pdmode = tk.StringVar(value="3")
        self.grid_spacing = tk.StringVar(value="10")
        self.use_outliers = tk.BooleanVar(value=False)
        self.use_boundary = tk.BooleanVar(value=True)
        self.prevent_overshoot = tk.BooleanVar(value=True)
        self.method = tk.StringVar(value="rbf")
        self.neighbors = tk.StringVar(value="50")
        self.smoothing = tk.StringVar(value="0.1")
        self.kernel = tk.StringVar(value="thin_plate_spline")
        self.voxel_size = tk.StringVar(value="auto")
        self.sigma = tk.StringVar(value="1.0")
        self.status = tk.StringVar(value="Status: Idle")

        outer = ttk.Frame(root, padding=12, style="App.TFrame")
        outer.grid(sticky="nsew")
        for column, minimum in enumerate((330, 410, 340)):
            outer.columnconfigure(column, weight=1, minsize=minimum)
        outer.rowconfigure(2, weight=1)
        source = ttk.LabelFrame(outer, text="Source file", padding=(12, 6), style="Card.TLabelframe")
        source.grid(row=0, column=0, columnspan=3, sticky="ew", pady=(0, 8))
        source.columnconfigure(1, weight=1)
        ttk.Label(source, text="Input XYZ file:", style="Bold.TLabel").grid(row=0, column=0, padx=(0, 12))
        file_entry = ttk.Entry(source, textvariable=self.path, width=80)
        file_entry.grid(row=0, column=1, sticky="ew")
        self.input_widgets.append(file_entry)
        browse = ttk.Button(source, text="Browse...", command=self.browse)
        browse.grid(row=0, column=2, padx=(12, 0))
        self.input_widgets.append(browse)
        ttk.Label(source, text="The five output files are saved beside the source XYZ file.",
                  style="Hint.TLabel").grid(row=1, column=1, sticky="w", pady=(5, 0))

        parameters = ttk.LabelFrame(outer, text="Processing parameters", padding=(12, 4), style="Card.TLabelframe")
        parameters.grid(row=1, column=0, sticky="nsew", padx=(0, 12))
        parameters.columnconfigure(1, weight=1)
        for row, label, variable, hint in (
            (0, "Min Dist:", self.min_dist, "XY spacing; 0 keeps distinct XYZ points."),
            (1, "Precision:", self.precision, "Output decimal places, from 0 to 15."),
            (2, "PDMODE:", self.pdmode, "DXF point style. PDSIZE = 0.5."),
            (3, "Grid Spacing:", self.grid_spacing, "Grid increment in coordinate units."),
        ):
            ttk.Label(parameters, text=label, style="Bold.TLabel").grid(
                row=row * 2, column=0, sticky="w", padx=(0, 12))
            entry = ttk.Entry(parameters, textvariable=variable, width=10)
            entry.grid(row=row * 2, column=1, sticky="ew")
            self.input_widgets.append(entry)
            ttk.Label(parameters, text=hint, style="Hint.TLabel", wraplength=300).grid(
                row=row * 2 + 1, column=0, columnspan=2, sticky="w", pady=(2, 4))
        for row, label, variable, hint in (
            (8, "Use robust outlier removal", self.use_outliers,
             "Screen the original cloud before thinning."),
            (10, "Use outer boundary", self.use_boundary,
             "Clip the grid to the survey footprint and export its boundary."),
        ):
            widget = ttk.Checkbutton(parameters, text=label, variable=variable)
            widget.grid(row=row, column=0, columnspan=2, sticky="w", pady=(3, 0))
            self.input_widgets.append(widget)
            ttk.Label(parameters, text=hint, style="Hint.TLabel", wraplength=285).grid(
                row=row + 1, column=0, columnspan=2, sticky="w", padx=(20, 0), pady=(2, 4))

        methods = ttk.LabelFrame(outer, text="Interpolation method", padding=(12, 4), style="Card.TLabelframe")
        methods.grid(row=1, column=1, sticky="nsew", padx=(0, 12))
        hints = {
            "rbf": "Smooth local RBF with bounded sparse-area safeguards.",
            "cubic": "Smooth cubic reconstruction inside the TIN.",
            "linear": "Planar triangles; values stay within vertex limits.",
            "voxel": "Reduce dense cells, fit cubic, then smooth the grid.",
        }
        for index, (key, label) in enumerate(METHODS.items()):
            widget = ttk.Radiobutton(methods, text=label, variable=self.method, value=key,
                                      command=self.update_options)
            widget.grid(row=2 * index, column=0, sticky="w", pady=(3, 0))
            self.input_widgets.append(widget)
            ttk.Label(methods, text=hints[key], style="Hint.TLabel", wraplength=365).grid(
                row=2 * index + 1, column=0, sticky="w", padx=(20, 0), pady=(1, 8))
        guard = ttk.Checkbutton(methods, text="Prevent local overshoot", variable=self.prevent_overshoot)
        guard.grid(row=10, column=0, sticky="w", pady=(6, 2))
        self.input_widgets.append(guard)
        ttk.Label(methods, text="Use bounded fits in sparse areas. Corrections are flagged "
                  "in the confidence file.", style="Hint.TLabel",
                  wraplength=385).grid(row=11, column=0, sticky="w", padx=(20, 0))

        options = ttk.LabelFrame(outer, text="Interpolation settings", padding=(12, 4), style="Card.TLabelframe")
        options.grid(row=1, column=2, sticky="nsew")
        options.columnconfigure(1, weight=1)
        ttk.Label(options, text="Local RBF", style="Bold.TLabel").grid(
            row=0, column=0, columnspan=2, sticky="w", pady=(0, 9))
        ttk.Label(options, text="Kernel:").grid(row=1, column=0, sticky="w", padx=(0, 10))
        kernel_widget = ttk.Combobox(options, textvariable=self.kernel,
                                     values=("thin_plate_spline", "cubic"), state="readonly", width=19)
        kernel_widget.grid(row=1, column=1, sticky="ew", pady=(0, 8))
        ttk.Label(options, text="Neighbours:").grid(row=2, column=0, sticky="w", padx=(0, 10))
        neighbours_widget = ttk.Entry(options, textvariable=self.neighbors, width=12)
        neighbours_widget.grid(row=2, column=1, sticky="ew", pady=(0, 8))
        ttk.Label(options, text="Smoothing:").grid(row=3, column=0, sticky="w", padx=(0, 10))
        smoothing_widget = ttk.Entry(options, textvariable=self.smoothing, width=12)
        smoothing_widget.grid(row=3, column=1, sticky="ew")
        ttk.Label(options, text="0 = exact interpolation. Use auto for tuning.", style="Hint.TLabel",
                  wraplength=315).grid(row=4, column=0, columnspan=2, sticky="w", pady=(4, 0))
        self.rbf_widgets = [kernel_widget, neighbours_widget, smoothing_widget]
        ttk.Separator(options).grid(row=5, column=0, columnspan=2, sticky="ew", pady=12)
        ttk.Label(options, text="Voxel smoothing", style="Bold.TLabel").grid(
            row=6, column=0, columnspan=2, sticky="w", pady=(0, 9))
        ttk.Label(options, text="Cell size:").grid(row=7, column=0, sticky="w", padx=(0, 10))
        voxel_widget = ttk.Entry(options, textvariable=self.voxel_size, width=12)
        voxel_widget.grid(row=7, column=1, sticky="ew", pady=(0, 8))
        ttk.Label(options, text="Gaussian sigma:").grid(row=8, column=0, sticky="w", padx=(0, 10))
        sigma_widget = ttk.Entry(options, textvariable=self.sigma, width=12)
        sigma_widget.grid(row=8, column=1, sticky="ew")
        ttk.Label(options, text="Sigma is in grid cells; 0 disables smoothing. "
                  "Auto cell size = Grid Spacing / 2.", style="Hint.TLabel", wraplength=315).grid(
            row=9, column=0, columnspan=2, sticky="w", pady=(4, 0))
        self.voxel_widgets = [voxel_widget, sigma_widget]
        self.input_widgets.extend(self.rbf_widgets + self.voxel_widgets)

        execution = ttk.LabelFrame(outer, text="Execution", padding=(12, 6), style="Card.TLabelframe")
        execution.grid(row=2, column=0, columnspan=3, sticky="nsew", pady=(8, 0))
        execution.columnconfigure(0, weight=1)
        execution.rowconfigure(1, weight=1)
        ttk.Label(execution, text="Processing log", style="Bold.TLabel").grid(row=0, column=0, sticky="w")
        ttk.Label(execution, text="Numba CPU acceleration", style="Hint.TLabel").grid(
            row=0, column=1, sticky="e", padx=(10, 14))
        self.cancel_button = ttk.Button(execution, text="Cancel", command=self.cancel, state="disabled")
        self.cancel_button.grid(row=0, column=2, padx=(0, 10))
        self.run_button = ttk.Button(execution, text="Run Conversion", command=self.run, style="Accent.TButton")
        self.run_button.grid(row=0, column=3, sticky="e")
        self.log_widget = ScrolledText(
            execution, height=5, width=134, wrap="word", foreground=text, background=panel,
            insertbackground=text, selectforeground=panel, selectbackground=accent,
            font=("Consolas" if sys.platform == "win32" else "DejaVu Sans Mono", 11),
            state="disabled", relief="solid", borderwidth=1, highlightthickness=1,
            highlightbackground=border, highlightcolor=accent, padx=6, pady=5)
        self.log_widget.grid(row=1, column=0, columnspan=4, sticky="nsew", pady=(8, 7))
        self.progressbar = ttk.Progressbar(execution, mode="determinate", maximum=100)
        self.progressbar.grid(row=2, column=0, columnspan=4, sticky="ew", pady=(0, 6))
        ttk.Label(execution, textvariable=self.status, foreground="#101827", wraplength=1190).grid(
            row=3, column=0, columnspan=4, sticky="w")
        self.update_options()
        root.bind("<Return>", lambda event: self.run() if not self.busy else None)
        root.after(100, self.poll)
        root.update_idletasks()
        width = max(1400, root.winfo_reqwidth())
        height = max(650, root.winfo_reqheight())
        root.geometry(f"{width}x{height}")
        file_entry.focus_set()


    def update_options(self) -> None:
        for widget in self.rbf_widgets:
            if self.busy or self.method.get() != "rbf":
                widget.configure(state="disabled")
            else:
                widget.configure(state="readonly" if isinstance(widget, self.ttk.Combobox) else "normal")
        for widget in self.voxel_widgets:
            widget.configure(state="normal" if not self.busy and self.method.get() == "voxel" else "disabled")

    def browse(self) -> None:
        from tkinter import filedialog
        path = filedialog.askopenfilename(title="Select source XYZ point cloud", filetypes=[
            ("XYZ point clouds", "*.xyz"), ("Text files", "*.txt"), ("All files", "*.*")])
        if path:
            self.path.set(path)

    @staticmethod
    def _number(text: str) -> float:
        return float(text.strip().replace(",", "."))

    def get_settings(self) -> Settings:
        method = self.method.get()
        smoothing_text = self.smoothing.get().strip().lower()
        voxel_text = self.voxel_size.get().strip().lower()
        try:
            settings = Settings(
                input_file=Path(self.path.get().strip().strip('"')).expanduser(),
                min_dist=self._number(self.min_dist.get()),
                precision=int(self.precision.get().strip()), pdmode=int(self.pdmode.get().strip()),
                grid_spacing=self._number(self.grid_spacing.get()),
                use_outliers=self.use_outliers.get(), use_boundary=self.use_boundary.get(), method=method,
                prevent_overshoot=self.prevent_overshoot.get(),
                neighbors=int(self.neighbors.get().strip()) if method == "rbf" else 50,
                smoothing=("auto" if smoothing_text == "auto" else self._number(smoothing_text)) if method == "rbf" else 0.1,
                kernel=self.kernel.get() if method == "rbf" else "thin_plate_spline",
                voxel_size=("auto" if voxel_text == "auto" else self._number(voxel_text)) if method == "voxel" else "auto",
                gaussian_sigma=self._number(self.sigma.get()) if method == "voxel" else 1.0)
        except ValueError as exc:
            raise ValueError("One or more numeric inputs are invalid. Use whole numbers for "
                             "Precision, PDMODE and Neighbours.") from exc
        settings.validate()
        return settings

    def append_log(self, message: str) -> None:
        self.log_widget.configure(state="normal")
        self.log_widget.insert("end", message + "\n")
        self.log_widget.see("end")
        self.log_widget.configure(state="disabled")

    def set_busy(self, busy: bool) -> None:
        self.busy = busy
        for widget in self.input_widgets:
            widget.configure(state="disabled" if busy else "normal")
        self.run_button.configure(state="disabled" if busy else "normal")
        self.cancel_button.configure(state="normal" if busy else "disabled")
        self.update_options()

    def run(self) -> None:
        from tkinter import messagebox
        if self.busy:
            return
        try:
            settings = self.get_settings()
        except (ValueError, OSError) as exc:
            messagebox.showerror(APP_TITLE, str(exc), parent=self.root)
            return
        self.log_widget.configure(state="normal")
        self.log_widget.delete("1.0", "end")
        self.log_widget.configure(state="disabled")
        self.cancel_event.clear()
        self.set_busy(True)
        self.progressbar["value"] = 0
        self.status.set("Starting conversion…")

        def work() -> None:
            try:
                result = process_xyz_to_dxf(settings, lambda message, fraction:
                    self.events.put(("progress", message, fraction)), self.cancel_event)
                self.events.put(("done", result, None))
            except Exception as exc:
                self.events.put(("cancelled" if isinstance(exc, CancelledError) else "error", str(exc), None))

        self.worker = threading.Thread(target=work, name="xyz2dxf-worker", daemon=False)
        self.worker.start()

    def cancel(self) -> None:
        if self.busy:
            self.cancel_event.set()
            self.cancel_button.configure(state="disabled")
            self.status.set("Cancelling after the current processing batch…")

    def poll(self) -> None:
        from tkinter import messagebox
        try:
            # Bound the work per poll so a backlog cannot starve Tk events.
            for _ in range(200):
                try:
                    kind, payload, fraction = self.events.get_nowait()
                except queue.Empty:
                    break
                if kind == "progress":
                    self.append_log(payload)
                    if not self.cancel_event.is_set():
                        self.status.set(payload)
                    if fraction is not None:
                        self.progressbar["value"] = 100 * fraction
                else:
                    if self.worker:
                        self.worker.join()
                        self.worker = None
                    self.set_busy(False)
                    if self.closing:
                        self.destroyed = True
                        self.root.destroy()
                        return
                    if kind == "done":
                        self.progressbar["value"] = 100
                        self.status.set(f"Done. {payload.filtered_points:,} retained points; "
                                        f"{payload.grid_points:,} grid nodes; {payload.elapsed:.2f} sec.")
                        self.append_log("Files created:\n" + "\n".join(str(p) for p in payload.outputs))
                        self.append_log(f"Finished in {payload.elapsed / 60:.2f} minutes.")
                    else:
                        self.status.set(payload)
                        self.append_log(payload)
                        if kind == "error":
                            messagebox.showerror(APP_TITLE, payload, parent=self.root)
        finally:
            if not self.destroyed:
                self.root.after(100, self.poll)

    def close(self) -> None:
        from tkinter import messagebox
        if not self.busy:
            self.destroyed = True
            self.root.destroy()
        elif messagebox.askyesno(APP_TITLE, "A conversion is running.\n\nCancel it and close?",
                                  default="no", parent=self.root):
            self.closing = True
            self.cancel()


def launch_gui() -> int:
    try:
        import tkinter as tk
    except ImportError:
        print("Tkinter is unavailable. Install Python's Tk support (python3-tk on Debian/Ubuntu).", file=sys.stderr)
        return 1
    if sys.platform == "win32":
        try:
            import ctypes
            ctypes.windll.shcore.SetProcessDpiAwareness(1)
        except (AttributeError, OSError):
            pass
    try:
        root = tk.Tk()
    except tk.TclError as exc:
        print(f"Cannot open the Tkinter display: {exc}\nUse --input FILE for command-line conversion.", file=sys.stderr)
        return 1
    ConverterApp(root)
    root.mainloop()
    return 0


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description="XYZ conditioning, SciPy interpolation and DXF export. No arguments opens the Tkinter GUI.")
    parser.add_argument("--input", type=Path, help="input XYZ file; enables command-line processing")
    parser.add_argument("--min-dist", type=float, default=1.0, help="minimum XY distance (default: 1)")
    parser.add_argument("--precision", type=int, default=2, help="output decimal places, 0–15 (default: 2)")
    parser.add_argument("--pdmode", type=int, default=3, help="DXF point display mode (default: 3)")
    parser.add_argument("--grid-spacing", type=float, default=10.0, help="grid increment (default: 10)")
    parser.add_argument("--method", choices=METHODS, default="rbf", help="interpolation method (default: rbf)")
    parser.add_argument("--outliers", action=argparse.BooleanOptionalAction, default=False,
                        help="enable robust outlier removal (default: disabled)")
    parser.add_argument("--no-boundary", action="store_true", help="keep the complete rectangular grid")
    parser.add_argument("--allow-overshoot", action="store_true", help="disable local range and sparse-support safeguards")
    parser.add_argument("--neighbors", type=int, default=50, help="local RBF neighbour count, 3–512 (default: 50)")
    parser.add_argument("--smoothing", default="0.1", help="nonnegative RBF smoothing or auto (default: 0.1)")
    parser.add_argument("--kernel", choices=("thin_plate_spline", "cubic"), default="thin_plate_spline")
    parser.add_argument("--voxel-size", default="auto", help="positive voxel cell size or auto (grid spacing / 2)")
    parser.add_argument("--gaussian-sigma", type=float, default=1.0, help="voxel-mode smoothing in grid cells, 0–4")
    parser.add_argument("--batch-size", type=int, default=4096, help="maximum base evaluation batch (default: 4096)")
    parser.add_argument("--validation-samples", type=int, default=256, help="maximum holdout count; 0 disables validation")
    arguments = sys.argv[1:] if argv is None else argv
    args = parser.parse_args(arguments)
    if args.input is None:
        if arguments:
            parser.error("Specify --input FILE for command-line processing, or run without arguments for the GUI.")
        return launch_gui()
    try:
        settings = Settings(
            input_file=args.input, min_dist=args.min_dist, precision=args.precision,
            pdmode=args.pdmode, grid_spacing=args.grid_spacing, use_outliers=args.outliers,
            use_boundary=not args.no_boundary, prevent_overshoot=not args.allow_overshoot,
            method=args.method, neighbors=args.neighbors,
            smoothing="auto" if args.smoothing.lower() == "auto" else float(args.smoothing),
            kernel=args.kernel, voxel_size="auto" if args.voxel_size.lower() == "auto" else float(args.voxel_size),
            gaussian_sigma=args.gaussian_sigma, batch_size=args.batch_size,
            validation_samples=args.validation_samples)
        result = process_xyz_to_dxf(settings, lambda message, fraction: print(message, flush=True))
        print(f"Finished in {result.elapsed / 60:.2f} minutes.", flush=True)
    except KeyboardInterrupt:
        print("Conversion interrupted; incomplete new data products were discarded.", file=sys.stderr)
        return 130
    except CancelledError as exc:
        print(str(exc), file=sys.stderr)
        return 130
    except (ValueError, OSError, MemoryError, RuntimeError) as exc:
        print(f"Error: {exc}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())


