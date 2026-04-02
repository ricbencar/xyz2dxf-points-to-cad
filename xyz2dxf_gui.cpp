/*****************************************************************************
 * XYZ to DXF Converter GUI
 * ---------------------------------------------------------------------------
 * Purpose
 *   Windows desktop application for the ingestion, conditioning, interpolation,
 *   quality assessment, and CAD-oriented export of large irregular XYZ point
 *   clouds. The program is intended for scientific and engineering workflows in
 *   which raw surveyed or computed elevation data must be converted into a
 *   cleaned, spatially regular, diagnostically traceable surface representation
 *   suitable for inspection, further numerical use, and DXF-based exchange.
 *
 * General description
 *   The application is not limited to geometric format conversion. It performs
 *   a complete surface-processing workflow beginning with raw XYZ input and
 *   ending with multiple engineering output products. The internal pipeline is
 *   designed to preserve numerical transparency and reproducibility while
 *   providing practical controls through a compact Windows GUI.
 *
 *   The code supports the treatment of irregularly distributed planar samples
 *   (x, y, z), where z typically represents elevation, depth, or another scalar
 *   field defined over a horizontal domain. The processed result is a regular
 *   grid together with confidence and diagnostic information that help the user
 *   interpret where interpolation is well supported and where extrapolation,
 *   boundary effects, or sparse neighborhoods may reduce reliability.
 *
 * Engineering scope
 *   In operational terms, the processing chain may include the following
 * stages, depending on the options selected by the user:
 *
 *     1. XYZ parsing and validation
 *        - Reads the input file line by line.
 *        - Accepts valid numeric triplets representing x, y, z coordinates.
 *        - Ignores malformed or non-data lines when appropriate.
 *        - Tracks total accepted points and reports progress during loading.
 *
 *     2. Optional robust local residual-based outlier filtering
 *        - Evaluates each point relative to its local neighborhood.
 *        - Estimates a local reference surface and local residual structure.
 *        - Rejects isolated or strongly inconsistent z-values when enabled.
 *        - Preserves points in poorly supported areas rather than forcing
 *          aggressive rejection where evidence is insufficient.
 *
 *     3. Deterministic minimum-distance thinning / duplicate suppression
 *        - Reduces redundant local point density in a repeatable manner.
 *        - Prevents extreme clustering from dominating interpolation support.
 *        - Falls back to exact duplicate removal when no positive thinning
 *          distance is requested.
 *        - Produces a more uniform working cloud for surface fitting.
 *
 *     4. Optional outer-boundary construction and enforcement
 *        - Builds a working outer-domain representation from the retained data.
 *        - Uses rasterized occupancy logic and loop extraction to derive a
 *          practical boundary envelope for the cloud.
 *        - Distinguishes interior supported regions from outside-domain zones.
 *        - Penalizes or suppresses predictions that lie outside the supported
 *          data footprint when boundary usage is enabled.
 *
 *     5. Interpolation-context assembly and automatic tuning support
 *        - Builds spatial search structures for neighborhood retrieval.
 *        - Estimates characteristic spacing and support scales.
 *        - Evaluates local sample density and support quality.
 *        - Supports internal parameter selection for more stable local fitting.
 *
 *     6. Regular-grid interpolation by the selected numerical method
 *        - Generates a structured output grid over the processed domain.
 *        - Computes z-values at regular grid nodes using the selected solver.
 *        - Adapts local support where required by the interpolation strategy.
 *        - Balances smoothness, locality, and robustness according to the
 *          chosen method.
 *
 *     7. Confidence, residual, fallback, and boundary diagnostics
 *        - Computes per-node indicators that summarize interpolation quality.
 *        - Identifies local residual levels and neighborhood support counts.
 *        - Marks cases where fallback logic or conservative behaviour was used.
 *        - Flags outside-boundary conditions and local clamping where relevant.
 *
 *     8. Streamed XYZ / confidence / DXF output generation
 *        - Writes cleaned point data for traceability.
 *        - Writes interpolated regular-grid values for downstream use.
 *        - Writes confidence and diagnostic fields for quality review.
 *        - Writes DXF output for CAD visualization and exchange.
 *        - Writes a report file summarizing settings, statistics, and status.
 *
 * Implemented interpolation modes
 *
 *   1. Bicubic Hermite + robust anisotropic MLS
 *      Primary smooth-surface mode intended for high-quality continuous surface
 *      reconstruction on a regular grid.
 *
 *      Conceptually, this mode proceeds in two stages:
 *        - local moving least-squares (MLS) fits are used to estimate nodal
 *          elevations and derivatives from irregular surrounding samples;
 *        - the estimated nodal quantities are then reconstructed over each grid
 *          cell using bicubic Hermite patches.
 *
 *      Main characteristics:
 *        - produces a smooth surface with continuous variation across cells;
 *        - uses local anisotropic support so the fit can adapt to directional
 *          sample distribution;
 *        - suitable where surface regularity and derivative consistency are
 *          important;
 *        - generally preferred for well-sampled terrain or bathymetric clouds
 *          when a smooth engineering surface is desired.
 *
 *   2. Local Thin Plate Spline (TPS)
 *      Flexible local interpolation mode intended for irregular, scattered, or
 *      geometrically complex point distributions.
 *
 *      This mode solves local thin-plate-spline systems over adaptive
 *      neighborhoods and blends the local solutions into the final prediction.
 *      Thin plate splines are attractive when stronger local flexibility is
 *      required to follow curved or spatially varying morphology.
 *
 *      Main characteristics:
 *        - strong ability to follow irregular surface variation;
 *        - adaptive neighborhood support in sparse or dense regions;
 *        - local regularization to reduce instability in difficult point sets;
 *        - appropriate when local curvature tracking is more important than
 *          obtaining the smooth patch structure of the bicubic method.
 *
 *   3. Clamped local MLS (bounded)
 *      Conservative interpolation mode intended for robust local estimation
 * with controlled overshoot behaviour.
 *
 *      This mode evaluates a local MLS surface and then constrains the final
 *      prediction to the local data envelope. It is therefore useful where
 *      conservative interpolation is preferred and where artificial peaks or
 *      depressions generated by local fitting should be limited.
 *
 *      Main characteristics:
 *        - local approximation with bounded output behaviour;
 *        - reduced overshoot risk in sparse or noisy neighborhoods;
 *        - practical for users who prefer conservative predictions;
 *        - especially useful when data quality is uneven or when the surface
 *          should remain close to directly observed local values.
 *
 * Main GUI inputs
 *   - Input File
 *       Source .xyz file chosen through the standard Windows file dialog.
 *
 *   - Min Dist
 *       Minimum horizontal spacing used for deterministic thinning of the input
 *       cloud. This helps suppress redundant near-duplicate samples and reduces
 *       excessive local clustering.
 *
 *   - Precision
 *       Number of decimal places written to the exported text-based outputs.
 *       This affects formatting only and does not change internal arithmetic.
 *
 *   - PDMODE
 *       DXF point-display mode written to the DXF header for point entities.
 *       This affects CAD display behaviour rather than the numerical solution.
 *
 *   - Grid Spacing
 *       Horizontal spacing of the regular output grid on which interpolated
 *       values and diagnostics are computed.
 *
 *   - Method
 *       Choice of interpolation scheme:
 *         bicubic MLS, local TPS, or clamped MLS.
 *
 *   - Options
 *       User controls for enabling or disabling:
 *         outlier filtering and outer-boundary enforcement.
 *
 * Primary output files
 *   - input.xyz.filtered.xyz
 *       Cleaned and thinned working point cloud used as interpolation support.
 *
 *   - input.xyz.grid.xyz
 *       Interpolated regular-grid node coordinates and scalar values.
 *
 *   - input.xyz.confidence.xyz
 *       Grid with confidence and auxiliary diagnostic fields used to assess the
 *       local quality and support of the interpolation.
 *
 *   - input.xyz.dxf
 *       DXF export containing processed point and/or grid representations for
 *       CAD viewing and exchange.
 *
 *   - input.xyz.rpt.txt
 *       Human-readable report containing status messages, summary statistics,
 *       processing settings, and numerical diagnostics.
 *
 * Numerical / implementation notes
 *   - All core numerical calculations are performed in double precision.
 *   - Output precision affects file formatting only.
 *   - If Min Dist <= 0, thinning falls back to exact duplicate suppression.
 *   - The interpolation routines explicitly check for insufficient local
 *     support, unstable local systems, singular or near-singular solves, and
 *     other local failure modes.
 *   - Boundary-aware logic is used to distinguish supported interior regions
 *     from outside-domain locations when boundary enforcement is enabled.
 *   - Confidence values are intended as engineering diagnostics rather than
 *     formal probabilistic uncertainty bounds.
 *   - Residual, fallback, clamping, and neighborhood indicators are exported so
 *     the user can evaluate where the computed surface is strongly supported
 * and where caution may be required.
 *   - The processing flow is designed to be deterministic for a given set of
 *     inputs and options, supporting reproducible engineering review.
 *
 * Build (MinGW / g++)
 *   g++ -O3 -fopenmp -march=native -std=c++17 -Wall -Wextra -pedantic \
 *       -Wconversion -Wsign-conversion -static -static-libgcc \
 *       -static-libstdc++ -mwindows -o xyz2dxf_gui.exe xyz2dxf_gui.cpp \
 *       -lkernel32 -lopengl32 -luuid -lcomdlg32 -lm
 ******************************************************************************/

#ifndef WIN32_LEAN_AND_MEAN
#define WIN32_LEAN_AND_MEAN
#endif

#ifndef NOMINMAX
#define NOMINMAX
#endif

#include <windows.h>
#include <commdlg.h>

#include <algorithm>
#include <array>
#include <atomic>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <fstream>
#include <functional>
#include <iomanip>
#include <limits>
#include <memory>
#include <new>
#include <numeric>
#include <sstream>
#include <stdexcept>
#include <string>
#include <thread>
#include <unordered_map>
#include <unordered_set>
#include <vector>

/*****************************************************************************
 * GUI Component IDs, etc.
 ****************************************************************************/
constexpr int IDC_INPUT_FILE = 1001;
constexpr int IDC_MIN_DIST = 1002;
constexpr int IDC_PRECISION = 1003;
constexpr int IDC_PDMODE = 1004;
constexpr int IDC_GRID_SPACING = 1005;
constexpr int IDC_BROWSE_BUTTON = 1007;
constexpr int IDC_RUN_BUTTON = 1008;
constexpr int IDC_STATUS_STATIC = 1009;
constexpr int IDC_RADIO_BICUBIC = 1010;
constexpr int IDC_RADIO_TPS = 1011;
constexpr int IDC_RADIO_CLAMPED_MLS = 1012;
constexpr int IDC_CHECK_OUTLIERS = 1013;
constexpr int IDC_CHECK_BOUNDARY = 1014;

enum InterpolationMethod {
  METHOD_BICUBIC = 0,
  METHOD_TPS = 1,
  METHOD_CLAMPED_MLS = 2
};

/*****************************************************************************
 * Data Structures
 ****************************************************************************/
struct Point3D {
  double x, y, z; // 3D coordinates

  auto operator==(const Point3D &other) const noexcept -> bool {
    return (x == other.x && y == other.y && z == other.z);
  }
};

// Hash functor for Point3D to enable usage in std::unordered_set
struct Point3DHash {
  auto operator()(const Point3D &p) const noexcept -> size_t {
    const auto h1 = std::hash<double>{}(p.x);
    const auto h2 = std::hash<double>{}(p.y);
    const auto h3 = std::hash<double>{}(p.z);

    size_t seed = 0;
    seed ^= h1 + 0x9e3779b97f4a7c15ULL + (seed << 6) + (seed >> 2);
    seed ^= h2 + 0x9e3779b97f4a7c15ULL + (seed << 6) + (seed >> 2);
    seed ^= h3 + 0x9e3779b97f4a7c15ULL + (seed << 6) + (seed >> 2);
    return seed;
  }
};

struct Cell2D {
  std::int64_t ix = 0;
  std::int64_t iy = 0;

  auto operator==(const Cell2D &other) const noexcept -> bool {
    return ix == other.ix && iy == other.iy;
  }
};

struct Cell2DHash {
  auto operator()(const Cell2D &c) const noexcept -> size_t {
    const auto h1 = std::hash<std::int64_t>{}(c.ix);
    const auto h2 = std::hash<std::int64_t>{}(c.iy);
    size_t seed = 0;
    seed ^= h1 + 0x9e3779b97f4a7c15ULL + (seed << 6) + (seed >> 2);
    seed ^= h2 + 0x9e3779b97f4a7c15ULL + (seed << 6) + (seed >> 2);
    return seed;
  }
};

/*****************************************************************************
 * cellCoord
 * ---------------------------------------------------------------------------
 * Purpose
 *   Converts a floating-point coordinate into an integer spatial-cell index
 * using floor-based binning.
 *****************************************************************************/
static auto cellCoord(double value, double cellSize) noexcept -> std::int64_t {
  return static_cast<std::int64_t>(std::floor(value / cellSize));
}

struct Grid2D {
  size_t nx = 0;
  size_t ny = 0;
  std::vector<double> values;

  Grid2D() = default;
  Grid2D(size_t nx_, size_t ny_, double init = 0.0)
      : nx(nx_), ny(ny_), values(nx_ * ny_, init) {}

  auto operator()(size_t i, size_t j) -> double & { return values[i * ny + j]; }

  auto operator()(size_t i, size_t j) const -> const double & {
    return values[i * ny + j];
  }
};

/*****************************************************************************
 * packSignedIntsToKey
 * ---------------------------------------------------------------------------
 * Purpose
 *   Packs two signed 32-bit cell indices into a single 64-bit key without
 *   invoking undefined behaviour from left-shifting negative signed values.
 *****************************************************************************/
static auto packSignedIntsToKey(int ix, int iy) noexcept -> std::uint64_t {
  const auto ux = static_cast<std::uint32_t>(ix);
  const auto uy = static_cast<std::uint32_t>(iy);
  return (static_cast<std::uint64_t>(ux) << 32U) |
         static_cast<std::uint64_t>(uy);
}

struct SpatialHash2D {
  double cellSize = 1.0;
  double xMin = 0.0;
  double yMin = 0.0;
  std::unordered_map<std::uint64_t, std::vector<size_t>> buckets;

  static auto cellKey(int ix, int iy) noexcept -> std::uint64_t {
    return packSignedIntsToKey(ix, iy);
  }

  void build(const std::vector<Point3D> &points, double requestedCellSize) {
    buckets.clear();
    cellSize = std::max(requestedCellSize, 1e-9);
    xMin = std::numeric_limits<double>::max();
    yMin = std::numeric_limits<double>::max();
    for (const auto &p : points) {
      if (p.x < xMin)
        xMin = p.x;
      if (p.y < yMin)
        yMin = p.y;
    }
    if (!std::isfinite(xMin))
      xMin = 0.0;
    if (!std::isfinite(yMin))
      yMin = 0.0;
    for (size_t i = 0; i < points.size(); ++i) {
      const int ix =
          static_cast<int>(std::floor((points[i].x - xMin) / cellSize));
      const int iy =
          static_cast<int>(std::floor((points[i].y - yMin) / cellSize));
      buckets[cellKey(ix, iy)].push_back(i);
    }
  }

  auto gatherNearest(const std::vector<Point3D> &points, double x, double y,
                     size_t desiredCount, size_t minCount) const
      -> std::vector<size_t> {
    std::vector<size_t> candidates;
    if (points.empty()) {
      return candidates;
    }

    desiredCount = std::min(desiredCount, points.size());
    minCount = std::min(minCount, points.size());

    const int ix0 = static_cast<int>(std::floor((x - xMin) / cellSize));
    const int iy0 = static_cast<int>(std::floor((y - yMin) / cellSize));

    const int maxRing = 64;
    for (int ring = 0; ring <= maxRing; ++ring) {
      for (int dx = -ring; dx <= ring; ++dx) {
        for (int dy = -ring; dy <= ring; ++dy) {
          if (std::max(std::abs(dx), std::abs(dy)) != ring) {
            continue;
          }
          const auto it = buckets.find(cellKey(ix0 + dx, iy0 + dy));
          if (it != buckets.end()) {
            candidates.insert(candidates.end(), it->second.begin(),
                              it->second.end());
          }
        }
      }
      if (candidates.size() >= minCount &&
          (candidates.size() >= desiredCount * 2 || ring >= 4)) {
        break;
      }
    }

    if (candidates.size() < minCount) {
      candidates.resize(points.size());
      std::iota(candidates.begin(), candidates.end(), static_cast<size_t>(0));
    }

    std::sort(candidates.begin(), candidates.end(),
              [&](size_t a, size_t b) -> bool {
                const double da = (points[a].x - x) * (points[a].x - x) +
                                  (points[a].y - y) * (points[a].y - y);
                const double db = (points[b].x - x) * (points[b].x - x) +
                                  (points[b].y - y) * (points[b].y - y);
                if (da != db) {
                  return da < db;
                }
                return a < b;
              });

    if (candidates.size() > desiredCount) {
      candidates.resize(desiredCount);
    }
    return candidates;
  }
};

struct ValidationMetrics {
  double rmse = std::numeric_limits<double>::infinity();
  double mae = std::numeric_limits<double>::infinity();
  double p95 = std::numeric_limits<double>::infinity();
  size_t count = 0;
  bool valid = false;
};

struct InterpolationDiagnostics {
  double confidence = 0.0;
  double residual = 0.0;
  double lambda = 0.0;
  double boundaryPenalty = 1.0;
  double outsideDistance = 0.0;
  size_t neighborCount = 0;
  size_t patchCount = 0;
  bool fallbackUsed = false;
  bool outsideHull = false;
  bool clampedLocal = false;
  std::string summary;
};

struct BicubicNodeData {
  double z = 0.0;
  double fx = 0.0;
  double fy = 0.0;
  double fxy = 0.0;
  double localResidual = 0.0;
  double localConfidence = 0.0;
  double supportRadius = 0.0;
  size_t neighborCount = 0;
  bool fallbackUsed = false;
  bool valid = false;
};

struct ConvexHull2D {
  std::vector<Point3D> vertices;
  std::vector<std::vector<Point3D>> loops;
  double cellSize = 1.0;
  double xMin = 0.0;
  double yMin = 0.0;
  size_t nx = 0;
  size_t ny = 0;
  std::vector<unsigned char> occupied;
  bool valid = false;
};

struct TPSModel {
  std::vector<Point3D> controlPoints;
  SpatialHash2D spatialIndex;
  ConvexHull2D hull;
  size_t neighborhoodSize = 48;
  double suggestedLambda = 1e-10;
  ValidationMetrics validation;
  bool valid = false;
};

struct GridDefinition {
  double xMin = 0.0;
  double yMin = 0.0;
  double xMax = 0.0;
  double yMax = 0.0;
  double gridSpacing = 1.0;
  size_t nx = 0;
  size_t ny = 0;
  size_t totalPoints = 0;
  bool valid = false;
};

struct TuningChoice {
  size_t neighborhoodTarget = 0;
  double supportMultiplier = 1.75;
  double ridgeFactor = 1.0;
  int basisOrder = 3;
  ValidationMetrics metrics;
};

struct ConfidenceSummary {
  size_t count = 0;
  size_t fallbackCount = 0;
  size_t outsideHullCount = 0;
  size_t clampedLocalCount = 0;
  double meanConfidence = 0.0;
  double meanResidual = 0.0;
};

constexpr UINT WM_APP_STATUS_UPDATE = WM_APP + 1;
constexpr UINT WM_APP_PROCESSING_DONE = WM_APP + 2;

/**
 * \brief Holds the 16 coefficients for a bicubic polynomial patch:
 *        z = sum_{m=0..3} sum_{n=0..3} a_{mn} * (x^m)(y^n)
 */
struct BicubicSpline {
  double a00, a01, a02, a03;
  double a10, a11, a12, a13;
  double a20, a21, a22, a23;
  double a30, a31, a32, a33;
};

/*****************************************************************************
 * Forward Declarations
 ****************************************************************************/
struct LocalFrame2D;
static auto solveDenseLinearSystem(std::vector<double> A, std::vector<double> b,
                                   size_t n, std::vector<double> &x) -> bool;
static auto averagePlanarSpacing(const std::vector<Point3D> &points) -> double;
static auto fitRobustLocalMLS(const std::vector<Point3D> &points,
                              const std::vector<size_t> &fitIndices,
                              const LocalFrame2D &frame, double supportRadius,
                              int basisOrder, double ridgeFactor,
                              std::vector<double> &coeffs,
                              ValidationMetrics *fitMetrics) -> bool;
static auto solveLocalTPSAtPoint(double x, double y, const TPSModel &model,
                                 double &z, std::string *diagnostics,
                                 InterpolationDiagnostics *interpDiagnostics)
    -> bool;
static void updateStatus(HWND hwndStatus, const std::string &message);
static auto openFileDialog(HWND hwnd) -> std::string;
static void computeBoundingBox(const std::vector<Point3D> &points, double &xMin,
                               double &xMax, double &yMin, double &yMax);

/*****************************************************************************
 * makeControlMenuHandle
 * ---------------------------------------------------------------------------
 * Purpose
 *   Builds a Windows HMENU-compatible control identifier for child-window
 * creation.
 *****************************************************************************/
static auto makeControlMenuHandle(int controlId) noexcept -> HMENU {
  return reinterpret_cast<HMENU>(static_cast<INT_PTR>(controlId));
}

static auto filterPointsGrid(const std::vector<Point3D> &points, double minDist)
    -> std::vector<Point3D>;
static auto
readXYZFile(const std::string &inputFileName, std::vector<Point3D> &points,
            size_t &totalPointsRead,
            const std::function<void(const std::string &)> &progressCallback,
            std::string *errorMessage) -> bool;

/*****************************************************************************
 * Basic GUI / geometry helpers
 ****************************************************************************/
static void updateStatus(HWND hwndStatus, const std::string &message) {
  SetWindowTextA(hwndStatus, message.c_str());
}

/*****************************************************************************
 * openFileDialog
 * ---------------------------------------------------------------------------
 * Purpose
 *   Opens the standard Windows file-selection dialog and returns the chosen
 * input path.
 *****************************************************************************/
static auto openFileDialog(HWND hwnd) -> std::string {
  OPENFILENAMEA ofn{};
  std::array<char, MAX_PATH> szFile{};

  ofn.lStructSize = sizeof(ofn);
  ofn.hwndOwner = hwnd;
  ofn.lpstrFile = szFile.data();
  ofn.nMaxFile = static_cast<DWORD>(szFile.size());
  ofn.lpstrFilter = "XYZ Files (*.xyz)\0*.xyz\0All Files (*.*)\0*.*\0";
  ofn.nFilterIndex = 1;
  ofn.Flags = OFN_PATHMUSTEXIST | OFN_FILEMUSTEXIST;

  if (GetOpenFileNameA(&ofn)) {
    return std::string(ofn.lpstrFile);
  }
  return "";
}

/*****************************************************************************
 * computeBoundingBox
 * ---------------------------------------------------------------------------
 * Purpose
 *   Computes the planar bounding box of a point set in x and y.
 *****************************************************************************/
static void computeBoundingBox(const std::vector<Point3D> &points, double &xMin,
                               double &xMax, double &yMin, double &yMax) {
  if (points.empty()) {
    xMin = 0.0;
    xMax = 0.0;
    yMin = 0.0;
    yMax = 0.0;
    return;
  }

  xMin = points.front().x;
  xMax = points.front().x;
  yMin = points.front().y;
  yMax = points.front().y;

  for (const auto &p : points) {
    xMin = std::min(xMin, p.x);
    xMax = std::max(xMax, p.x);
    yMin = std::min(yMin, p.y);
    yMax = std::max(yMax, p.y);
  }
}

/*****************************************************************************
 * filterPointsGrid
 * ---------------------------------------------------------------------------
 * Purpose
 *   Applies deterministic horizontal thinning. If the minimum distance is
 * non-positive, it falls back to exact duplicate removal.
 *****************************************************************************/
static auto filterPointsGrid(const std::vector<Point3D> &points, double minDist)
    -> std::vector<Point3D> {
  if (points.empty()) {
    return {};
  }

  if (!(minDist > 0.0)) {
    std::unordered_set<Point3D, Point3DHash> seen;
    std::vector<Point3D> uniquePoints;
    uniquePoints.reserve(points.size());
    for (const auto &p : points) {
      if (seen.emplace(p).second) {
        uniquePoints.push_back(p);
      }
    }
    std::sort(uniquePoints.begin(), uniquePoints.end(),
              [](const Point3D &a, const Point3D &b) -> bool {
                if (a.x != b.x)
                  return a.x < b.x;
                if (a.y != b.y)
                  return a.y < b.y;
                return a.z < b.z;
              });
    return uniquePoints;
  }

  struct CandidatePoint {
    Point3D p;
    Cell2D cell;
    double score = 0.0;
  };

  std::vector<CandidatePoint> candidates;
  candidates.reserve(points.size());
  for (const auto &p : points) {
    const Cell2D cell{cellCoord(p.x, minDist), cellCoord(p.y, minDist)};
    const double cx = (static_cast<double>(cell.ix) + 0.5) * minDist;
    const double cy = (static_cast<double>(cell.iy) + 0.5) * minDist;
    const double dx = p.x - cx;
    const double dy = p.y - cy;
    candidates.push_back({p, cell, dx * dx + dy * dy});
  }
  std::stable_sort(
      candidates.begin(), candidates.end(),
      [](const CandidatePoint &a, const CandidatePoint &b) -> bool {
        if (a.cell.ix != b.cell.ix)
          return a.cell.ix < b.cell.ix;
        if (a.cell.iy != b.cell.iy)
          return a.cell.iy < b.cell.iy;
        if (a.score != b.score)
          return a.score < b.score;
        if (a.p.x != b.p.x)
          return a.p.x < b.p.x;
        if (a.p.y != b.p.y)
          return a.p.y < b.p.y;
        return a.p.z < b.p.z;
      });

  const double minDistSq = minDist * minDist;
  std::unordered_map<Cell2D, std::vector<Point3D>, Cell2DHash> sparseGrid;
  sparseGrid.reserve(std::min<size_t>(candidates.size(), 1U << 20));

  std::vector<Point3D> accepted;
  accepted.reserve(points.size());

  for (const auto &cand : candidates) {
    bool tooClose = false;
    for (std::int64_t dx = -1; dx <= 1 && !tooClose; ++dx) {
      for (std::int64_t dy = -1; dy <= 1 && !tooClose; ++dy) {
        const auto it =
            sparseGrid.find(Cell2D{cand.cell.ix + dx, cand.cell.iy + dy});
        if (it == sparseGrid.end()) {
          continue;
        }
        for (const auto &q : it->second) {
          const double ddx = cand.p.x - q.x;
          const double ddy = cand.p.y - q.y;
          if (ddx * ddx + ddy * ddy < minDistSq) {
            tooClose = true;
            break;
          }
        }
      }
    }
    if (!tooClose) {
      sparseGrid[cand.cell].push_back(cand.p);
      accepted.push_back(cand.p);
    }
  }

  return accepted;
}

/*****************************************************************************
 * readXYZFile
 * ---------------------------------------------------------------------------
 * Purpose
 *   Streams an XYZ file from disk, validates parseable rows, and reports
 * progress while loading points.
 *****************************************************************************/
static auto
readXYZFile(const std::string &inputFileName, std::vector<Point3D> &points,
            size_t &totalPointsRead,
            const std::function<void(const std::string &)> &progressCallback,
            std::string *errorMessage) -> bool {
  points.clear();
  totalPointsRead = 0U;

  std::ifstream inFile(inputFileName);
  if (!inFile.is_open()) {
    if (errorMessage) {
      *errorMessage = "Error: Unable to open input file.";
    }
    return false;
  }

  constexpr size_t progressStep = 1000000U;
  std::string line;
  points.reserve(1U << 20);

  while (std::getline(inFile, line)) {
    if (line.empty() || line[0] == '#') {
      continue;
    }
    std::replace(line.begin(), line.end(), ',', ' ');
    std::istringstream ss(line);
    Point3D p{};
    if (!(ss >> p.x >> p.y >> p.z)) {
      continue;
    }
    ++totalPointsRead;
    points.push_back(p);
    if (progressCallback && (totalPointsRead % progressStep) == 0U) {
      std::ostringstream msg;
      msg << "Read " << totalPointsRead << " XYZ points ...";
      progressCallback(msg.str());
    }
  }

  if (points.empty()) {
    if (errorMessage) {
      *errorMessage = "Error: No valid points found in file.";
    }
    return false;
  }

  if (progressCallback) {
    std::ostringstream msg;
    msg << "Input read complete: " << points.size() << " valid XYZ points.";
    progressCallback(msg.str());
  }
  return true;
}

/*****************************************************************************
 * computeMedianUnchecked
 * ---------------------------------------------------------------------------
 * Purpose
 *   Computes the median of a numeric sample without requiring prior sorting by
 * the caller.
 *****************************************************************************/
static auto computeMedianUnchecked(std::vector<double> values) -> double {
  if (values.empty()) {
    return 0.0;
  }
  const size_t mid = values.size() / 2U;
  std::nth_element(values.begin(),
                   values.begin() + static_cast<std::ptrdiff_t>(mid),
                   values.end());
  double med = values[mid];
  if ((values.size() % 2U) == 0U) {
    std::nth_element(values.begin(),
                     values.begin() + static_cast<std::ptrdiff_t>(mid - 1U),
                     values.end());
    med = 0.5 * (med + values[mid - 1U]);
  }
  return med;
}

/*****************************************************************************
 * makeValidationMetricsFromErrors
 * ---------------------------------------------------------------------------
 * Purpose
 *   Builds RMSE, MAE, and percentile error metrics from a residual/error
 * sample.
 *****************************************************************************/
static auto makeValidationMetricsFromErrors(std::vector<double> errors)
    -> ValidationMetrics {
  ValidationMetrics out;
  if (errors.empty()) {
    return out;
  }
  double sumSq = 0.0;
  double sumAbs = 0.0;
  for (double e : errors) {
    const double ae = std::fabs(e);
    sumSq += ae * ae;
    sumAbs += ae;
  }
  std::sort(errors.begin(), errors.end(), [](double a, double b) -> bool {
    return std::fabs(a) < std::fabs(b);
  });
  const size_t q95 = std::min(errors.size() - 1U, (95U * errors.size()) / 100U);
  out.count = errors.size();
  out.rmse = std::sqrt(sumSq / static_cast<double>(errors.size()));
  out.mae = sumAbs / static_cast<double>(errors.size());
  out.p95 = std::fabs(errors[q95]);
  out.valid = true;
  return out;
}

/*****************************************************************************
 * validationScore
 * ---------------------------------------------------------------------------
 * Purpose
 *   Combines validation metrics into a single scalar score used for tuning
 * comparisons.
 *****************************************************************************/
static auto validationScore(const ValidationMetrics &m) -> double {
  if (!m.valid) {
    return std::numeric_limits<double>::infinity();
  }
  return m.rmse + 0.35 * m.p95 + 0.15 * m.mae;
}

/*****************************************************************************
 * calibratedConfidenceFromObservedError
 * ---------------------------------------------------------------------------
 * Purpose
 *   Converts a raw confidence estimate into a calibrated confidence using
 * observed validation error levels.
 *****************************************************************************/
static auto calibratedConfidenceFromObservedError(
    double rawConfidence, double localResidual,
    const ValidationMetrics &observedMetrics) -> double {
  if (!std::isfinite(rawConfidence) || !(rawConfidence > 0.0)) {
    rawConfidence = 0.05;
  }
  if (!observedMetrics.valid) {
    return std::max(0.05, std::min(1.0, rawConfidence));
  }
  const double scale =
      std::max({observedMetrics.rmse, 0.5 * observedMetrics.p95, 1e-8});
  const double penalty = 1.0 / (1.0 + std::max(localResidual, 0.0) / scale);
  return std::max(0.02, std::min(1.0, rawConfidence * penalty));
}

/*****************************************************************************
 * makeDeterministicHoldoutIndices
 * ---------------------------------------------------------------------------
 * Purpose
 *   Builds deterministic holdout indices for repeatable validation and tuning
 * passes.
 *****************************************************************************/
static auto makeDeterministicHoldoutIndices(size_t n, size_t maxCount)
    -> std::vector<size_t> {
  std::vector<size_t> out;
  if (n == 0U || maxCount == 0U) {
    return out;
  }
  const size_t target = std::min(maxCount, n);
  const size_t step = std::max<size_t>(1U, n / target);
  for (size_t i = step / 2U; i < n && out.size() < target; i += step) {
    out.push_back(i);
  }
  if (out.empty()) {
    out.push_back(0U);
  }
  return out;
}

/*****************************************************************************
 * polygonArea2D
 * ---------------------------------------------------------------------------
 * Purpose
 *   Computes the signed area of a 2D polygon represented by Point3D vertices.
 *****************************************************************************/
static auto polygonArea2D(const std::vector<Point3D> &poly) -> double {
  if (poly.size() < 3U) {
    return 0.0;
  }
  double area = 0.0;
  for (size_t i = 0; i < poly.size(); ++i) {
    const Point3D &a = poly[i];
    const Point3D &b = poly[(i + 1U) % poly.size()];
    area += a.x * b.y - b.x * a.y;
  }
  return 0.5 * area;
}

/*****************************************************************************
 * boundaryMaskAt
 * ---------------------------------------------------------------------------
 * Purpose
 *   Queries whether a rasterized boundary mask cell is occupied.
 *****************************************************************************/
static auto boundaryMaskAt(const ConvexHull2D &hull, int ix, int iy) -> bool {
  if (ix < 0 || iy < 0) {
    return false;
  }
  const auto sx = static_cast<size_t>(ix);
  const auto sy = static_cast<size_t>(iy);
  if (sx >= hull.nx || sy >= hull.ny) {
    return false;
  }
  const size_t idx = sx * hull.ny + sy;
  return idx < hull.occupied.size() && hull.occupied[idx] != 0U;
}

/*****************************************************************************
 * simplifyAxisAlignedLoop
 * ---------------------------------------------------------------------------
 * Purpose
 *   Removes redundant collinear vertices from a raster-derived boundary loop.
 *****************************************************************************/
static void simplifyAxisAlignedLoop(std::vector<Point3D> &loop) {
  if (loop.size() < 3U) {
    return;
  }

  bool changed = true;
  while (changed && loop.size() >= 3U) {
    changed = false;
    std::vector<Point3D> simplified;
    simplified.reserve(loop.size());
    const size_t n = loop.size();
    for (size_t i = 0; i < n; ++i) {
      const Point3D &prev = loop[(i + n - 1U) % n];
      const Point3D &cur = loop[i];
      const Point3D &next = loop[(i + 1U) % n];
      const double v1x = cur.x - prev.x;
      const double v1y = cur.y - prev.y;
      const double v2x = next.x - cur.x;
      const double v2y = next.y - cur.y;
      const double cross = v1x * v2y - v1y * v2x;
      const double dot = v1x * v2x + v1y * v2y;
      if (std::fabs(cross) <= 1e-12 && dot >= 0.0) {
        changed = true;
        continue;
      }
      simplified.push_back(cur);
    }
    if (!simplified.empty()) {
      loop.swap(simplified);
    } else {
      break;
    }
  }
}

/*****************************************************************************
 * dilateBinaryMask
 * ---------------------------------------------------------------------------
 * Purpose
 *   Performs binary-mask dilation on the raster boundary representation.
 *****************************************************************************/
static auto dilateBinaryMask(const std::vector<unsigned char> &mask, size_t nx,
                             size_t ny, int radius)
    -> std::vector<unsigned char> {
  if (mask.empty() || nx == 0U || ny == 0U || radius <= 0) {
    return mask;
  }

  std::vector<unsigned char> out(nx * ny, 0U);
  for (size_t ix = 0; ix < nx; ++ix) {
    for (size_t iy = 0; iy < ny; ++iy) {
      if (mask[ix * ny + iy] == 0U) {
        continue;
      }
      const int ix0 = static_cast<int>(ix);
      const int iy0 = static_cast<int>(iy);
      for (int dx = -radius; dx <= radius; ++dx) {
        for (int dy = -radius; dy <= radius; ++dy) {
          const int jx = ix0 + dx;
          const int jy = iy0 + dy;
          if (jx < 0 || jy < 0) {
            continue;
          }
          const auto sx = static_cast<size_t>(jx);
          const auto sy = static_cast<size_t>(jy);
          if (sx < nx && sy < ny) {
            out[sx * ny + sy] = 1U;
          }
        }
      }
    }
  }
  return out;
}

/*****************************************************************************
 * erodeBinaryMask
 * ---------------------------------------------------------------------------
 * Purpose
 *   Performs binary-mask erosion on the raster boundary representation.
 *****************************************************************************/
static auto erodeBinaryMask(const std::vector<unsigned char> &mask, size_t nx,
                            size_t ny, int radius)
    -> std::vector<unsigned char> {
  if (mask.empty() || nx == 0U || ny == 0U || radius <= 0) {
    return mask;
  }

  std::vector<unsigned char> out(nx * ny, 0U);
  for (size_t ix = 0; ix < nx; ++ix) {
    for (size_t iy = 0; iy < ny; ++iy) {
      bool keep = true;
      const int ix0 = static_cast<int>(ix);
      const int iy0 = static_cast<int>(iy);
      for (int dx = -radius; dx <= radius && keep; ++dx) {
        for (int dy = -radius; dy <= radius; ++dy) {
          const int jx = ix0 + dx;
          const int jy = iy0 + dy;
          if (jx < 0 || jy < 0) {
            keep = false;
            break;
          }
          const auto sx = static_cast<size_t>(jx);
          const auto sy = static_cast<size_t>(jy);
          if (sx >= nx || sy >= ny || mask[sx * ny + sy] == 0U) {
            keep = false;
            break;
          }
        }
      }
      if (keep) {
        out[ix * ny + iy] = 1U;
      }
    }
  }
  return out;
}

/*****************************************************************************
 * closeBinaryMask
 * ---------------------------------------------------------------------------
 * Purpose
 *   Applies morphological closing to connect nearby occupied raster cells.
 *****************************************************************************/
static auto closeBinaryMask(const std::vector<unsigned char> &mask, size_t nx,
                            size_t ny, int radius)
    -> std::vector<unsigned char> {
  if (mask.empty() || nx == 0U || ny == 0U || radius <= 0) {
    return mask;
  }
  return erodeBinaryMask(dilateBinaryMask(mask, nx, ny, radius), nx, ny,
                         radius);
}

/*****************************************************************************
 * countConnectedComponents
 * ---------------------------------------------------------------------------
 * Purpose
 *   Counts connected occupied components in a raster mask.
 *****************************************************************************/
static auto countConnectedComponents(const std::vector<unsigned char> &mask,
                                     size_t nx, size_t ny) -> size_t {
  if (mask.empty() || nx == 0U || ny == 0U) {
    return 0U;
  }

  std::vector<unsigned char> visited(nx * ny, 0U);
  std::vector<size_t> stack;
  stack.reserve(1024U);
  size_t componentCount = 0U;

  for (size_t ix = 0; ix < nx; ++ix) {
    for (size_t iy = 0; iy < ny; ++iy) {
      const size_t start = ix * ny + iy;
      if (mask[start] == 0U || visited[start] != 0U) {
        continue;
      }
      ++componentCount;
      visited[start] = 1U;
      stack.clear();
      stack.push_back(start);
      while (!stack.empty()) {
        const size_t idx = stack.back();
        stack.pop_back();

        const size_t cx = idx / ny;
        const size_t cy = idx % ny;
        const int x0 = static_cast<int>(cx);
        const int y0 = static_cast<int>(cy);
        const int offsets[4][2] = {{-1, 0}, {1, 0}, {0, -1}, {0, 1}};
        for (const auto &off : offsets) {
          const int nxCell = x0 + off[0];
          const int nyCell = y0 + off[1];
          if (nxCell < 0 || nyCell < 0) {
            continue;
          }
          const auto sx = static_cast<size_t>(nxCell);
          const auto sy = static_cast<size_t>(nyCell);
          if (sx >= nx || sy >= ny) {
            continue;
          }
          const size_t nidx = sx * ny + sy;
          if (mask[nidx] != 0U && visited[nidx] == 0U) {
            visited[nidx] = 1U;
            stack.push_back(nidx);
          }
        }
      }
    }
  }
  return componentCount;
}

/*****************************************************************************
 * fillInteriorHoles
 * ---------------------------------------------------------------------------
 * Purpose
 *   Fills enclosed holes inside a raster mask while preserving the exterior
 * background.
 *****************************************************************************/
static auto fillInteriorHoles(const std::vector<unsigned char> &mask, size_t nx,
                              size_t ny) -> std::vector<unsigned char> {
  if (mask.empty() || nx == 0U || ny == 0U) {
    return mask;
  }

  std::vector<unsigned char> outside(nx * ny, 0U);
  std::vector<size_t> stack;
  stack.reserve(2U * (nx + ny) + 16U);

  auto tryPush = [&](size_t ix, size_t iy) -> void {
    const size_t idx = ix * ny + iy;
    if (mask[idx] == 0U && outside[idx] == 0U) {
      outside[idx] = 1U;
      stack.push_back(idx);
    }
  };

  for (size_t ix = 0; ix < nx; ++ix) {
    tryPush(ix, 0U);
    tryPush(ix, ny - 1U);
  }
  for (size_t iy = 0; iy < ny; ++iy) {
    tryPush(0U, iy);
    tryPush(nx - 1U, iy);
  }

  while (!stack.empty()) {
    const size_t idx = stack.back();
    stack.pop_back();

    const size_t cx = idx / ny;
    const size_t cy = idx % ny;
    const int x0 = static_cast<int>(cx);
    const int y0 = static_cast<int>(cy);

    const int offsets[4][2] = {{-1, 0}, {1, 0}, {0, -1}, {0, 1}};
    for (const auto &off : offsets) {
      const int nxCell = x0 + off[0];
      const int nyCell = y0 + off[1];
      if (nxCell < 0 || nyCell < 0) {
        continue;
      }
      const auto sx = static_cast<size_t>(nxCell);
      const auto sy = static_cast<size_t>(nyCell);
      if (sx >= nx || sy >= ny) {
        continue;
      }
      const size_t nidx = sx * ny + sy;
      if (mask[nidx] == 0U && outside[nidx] == 0U) {
        outside[nidx] = 1U;
        stack.push_back(nidx);
      }
    }
  }

  std::vector<unsigned char> filled = mask;
  for (size_t idx = 0; idx < filled.size(); ++idx) {
    if (filled[idx] == 0U && outside[idx] == 0U) {
      filled[idx] = 1U;
    }
  }
  return filled;
}

/*****************************************************************************
 * buildConvexHull2D
 * ---------------------------------------------------------------------------
 * Purpose
 *   Builds the working outer boundary representation from the filtered point
 * set using rasterization and loop extraction.
 *****************************************************************************/
static auto buildConvexHull2D(const std::vector<Point3D> &points,
                              double requestedCellSize = 0.0) -> ConvexHull2D {
  ConvexHull2D hull;
  if (points.empty()) {
    return hull;
  }

  double dataXMin = 0.0;
  double dataXMax = 0.0;
  double dataYMin = 0.0;
  double dataYMax = 0.0;
  computeBoundingBox(points, dataXMin, dataXMax, dataYMin, dataYMax);

  const double spacing = averagePlanarSpacing(points);
  hull.cellSize = (requestedCellSize > 0.0) ? std::max(requestedCellSize, 1e-6)
                                            : std::max(2.5 * spacing, 1e-6);

  const int marginCells = 2;
  hull.xMin = dataXMin - static_cast<double>(marginCells) * hull.cellSize;
  hull.yMin = dataYMin - static_cast<double>(marginCells) * hull.cellSize;

  hull.nx =
      static_cast<size_t>(std::floor((dataXMax - hull.xMin) / hull.cellSize)) +
      static_cast<size_t>(2 * marginCells + 2);
  hull.ny =
      static_cast<size_t>(std::floor((dataYMax - hull.yMin) / hull.cellSize)) +
      static_cast<size_t>(2 * marginCells + 2);
  hull.nx = std::max<size_t>(hull.nx, 3U);
  hull.ny = std::max<size_t>(hull.ny, 3U);

  std::vector<unsigned char> raster(hull.nx * hull.ny, 0U);
  for (const auto &p : points) {
    const int ix =
        static_cast<int>(std::floor((p.x - hull.xMin) / hull.cellSize));
    const int iy =
        static_cast<int>(std::floor((p.y - hull.yMin) / hull.cellSize));
    if (ix >= 0 && iy >= 0) {
      const auto sx = static_cast<size_t>(ix);
      const auto sy = static_cast<size_t>(iy);
      if (sx < hull.nx && sy < hull.ny) {
        raster[sx * hull.ny + sy] = 1U;
      }
    }
  }

  std::vector<unsigned char> chosenMask = raster;
  size_t chosenComponents =
      countConnectedComponents(chosenMask, hull.nx, hull.ny);

  constexpr int maxClosingRadius = 12;
  if (chosenComponents > 1U) {
    for (int radius = 1; radius <= maxClosingRadius; ++radius) {
      std::vector<unsigned char> candidate =
          closeBinaryMask(raster, hull.nx, hull.ny, radius);
      candidate = fillInteriorHoles(candidate, hull.nx, hull.ny);
      const size_t componentCount =
          countConnectedComponents(candidate, hull.nx, hull.ny);
      chosenMask = std::move(candidate);
      chosenComponents = componentCount;
      if (componentCount <= 1U) {
        break;
      }
    }
  }

  if (chosenComponents > 1U) {
    int radius = maxClosingRadius + 1;
    const int hardMaxRadius = 32;
    while (chosenComponents > 1U && radius <= hardMaxRadius) {
      std::vector<unsigned char> candidate =
          closeBinaryMask(raster, hull.nx, hull.ny, radius);
      candidate = fillInteriorHoles(candidate, hull.nx, hull.ny);
      const size_t componentCount =
          countConnectedComponents(candidate, hull.nx, hull.ny);
      chosenMask = std::move(candidate);
      chosenComponents = componentCount;
      ++radius;
    }
  }

  hull.occupied = fillInteriorHoles(chosenMask, hull.nx, hull.ny);
  for (const auto &p : points) {
    const int ix =
        static_cast<int>(std::floor((p.x - hull.xMin) / hull.cellSize));
    const int iy =
        static_cast<int>(std::floor((p.y - hull.yMin) / hull.cellSize));
    if (ix >= 0 && iy >= 0) {
      const auto sx = static_cast<size_t>(ix);
      const auto sy = static_cast<size_t>(iy);
      if (sx < hull.nx && sy < hull.ny) {
        hull.occupied[sx * hull.ny + sy] = 1U;
      }
    }
  }
  hull.occupied = fillInteriorHoles(hull.occupied, hull.nx, hull.ny);
  hull.valid = std::any_of(hull.occupied.begin(), hull.occupied.end(),
                           [](unsigned char v) -> bool { return v != 0U; });
  if (!hull.valid) {
    return hull;
  }

  struct Edge2D {
    int x0 = 0;
    int y0 = 0;
    int x1 = 0;
    int y1 = 0;
    bool used = false;
  };

  auto vertexKey = [](int x, int y) noexcept -> std::uint64_t {
    return packSignedIntsToKey(x, y);
  };

  std::vector<Edge2D> edges;
  edges.reserve(hull.nx * hull.ny);
  std::unordered_map<std::uint64_t, std::vector<size_t>> outgoing;

  auto addEdge = [&](int x0, int y0, int x1, int y1) -> void {
    const size_t idx = edges.size();
    edges.push_back(Edge2D{x0, y0, x1, y1, false});
    outgoing[vertexKey(x0, y0)].push_back(idx);
  };

  for (size_t ix = 0; ix < hull.nx; ++ix) {
    for (size_t iy = 0; iy < hull.ny; ++iy) {
      if (hull.occupied[ix * hull.ny + iy] == 0U) {
        continue;
      }
      const int x = static_cast<int>(ix);
      const int y = static_cast<int>(iy);
      if (!boundaryMaskAt(hull, x, y - 1))
        addEdge(x, y, x + 1, y);
      if (!boundaryMaskAt(hull, x + 1, y))
        addEdge(x + 1, y, x + 1, y + 1);
      if (!boundaryMaskAt(hull, x, y + 1))
        addEdge(x + 1, y + 1, x, y + 1);
      if (!boundaryMaskAt(hull, x - 1, y))
        addEdge(x, y + 1, x, y);
    }
  }

  for (size_t e = 0; e < edges.size(); ++e) {
    if (edges[e].used) {
      continue;
    }

    std::vector<Point3D> loop;
    const int startX = edges[e].x0;
    const int startY = edges[e].y0;
    int curX = startX;
    int curY = startY;
    size_t currentEdge = e;
    bool closed = false;

    while (true) {
      Edge2D &edge = edges[currentEdge];
      if (edge.used) {
        break;
      }
      edge.used = true;
      if (loop.empty()) {
        loop.push_back(
            {hull.xMin + static_cast<double>(edge.x0) * hull.cellSize,
             hull.yMin + static_cast<double>(edge.y0) * hull.cellSize, 0.0});
      }
      loop.push_back({hull.xMin + static_cast<double>(edge.x1) * hull.cellSize,
                      hull.yMin + static_cast<double>(edge.y1) * hull.cellSize,
                      0.0});
      curX = edge.x1;
      curY = edge.y1;
      if (curX == startX && curY == startY) {
        closed = true;
        break;
      }
      const auto it = outgoing.find(vertexKey(curX, curY));
      if (it == outgoing.end()) {
        break;
      }
      bool foundNext = false;
      for (size_t nextEdge : it->second) {
        if (!edges[nextEdge].used) {
          currentEdge = nextEdge;
          foundNext = true;
          break;
        }
      }
      if (!foundNext) {
        break;
      }
    }

    if (closed && loop.size() >= 4U) {
      loop.pop_back();
      simplifyAxisAlignedLoop(loop);
      if (loop.size() >= 3U) {
        hull.loops.push_back(std::move(loop));
      }
    }
  }

  if (!hull.loops.empty()) {
    size_t bestIdx = 0U;
    double bestArea = 0.0;
    for (size_t i = 0; i < hull.loops.size(); ++i) {
      const double area = std::fabs(polygonArea2D(hull.loops[i]));
      if (area > bestArea) {
        bestArea = area;
        bestIdx = i;
      }
    }
    hull.vertices = hull.loops[bestIdx];
    std::vector<std::vector<Point3D>> singleLoop;
    singleLoop.push_back(hull.vertices);
    hull.loops.swap(singleLoop);
  }

  return hull;
}

/*****************************************************************************
 * distancePointToSegment2D
 * ---------------------------------------------------------------------------
 * Purpose
 *   Computes the minimum planar distance from a point to a line segment.
 *****************************************************************************/
static auto distancePointToSegment2D(double x, double y, double x0, double y0,
                                     double x1, double y1) -> double {
  const double dx = x1 - x0;
  const double dy = y1 - y0;
  const double l2 = dx * dx + dy * dy;
  if (!(l2 > 0.0)) {
    const double dx0 = x - x0;
    const double dy0 = y - y0;
    return std::sqrt(dx0 * dx0 + dy0 * dy0);
  }
  const double t =
      std::min(1.0, std::max(0.0, ((x - x0) * dx + (y - y0) * dy) / l2));
  const double px = x0 + t * dx;
  const double py = y0 + t * dy;
  const double ddx = x - px;
  const double ddy = y - py;
  return std::sqrt(ddx * ddx + ddy * ddy);
}

/*****************************************************************************
 * pointInPolygonLoop2D
 * ---------------------------------------------------------------------------
 * Purpose
 *   Evaluates whether a planar point lies inside or on a polygon loop.
 *****************************************************************************/
static auto pointInPolygonLoop2D(const std::vector<Point3D> &loop, double x,
                                 double y) -> bool {
  if (loop.size() < 3U) {
    return false;
  }

  bool inside = false;
  for (size_t i = 0, j = loop.size() - 1U; i < loop.size(); j = i++) {
    const double xi = loop[i].x;
    const double yi = loop[i].y;
    const double xj = loop[j].x;
    const double yj = loop[j].y;
    const double dist = distancePointToSegment2D(x, y, xi, yi, xj, yj);
    if (dist <= 1e-10) {
      return true;
    }
    const bool intersects =
        ((yi > y) != (yj > y)) &&
        (x < (xj - xi) * (y - yi) / std::max(yj - yi, 1e-30) + xi);
    if (intersects) {
      inside = !inside;
    }
  }
  return inside;
}

/*****************************************************************************
 * pointInConvexHull2D
 * ---------------------------------------------------------------------------
 * Purpose
 *   Tests whether a point lies inside the active outer-boundary representation.
 *****************************************************************************/
static auto pointInConvexHull2D(const ConvexHull2D &hull, double x, double y)
    -> bool {
  if (!hull.valid) {
    return false;
  }

  const double cs = std::max(hull.cellSize, 1e-6);
  const int ix = static_cast<int>(std::floor((x - hull.xMin) / cs));
  const int iy = static_cast<int>(std::floor((y - hull.yMin) / cs));
  if (boundaryMaskAt(hull, ix, iy)) {
    return true;
  }

  bool inside = false;
  for (const auto &loop : hull.loops) {
    if (pointInPolygonLoop2D(loop, x, y)) {
      inside = !inside;
    }
  }
  return inside;
}

/*****************************************************************************
 * distanceToConvexHull2D
 * ---------------------------------------------------------------------------
 * Purpose
 *   Computes the planar distance from a point to the active outer boundary.
 *****************************************************************************/
static auto distanceToConvexHull2D(const ConvexHull2D &hull, double x, double y)
    -> double {
  if (!hull.valid) {
    return 0.0;
  }
  if (pointInConvexHull2D(hull, x, y)) {
    return 0.0;
  }

  double best = std::numeric_limits<double>::infinity();
  for (const auto &loop : hull.loops) {
    if (loop.size() < 2U) {
      continue;
    }
    for (size_t i = 0; i < loop.size(); ++i) {
      const Point3D &a = loop[i];
      const Point3D &b = loop[(i + 1U) % loop.size()];
      best = std::min(best, distancePointToSegment2D(x, y, a.x, a.y, b.x, b.y));
    }
  }
  return std::isfinite(best) ? best : 0.0;
}

/*****************************************************************************
 * solveWeightedPlaneAtPoint
 * ---------------------------------------------------------------------------
 * Purpose
 *   Fits a locally weighted planar model and optionally reports residual
 * metrics.
 *****************************************************************************/
static auto solveWeightedPlaneAtPoint(
    const std::vector<Point3D> &points, const std::vector<size_t> &neighbors,
    double x0, double y0, double &a, double &bx, double &by,
    ValidationMetrics *residualMetrics,
    size_t excludeIndex = std::numeric_limits<size_t>::max()) -> bool {
  a = 0.0;
  bx = 0.0;
  by = 0.0;
  if (neighbors.size() < 3U) {
    return false;
  }

  double meanR2 = 0.0;
  size_t usedCount = 0U;
  for (size_t idx : neighbors) {
    const double dx = points[idx].x - x0;
    const double dy = points[idx].y - y0;
    meanR2 += dx * dx + dy * dy;
    ++usedCount;
  }
  if (usedCount < 3U) {
    return false;
  }
  meanR2 /= static_cast<double>(usedCount);
  const double h = std::max(std::sqrt(std::max(meanR2, 1e-12)), 1e-6);

  std::vector<double> ATA(9U, 0.0);
  std::vector<double> ATb(3U, 0.0);
  for (size_t idx : neighbors) {
    const double dx = points[idx].x - x0;
    const double dy = points[idx].y - y0;
    const double r2 = dx * dx + dy * dy;
    const double w = 1.0 / (1.0 + r2 / (h * h));
    const double phi[3] = {1.0, dx, dy};
    for (size_t r = 0; r < 3U; ++r) {
      ATb[r] += w * phi[r] * points[idx].z;
      for (size_t c = 0; c < 3U; ++c) {
        ATA[r * 3U + c] += w * phi[r] * phi[c];
      }
    }
  }
  for (size_t d = 0; d < 3U; ++d) {
    ATA[d * 3U + d] += 1e-12;
  }

  std::vector<double> coeffs;
  if (!solveDenseLinearSystem(ATA, ATb, 3U, coeffs) || coeffs.size() < 3U) {
    return false;
  }
  a = coeffs[0];
  bx = coeffs[1];
  by = coeffs[2];

  if (residualMetrics) {
    std::vector<double> residuals;
    residuals.reserve(neighbors.size());
    for (size_t idx : neighbors) {
      if (idx == excludeIndex) {
        continue;
      }
      const double dx = points[idx].x - x0;
      const double dy = points[idx].y - y0;
      residuals.push_back(points[idx].z - (a + bx * dx + by * dy));
    }
    *residualMetrics = makeValidationMetricsFromErrors(std::move(residuals));
  }
  return std::isfinite(a) && std::isfinite(bx) && std::isfinite(by);
}

/*****************************************************************************
 * removeZOutliersInPlace
 * ---------------------------------------------------------------------------
 * Purpose
 *   Performs robust local residual-based outlier rejection directly on the
 * working point vector.
 *****************************************************************************/
static auto removeZOutliersInPlace(
    std::vector<Point3D> &points, double neighborDist, double zThresholdFactor,
    const std::function<void(const std::string &)> &progressCallback,
    std::string *errorMessage) -> bool {
  if (points.empty()) {
    return true;
  }
  if (!(zThresholdFactor > 0.0) || !(neighborDist > 0.0)) {
    return true;
  }

  const double neighborDistSq = neighborDist * neighborDist;
  SpatialHash2D spatialIndex;
  spatialIndex.build(points, std::max(neighborDist, 1e-6));
  std::vector<unsigned char> keep(points.size(), 1U);

  const size_t reportStep =
      std::max<size_t>(250000U, points.size() / 100U + 1U);
  size_t kept = 0U;
  size_t rejected = 0U;

  for (size_t i = 0; i < points.size(); ++i) {
    const Point3D &pi = points[i];
    const auto candidates = spatialIndex.gatherNearest(
        points, pi.x, pi.y, std::min<size_t>(points.size(), 96U),
        std::min<size_t>(points.size(), 12U));
    std::vector<size_t> neighbors;
    neighbors.reserve(candidates.size());
    for (size_t idx : candidates) {
      const double dx = points[idx].x - pi.x;
      const double dy = points[idx].y - pi.y;
      if (idx != i && dx * dx + dy * dy <= neighborDistSq) {
        neighbors.push_back(idx);
      }
    }
    if (neighbors.size() < 6U) {
      keep[i] = 1U;
      ++kept;
      continue;
    }

    double a = 0.0;
    double bx = 0.0;
    double by = 0.0;
    ValidationMetrics residualMetrics;
    if (!solveWeightedPlaneAtPoint(points, neighbors, pi.x, pi.y, a, bx, by,
                                   &residualMetrics, i)) {
      keep[i] = 1U;
      ++kept;
      continue;
    }

    std::vector<double> residuals;
    residuals.reserve(neighbors.size());
    for (size_t idx : neighbors) {
      const double dx = points[idx].x - pi.x;
      const double dy = points[idx].y - pi.y;
      residuals.push_back(points[idx].z - (a + bx * dx + by * dy));
    }
    const double medianResidual = computeMedianUnchecked(residuals);
    for (double &r : residuals) {
      r = std::fabs(r - medianResidual);
    }
    const double mad = computeMedianUnchecked(residuals);
    const double robustSigma = std::max(1.4826 * mad, 1e-8);
    const double selfResidual = std::fabs(pi.z - a - medianResidual);
    const bool keepPoint = (selfResidual <= zThresholdFactor * robustSigma);
    keep[i] = keepPoint ? 1U : 0U;
    if (keepPoint) {
      ++kept;
    } else {
      ++rejected;
    }

    if (progressCallback && ((i + 1U) % reportStep) == 0U) {
      std::ostringstream msg;
      const double pct = 100.0 * static_cast<double>(i + 1U) /
                         static_cast<double>(points.size());
      msg << "Robust residual outlier scan " << (i + 1U) << "/" << points.size()
          << ", kept=" << kept << ", rejected=" << rejected << " ("
          << std::fixed << std::setprecision(1) << pct << "%)";
      progressCallback(msg.str());
    }
  }

  size_t writeIdx = 0U;
  for (size_t i = 0; i < points.size(); ++i) {
    if (keep[i] != 0U) {
      if (writeIdx != i) {
        points[writeIdx] = points[i];
      }
      ++writeIdx;
    }
  }
  points.resize(writeIdx);
  points.shrink_to_fit();

  if (points.empty()) {
    if (errorMessage) {
      *errorMessage =
          "Error: All points were removed during robust residual outlier "
          "filtering.";
    }
    return false;
  }
  if (progressCallback) {
    std::ostringstream msg;
    msg << "Robust residual outlier filtering complete: retained "
        << points.size() << " points, removed " << rejected << ".";
    progressCallback(msg.str());
  }
  return true;
} /****************************************************************************
   * 7) High-accuracy interpolation helpers
   ****************************************************************************/
static auto averagePlanarSpacing(const std::vector<Point3D> &points) -> double {
  if (points.size() < 2) {
    return 1.0;
  }

  double xMin, xMax, yMin, yMax;
  computeBoundingBox(points, xMin, xMax, yMin, yMax);
  const double area = std::max((xMax - xMin) * (yMax - yMin), 1e-12);
  return std::sqrt(area / static_cast<double>(points.size()));
}

/*****************************************************************************
 * recommendedBoundaryCellSize
 * ---------------------------------------------------------------------------
 * Purpose
 *   Chooses a practical raster cell size for boundary construction from point
 * spacing and thinning hints.
 *****************************************************************************/
static auto recommendedBoundaryCellSize(const std::vector<Point3D> &points,
                                        double minDistHint) -> double {
  const double spacing = averagePlanarSpacing(points);
  if (minDistHint > 0.0 && std::isfinite(minDistHint)) {
    return std::max(minDistHint, 1e-6);
  }
  return std::max(spacing, 1e-6);
}

/*****************************************************************************
 * solveDenseLinearSystem
 * ---------------------------------------------------------------------------
 * Purpose
 *   Solves a small dense linear system with the internal direct solver used
 * throughout the interpolation routines.
 *****************************************************************************/
static auto solveDenseLinearSystem(std::vector<double> A, std::vector<double> b,
                                   size_t n, std::vector<double> &x) -> bool {
  auto idx = [n](size_t r, size_t c) -> size_t { return r * n + c; };

  x.assign(n, 0.0);
  for (size_t k = 0; k < n; ++k) {
    size_t pivot = k;
    double pivotAbs = std::fabs(A[idx(k, k)]);
    for (size_t r = k + 1; r < n; ++r) {
      const double cand = std::fabs(A[idx(r, k)]);
      if (cand > pivotAbs) {
        pivotAbs = cand;
        pivot = r;
      }
    }
    if (!(pivotAbs > 1e-18) || !std::isfinite(pivotAbs)) {
      return false;
    }
    if (pivot != k) {
      for (size_t c = k; c < n; ++c) {
        std::swap(A[idx(k, c)], A[idx(pivot, c)]);
      }
      std::swap(b[k], b[pivot]);
    }

    const double diag = A[idx(k, k)];
    for (size_t r = k + 1; r < n; ++r) {
      const double factor = A[idx(r, k)] / diag;
      A[idx(r, k)] = 0.0;
      if (factor == 0.0) {
        continue;
      }
      for (size_t c = k + 1; c < n; ++c) {
        A[idx(r, c)] -= factor * A[idx(k, c)];
      }
      b[r] -= factor * b[k];
    }
  }

  for (size_t ii = n; ii-- > 0;) {
    double sum = b[ii];
    for (size_t c = ii + 1; c < n; ++c) {
      sum -= A[idx(ii, c)] * x[c];
    }
    const double diag = A[idx(ii, ii)];
    if (!(std::fabs(diag) > 1e-18) || !std::isfinite(diag)) {
      return false;
    }
    x[ii] = sum / diag;
  }
  return true;
}

/*****************************************************************************
 * compactWeightWendland
 * ---------------------------------------------------------------------------
 * Purpose
 *   Evaluates the compact Wendland weight function used by local MLS fits.
 *****************************************************************************/
static auto compactWeightWendland(double distance, double supportRadius)
    -> double {
  if (!(supportRadius > 0.0)) {
    return 1.0;
  }
  const double q = distance / supportRadius;
  if (q >= 1.0) {
    return 0.0;
  }
  const double oneMinus = 1.0 - q;
  return oneMinus * oneMinus * oneMinus * oneMinus * (4.0 * q + 1.0);
}

/*****************************************************************************
 * inverseDistanceFallback
 * ---------------------------------------------------------------------------
 * Purpose
 *   Computes a conservative inverse-distance prediction when the preferred
 * local model is unavailable.
 *****************************************************************************/
static auto inverseDistanceFallback(const std::vector<Point3D> &points,
                                    const SpatialHash2D &spatialIndex, double x,
                                    double y) -> double {
  const auto neighbors = spatialIndex.gatherNearest(
      points, x, y, std::min<size_t>(16, points.size()), 1);
  double sumW = 0.0;
  double sumZ = 0.0;
  for (size_t idx : neighbors) {
    const double dx = points[idx].x - x;
    const double dy = points[idx].y - y;
    const double w = 1.0 / (1e-12 + dx * dx + dy * dy);
    sumW += w;
    sumZ += w * points[idx].z;
  }
  if (sumW <= 0.0) {
    return points.empty() ? 0.0 : points.front().z;
  }
  return sumZ / sumW;
}

struct LocalFrame2D {
  double x0 = 0.0;
  double y0 = 0.0;
  double e1x = 1.0;
  double e1y = 0.0;
  double e2x = 0.0;
  double e2y = 1.0;
  double s1 = 1.0;
  double s2 = 1.0;
};

/*****************************************************************************
 * computeLocalPrincipalFrame
 * ---------------------------------------------------------------------------
 * Purpose
 *   Builds a local principal frame used to orient anisotropic neighborhoods.
 *****************************************************************************/
static void computeLocalPrincipalFrame(const std::vector<Point3D> &points,
                                       const std::vector<size_t> &neighbors,
                                       double x0, double y0,
                                       double referenceScale,
                                       LocalFrame2D &frame) {
  frame = LocalFrame2D{};
  frame.x0 = x0;
  frame.y0 = y0;

  if (neighbors.empty()) {
    frame.s1 = std::max(referenceScale, 1e-6);
    frame.s2 = std::max(referenceScale, 1e-6);
    return;
  }

  double meanR2 = 0.0;
  for (size_t idx : neighbors) {
    const double dx = points[idx].x - x0;
    const double dy = points[idx].y - y0;
    meanR2 += dx * dx + dy * dy;
  }
  meanR2 /= static_cast<double>(neighbors.size());
  const double h =
      std::max({std::sqrt(std::max(meanR2, 1e-12)), referenceScale, 1e-6});

  double sxx = 0.0;
  double syy = 0.0;
  double sxy = 0.0;
  double wsum = 0.0;
  for (size_t idx : neighbors) {
    const double dx = points[idx].x - x0;
    const double dy = points[idx].y - y0;
    const double r2 = dx * dx + dy * dy;
    const double w = 1.0 / (1.0 + r2 / (h * h));
    sxx += w * dx * dx;
    syy += w * dy * dy;
    sxy += w * dx * dy;
    wsum += w;
  }

  if (!(wsum > 0.0)) {
    frame.s1 = std::max(referenceScale, 1e-6);
    frame.s2 = std::max(referenceScale, 1e-6);
    return;
  }

  sxx /= wsum;
  syy /= wsum;
  sxy /= wsum;

  const double tr = sxx + syy;
  const double det_term =
      std::sqrt(std::max((sxx - syy) * (sxx - syy) + 4.0 * sxy * sxy, 0.0));
  const double lambda1 = std::max(0.5 * (tr + det_term), 0.0);
  const double lambda2 = std::max(0.5 * (tr - det_term), 0.0);
  const double theta = 0.5 * std::atan2(2.0 * sxy, sxx - syy);

  frame.e1x = std::cos(theta);
  frame.e1y = std::sin(theta);
  frame.e2x = -frame.e1y;
  frame.e2y = frame.e1x;
  frame.s1 =
      std::max(std::sqrt(lambda1), 0.75 * std::max(referenceScale, 1e-6));
  frame.s2 =
      std::max(std::sqrt(lambda2), 0.35 * std::max(referenceScale, 1e-6));
}

/*****************************************************************************
 * toLocalUV
 * ---------------------------------------------------------------------------
 * Purpose
 *   Transforms a global point into the local principal-frame coordinates.
 *****************************************************************************/
static void toLocalUV(const LocalFrame2D &frame, double x, double y, double &u,
                      double &v) {
  const double dx = x - frame.x0;
  const double dy = y - frame.y0;
  u = (dx * frame.e1x + dy * frame.e1y) / frame.s1;
  v = (dx * frame.e2x + dy * frame.e2y) / frame.s2;
}

/*****************************************************************************
 * robustSupportRadiusFromLocalUV
 * ---------------------------------------------------------------------------
 * Purpose
 *   Computes a robust support radius for local fitting in the transformed
 * neighborhood space.
 *****************************************************************************/
static auto robustSupportRadiusFromLocalUV(const std::vector<Point3D> &points,
                                           const std::vector<size_t> &neighbors,
                                           const LocalFrame2D &frame,
                                           double minimumRadius) -> double {
  std::vector<double> radii;
  radii.reserve(neighbors.size());
  for (size_t idx : neighbors) {
    double u = 0.0;
    double v = 0.0;
    toLocalUV(frame, points[idx].x, points[idx].y, u, v);
    radii.push_back(std::sqrt(u * u + v * v));
  }
  if (radii.empty()) {
    return std::max(1.5, minimumRadius);
  }
  std::sort(radii.begin(), radii.end());
  const size_t q80 = std::min(radii.size() - 1U, (radii.size() * 4U) / 5U);
  const size_t q90 = std::min(radii.size() - 1U, (radii.size() * 9U) / 10U);
  return std::max({minimumRadius, 1.5, radii[q80], 0.75 * radii[q90]});
}

/*****************************************************************************
 * selectAdaptiveNeighborhood
 * ---------------------------------------------------------------------------
 * Purpose
 *   Selects an adaptive neighborhood around a query location using the spatial
 * hash.
 *****************************************************************************/
static auto selectAdaptiveNeighborhood(const std::vector<Point3D> &points,
                                       const std::vector<size_t> &candidates,
                                       const LocalFrame2D &frame,
                                       size_t minKeep, size_t maxKeep,
                                       double supportRadius)
    -> std::vector<size_t> {
  struct Entry {
    size_t idx = 0;
    double rho = 0.0;
  };

  std::vector<Entry> ranked;
  ranked.reserve(candidates.size());
  for (size_t idx : candidates) {
    double u = 0.0;
    double v = 0.0;
    toLocalUV(frame, points[idx].x, points[idx].y, u, v);
    ranked.push_back({idx, std::sqrt(u * u + v * v)});
  }
  std::sort(ranked.begin(), ranked.end(),
            [](const Entry &a, const Entry &b) -> bool {
              if (a.rho != b.rho) {
                return a.rho < b.rho;
              }
              return a.idx < b.idx;
            });

  std::vector<size_t> selected;
  selected.reserve(std::min(maxKeep, ranked.size()));
  for (const auto &entry : ranked) {
    if (entry.rho <= supportRadius || selected.size() < minKeep) {
      selected.push_back(entry.idx);
      if (selected.size() >= maxKeep) {
        break;
      }
    } else {
      break;
    }
  }
  if (selected.size() < minKeep) {
    selected.clear();
    for (size_t k = 0; k < ranked.size() && k < maxKeep; ++k) {
      selected.push_back(ranked[k].idx);
    }
  }
  return selected;
}

/*****************************************************************************
 * pointDistanceXY
 * ---------------------------------------------------------------------------
 * Purpose
 *   Computes planar distance between two points.
 *****************************************************************************/
static auto pointDistanceXY(double x0, double y0, double x1, double y1)
    -> double {
  const double dx = x0 - x1;
  const double dy = y0 - y1;
  return std::sqrt(dx * dx + dy * dy);
}

/*****************************************************************************
 * polynomialBasisCount
 * ---------------------------------------------------------------------------
 * Purpose
 *   Returns the number of polynomial terms for the requested basis order.
 *****************************************************************************/
static auto polynomialBasisCount(int basisOrder) -> size_t {
  return (basisOrder >= 3) ? 10U : 6U;
}

/*****************************************************************************
 * fillPolynomialBasis
 * ---------------------------------------------------------------------------
 * Purpose
 *   Fills the polynomial basis vector used by the MLS solver.
 *****************************************************************************/
static void fillPolynomialBasis(double u, double v, int basisOrder,
                                double *phi) {
  phi[0] = 1.0;
  phi[1] = u;
  phi[2] = v;
  phi[3] = 0.5 * u * u;
  phi[4] = u * v;
  phi[5] = 0.5 * v * v;
  if (basisOrder >= 3) {
    phi[6] = (u * u * u) / 6.0;
    phi[7] = 0.5 * u * u * v;
    phi[8] = 0.5 * u * v * v;
    phi[9] = (v * v * v) / 6.0;
  }
}

/*****************************************************************************
 * evaluatePolynomialValueFromCoeffs
 * ---------------------------------------------------------------------------
 * Purpose
 *   Evaluates the fitted polynomial at a local query position.
 *****************************************************************************/
static auto evaluatePolynomialValueFromCoeffs(const std::vector<double> &coeffs,
                                              int basisOrder, double u,
                                              double v) -> double {
  double phi[10] = {};
  fillPolynomialBasis(u, v, basisOrder, phi);
  const size_t basisCount = polynomialBasisCount(basisOrder);
  double value = 0.0;
  for (size_t i = 0; i < basisCount; ++i) {
    value += coeffs[i] * phi[i];
  }
  return value;
}

/*****************************************************************************
 * evaluateMLSValidation
 * ---------------------------------------------------------------------------
 * Purpose
 *   Evaluates MLS predictive quality on the supplied validation subset.
 *****************************************************************************/
static auto evaluateMLSValidation(const std::vector<Point3D> &points,
                                  const std::vector<size_t> &indices,
                                  const LocalFrame2D &frame,
                                  const std::vector<double> &coeffs,
                                  int basisOrder) -> ValidationMetrics {
  std::vector<double> errors;
  errors.reserve(indices.size());
  for (size_t idx : indices) {
    double u = 0.0;
    double v = 0.0;
    toLocalUV(frame, points[idx].x, points[idx].y, u, v);
    const double pred =
        evaluatePolynomialValueFromCoeffs(coeffs, basisOrder, u, v);
    errors.push_back(pred - points[idx].z);
  }
  return makeValidationMetricsFromErrors(std::move(errors));
}

/*****************************************************************************
 * evaluateLocalMLSPredictiveValidation
 * ---------------------------------------------------------------------------
 * Purpose
 *   Performs localized predictive validation for a candidate MLS fit.
 *****************************************************************************/
static auto evaluateLocalMLSPredictiveValidation(
    const std::vector<Point3D> &points, const std::vector<size_t> &neighbors,
    const LocalFrame2D &frame, double supportRadius, int basisOrder,
    double ridgeFactor, size_t maxChecks = 8U) -> ValidationMetrics {
  std::vector<double> errors;
  if (neighbors.size() < polynomialBasisCount(basisOrder) + 2U) {
    return ValidationMetrics{};
  }

  const size_t step =
      std::max<size_t>(1U, neighbors.size() / std::max<size_t>(1U, maxChecks));
  for (size_t k = step / 2U; k < neighbors.size() && errors.size() < maxChecks;
       k += step) {
    const size_t holdoutIdx = neighbors[k];
    std::vector<size_t> fitIndices;
    fitIndices.reserve(neighbors.size() - 1U);
    for (size_t idx : neighbors) {
      if (idx != holdoutIdx) {
        fitIndices.push_back(idx);
      }
    }
    if (fitIndices.size() < polynomialBasisCount(basisOrder)) {
      continue;
    }
    std::vector<double> coeffs;
    if (!fitRobustLocalMLS(points, fitIndices, frame, supportRadius, basisOrder,
                           ridgeFactor, coeffs, nullptr)) {
      continue;
    }
    double u = 0.0;
    double v = 0.0;
    toLocalUV(frame, points[holdoutIdx].x, points[holdoutIdx].y, u, v);
    const double pred =
        evaluatePolynomialValueFromCoeffs(coeffs, basisOrder, u, v);
    if (std::isfinite(pred)) {
      errors.push_back(pred - points[holdoutIdx].z);
    }
  }
  return makeValidationMetricsFromErrors(std::move(errors));
}

/*****************************************************************************
 * fitRobustLocalMLS
 * ---------------------------------------------------------------------------
 * Purpose
 *   Fits a robust local MLS polynomial model with optional validation
 * diagnostics.
 *****************************************************************************/
static auto fitRobustLocalMLS(const std::vector<Point3D> &points,
                              const std::vector<size_t> &fitIndices,
                              const LocalFrame2D &frame, double supportRadius,
                              int basisOrder, double ridgeFactor,
                              std::vector<double> &coeffs,
                              ValidationMetrics *fitMetrics) -> bool {
  const size_t basisCount = polynomialBasisCount(basisOrder);
  coeffs.assign(basisCount, 0.0);
  if (fitIndices.size() < basisCount) {
    return false;
  }

  std::vector<double> robustWeights(fitIndices.size(), 1.0);
  ValidationMetrics localFitMetrics;

  for (int iter = 0; iter < 3; ++iter) {
    std::vector<double> ATA(basisCount * basisCount, 0.0);
    std::vector<double> ATb(basisCount, 0.0);
    double baseTrace = 0.0;

    for (size_t k = 0; k < fitIndices.size(); ++k) {
      const size_t idx = fitIndices[k];
      double u = 0.0;
      double v = 0.0;
      toLocalUV(frame, points[idx].x, points[idx].y, u, v);
      const double rho = std::sqrt(u * u + v * v);
      const double wCompact = compactWeightWendland(rho, supportRadius * 1.001);
      const double wGeom = std::max(wCompact / (1.0 + rho * rho), 1e-14);
      const double w = std::max(wGeom * robustWeights[k], 1e-14);
      double phi[10] = {};
      fillPolynomialBasis(u, v, basisOrder, phi);
      for (size_t r = 0; r < basisCount; ++r) {
        ATb[r] += w * phi[r] * points[idx].z;
        for (size_t c = 0; c < basisCount; ++c) {
          ATA[r * basisCount + c] += w * phi[r] * phi[c];
        }
      }
    }
    for (size_t d = 0; d < basisCount; ++d) {
      baseTrace += ATA[d * basisCount + d];
    }
    const double ridge = std::max(1e-12, ridgeFactor * 1e-11 * baseTrace /
                                             static_cast<double>(basisCount));
    for (size_t d = 0; d < basisCount; ++d) {
      ATA[d * basisCount + d] += ridge;
    }
    if (!solveDenseLinearSystem(ATA, ATb, basisCount, coeffs)) {
      return false;
    }

    localFitMetrics =
        evaluateMLSValidation(points, fitIndices, frame, coeffs, basisOrder);
    std::vector<double> residuals;
    residuals.reserve(fitIndices.size());
    for (size_t idx : fitIndices) {
      double u = 0.0;
      double v = 0.0;
      toLocalUV(frame, points[idx].x, points[idx].y, u, v);
      const double pred =
          evaluatePolynomialValueFromCoeffs(coeffs, basisOrder, u, v);
      residuals.push_back(points[idx].z - pred);
    }
    const double med = computeMedianUnchecked(residuals);
    for (double &r : residuals) {
      r = std::fabs(r - med);
    }
    const double mad = computeMedianUnchecked(residuals);
    const double sigma = std::max(1.4826 * mad, 1e-8);
    for (size_t k = 0; k < fitIndices.size(); ++k) {
      const size_t idx = fitIndices[k];
      double u = 0.0;
      double v = 0.0;
      toLocalUV(frame, points[idx].x, points[idx].y, u, v);
      const double pred =
          evaluatePolynomialValueFromCoeffs(coeffs, basisOrder, u, v);
      const double t = std::fabs(points[idx].z - pred) / (4.685 * sigma);
      if (t >= 1.0) {
        robustWeights[k] = 1e-4;
      } else {
        const double oneMinus = 1.0 - t * t;
        robustWeights[k] = std::max(oneMinus * oneMinus, 1e-4);
      }
    }
  }

  if (fitMetrics) {
    *fitMetrics = localFitMetrics;
  }
  return std::all_of(coeffs.begin(), coeffs.end(),
                     [](double v) -> bool { return std::isfinite(v); });
}

/*****************************************************************************
 * defaultTuningChoice
 * ---------------------------------------------------------------------------
 * Purpose
 *   Builds the default interpolation tuning choice used as an initial
 * candidate.
 *****************************************************************************/
static auto defaultTuningChoice() -> TuningChoice {
  TuningChoice t;
  t.neighborhoodTarget = 64U;
  t.supportMultiplier = 1.75;
  t.ridgeFactor = 1.0;
  t.basisOrder = 3;
  return t;
}

/*****************************************************************************
 * limitHermiteNodeDerivatives
 * ---------------------------------------------------------------------------
 * Purpose
 *   Applies slope limiting to bicubic Hermite nodal derivatives.
 *****************************************************************************/
static void limitHermiteNodeDerivatives(BicubicNodeData &node,
                                        double gridSpacing, double localZMin,
                                        double localZMax) {
  const double zRange = std::max(localZMax - localZMin, 1e-8);
  const double slopeLimit = 3.0 * zRange / std::max(gridSpacing, 1e-6);
  const double twistLimit =
      4.0 * zRange / std::max(gridSpacing * gridSpacing, 1e-6);
  node.fx = std::max(-slopeLimit, std::min(slopeLimit, node.fx));
  node.fy = std::max(-slopeLimit, std::min(slopeLimit, node.fy));
  node.fxy = std::max(-twistLimit, std::min(twistLimit, node.fxy));
}

/*****************************************************************************
 * minmod2
 * ---------------------------------------------------------------------------
 * Purpose
 *   Evaluates a two-argument minmod limiter used by derivative limiting.
 *****************************************************************************/
static auto minmod2(double a, double b) -> double {
  if (a * b <= 0.0) {
    return 0.0;
  }
  return (std::fabs(a) < std::fabs(b)) ? a : b;
}

/*****************************************************************************
 * limitHermiteCellDerivatives
 * ---------------------------------------------------------------------------
 * Purpose
 *   Applies consistency limiting to the derivatives used within a bicubic cell.
 *****************************************************************************/
static void limitHermiteCellDerivatives(BicubicNodeData &n00,
                                        BicubicNodeData &n10,
                                        BicubicNodeData &n01,
                                        BicubicNodeData &n11, double hx,
                                        double hy) {
  const double sx0 = (n10.z - n00.z) / std::max(hx, 1e-6);
  const double sx1 = (n11.z - n01.z) / std::max(hx, 1e-6);
  const double sy0 = (n01.z - n00.z) / std::max(hy, 1e-6);
  const double sy1 = (n11.z - n10.z) / std::max(hy, 1e-6);

  const double limFx00 = minmod2(n00.fx, sx0);
  const double limFx10 = minmod2(n10.fx, sx0);
  const double limFx01 = minmod2(n01.fx, sx1);
  const double limFx11 = minmod2(n11.fx, sx1);
  const double limFy00 = minmod2(n00.fy, sy0);
  const double limFy01 = minmod2(n01.fy, sy0);
  const double limFy10 = minmod2(n10.fy, sy1);
  const double limFy11 = minmod2(n11.fy, sy1);

  n00.fx = limFx00;
  n10.fx = limFx10;
  n01.fx = limFx01;
  n11.fx = limFx11;
  n00.fy = limFy00;
  n01.fy = limFy01;
  n10.fy = limFy10;
  n11.fy = limFy11;

  const double crossBase = 0.5 * ((sy1 - sy0) / std::max(hx, 1e-6) +
                                  (sx1 - sx0) / std::max(hy, 1e-6));
  const double crossLimit = 2.0 * std::fabs(crossBase) + 1e-8;
  n00.fxy =
      std::max(-crossLimit, std::min(crossLimit, minmod2(n00.fxy, crossBase)));
  n10.fxy =
      std::max(-crossLimit, std::min(crossLimit, minmod2(n10.fxy, crossBase)));
  n01.fxy =
      std::max(-crossLimit, std::min(crossLimit, minmod2(n01.fxy, crossBase)));
  n11.fxy =
      std::max(-crossLimit, std::min(crossLimit, minmod2(n11.fxy, crossBase)));
}

/*****************************************************************************
 * fitLocalPolynomialNode
 * ---------------------------------------------------------------------------
 * Purpose
 *   Estimates bicubic-node elevation and derivatives at a grid node from nearby
 * scattered data.
 *****************************************************************************/
static auto fitLocalPolynomialNode(
    const std::vector<Point3D> &points, const SpatialHash2D &spatialIndex,
    double x0, double y0, double gridSpacing, BicubicNodeData &node,
    const TuningChoice *tuningChoice = nullptr,
    InterpolationDiagnostics *diagnostics = nullptr,
    size_t excludeIndex = std::numeric_limits<size_t>::max()) -> bool {
  node = BicubicNodeData{};
  if (diagnostics) {
    *diagnostics = InterpolationDiagnostics{};
  }
  if (points.empty()) {
    return false;
  }

  const TuningChoice tuning =
      tuningChoice ? *tuningChoice : defaultTuningChoice();
  const size_t probeCount = std::min<size_t>(
      points.size(), std::max<size_t>(tuning.neighborhoodTarget * 3U, 192U));
  const size_t minCount = std::min<size_t>(points.size(), 24U);
  auto candidates =
      spatialIndex.gatherNearest(points, x0, y0, probeCount, minCount);
  if (excludeIndex != std::numeric_limits<size_t>::max()) {
    candidates.erase(
        std::remove(candidates.begin(), candidates.end(), excludeIndex),
        candidates.end());
  }
  if (candidates.empty()) {
    return false;
  }

  LocalFrame2D probeFrame;
  computeLocalPrincipalFrame(points, candidates, x0, y0, gridSpacing,
                             probeFrame);
  const double baseSupport =
      robustSupportRadiusFromLocalUV(points, candidates, probeFrame, 1.6);
  const size_t maxKeep = std::min<size_t>(
      points.size(), std::max<size_t>(tuning.neighborhoodTarget + 24U, 96U));
  auto neighbors = selectAdaptiveNeighborhood(
      points, candidates, probeFrame, minCount, maxKeep,
      std::max(1.5, baseSupport * tuning.supportMultiplier));
  if (neighbors.size() <
      std::max<size_t>(12U, polynomialBasisCount(tuning.basisOrder) + 4U)) {
    node.z = inverseDistanceFallback(points, spatialIndex, x0, y0);
    node.localConfidence = 0.1;
    node.fallbackUsed = true;
    node.valid = true;
    if (diagnostics) {
      diagnostics->fallbackUsed = true;
      diagnostics->confidence = 0.1;
      diagnostics->neighborCount = neighbors.size();
      diagnostics->summary =
          "MLS fallback to IDW due to insufficient neighborhood.";
    }
    return true;
  }

  LocalFrame2D frame;
  computeLocalPrincipalFrame(points, neighbors, x0, y0, gridSpacing, frame);
  const double supportRadius = std::max(
      1.5, robustSupportRadiusFromLocalUV(points, neighbors, frame, 1.5) *
               tuning.supportMultiplier);

  std::vector<size_t> fitIndices = neighbors;

  std::vector<double> coeffs;
  ValidationMetrics fitMetrics;
  if (!fitRobustLocalMLS(points, fitIndices, frame, supportRadius,
                         tuning.basisOrder, tuning.ridgeFactor, coeffs,
                         &fitMetrics)) {
    node.z = inverseDistanceFallback(points, spatialIndex, x0, y0);
    node.localConfidence = 0.1;
    node.fallbackUsed = true;
    node.valid = true;
    if (diagnostics) {
      diagnostics->fallbackUsed = true;
      diagnostics->confidence = 0.1;
      diagnostics->neighborCount = neighbors.size();
      diagnostics->summary =
          "MLS fallback to IDW because the robust solve failed.";
    }
    return true;
  }

  ValidationMetrics validationMetrics = evaluateLocalMLSPredictiveValidation(
      points, neighbors, frame, supportRadius, tuning.basisOrder,
      tuning.ridgeFactor);
  if (!validationMetrics.valid) {
    validationMetrics = fitMetrics;
  }

  const double du_dx = frame.e1x / frame.s1;
  const double du_dy = frame.e1y / frame.s1;
  const double dv_dx = frame.e2x / frame.s2;
  const double dv_dy = frame.e2y / frame.s2;

  const double fu = coeffs[1];
  const double fv = coeffs[2];
  const double fuu = coeffs[3];
  const double fuv = coeffs[4];
  const double fvv = coeffs[5];

  node.z = coeffs[0];
  node.fx = fu * du_dx + fv * dv_dx;
  node.fy = fu * du_dy + fv * dv_dy;
  node.fxy = fuu * du_dx * du_dy + fuv * (du_dx * dv_dy + dv_dx * du_dy) +
             fvv * dv_dx * dv_dy;
  node.localResidual =
      validationMetrics.valid ? validationMetrics.rmse : fitMetrics.rmse;
  const double rawConfidence =
      validationMetrics.valid
          ? (1.0 / (1.0 + 1.5 * validationScore(validationMetrics)))
          : 0.35;
  node.localConfidence = std::clamp(rawConfidence, 0.05, 1.0);
  node.supportRadius = supportRadius;
  node.neighborCount = neighbors.size();
  node.valid = std::isfinite(node.z) && std::isfinite(node.fx) &&
               std::isfinite(node.fy) && std::isfinite(node.fxy);

  double localZMin = std::numeric_limits<double>::infinity();
  double localZMax = -std::numeric_limits<double>::infinity();
  for (size_t idx : neighbors) {
    localZMin = std::min(localZMin, points[idx].z);
    localZMax = std::max(localZMax, points[idx].z);
  }
  if (std::isfinite(localZMin) && std::isfinite(localZMax)) {
    limitHermiteNodeDerivatives(node, gridSpacing, localZMin, localZMax);
  }

  if (!node.valid) {
    node.z = inverseDistanceFallback(points, spatialIndex, x0, y0);
    node.fx = 0.0;
    node.fy = 0.0;
    node.fxy = 0.0;
    node.localConfidence = 0.1;
    node.fallbackUsed = true;
    node.valid = true;
  }

  if (diagnostics) {
    diagnostics->confidence = node.localConfidence;
    diagnostics->residual = node.localResidual;
    diagnostics->neighborCount = node.neighborCount;
    diagnostics->patchCount = 1U;
    diagnostics->fallbackUsed = node.fallbackUsed;
    diagnostics->summary = validationMetrics.valid
                               ? "robust anisotropic MLS"
                               : "robust anisotropic MLS (fit-metric only)";
  }
  return true;
}

/*****************************************************************************
 * evaluateBicubicHermite
 * ---------------------------------------------------------------------------
 * Purpose
 *   Evaluates a bicubic Hermite patch at local cell coordinates.
 *****************************************************************************/
static auto evaluateBicubicHermite(const BicubicNodeData &n00,
                                   const BicubicNodeData &n10,
                                   const BicubicNodeData &n01,
                                   const BicubicNodeData &n11, double u,
                                   double v, double hx, double hy) -> double {
  const double h00u = 2.0 * u * u * u - 3.0 * u * u + 1.0;
  const double h10u = u * u * u - 2.0 * u * u + u;
  const double h01u = -2.0 * u * u * u + 3.0 * u * u;
  const double h11u = u * u * u - u * u;

  const double h00v = 2.0 * v * v * v - 3.0 * v * v + 1.0;
  const double h10v = v * v * v - 2.0 * v * v + v;
  const double h01v = -2.0 * v * v * v + 3.0 * v * v;
  const double h11v = v * v * v - v * v;

  return n00.z * h00u * h00v + n10.z * h01u * h00v + n01.z * h00u * h01v +
         n11.z * h01u * h01v +
         hx * (n00.fx * h10u * h00v + n10.fx * h11u * h00v +
               n01.fx * h10u * h01v + n11.fx * h11u * h01v) +
         hy * (n00.fy * h00u * h10v + n10.fy * h01u * h10v +
               n01.fy * h00u * h11v + n11.fy * h01u * h11v) +
         hx * hy *
             (n00.fxy * h10u * h10v + n10.fxy * h11u * h10v +
              n01.fxy * h10u * h11v + n11.fxy * h11u * h11v);
}

/*****************************************************************************
 * buildTPSModel
 * ---------------------------------------------------------------------------
 * Purpose
 *   Builds the TPS working model and its validation diagnostics.
 *****************************************************************************/
static auto buildTPSModel(std::vector<Point3D> pts, TPSModel &model,
                          std::string *diagnostics) -> bool {
  model = TPSModel{};
  if (pts.empty()) {
    if (diagnostics) {
      *diagnostics = "TPS solver received no control points.";
    }
    return false;
  }

  const size_t pointCount = pts.size();
  const double spacing = averagePlanarSpacing(pts);
  model.controlPoints = std::move(pts);
  model.spatialIndex.build(model.controlPoints, std::max(3.0 * spacing, 1e-6));
  model.neighborhoodSize = std::min<size_t>(
      {std::max<size_t>(64U,
                        static_cast<size_t>(
                            3.0 * std::sqrt(static_cast<double>(pointCount)))),
       static_cast<unsigned long long>(160U), pointCount});
  model.suggestedLambda = 1e-9;
  model.hull = buildConvexHull2D(model.controlPoints);

  std::vector<double> learnedLambdas;
  std::vector<double> errors;
  const auto holdout =
      makeDeterministicHoldoutIndices(model.controlPoints.size(), 32U);
  if (holdout.size() >= 8U) {
    std::vector<Point3D> training;
    training.reserve(model.controlPoints.size() - holdout.size());
    std::vector<unsigned char> isHoldout(model.controlPoints.size(), 0U);
    for (size_t idx : holdout)
      isHoldout[idx] = 1U;
    for (size_t i = 0; i < model.controlPoints.size(); ++i) {
      if (isHoldout[i] == 0U) {
        training.push_back(model.controlPoints[i]);
      }
    }
    TPSModel trainingModel;
    trainingModel.controlPoints = training;
    trainingModel.spatialIndex.build(
        trainingModel.controlPoints,
        std::max(3.0 * averagePlanarSpacing(trainingModel.controlPoints),
                 1e-6));
    trainingModel.neighborhoodSize = std::min<size_t>(
        {std::max<size_t>(64U,
                          static_cast<size_t>(
                              3.0 * std::sqrt(static_cast<double>(
                                        trainingModel.controlPoints.size())))),
         static_cast<unsigned long long>(160U),
         trainingModel.controlPoints.size()});
    trainingModel.suggestedLambda = model.suggestedLambda;
    trainingModel.hull = buildConvexHull2D(trainingModel.controlPoints);
    trainingModel.valid = true;
    for (size_t idx : holdout) {
      double pred = 0.0;
      InterpolationDiagnostics d;
      if (solveLocalTPSAtPoint(model.controlPoints[idx].x,
                               model.controlPoints[idx].y, trainingModel, pred,
                               nullptr, &d) &&
          std::isfinite(pred)) {
        errors.push_back(pred - model.controlPoints[idx].z);
        if (d.lambda > 0.0 && std::isfinite(d.lambda)) {
          learnedLambdas.push_back(d.lambda);
        }
      }
    }
  }
  model.validation = makeValidationMetricsFromErrors(std::move(errors));
  if (!learnedLambdas.empty()) {
    std::sort(learnedLambdas.begin(), learnedLambdas.end());
    model.suggestedLambda = learnedLambdas[learnedLambdas.size() / 2U];
  }
  model.valid = true;

  if (diagnostics) {
    std::ostringstream oss;
    oss << "High-accuracy local TPS controls=" << model.controlPoints.size()
        << ", target neighborhood=" << model.neighborhoodSize
        << ", support-mask cells="
        << std::count(model.hull.occupied.begin(), model.hull.occupied.end(),
                      static_cast<unsigned char>(1))
        << ", learned lambda=" << model.suggestedLambda;
    if (model.validation.valid) {
      oss << ", holdout RMSE=" << model.validation.rmse;
    }
    *diagnostics = oss.str();
  }
  return true;
}

/*****************************************************************************
 * solveLocalTPSAtPoint
 * ---------------------------------------------------------------------------
 * Purpose
 *   Evaluates the TPS interpolant at a single query location with diagnostics.
 *****************************************************************************/
static auto
solveLocalTPSAtPoint(double x, double y, const TPSModel &model, double &z,
                     std::string *diagnostics,
                     InterpolationDiagnostics *interpDiagnostics = nullptr)
    -> bool {
  z = 0.0;
  if (interpDiagnostics) {
    *interpDiagnostics = InterpolationDiagnostics{};
  }
  if (!model.valid || model.controlPoints.empty()) {
    return false;
  }

  auto phi = [](double r2) -> double {
    if (r2 <= 1e-24) {
      return 0.0;
    }
    return r2 * std::log(r2);
  };

  auto solvePatchSystem =
      [&](const std::vector<size_t> &neighbors, const LocalFrame2D &frame,
          const std::vector<size_t> &fitSelection, double lambda,
          std::vector<double> &solution, std::vector<double> &un,
          std::vector<double> &vn) -> bool {
    const size_t m = fitSelection.size();
    if (m < 8U) {
      return false;
    }
    const size_t N = m + 3U;
    un.assign(m, 0.0);
    vn.assign(m, 0.0);
    std::vector<double> rhs(N, 0.0);
    for (size_t i = 0; i < m; ++i) {
      const Point3D &p = model.controlPoints[neighbors[fitSelection[i]]];
      toLocalUV(frame, p.x, p.y, un[i], vn[i]);
      rhs[i] = p.z;
    }

    std::vector<double> A(N * N, 0.0);
    for (size_t i = 0; i < m; ++i) {
      A[i * N + i] = lambda;
      A[i * N + m] = 1.0;
      A[i * N + (m + 1U)] = un[i];
      A[i * N + (m + 2U)] = vn[i];
      A[m * N + i] = 1.0;
      A[(m + 1U) * N + i] = un[i];
      A[(m + 2U) * N + i] = vn[i];
    }
    for (size_t i = 0; i < m; ++i) {
      for (size_t j = i + 1U; j < m; ++j) {
        const double du = un[i] - un[j];
        const double dv = vn[i] - vn[j];
        const double value = phi(du * du + dv * dv);
        A[i * N + j] = value;
        A[j * N + i] = value;
      }
    }
    return solveDenseLinearSystem(A, rhs, N, solution) &&
           std::all_of(solution.begin(), solution.end(),
                       [](double v) -> bool { return std::isfinite(v); });
  };

  auto evaluateTPSPatch =
      [&](const std::vector<double> &solution, const std::vector<double> &un,
          const std::vector<double> &vn, const LocalFrame2D &frame, double qx,
          double qy) -> double {
    const size_t m = un.size();
    const size_t N = m + 3U;
    if (solution.size() < N) {
      return std::numeric_limits<double>::quiet_NaN();
    }
    double uq = 0.0;
    double vq = 0.0;
    toLocalUV(frame, qx, qy, uq, vq);
    double value = solution[m] + solution[m + 1U] * uq + solution[m + 2U] * vq;
    for (size_t i = 0; i < m; ++i) {
      const double du = un[i] - uq;
      const double dv = vn[i] - vq;
      value += solution[i] * phi(du * du + dv * dv);
    }
    return value;
  };

  auto solveLocalPatch = [&](double anchorX, double anchorY, double &zPatch,
                             double &patchResidual, double &patchSupport,
                             double &chosenLambda,
                             size_t &usedNeighbors) -> bool {
    const size_t probeCount =
        std::min<size_t>(model.controlPoints.size(),
                         std::max<size_t>(model.neighborhoodSize * 2U, 96U));
    const size_t minCount = std::min<size_t>(model.controlPoints.size(), 24U);
    auto candidates = model.spatialIndex.gatherNearest(
        model.controlPoints, anchorX, anchorY, probeCount, minCount);
    if (candidates.size() < 8U) {
      return false;
    }

    LocalFrame2D probeFrame;
    computeLocalPrincipalFrame(model.controlPoints, candidates, anchorX,
                               anchorY, model.spatialIndex.cellSize,
                               probeFrame);
    const double probeSupport = robustSupportRadiusFromLocalUV(
        model.controlPoints, candidates, probeFrame, 1.75);
    const size_t maxKeep =
        std::min<size_t>(model.controlPoints.size(), model.neighborhoodSize);
    auto neighbors =
        selectAdaptiveNeighborhood(model.controlPoints, candidates, probeFrame,
                                   minCount, maxKeep, probeSupport);
    if (neighbors.size() < 8U) {
      return false;
    }

    LocalFrame2D frame;
    computeLocalPrincipalFrame(model.controlPoints, neighbors, anchorX, anchorY,
                               model.spatialIndex.cellSize, frame);
    patchSupport = robustSupportRadiusFromLocalUV(model.controlPoints,
                                                  neighbors, frame, 1.5);
    usedNeighbors = neighbors.size();

    std::vector<size_t> fitSelection;
    std::vector<size_t> validationSelection;
    for (size_t k = 0; k < neighbors.size(); ++k) {
      if (neighbors.size() >= 20U && (k % 5U) == 0U) {
        validationSelection.push_back(k);
      } else {
        fitSelection.push_back(k);
      }
    }
    if (fitSelection.size() < 8U) {
      fitSelection.clear();
      for (size_t k = 0; k < neighbors.size(); ++k) {
        fitSelection.push_back(k);
      }
      validationSelection.clear();
    }

    std::vector<double> bestX;
    std::vector<double> bestUn;
    std::vector<double> bestVn;
    ValidationMetrics bestMetrics;
    double bestScore = std::numeric_limits<double>::infinity();
    bool solved = false;
    for (double lambda : {1e-12, 1e-11, 1e-10, 1e-9, 1e-8, 1e-6}) {
      std::vector<double> xvec;
      std::vector<double> un;
      std::vector<double> vn;
      if (!solvePatchSystem(neighbors, frame, fitSelection, lambda, xvec, un,
                            vn)) {
        continue;
      }

      ValidationMetrics metrics;
      if (!validationSelection.empty()) {
        std::vector<double> errors;
        errors.reserve(validationSelection.size());
        for (size_t localIdx : validationSelection) {
          const Point3D &vp = model.controlPoints[neighbors[localIdx]];
          const double pred = evaluateTPSPatch(xvec, un, vn, frame, vp.x, vp.y);
          if (std::isfinite(pred)) {
            errors.push_back(pred - vp.z);
          }
        }
        metrics = makeValidationMetricsFromErrors(std::move(errors));
      }
      if (!metrics.valid) {
        std::vector<double> errors;
        errors.reserve(fitSelection.size());
        for (size_t localIdx : fitSelection) {
          const Point3D &fp = model.controlPoints[neighbors[localIdx]];
          const double pred = evaluateTPSPatch(xvec, un, vn, frame, fp.x, fp.y);
          if (std::isfinite(pred)) {
            errors.push_back(pred - fp.z);
          }
        }
        metrics = makeValidationMetricsFromErrors(std::move(errors));
      }
      const double score = validationScore(metrics);
      if (score < bestScore) {
        bestScore = score;
        bestMetrics = metrics;
        bestX = std::move(xvec);
        bestUn = std::move(un);
        bestVn = std::move(vn);
        chosenLambda = lambda;
        solved = true;
      }
    }

    if (!solved) {
      return false;
    }

    zPatch = evaluateTPSPatch(bestX, bestUn, bestVn, frame, x, y);
    patchResidual = bestMetrics.valid ? bestMetrics.rmse : 1.0;
    return std::isfinite(zPatch);
  };

  const size_t anchorProbe = std::min<size_t>(model.controlPoints.size(), 8U);
  const auto anchorCandidates = model.spatialIndex.gatherNearest(
      model.controlPoints, x, y, anchorProbe, 1U);
  if (anchorCandidates.empty()) {
    return false;
  }

  double weightedZ = 0.0;
  double weightSum = 0.0;
  double bestSingleZ = 0.0;
  double bestSingleResidual = std::numeric_limits<double>::infinity();
  double bestLambda = model.suggestedLambda;
  size_t successes = 0U;
  size_t meanNeighbors = 0U;

  const size_t blendCount = std::min<size_t>(4U, anchorCandidates.size());
  for (size_t a = 0; a < blendCount; ++a) {
    const Point3D &anchor = model.controlPoints[anchorCandidates[a]];
    double zPatch = 0.0;
    double residual = std::numeric_limits<double>::infinity();
    double supportRadius = 1.5;
    double lambda = model.suggestedLambda;
    size_t usedNeighbors = 0U;
    if (!solveLocalPatch(anchor.x, anchor.y, zPatch, residual, supportRadius,
                         lambda, usedNeighbors)) {
      continue;
    }

    const double d = pointDistanceXY(x, y, anchor.x, anchor.y);
    const double supportXY = std::max(
        1.0, supportRadius * std::max(model.spatialIndex.cellSize, 1e-6));
    const double w = std::max(compactWeightWendland(d, supportXY * 1.001) /
                                  std::sqrt(1e-8 + residual),
                              1e-12);
    weightedZ += w * zPatch;
    weightSum += w;
    ++successes;
    meanNeighbors += usedNeighbors;

    if (residual < bestSingleResidual) {
      bestSingleResidual = residual;
      bestSingleZ = zPatch;
      bestLambda = lambda;
    }
  }

  bool fallbackUsed = false;
  if (weightSum > 0.0) {
    z = weightedZ / weightSum;
  } else if (std::isfinite(bestSingleZ)) {
    z = bestSingleZ;
  } else {
    z = inverseDistanceFallback(model.controlPoints, model.spatialIndex, x, y);
    fallbackUsed = true;
  }

  const bool outsideHull =
      model.hull.valid ? !pointInConvexHull2D(model.hull, x, y) : false;
  const double outsideDistance =
      outsideHull ? distanceToConvexHull2D(model.hull, x, y) : 0.0;
  const double boundaryPenalty =
      outsideHull
          ? (1.0 / (1.0 + outsideDistance /
                              std::max(model.spatialIndex.cellSize, 1e-6)))
          : 1.0;
  z *= 1.0;

  if (interpDiagnostics) {
    interpDiagnostics->confidence = calibratedConfidenceFromObservedError(
        boundaryPenalty / (1.0 + 4.0 * std::max(bestSingleResidual, 0.0)),
        std::max(bestSingleResidual, 0.0), model.validation);
    interpDiagnostics->residual =
        std::isfinite(bestSingleResidual) ? bestSingleResidual : 1.0;
    interpDiagnostics->lambda = bestLambda;
    interpDiagnostics->boundaryPenalty = boundaryPenalty;
    interpDiagnostics->outsideDistance = outsideDistance;
    interpDiagnostics->neighborCount =
        (successes > 0U) ? (meanNeighbors / successes) : 0U;
    interpDiagnostics->patchCount = successes;
    interpDiagnostics->fallbackUsed = fallbackUsed;
    interpDiagnostics->outsideHull = outsideHull;
    std::ostringstream oss;
    oss << "TPS patches=" << successes << ", lambda=" << bestLambda
        << ", residual=" << bestSingleResidual
        << (outsideHull ? ", outside-hull" : "");
    interpDiagnostics->summary = oss.str();
  }
  if (diagnostics) {
    if (interpDiagnostics) {
      *diagnostics = interpDiagnostics->summary;
    } else {
      std::ostringstream oss;
      oss << "TPS patches=" << successes << ", residual=" << bestSingleResidual;
      *diagnostics = oss.str();
    }
  }
  return std::isfinite(z);
}

/****************************************************************************
 * 8) Bicubic interpolation via local polynomial node estimation
 ****************************************************************************/
static auto formatPointXYZLine(const Point3D &p, int precision) -> std::string {
  std::ostringstream oss;
  oss << std::fixed << std::setprecision(precision) << p.x << ' ' << p.y << ' '
      << p.z << "\n";
  return oss.str();
}

/*****************************************************************************
 * checkedMultiplySizeT
 * ---------------------------------------------------------------------------
 * Purpose
 *   Safely multiplies size_t dimensions and detects overflow before allocation.
 *****************************************************************************/
static auto checkedMultiplySizeT(size_t a, size_t b, size_t &out) -> bool {
  if (a == 0 || b == 0) {
    out = 0;
    return true;
  }
  if (a > (std::numeric_limits<size_t>::max() / b)) {
    out = 0;
    return false;
  }
  out = a * b;
  return true;
}

/*****************************************************************************
 * buildGridDefinition
 * ---------------------------------------------------------------------------
 * Purpose
 *   Builds the regular output-grid definition from bounds and spacing.
 *****************************************************************************/
static auto buildGridDefinition(const std::vector<Point3D> &points,
                                double gridSpacing, GridDefinition &grid,
                                std::string *errorMessage) -> bool {
  grid = GridDefinition{};
  if (points.empty()) {
    if (errorMessage) {
      *errorMessage =
          "Error: no points available to define interpolation grid.";
    }
    return false;
  }
  if (!(gridSpacing > 0.0) || !std::isfinite(gridSpacing)) {
    if (errorMessage) {
      *errorMessage = "Error: gridSpacing must be a positive finite value.";
    }
    return false;
  }

  double dataXMin, dataXMax, dataYMin, dataYMax;
  computeBoundingBox(points, dataXMin, dataXMax, dataYMin, dataYMax);

  const double margin = 1.5 * gridSpacing;
  grid.xMin = dataXMin - margin;
  grid.xMax = dataXMax + margin;
  grid.yMin = dataYMin - margin;
  grid.yMax = dataYMax + margin;
  grid.gridSpacing = gridSpacing;

  const double width = std::max(grid.xMax - grid.xMin, 1.0);
  const double height = std::max(grid.yMax - grid.yMin, 1.0);
  grid.nx = static_cast<size_t>(std::ceil(width / gridSpacing)) + 1U;
  grid.ny = static_cast<size_t>(std::ceil(height / gridSpacing)) + 1U;
  if (grid.nx < 2U || grid.ny < 2U) {
    if (errorMessage) {
      *errorMessage = "Error: interpolation grid is degenerate.";
    }
    return false;
  }
  if (!checkedMultiplySizeT(grid.nx, grid.ny, grid.totalPoints)) {
    if (errorMessage) {
      *errorMessage = "Error: interpolation grid dimensions overflow size_t.";
    }
    return false;
  }

  grid.valid = true;
  return true;
}

struct BicubicStreamingContext {
  const std::vector<Point3D> *points = nullptr;
  SpatialHash2D spatialIndex;
  ConvexHull2D hull;
  GridDefinition grid;
  TuningChoice tuning;
  ValidationMetrics validation;
};

struct ClampedMLSStreamingContext {
  const std::vector<Point3D> *points = nullptr;
  SpatialHash2D spatialIndex;
  ConvexHull2D hull;
  GridDefinition grid;
  TuningChoice tuning;
  ValidationMetrics validation;
};

struct TPSStreamingContext {
  TPSModel model;
  GridDefinition grid;
};

/*****************************************************************************
 * evaluateGlobalMLSHoldout
 * ---------------------------------------------------------------------------
 * Purpose
 *   Evaluates holdout accuracy for a candidate global MLS/bicubic tuning
 * choice.
 *****************************************************************************/
static auto evaluateGlobalMLSHoldout(const std::vector<Point3D> &points,
                                     const SpatialHash2D &spatialIndex,
                                     double gridSpacing,
                                     const TuningChoice &tuning,
                                     bool clampedMode) -> ValidationMetrics {
  std::vector<double> errors;
  const auto holdout = makeDeterministicHoldoutIndices(points.size(), 48U);
  errors.reserve(holdout.size());
  for (size_t idx : holdout) {
    BicubicNodeData node;
    InterpolationDiagnostics diag;
    if (!fitLocalPolynomialNode(points, spatialIndex, points[idx].x,
                                points[idx].y, gridSpacing, node, &tuning,
                                &diag, idx)) {
      continue;
    }
    double pred = node.z;
    if (clampedMode) {
      auto neighbors =
          spatialIndex.gatherNearest(points, points[idx].x, points[idx].y,
                                     std::min<size_t>(points.size(), 32U), 1U);
      double zMin = std::numeric_limits<double>::infinity();
      double zMax = -std::numeric_limits<double>::infinity();
      for (size_t j : neighbors) {
        if (j == idx)
          continue;
        zMin = std::min(zMin, points[j].z);
        zMax = std::max(zMax, points[j].z);
      }
      if (std::isfinite(zMin) && std::isfinite(zMax)) {
        pred = std::max(zMin, std::min(zMax, pred));
      }
    }
    if (std::isfinite(pred)) {
      errors.push_back(pred - points[idx].z);
    }
  }
  return makeValidationMetricsFromErrors(std::move(errors));
}

/*****************************************************************************
 * tuneBicubicChoice
 * ---------------------------------------------------------------------------
 * Purpose
 *   Selects bicubic/MLS tuning parameters from deterministic candidate testing.
 *****************************************************************************/
static auto tuneBicubicChoice(const std::vector<Point3D> &points,
                              const SpatialHash2D &spatialIndex,
                              double gridSpacing) -> TuningChoice {
  TuningChoice best = defaultTuningChoice();
  if (points.empty()) {
    return best;
  }

  std::vector<size_t> neighborhoodOptions = {48U, 64U, 96U, 128U};
  std::vector<double> supportOptions = {1.75, 2.25, 2.75};
  std::vector<double> ridgeOptions = {1.0, 10.0, 50.0};
  std::vector<int> orderOptions = {2, 3};

  double bestScore = std::numeric_limits<double>::infinity();
  for (size_t neigh : neighborhoodOptions) {
    neigh = std::min(neigh, points.size());
    for (double support : supportOptions) {
      for (double ridge : ridgeOptions) {
        for (int order : orderOptions) {
          if (order >= 3 && points.size() < 18U) {
            continue;
          }
          TuningChoice candidate;
          candidate.neighborhoodTarget = std::max<size_t>(neigh, 24U);
          candidate.supportMultiplier = support;
          candidate.ridgeFactor = ridge;
          candidate.basisOrder = order;
          candidate.metrics = evaluateGlobalMLSHoldout(
              points, spatialIndex, gridSpacing, candidate, false);
          if (!candidate.metrics.valid) {
            continue;
          }
          const double score = validationScore(candidate.metrics);
          if (score < bestScore) {
            bestScore = score;
            best = candidate;
          }
        }
      }
    }
  }
  return best;
}

/*****************************************************************************
 * buildBicubicStreamingContext
 * ---------------------------------------------------------------------------
 * Purpose
 *   Builds the bicubic streaming context used during regular-grid output
 * generation.
 *****************************************************************************/
static auto buildBicubicStreamingContext(const std::vector<Point3D> &points,
                                         double gridSpacing,
                                         double boundaryCellSize,
                                         BicubicStreamingContext &ctx,
                                         std::string *errorMessage) -> bool {
  ctx = BicubicStreamingContext{};
  if (!buildGridDefinition(points, gridSpacing, ctx.grid, errorMessage)) {
    return false;
  }
  const double spacing = averagePlanarSpacing(points);
  ctx.points = &points;
  ctx.spatialIndex.build(points, std::max(2.5 * spacing, 1.5 * gridSpacing));
  ctx.hull = buildConvexHull2D(points, boundaryCellSize);
  ctx.tuning = tuneBicubicChoice(points, ctx.spatialIndex, gridSpacing);
  ctx.validation = evaluateGlobalMLSHoldout(points, ctx.spatialIndex,
                                            gridSpacing, ctx.tuning, false);
  return true;
}

/*****************************************************************************
 * buildClampedMLSStreamingContext
 * ---------------------------------------------------------------------------
 * Purpose
 *   Builds the clamped-MLS streaming context used during regular-grid output
 * generation.
 *****************************************************************************/
static auto buildClampedMLSStreamingContext(const std::vector<Point3D> &points,
                                            double gridSpacing,
                                            double boundaryCellSize,
                                            ClampedMLSStreamingContext &ctx,
                                            std::string *errorMessage) -> bool {
  ctx = ClampedMLSStreamingContext{};
  if (!buildGridDefinition(points, gridSpacing, ctx.grid, errorMessage)) {
    return false;
  }
  const double spacing = averagePlanarSpacing(points);
  ctx.points = &points;
  ctx.spatialIndex.build(points, std::max(2.5 * spacing, 1.5 * gridSpacing));
  ctx.hull = buildConvexHull2D(points, boundaryCellSize);
  ctx.tuning = tuneBicubicChoice(points, ctx.spatialIndex, gridSpacing);
  ctx.validation = evaluateGlobalMLSHoldout(points, ctx.spatialIndex,
                                            gridSpacing, ctx.tuning, true);
  return true;
}

/*****************************************************************************
 * buildTPSStreamingContext
 * ---------------------------------------------------------------------------
 * Purpose
 *   Builds the TPS streaming context used during regular-grid output
 * generation.
 *****************************************************************************/
static auto buildTPSStreamingContext(std::vector<Point3D> controlPoints,
                                     const GridDefinition &extentGrid,
                                     double boundaryCellSize,
                                     TPSStreamingContext &ctx,
                                     std::string *errorMessage) -> bool {
  ctx = TPSStreamingContext{};
  ctx.grid = extentGrid;
  if (!ctx.grid.valid) {
    if (errorMessage) {
      *errorMessage = "Error: invalid TPS grid definition.";
    }
    return false;
  }
  std::string diagnostics;
  if (!buildTPSModel(std::move(controlPoints), ctx.model, &diagnostics)) {
    if (errorMessage) {
      *errorMessage = diagnostics.empty() ? "Error: failed to build TPS model."
                                          : diagnostics;
    }
    return false;
  }
  ctx.model.hull = buildConvexHull2D(ctx.model.controlPoints, boundaryCellSize);
  return true;
}

/*****************************************************************************
 * evaluateBicubicStreamPointDetailed
 * ---------------------------------------------------------------------------
 * Purpose
 *   Evaluates one bicubic grid point together with detailed interpolation
 * diagnostics.
 *****************************************************************************/
static auto evaluateBicubicStreamPointDetailed(
    const BicubicStreamingContext &ctx, size_t i, size_t j, Point3D &outPoint,
    InterpolationDiagnostics *diagnostics, std::string *errorMessage) -> bool {
  if (diagnostics) {
    *diagnostics = InterpolationDiagnostics{};
  }
  if (!ctx.points || !ctx.grid.valid || ctx.grid.nx < 2U || ctx.grid.ny < 2U) {
    if (errorMessage) {
      *errorMessage = "Error: bicubic context is invalid.";
    }
    return false;
  }

  const double gs = ctx.grid.gridSpacing;
  const double gx = ctx.grid.xMin + static_cast<double>(i) * gs;
  const double gy = ctx.grid.yMin + static_cast<double>(j) * gs;

  const size_t ci = (i + 1U < ctx.grid.nx) ? i : (ctx.grid.nx - 2U);
  const size_t cj = (j + 1U < ctx.grid.ny) ? j : (ctx.grid.ny - 2U);

  const double xCell0 = ctx.grid.xMin + static_cast<double>(ci) * gs;
  const double yCell0 = ctx.grid.yMin + static_cast<double>(cj) * gs;
  const double u = std::min(std::max((gx - xCell0) / gs, 0.0), 1.0);
  const double v = std::min(std::max((gy - yCell0) / gs, 0.0), 1.0);

  BicubicNodeData n00, n10, n01, n11;
  InterpolationDiagnostics d00, d10, d01, d11;
  const bool ok00 =
      fitLocalPolynomialNode(*ctx.points, ctx.spatialIndex, xCell0, yCell0, gs,
                             n00, &ctx.tuning, &d00);
  const bool ok10 =
      fitLocalPolynomialNode(*ctx.points, ctx.spatialIndex, xCell0 + gs, yCell0,
                             gs, n10, &ctx.tuning, &d10);
  const bool ok01 =
      fitLocalPolynomialNode(*ctx.points, ctx.spatialIndex, xCell0, yCell0 + gs,
                             gs, n01, &ctx.tuning, &d01);
  const bool ok11 =
      fitLocalPolynomialNode(*ctx.points, ctx.spatialIndex, xCell0 + gs,
                             yCell0 + gs, gs, n11, &ctx.tuning, &d11);

  double gz = 0.0;
  bool fallback = false;
  if (ok00 && ok10 && ok01 && ok11 && n00.valid && n10.valid && n01.valid &&
      n11.valid) {
    limitHermiteCellDerivatives(n00, n10, n01, n11, gs, gs);
    gz = evaluateBicubicHermite(n00, n10, n01, n11, u, v, gs, gs);
    const double zMin = std::min({n00.z, n10.z, n01.z, n11.z});
    const double zMax = std::max({n00.z, n10.z, n01.z, n11.z});
    const double pad = 0.05 * std::max(zMax - zMin, 1e-8);
    gz = std::max(zMin - pad, std::min(zMax + pad, gz));
  } else {
    gz = inverseDistanceFallback(*ctx.points, ctx.spatialIndex, gx, gy);
    fallback = true;
  }

  if (!std::isfinite(gz)) {
    if (errorMessage) {
      *errorMessage =
          "Error: bicubic Hermite evaluation produced a non-finite value.";
    }
    return false;
  }

  const bool outsideHull =
      ctx.hull.valid ? !pointInConvexHull2D(ctx.hull, gx, gy) : false;
  const double outsideDistance =
      outsideHull ? distanceToConvexHull2D(ctx.hull, gx, gy) : 0.0;
  const double boundaryPenalty =
      outsideHull ? (1.0 / (1.0 + outsideDistance / std::max(gs, 1e-6))) : 1.0;

  outPoint = {gx, gy, gz};
  if (diagnostics) {
    diagnostics->residual =
        0.25 * (d00.residual + d10.residual + d01.residual + d11.residual);
    const double meanNodeConfidence = 0.25 * (d00.confidence + d10.confidence +
                                              d01.confidence + d11.confidence);
    diagnostics->confidence =
        std::clamp(meanNodeConfidence * boundaryPenalty, 0.02, 1.0);
    diagnostics->neighborCount = (d00.neighborCount + d10.neighborCount +
                                  d01.neighborCount + d11.neighborCount) /
                                 4U;
    diagnostics->patchCount = 4U;
    diagnostics->fallbackUsed = fallback;
    diagnostics->outsideHull = outsideHull;
    diagnostics->outsideDistance = outsideDistance;
    diagnostics->boundaryPenalty = boundaryPenalty;
    diagnostics->summary =
        fallback ? "bicubic fallback" : "bicubic Hermite with robust MLS";
  }
  return std::isfinite(outPoint.x) && std::isfinite(outPoint.y) &&
         std::isfinite(outPoint.z);
}

/*****************************************************************************
 * evaluateBicubicStreamPoint
 * ---------------------------------------------------------------------------
 * Purpose
 *   Evaluates one bicubic grid point for the regular output stream.
 *****************************************************************************/
static auto evaluateBicubicStreamPoint(const BicubicStreamingContext &ctx,
                                       size_t i, size_t j, Point3D &outPoint,
                                       std::string *errorMessage) -> bool {
  return evaluateBicubicStreamPointDetailed(ctx, i, j, outPoint, nullptr,
                                            errorMessage);
}

/*****************************************************************************
 * evaluateClampedMLSStreamPointDetailed
 * ---------------------------------------------------------------------------
 * Purpose
 *   Evaluates one clamped-MLS grid point together with detailed diagnostics.
 *****************************************************************************/
static auto
evaluateClampedMLSStreamPointDetailed(const ClampedMLSStreamingContext &ctx,
                                      size_t i, size_t j, Point3D &outPoint,
                                      InterpolationDiagnostics *diagnostics,
                                      std::string *errorMessage) -> bool {
  if (diagnostics) {
    *diagnostics = InterpolationDiagnostics{};
  }
  if (!ctx.points || !ctx.grid.valid) {
    if (errorMessage) {
      *errorMessage = "Error: clamped local MLS context is invalid.";
    }
    return false;
  }
  const double gx =
      ctx.grid.xMin + static_cast<double>(i) * ctx.grid.gridSpacing;
  const double gy =
      ctx.grid.yMin + static_cast<double>(j) * ctx.grid.gridSpacing;

  BicubicNodeData node;
  InterpolationDiagnostics localDiag;
  if (!fitLocalPolynomialNode(*ctx.points, ctx.spatialIndex, gx, gy,
                              ctx.grid.gridSpacing, node, &ctx.tuning,
                              &localDiag)) {
    if (errorMessage) {
      *errorMessage = "Error: local monotone MLS evaluation failed.";
    }
    return false;
  }

  auto neighbors = ctx.spatialIndex.gatherNearest(
      *ctx.points, gx, gy, std::min<size_t>(ctx.points->size(), 32U), 1U);
  double zMin = std::numeric_limits<double>::infinity();
  double zMax = -std::numeric_limits<double>::infinity();
  for (size_t idx : neighbors) {
    zMin = std::min(zMin, (*ctx.points)[idx].z);
    zMax = std::max(zMax, (*ctx.points)[idx].z);
  }
  double gz = node.z;
  bool clamped = false;
  if (std::isfinite(zMin) && std::isfinite(zMax)) {
    const double clampedZ = std::max(zMin, std::min(zMax, gz));
    clamped = (clampedZ != gz);
    gz = clampedZ;
  }

  const bool outsideHull =
      ctx.hull.valid ? !pointInConvexHull2D(ctx.hull, gx, gy) : false;
  const double outsideDistance =
      outsideHull ? distanceToConvexHull2D(ctx.hull, gx, gy) : 0.0;
  const double boundaryPenalty =
      outsideHull ? (1.0 / (1.0 + outsideDistance /
                                      std::max(ctx.grid.gridSpacing, 1e-6)))
                  : 1.0;

  outPoint = {gx, gy, gz};
  if (diagnostics) {
    *diagnostics = localDiag;
    diagnostics->outsideHull = outsideHull;
    diagnostics->outsideDistance = outsideDistance;
    diagnostics->boundaryPenalty = boundaryPenalty;
    diagnostics->confidence = calibratedConfidenceFromObservedError(
        diagnostics->confidence * boundaryPenalty, diagnostics->residual,
        ctx.validation);
    diagnostics->clampedLocal = clamped;
    diagnostics->summary =
        clamped ? "clamped local MLS (clamped)" : "clamped local MLS";
  }
  return std::isfinite(outPoint.z);
}

/*****************************************************************************
 * evaluateClampedMLSStreamPoint
 * ---------------------------------------------------------------------------
 * Purpose
 *   Evaluates one clamped-MLS grid point for the regular output stream.
 *****************************************************************************/
static auto evaluateClampedMLSStreamPoint(const ClampedMLSStreamingContext &ctx,
                                          size_t i, size_t j, Point3D &outPoint,
                                          std::string *errorMessage) -> bool {
  return evaluateClampedMLSStreamPointDetailed(ctx, i, j, outPoint, nullptr,
                                               errorMessage);
}

/*****************************************************************************
 * evaluateTPSStreamPointDetailed
 * ---------------------------------------------------------------------------
 * Purpose
 *   Evaluates one TPS grid point together with detailed diagnostics.
 *****************************************************************************/
static auto evaluateTPSStreamPointDetailed(
    const TPSStreamingContext &ctx, size_t i, size_t j, Point3D &outPoint,
    InterpolationDiagnostics *diagnostics, std::string *errorMessage) -> bool {
  if (diagnostics) {
    *diagnostics = InterpolationDiagnostics{};
  }
  if (!ctx.grid.valid || !ctx.model.valid) {
    if (errorMessage) {
      *errorMessage = "Error: TPS context is invalid.";
    }
    return false;
  }

  const double gx =
      ctx.grid.xMin + static_cast<double>(i) * ctx.grid.gridSpacing;
  const double gy =
      ctx.grid.yMin + static_cast<double>(j) * ctx.grid.gridSpacing;
  double gz = 0.0;
  InterpolationDiagnostics localDiag;
  if (!solveLocalTPSAtPoint(gx, gy, ctx.model, gz, nullptr, &localDiag) ||
      !std::isfinite(gz)) {
    gz = inverseDistanceFallback(ctx.model.controlPoints,
                                 ctx.model.spatialIndex, gx, gy);
    localDiag.fallbackUsed = true;
    localDiag.confidence = 0.05;
  }
  outPoint = {gx, gy, gz};
  if (diagnostics) {
    *diagnostics = localDiag;
  }
  return std::isfinite(outPoint.x) && std::isfinite(outPoint.y) &&
         std::isfinite(outPoint.z);
}

/*****************************************************************************
 * evaluateTPSStreamPoint
 * ---------------------------------------------------------------------------
 * Purpose
 *   Evaluates one TPS grid point for the regular output stream.
 *****************************************************************************/
static auto evaluateTPSStreamPoint(const TPSStreamingContext &ctx, size_t i,
                                   size_t j, Point3D &outPoint,
                                   std::string *errorMessage) -> bool {
  return evaluateTPSStreamPointDetailed(ctx, i, j, outPoint, nullptr,
                                        errorMessage);
}

/*****************************************************************************
 * maybeReportRowProgress
 * ---------------------------------------------------------------------------
 * Purpose
 *   Reports periodic row-wise progress during streamed output generation.
 *****************************************************************************/
static void maybeReportRowProgress(
    const std::function<void(const std::string &)> &progressCallback,
    const char *phaseLabel, size_t rowIndex, size_t rowCount) {
  if (!progressCallback || rowCount == 0U) {
    return;
  }
  const size_t step = std::max<size_t>(1U, rowCount / 50U);
  const bool shouldReport = (rowIndex == 0U) || ((rowIndex + 1U) == rowCount) ||
                            (((rowIndex + 1U) % step) == 0U);
  if (!shouldReport) {
    return;
  }

  const double pct = (100.0 * static_cast<double>(rowIndex + 1U)) /
                     static_cast<double>(rowCount);
  std::ostringstream oss;
  oss << phaseLabel << " row " << (rowIndex + 1U) << "/" << rowCount << " ("
      << std::fixed << std::setprecision(1) << pct << "%)";
  progressCallback(oss.str());
}

/*****************************************************************************
 * writeGridXYZStream
 * ---------------------------------------------------------------------------
 * Purpose
 *   Streams the interpolated regular grid to XYZ without storing the full grid
 * in memory.
 *****************************************************************************/
static auto writeGridXYZStream(
    const std::string &outputFileName, const GridDefinition &grid,
    const ConvexHull2D *boundary, int precision,
    const std::function<bool(size_t, size_t, Point3D &, std::string *)>
        &evaluator,
    const std::function<void(const std::string &)> &progressCallback,
    std::string *errorMessage) -> bool {
  std::ofstream outFile(outputFileName, std::ios::out | std::ios::binary);
  if (!outFile.is_open()) {
    if (errorMessage) {
      *errorMessage = "Error creating grid XYZ file: " + outputFileName;
    }
    return false;
  }
  std::vector<char> fileBuffer(8U * 1024U * 1024U);
  outFile.rdbuf()->pubsetbuf(fileBuffer.data(),
                             static_cast<std::streamsize>(fileBuffer.size()));
  outFile << std::fixed << std::setprecision(precision);

  if (progressCallback) {
    std::ostringstream msg;
    msg << "grid.xyz writer opened. Rows=" << grid.nx << ", cols=" << grid.ny
        << ", planned nodes=" << grid.totalPoints;
    progressCallback(msg.str());
  }

  std::string buffer;
  buffer.reserve(4U * 1024U * 1024U);
  size_t writtenCount = 0U;

  for (size_t i = 0; i < grid.nx; ++i) {
    for (size_t j = 0; j < grid.ny; ++j) {
      const double gx = grid.xMin + static_cast<double>(i) * grid.gridSpacing;
      const double gy = grid.yMin + static_cast<double>(j) * grid.gridSpacing;
      if (boundary && boundary->valid &&
          !pointInConvexHull2D(*boundary, gx, gy)) {
        continue;
      }

      Point3D p{};
      std::string localError;
      if (!evaluator(i, j, p, &localError)) {
        if (errorMessage) {
          *errorMessage = localError.empty()
                              ? ("Error evaluating interpolation point for: " +
                                 outputFileName)
                              : localError;
        }
        return false;
      }
      buffer += formatPointXYZLine(p, precision);
      ++writtenCount;
      if (buffer.size() >= (4U * 1024U * 1024U)) {
        outFile << buffer;
        if (!outFile.good()) {
          if (errorMessage) {
            *errorMessage = "Error writing grid XYZ file: " + outputFileName;
          }
          return false;
        }
        buffer.clear();
      }
    }
    maybeReportRowProgress(progressCallback, "Writing grid.xyz:", i, grid.nx);
  }

  outFile << buffer;
  const bool ok = outFile.good();
  outFile.close();
  if (ok && progressCallback) {
    std::ostringstream msg;
    msg << "grid.xyz written nodes inside boundary: " << writtenCount;
    progressCallback(msg.str());
  }
  if (!ok && errorMessage) {
    *errorMessage = "Error writing grid XYZ file: " + outputFileName;
  }
  return ok;
}

/*****************************************************************************
 * writeConfidenceXYZStream
 * ---------------------------------------------------------------------------
 * Purpose
 *   Streams confidence and diagnostic fields for each grid node to the
 * confidence XYZ output.
 *****************************************************************************/
static auto writeConfidenceXYZStream(
    const std::string &outputFileName, const GridDefinition &grid,
    const ConvexHull2D *boundary, int precision,
    const std::function<bool(size_t, size_t, Point3D &,
                             InterpolationDiagnostics &, std::string *)>
        &evaluator,
    const std::function<void(const std::string &)> &progressCallback,
    ConfidenceSummary &summary, std::string *errorMessage) -> bool {
  summary = ConfidenceSummary{};
  std::ofstream outFile(outputFileName, std::ios::out | std::ios::binary);
  if (!outFile.is_open()) {
    if (errorMessage) {
      *errorMessage = "Error creating confidence XYZ file: " + outputFileName;
    }
    return false;
  }
  outFile << std::fixed << std::setprecision(precision);
  outFile << "# x y z confidence residual lambda boundaryPenalty fallback "
             "outsideHull clampedLocal neighborCount patchCount\n";

  for (size_t i = 0; i < grid.nx; ++i) {
    for (size_t j = 0; j < grid.ny; ++j) {
      const double gx = grid.xMin + static_cast<double>(i) * grid.gridSpacing;
      const double gy = grid.yMin + static_cast<double>(j) * grid.gridSpacing;
      if (boundary && boundary->valid &&
          !pointInConvexHull2D(*boundary, gx, gy)) {
        continue;
      }

      Point3D p{};
      InterpolationDiagnostics d;
      std::string localError;
      if (!evaluator(i, j, p, d, &localError)) {
        if (errorMessage) {
          *errorMessage =
              localError.empty()
                  ? ("Error evaluating confidence point for: " + outputFileName)
                  : localError;
        }
        return false;
      }
      outFile << p.x << ' ' << p.y << ' ' << p.z << ' ' << d.confidence << ' '
              << d.residual << ' ' << d.lambda << ' ' << d.boundaryPenalty
              << ' ' << (d.fallbackUsed ? 1 : 0) << ' '
              << (d.outsideHull ? 1 : 0) << ' ' << (d.clampedLocal ? 1 : 0)
              << ' ' << d.neighborCount << ' ' << d.patchCount << "\n";
      ++summary.count;
      summary.meanConfidence += d.confidence;
      summary.meanResidual += d.residual;
      summary.fallbackCount += d.fallbackUsed ? 1U : 0U;
      summary.outsideHullCount += d.outsideHull ? 1U : 0U;
      summary.clampedLocalCount += d.clampedLocal ? 1U : 0U;
    }
    maybeReportRowProgress(progressCallback, "Writing confidence.xyz:", i,
                           grid.nx);
  }
  if (summary.count > 0U) {
    summary.meanConfidence /= static_cast<double>(summary.count);
    summary.meanResidual /= static_cast<double>(summary.count);
  }
  const bool ok = outFile.good();
  outFile.close();
  if (ok && progressCallback) {
    std::ostringstream msg;
    msg << "confidence.xyz written nodes inside boundary: " << summary.count;
    progressCallback(msg.str());
  }
  if (!ok && errorMessage) {
    *errorMessage = "Error writing confidence XYZ file: " + outputFileName;
  }
  return ok;
}

/*****************************************************************************
 * computeBoundingBoxFromXYZFile
 * ---------------------------------------------------------------------------
 * Purpose
 *   Computes the source-data bounding box directly from the XYZ file when
 * needed for streamed writing.
 *****************************************************************************/
static auto computeBoundingBoxFromXYZFile(const std::string &inputFileName,
                                          double &xmin, double &xmax,
                                          double &ymin, double &ymax) -> bool {
  std::ifstream inFile(inputFileName);
  if (!inFile.is_open()) {
    return false;
  }

  xmin = std::numeric_limits<double>::max();
  xmax = -std::numeric_limits<double>::max();
  ymin = std::numeric_limits<double>::max();
  ymax = -std::numeric_limits<double>::max();

  std::string line;
  bool found = false;
  while (std::getline(inFile, line)) {
    if (line.empty() || line[0] == '#') {
      continue;
    }
    std::replace(line.begin(), line.end(), ',', ' ');
    std::istringstream ss(line);
    Point3D p{};
    if (!(ss >> p.x >> p.y >> p.z)) {
      continue;
    }
    xmin = std::min(xmin, p.x);
    xmax = std::max(xmax, p.x);
    ymin = std::min(ymin, p.y);
    ymax = std::max(ymax, p.y);
    found = true;
  }
  return found;
}

/*****************************************************************************
 * appendBoundaryLoopToDXFBuffer
 * ---------------------------------------------------------------------------
 * Purpose
 *   Appends one boundary loop to the DXF output buffer.
 *****************************************************************************/
static void appendBoundaryLoopToDXFBuffer(std::string &buffer,
                                          const std::vector<Point3D> &loop,
                                          int precision) {
  if (loop.size() < 2U) {
    return;
  }

  std::ostringstream ss;
  ss << std::fixed << std::setprecision(precision);

  ss << "0\nPOLYLINE\n"
     << "8\nboundary\n"
     << "66\n1\n"
     << "70\n1\n"
     << "10\n0.0\n20\n0.0\n30\n0.0\n";

  for (const auto &p : loop) {
    ss << "0\nVERTEX\n"
       << "8\nboundary\n"
       << "10\n"
       << p.x << "\n"
       << "20\n"
       << p.y << "\n"
       << "30\n0.0\n"
       << "70\n0\n";
  }

  ss << "0\nSEQEND\n8\nboundary\n";
  buffer += ss.str();
}

/*****************************************************************************
 * writeDXFStreamFromXYZFile
 * ---------------------------------------------------------------------------
 * Purpose
 *   Writes the DXF output by streaming filtered points and grid nodes from
 * file/context data.
 *****************************************************************************/
static auto writeDXFStreamFromXYZFile(
    const std::string &outputFileName, const std::string &xyzInputFileName,
    const GridDefinition &grid, const ConvexHull2D *boundary, int precision,
    int pdmode, bool hasGrid,
    const std::function<bool(size_t, size_t, Point3D &, std::string *)>
        &evaluator,
    const std::function<void(const std::string &)> &progressCallback,
    std::string *errorMessage) -> bool {
  double xmin = 0.0;
  double xmax = 0.0;
  double ymin = 0.0;
  double ymax = 0.0;
  if (!computeBoundingBoxFromXYZFile(xyzInputFileName, xmin, xmax, ymin,
                                     ymax)) {
    if (errorMessage) {
      *errorMessage = "Error reading filtered XYZ file for DXF bounding box: " +
                      xyzInputFileName;
    }
    return false;
  }

  if (hasGrid && grid.valid) {
    const double gridXMax =
        grid.xMin + static_cast<double>(grid.nx - 1U) * grid.gridSpacing;
    const double gridYMax =
        grid.yMin + static_cast<double>(grid.ny - 1U) * grid.gridSpacing;
    xmin = std::min(xmin, grid.xMin);
    xmax = std::max(xmax, gridXMax);
    ymin = std::min(ymin, grid.yMin);
    ymax = std::max(ymax, gridYMax);
  }
  if (boundary && boundary->valid) {
    for (const auto &loop : boundary->loops) {
      for (const auto &p : loop) {
        xmin = std::min(xmin, p.x);
        xmax = std::max(xmax, p.x);
        ymin = std::min(ymin, p.y);
        ymax = std::max(ymax, p.y);
      }
    }
  }

  std::ofstream outFile(outputFileName, std::ios::out | std::ios::binary);
  if (!outFile.is_open()) {
    if (errorMessage) {
      *errorMessage = "Error creating DXF file: " + outputFileName;
    }
    return false;
  }
  std::vector<char> fileBuffer(8U * 1024U * 1024U);
  outFile.rdbuf()->pubsetbuf(fileBuffer.data(),
                             static_cast<std::streamsize>(fileBuffer.size()));
  outFile << std::fixed << std::setprecision(precision);

  if (progressCallback) {
    std::ostringstream msg;
    msg << "DXF writer opened. Filtered source=" << xyzInputFileName
        << ", grid nodes=" << grid.totalPoints;
    progressCallback(msg.str());
  }

  const double centerX = 0.5 * (xmin + xmax);
  const double centerY = 0.5 * (ymin + ymax);
  double viewSize = std::max(xmax - xmin, ymax - ymin) * 1.1;
  if (!(viewSize > 0.0) || !std::isfinite(viewSize)) {
    viewSize = 1.0;
  }

  outFile << "0\nSECTION\n2\nHEADER\n"
          << "9\n$PDMODE\n70\n"
          << pdmode << "\n"
          << "9\n$PDSIZE\n40\n0.5\n"
          << "0\nENDSEC\n";

  outFile << "0\nSECTION\n2\nTABLES\n"
          << "0\nTABLE\n2\nLAYER\n"
          << "0\nLAYER\n2\nxyz_points\n70\n0\n62\n7\n6\nCONTINUOUS\n"
          << "0\nLAYER\n2\nxyz_labels\n70\n0\n62\n3\n6\nCONTINUOUS\n"
          << "0\nLAYER\n2\nboundary\n70\n0\n62\n1\n6\nCONTINUOUS\n";
  if (hasGrid) {
    outFile << "0\nLAYER\n2\ngrid_points\n70\n0\n62\n5\n6\nCONTINUOUS\n"
            << "0\nLAYER\n2\ngrid_labels\n70\n0\n62\n4\n6\nCONTINUOUS\n";
  }
  outFile << "0\nENDTAB\n"
          << "0\nTABLE\n2\nVPORT\n"
          << "0\nVPORT\n2\n*ACTIVE\n10\n0.0\n20\n0.0\n11\n1.0\n21\n1.0\n12\n"
          << centerX << "\n22\n"
          << centerY << "\n40\n"
          << viewSize << "\n"
          << "0\nENDTAB\n0\nENDSEC\n";

  outFile << "0\nSECTION\n2\nENTITIES\n";

  auto appendPointAndLabel = [&](std::string &buffer, const Point3D &p,
                                 const std::string &layerPoints,
                                 const std::string &layerLabels) -> void {
    std::ostringstream ss;
    ss << std::fixed << std::setprecision(precision);
    ss << "0\nPOINT\n8\n"
       << layerPoints << "\n10\n"
       << p.x << "\n20\n"
       << p.y << "\n30\n"
       << p.z << "\n0\nTEXT\n8\n"
       << layerLabels << "\n10\n"
       << (p.x + 0.2) << "\n20\n"
       << (p.y + 0.2) << "\n30\n0.0\n40\n1.0\n1\n"
       << p.z << "\n";
    buffer += ss.str();
  };

  std::string buffer;
  buffer.reserve(256U * 1024U);

  if (boundary && boundary->valid) {
    for (const auto &loop : boundary->loops) {
      appendBoundaryLoopToDXFBuffer(buffer, loop, precision);
      if (buffer.size() >= (256U * 1024U)) {
        outFile << buffer;
        if (!outFile.good()) {
          if (errorMessage) {
            *errorMessage =
                "Error writing DXF boundary entities: " + outputFileName;
          }
          return false;
        }
        buffer.clear();
      }
    }
  }

  std::ifstream inFile(xyzInputFileName);
  if (!inFile.is_open()) {
    if (errorMessage) {
      *errorMessage =
          "Error reading filtered XYZ file for DXF: " + xyzInputFileName;
    }
    return false;
  }

  std::string line;
  while (std::getline(inFile, line)) {
    if (line.empty() || line[0] == '#') {
      continue;
    }
    std::replace(line.begin(), line.end(), ',', ' ');
    std::istringstream ss(line);
    Point3D p{};
    if (!(ss >> p.x >> p.y >> p.z)) {
      continue;
    }
    appendPointAndLabel(buffer, p, "xyz_points", "xyz_labels");
    if (buffer.size() >= (64U * 1024U)) {
      outFile << buffer;
      if (!outFile.good()) {
        if (errorMessage) {
          *errorMessage = "Error writing DXF file: " + outputFileName;
        }
        return false;
      }
      buffer.clear();
    }
  }

  size_t writtenGridCount = 0U;
  if (hasGrid) {
    for (size_t i = 0; i < grid.nx; ++i) {
      for (size_t j = 0; j < grid.ny; ++j) {
        const double gx = grid.xMin + static_cast<double>(i) * grid.gridSpacing;
        const double gy = grid.yMin + static_cast<double>(j) * grid.gridSpacing;
        if (boundary && boundary->valid &&
            !pointInConvexHull2D(*boundary, gx, gy)) {
          continue;
        }

        Point3D p{};
        std::string localError;
        if (!evaluator(i, j, p, &localError)) {
          if (errorMessage) {
            *errorMessage =
                localError.empty()
                    ? ("Error evaluating interpolation point for DXF: " +
                       outputFileName)
                    : localError;
          }
          return false;
        }
        appendPointAndLabel(buffer, p, "grid_points", "grid_labels");
        ++writtenGridCount;
        if (buffer.size() >= (4U * 1024U * 1024U)) {
          outFile << buffer;
          if (!outFile.good()) {
            if (errorMessage) {
              *errorMessage = "Error writing DXF file: " + outputFileName;
            }
            return false;
          }
          buffer.clear();
        }
      }
      maybeReportRowProgress(progressCallback, "Writing DXF:", i, grid.nx);
    }
  }

  outFile << buffer;
  outFile << "0\nENDSEC\n0\nEOF\n";
  const bool ok = outFile.good();
  outFile.close();
  if (ok && progressCallback) {
    std::ostringstream msg;
    msg << "DXF written grid nodes inside boundary: " << writtenGridCount;
    progressCallback(msg.str());
  }
  if (!ok && errorMessage) {
    *errorMessage = "Error writing DXF file: " + outputFileName;
  }
  return ok;
} /****************************************************************************
   * 11) Writers
   ****************************************************************************/
static auto writeXYZWithProgress(
    const std::string &outputFileName, const std::vector<Point3D> &points,
    int precision, const char *stageName, std::string *errorMessage,
    const std::function<void(const std::string &)> *statusUpdate = nullptr)
    -> bool {
  std::ofstream outFile(outputFileName, std::ios::out | std::ios::binary);
  if (!outFile.is_open()) {
    if (errorMessage) {
      *errorMessage = std::string("Error creating ") + stageName +
                      " file: " + outputFileName;
    }
    return false;
  }

  std::vector<char> fileBuffer(8U * 1024U * 1024U);
  outFile.rdbuf()->pubsetbuf(fileBuffer.data(),
                             static_cast<std::streamsize>(fileBuffer.size()));
  outFile << std::fixed << std::setprecision(precision);

  std::string chunk;
  chunk.reserve(4U * 1024U * 1024U);

  const size_t total = points.size();
  const size_t reportStep = std::max<size_t>(50000U, total / 100U + 1U);
  size_t nextReport = reportStep;

  for (size_t i = 0; i < total; ++i) {
    chunk += formatPointXYZLine(points[i], precision);
    if (chunk.size() >= 4U * 1024U * 1024U) {
      outFile.write(chunk.data(), static_cast<std::streamsize>(chunk.size()));
      chunk.clear();
      if (!outFile.good()) {
        if (errorMessage) {
          *errorMessage = std::string("Error writing ") + stageName +
                          " file: " + outputFileName;
        }
        return false;
      }
    }

    if (statusUpdate && i + 1U >= nextReport) {
      std::ostringstream msg;
      const double pct = total > 0U ? (100.0 * static_cast<double>(i + 1U) /
                                       static_cast<double>(total))
                                    : 100.0;
      msg << "Writing " << stageName << " ... " << (i + 1U) << "/" << total
          << " (" << std::fixed << std::setprecision(1) << pct << "%)";
      (*statusUpdate)(msg.str());
      nextReport += reportStep;
    }
  }

  if (!chunk.empty()) {
    outFile.write(chunk.data(), static_cast<std::streamsize>(chunk.size()));
  }

  const bool ok = outFile.good();
  outFile.close();
  if (!ok) {
    if (errorMessage) {
      *errorMessage = std::string("Error writing ") + stageName +
                      " file: " + outputFileName;
    }
    return false;
  }

  if (statusUpdate) {
    std::ostringstream msg;
    msg << "Writing " << stageName << " ... " << total << "/" << total
        << " (100.0%)";
    (*statusUpdate)(msg.str());
  }
  return true;
}

/*****************************************************************************
 * writeFilteredXYZ
 * ---------------------------------------------------------------------------
 * Purpose
 *   Writes the filtered/thinned point cloud to the filtered XYZ output file.
 *****************************************************************************/
static auto writeFilteredXYZ(
    const std::string &outputFileName,
    const std::vector<Point3D> &filteredPoints, int precision,
    std::string *errorMessage,
    const std::function<void(const std::string &)> *statusUpdate = nullptr)
    -> bool {
  return writeXYZWithProgress(outputFileName, filteredPoints, precision,
                              "filtered.xyz", errorMessage, statusUpdate);
}

/*****************************************************************************
 * buildBoundaryAndClipFilteredPoints
 * ---------------------------------------------------------------------------
 * Purpose
 *   Builds and iteratively refines the working outer boundary so the filtered
 * point set is fully contained.
 *****************************************************************************/
static auto buildBoundaryAndClipFilteredPoints(
    std::vector<Point3D> &filteredPoints, double minDistHint,
    ConvexHull2D &boundary,
    const std::function<void(const std::string &)> &progressCallback,
    std::string *errorMessage) -> bool {
  boundary = ConvexHull2D{};
  if (filteredPoints.empty()) {
    if (errorMessage) {
      *errorMessage =
          "Error: Cannot build boundary from an empty filtered point set.";
    }
    return false;
  }

  double boundaryCellSize =
      recommendedBoundaryCellSize(filteredPoints, minDistHint);
  constexpr size_t maxIterations = 6U;
  for (size_t iter = 0U; iter < maxIterations; ++iter) {
    ConvexHull2D candidate =
        buildConvexHull2D(filteredPoints, boundaryCellSize);
    if (!candidate.valid) {
      if (errorMessage) {
        *errorMessage =
            "Error: Failed to build a valid outer boundary from the filtered "
            "XYZ points.";
      }
      return false;
    }

    size_t outsideCount = 0U;
    for (const auto &p : filteredPoints) {
      if (!pointInConvexHull2D(candidate, p.x, p.y)) {
        ++outsideCount;
      }
    }

    if (progressCallback) {
      std::ostringstream msg;
      msg << "Boundary refinement iteration " << (iter + 1U)
          << ": filtered XYZ points outside boundary = " << outsideCount
          << ", loops=" << candidate.loops.size()
          << ", cellSize=" << boundaryCellSize;
      progressCallback(msg.str());
    }

    boundary = std::move(candidate);
    if (outsideCount == 0U) {
      return true;
    }

    boundaryCellSize *= 1.35;
  }

  if (errorMessage) {
    *errorMessage =
        "Error: Unable to build an outer boundary that contains all filtered "
        "XYZ points.";
  }
  return false;
}

/*****************************************************************************
 * 12) processXYZtoDXF
 ****************************************************************************/
static auto
processXYZtoDXF(const std::string &inputFileName, double minDist, int precision,
                int pdmode, double gridSpacing, bool useOutlierRemoval,
                bool useBoundary, InterpolationMethod methodUsed,
                std::function<void(const std::string &)> statusUpdate,
                std::string &finalMessage) -> bool {
  std::vector<std::string> reportLines;

  std::function<void(const std::string &)> reportStatus =
      [&](const std::string &msg) -> void {
    reportLines.push_back(msg);
    statusUpdate(msg);
  };

  auto writeReportFile = [&]() -> void {
    const std::string reportFileName = inputFileName + ".rpt.txt";
    std::ofstream rptFile(reportFileName);
    if (rptFile.is_open()) {
      for (const auto &line : reportLines) {
        rptFile << line << "\n";
      }
    }
  };

  auto finish = [&](bool success, const std::string &msg) -> bool {
    finalMessage = msg;
    if (!msg.empty()) {
      reportStatus(msg);
    }
    writeReportFile();
    return success;
  };

  try {
    const auto startTime = std::chrono::high_resolution_clock::now();

    if (minDist < 0.0) {
      return finish(false, "Error: minDist must be >= 0.");
    }
    if (precision < 0 || precision > 15) {
      return finish(false, "Error: precision must be between 0 and 15.");
    }
    if (!(gridSpacing > 0.0) || !std::isfinite(gridSpacing)) {
      return finish(false,
                    "Error: gridSpacing must be a positive finite value.");
    }

    reportStatus("Phase 1/7: reading original XYZ input ...");
    size_t totalPointsRead = 0U;
    std::vector<Point3D> filteredPoints;
    std::string readError;
    if (!readXYZFile(
            inputFileName, filteredPoints, totalPointsRead,
            [&](const std::string &msg) -> void { reportStatus(msg); },
            &readError)) {
      return finish(false, readError.empty()
                               ? "Error: No valid points found in file."
                               : readError);
    }
    {
      std::ostringstream msg;
      msg << "Total valid points read: " << totalPointsRead;
      reportStatus(msg.str());
    }

    if (useOutlierRemoval) {
      reportStatus(
          "Phase 2/7: robust local residual outlier filtering on original XYZ "
          "points ...");
      const double rawSpacing = averagePlanarSpacing(filteredPoints);
      const double neighborDist =
          std::max({5.0 * minDist, 3.0 * rawSpacing, 0.01});
      const double zThresholdFactor = 3.5;
      {
        std::ostringstream msg;
        msg << "Robust outlier parameters: neighborDist=" << neighborDist
            << ", zThresholdFactor=" << zThresholdFactor;
        reportStatus(msg.str());
      }
      std::string outlierError;
      if (!removeZOutliersInPlace(
              filteredPoints, neighborDist, zThresholdFactor,
              [&](const std::string &msg) -> void { reportStatus(msg); },
              &outlierError)) {
        return finish(false, outlierError.empty()
                                 ? "Error: All original XYZ points were "
                                   "removed during robust outlier filtering."
                                 : outlierError);
      }
      {
        std::ostringstream msg;
        msg << "Points after robust outlier removal on original XYZ data: "
            << filteredPoints.size();
        reportStatus(msg.str());
      }
    } else {
      reportStatus(
          "Phase 2/7: robust outlier removal disabled by user. Original XYZ "
          "points kept as read.");
    }

    reportStatus(
        "Phase 3/7: deterministic minDist thinning after optional outlier "
        "rejection ...");
    filteredPoints = filterPointsGrid(filteredPoints, minDist);
    {
      std::ostringstream msg;
      msg << "Points after minDist filter: " << filteredPoints.size();
      reportStatus(msg.str());
    }
    if (filteredPoints.empty()) {
      return finish(false,
                    "Error: All points were removed by the minDist filter.");
    }

    GridDefinition grid;
    ConvexHull2D outputBoundary;
    const ConvexHull2D *outputBoundaryPtr = nullptr;
    if (useBoundary) {
      reportStatus(
          "Phase 4/7: building a single tight outer boundary that contains all "
          "filtered XYZ points ...");
      std::string boundaryError;
      if (!buildBoundaryAndClipFilteredPoints(filteredPoints, minDist,
                                              outputBoundary, reportStatus,
                                              &boundaryError)) {
        return finish(false, boundaryError);
      }
      {
        std::ostringstream msg;
        const size_t occupiedCells = static_cast<size_t>(std::count(
            outputBoundary.occupied.begin(), outputBoundary.occupied.end(),
            static_cast<unsigned char>(1)));
        msg << "Boundary contains filtered XYZ points: "
            << filteredPoints.size()
            << ", loops=" << outputBoundary.loops.size()
            << ", occupiedCells=" << occupiedCells;
        reportStatus(msg.str());
      }
      outputBoundaryPtr = &outputBoundary;
    } else {
      reportStatus(
          "Phase 4/7: boundary definition/usage disabled by user. All filtered "
          "XYZ points and all interpolated nodes are kept.");
    }

    reportStatus(
        "Phase 5/7: writing filtered.xyz to disk before interpolation ...");
    {
      std::string ioError;
      if (!writeFilteredXYZ(inputFileName + ".filtered.xyz", filteredPoints,
                            precision, &ioError, &reportStatus)) {
        return finish(false, ioError);
      }
    }
    reportStatus("filtered.xyz successfully written.");

    std::function<bool(size_t, size_t, Point3D &, std::string *)> evaluator;
    std::function<bool(size_t, size_t, Point3D &, InterpolationDiagnostics &,
                       std::string *)>
        diagEvaluator;
    std::string methodLabel;
    ValidationMetrics methodValidation;

    if (methodUsed == METHOD_TPS) {
      methodLabel = "TPS";
      GridDefinition fullExtentGrid;
      std::string gridError;
      if (!buildGridDefinition(filteredPoints, gridSpacing, fullExtentGrid,
                               &gridError)) {
        return finish(false, gridError);
      }

      reportStatus("TPS control set: using all filtered points.");
      {
        std::ostringstream msg;
        msg << "TPS control-point count: " << filteredPoints.size();
        reportStatus(msg.str());
      }
      if (filteredPoints.empty()) {
        return finish(false,
                      "Error: TPS received an empty filtered point cloud.");
      }

      TPSStreamingContext tpsCtx;
      std::string ctxError;
      const double tpsBoundaryCellSize =
          recommendedBoundaryCellSize(filteredPoints, minDist);
      if (!buildTPSStreamingContext(std::move(filteredPoints), fullExtentGrid,
                                    tpsBoundaryCellSize, tpsCtx, &ctxError)) {
        return finish(false, ctxError);
      }
      if (useBoundary) {
        tpsCtx.model.hull = outputBoundary;
      }
      grid = tpsCtx.grid;
      methodValidation = tpsCtx.model.validation;
      evaluator = [tpsCtx](size_t i, size_t j, Point3D &p,
                           std::string *err) -> bool {
        return evaluateTPSStreamPoint(tpsCtx, i, j, p, err);
      };
      diagEvaluator = [tpsCtx](size_t i, size_t j, Point3D &p,
                               InterpolationDiagnostics &d,
                               std::string *err) -> bool {
        return evaluateTPSStreamPointDetailed(tpsCtx, i, j, p, &d, err);
      };
      reportStatus(
          "Phase 6/7: TPS context ready with cross-validated local lambda "
          "selection.");
    } else if (methodUsed == METHOD_CLAMPED_MLS) {
      methodLabel = "Clamped local MLS (bounded)";
      ClampedMLSStreamingContext monoCtx;
      std::string ctxError;
      if (!buildClampedMLSStreamingContext(
              filteredPoints, gridSpacing,
              recommendedBoundaryCellSize(filteredPoints, minDist), monoCtx,
              &ctxError)) {
        return finish(false, ctxError);
      }
      grid = monoCtx.grid;
      methodValidation = monoCtx.validation;
      evaluator = [monoCtx](size_t i, size_t j, Point3D &p,
                            std::string *err) -> bool {
        return evaluateClampedMLSStreamPoint(monoCtx, i, j, p, err);
      };
      diagEvaluator = [monoCtx](size_t i, size_t j, Point3D &p,
                                InterpolationDiagnostics &d,
                                std::string *err) -> bool {
        return evaluateClampedMLSStreamPointDetailed(monoCtx, i, j, p, &d, err);
      };
      reportStatus(
          "Phase 6/7: clamped local MLS context ready with predictive holdout "
          "validation.");
    } else {
      methodLabel = "Bicubic Hermite + robust MLS";
      BicubicStreamingContext bicCtx;
      std::string ctxError;
      if (!buildBicubicStreamingContext(
              filteredPoints, gridSpacing,
              recommendedBoundaryCellSize(filteredPoints, minDist), bicCtx,
              &ctxError)) {
        return finish(false, ctxError);
      }
      grid = bicCtx.grid;
      methodValidation = bicCtx.validation;
      {
        std::ostringstream msg;
        msg << "Selected bicubic tuning: neighbors="
            << bicCtx.tuning.neighborhoodTarget
            << ", supportMultiplier=" << bicCtx.tuning.supportMultiplier
            << ", ridgeFactor=" << bicCtx.tuning.ridgeFactor
            << ", basisOrder=" << bicCtx.tuning.basisOrder;
        reportStatus(msg.str());
      }
      evaluator = [bicCtx](size_t i, size_t j, Point3D &p,
                           std::string *err) -> bool {
        return evaluateBicubicStreamPoint(bicCtx, i, j, p, err);
      };
      diagEvaluator = [bicCtx](size_t i, size_t j, Point3D &p,
                               InterpolationDiagnostics &d,
                               std::string *err) -> bool {
        return evaluateBicubicStreamPointDetailed(bicCtx, i, j, p, &d, err);
      };
      reportStatus(
          "Phase 6/7: bicubic context ready with predictive robust anisotropic "
          "MLS and 2D derivative limiting.");
    }

    if (methodValidation.valid) {
      std::ostringstream msg;
      msg << "Validation summary for " << methodLabel
          << ": RMSE=" << methodValidation.rmse
          << ", MAE=" << methodValidation.mae
          << ", P95=" << methodValidation.p95
          << ", n=" << methodValidation.count;
      reportStatus(msg.str());
    }

    {
      std::ostringstream msg;
      msg << "Grid nodes planned: " << grid.totalPoints << " (" << grid.nx
          << " x " << grid.ny << ")"
          << ", spacing=" << grid.gridSpacing;
      reportStatus(msg.str());
    }
    if (useBoundary) {
      const size_t occupiedCells = static_cast<size_t>(std::count(
          outputBoundary.occupied.begin(), outputBoundary.occupied.end(),
          static_cast<unsigned char>(1)));
      std::ostringstream msg;
      msg << "Boundary ready from source XYZ cloud using grid/raster occupancy "
             "+ contour extraction: cellSize="
          << outputBoundary.cellSize << ", occupiedCells=" << occupiedCells
          << ", loops=" << outputBoundary.loops.size();
      reportStatus(msg.str());
    } else {
      reportStatus(
          "Boundary usage disabled: streamed outputs will not be clipped by an "
          "outer polygon.");
    }

    reportStatus("Phase 6/7: starting streamed grid.xyz writing ...");
    {
      std::string ioError;
      if (!writeGridXYZStream(
              inputFileName + ".grid.xyz", grid, outputBoundaryPtr, precision,
              evaluator,
              [&](const std::string &msg) -> void { reportStatus(msg); },
              &ioError)) {
        return finish(false, ioError);
      }
    }

    reportStatus("Phase 5b/6: writing confidence diagnostics grid ...");
    ConfidenceSummary confidenceSummary;
    {
      std::string ioError;
      if (!writeConfidenceXYZStream(
              inputFileName + ".confidence.xyz", grid, outputBoundaryPtr,
              precision, diagEvaluator,
              [&](const std::string &msg) -> void { reportStatus(msg); },
              confidenceSummary, &ioError)) {
        return finish(false, ioError);
      }
    }
    {
      std::ostringstream msg;
      msg << "Confidence summary: meanConfidence="
          << confidenceSummary.meanConfidence
          << ", meanResidual=" << confidenceSummary.meanResidual
          << ", fallbacks=" << confidenceSummary.fallbackCount
          << ", outsideHull=" << confidenceSummary.outsideHullCount
          << ", clampedLocal=" << confidenceSummary.clampedLocalCount;
      reportStatus(msg.str());
    }

    reportStatus("Phase 7/7: starting streamed DXF generation ...");
    {
      std::string ioError;
      if (!writeDXFStreamFromXYZFile(
              inputFileName + ".dxf", inputFileName + ".filtered.xyz", grid,
              outputBoundaryPtr, precision, pdmode,
              grid.valid && grid.totalPoints > 0U, evaluator,
              [&](const std::string &msg) -> void { reportStatus(msg); },
              &ioError)) {
        return finish(false, ioError);
      }
    }

    const auto endTime = std::chrono::high_resolution_clock::now();
    const double elapsedSec =
        std::chrono::duration<double>(endTime - startTime).count();
    const int totalSeconds = static_cast<int>(std::round(elapsedSec));

    std::ostringstream finishMsg;
    finishMsg << "Done. Total time: " << totalSeconds << " sec.";
    return finish(true, finishMsg.str());
  } catch (const std::bad_alloc &) {
    return finish(false, "Error: insufficient memory during interpolation. The "
                         "process now aborts safely instead of crashing.");
  } catch (const std::exception &e) {
    return finish(false, std::string("Error: ") + e.what());
  } catch (...) {
    return finish(false, "Error: unexpected failure during interpolation.");
  }
}

struct GuiStatusMessage {
  std::string message;
};

struct GuiCompletionMessage {
  bool success = false;
  std::string message;
};

struct GuiWorkerState {
  std::thread worker;
  std::atomic_bool running{false};

  ~GuiWorkerState() {
    if (worker.joinable()) {
      worker.join();
    }
  }
};

/****************************************************************************
 * Professional GUI layout helpers
 ****************************************************************************/
struct GuiLayoutHandles {
  HWND grpSource = nullptr;
  HWND lblInputFile = nullptr;
  HWND hInputFile = nullptr;
  HWND btnBrowse = nullptr;
  HWND txtSourceHint = nullptr;

  HWND grpParameters = nullptr;
  HWND lblMinDist = nullptr;
  HWND hMinDist = nullptr;
  HWND txtMinDist = nullptr;
  HWND lblPrecision = nullptr;
  HWND hPrecision = nullptr;
  HWND txtPrecision = nullptr;
  HWND lblPDMODE = nullptr;
  HWND hPDMODE = nullptr;
  HWND txtPDMODE = nullptr;
  HWND lblGridSpacing = nullptr;
  HWND hGridSpacing = nullptr;
  HWND txtGridSpacing = nullptr;
  HWND hCheckOutlier = nullptr;
  HWND txtCheckOutlier = nullptr;
  HWND hCheckBoundary = nullptr;
  HWND txtCheckBoundary = nullptr;

  HWND grpMethod = nullptr;
  HWND txtMethodIntro = nullptr;
  HWND hRadioBicubic = nullptr;
  HWND txtBicubic = nullptr;
  HWND hRadioTPS = nullptr;
  HWND txtTPS = nullptr;
  HWND hRadioClamped = nullptr;
  HWND txtClamped = nullptr;

  HWND grpExecution = nullptr;
  HWND btnRun = nullptr;
  HWND txtExecutionHint = nullptr;
  HWND hwndStatus = nullptr;
};

struct GuiWindowState {
  GuiLayoutHandles gui;
  HFONT hGuiFont = nullptr;
  HFONT hGuiFontBold = nullptr;
  std::shared_ptr<std::atomic_bool> windowAlive =
      std::make_shared<std::atomic_bool>(true);
  GuiWorkerState workerState;

  ~GuiWindowState() {
    if (hGuiFontBold != nullptr && hGuiFontBold != hGuiFont) {
      DeleteObject(hGuiFontBold);
      hGuiFontBold = nullptr;
    }
    if (hGuiFont != nullptr) {
      DeleteObject(hGuiFont);
      hGuiFont = nullptr;
    }
  }
};

static auto GetGuiWindowState(HWND hwnd) noexcept -> GuiWindowState * {
  return reinterpret_cast<GuiWindowState *>(
      GetWindowLongPtrA(hwnd, GWLP_USERDATA));
}

/*****************************************************************************
 * ApplyFontToWindowTree
 * ---------------------------------------------------------------------------
 * Purpose
 *   Recursively applies the chosen GUI font to a window hierarchy.
 *****************************************************************************/
static void ApplyFontToWindowTree(const GuiLayoutHandles &g, HFONT font,
                                  HFONT fontBold) {
  const HWND handles[] = {
      g.grpSource,        g.lblInputFile,    g.hInputFile,
      g.btnBrowse,        g.txtSourceHint,   g.grpParameters,
      g.lblMinDist,       g.hMinDist,        g.txtMinDist,
      g.lblPrecision,     g.hPrecision,      g.txtPrecision,
      g.lblPDMODE,        g.hPDMODE,         g.txtPDMODE,
      g.lblGridSpacing,   g.hGridSpacing,    g.txtGridSpacing,
      g.hCheckOutlier,    g.txtCheckOutlier, g.hCheckBoundary,
      g.txtCheckBoundary, g.grpMethod,       g.txtMethodIntro,
      g.hRadioBicubic,    g.txtBicubic,      g.hRadioTPS,
      g.txtTPS,           g.hRadioClamped,   g.txtClamped,
      g.grpExecution,     g.btnRun,          g.txtExecutionHint,
      g.hwndStatus};
  for (HWND h : handles) {
    if (h != nullptr) {
      SendMessageA(h, WM_SETFONT, reinterpret_cast<WPARAM>(font), TRUE);
    }
  }
  const HWND groupBoxes[] = {g.grpSource, g.grpParameters, g.grpMethod,
                             g.grpExecution};
  for (HWND h : groupBoxes) {
    if (h != nullptr) {
      SendMessageA(
          h, WM_SETFONT,
          reinterpret_cast<WPARAM>(fontBold != nullptr ? fontBold : font),
          TRUE);
    }
  }
  if (g.btnRun != nullptr && fontBold != nullptr) {
    SendMessageA(g.btnRun, WM_SETFONT, reinterpret_cast<WPARAM>(fontBold),
                 TRUE);
  }
}

/*****************************************************************************
 * ApplyProfessionalGuiLayout
 * ---------------------------------------------------------------------------
 * Purpose
 *   Creates and applies the final GUI sizing, font, and control layout.
 *****************************************************************************/
static void ApplyProfessionalGuiLayout(HWND hwnd, const GuiLayoutHandles &g) {
  RECT rc{};
  GetClientRect(hwnd, &rc);
  const int clientW = rc.right - rc.left;

  const int margin = 16;
  const int sectionGap = 12;
  const int labelW = 118;
  const int editW = 63;
  const int browseW = 102;
  const int topGroupH = 96;
  const int leftGroupH = 360;
  const int rightGroupH = 392;
  const int bottomGroupH = 126;
  const int midTop = margin + topGroupH + sectionGap;
  const int leftW = std::max(455, (clientW - 3 * margin) / 2);
  const int rightW = clientW - 3 * margin - leftW;

  const int xLeft = margin;
  const int xRight = margin + leftW + margin;
  const int yBottom = midTop + std::max(leftGroupH, rightGroupH) + sectionGap;

  MoveWindow(g.grpSource, margin, margin, clientW - 2 * margin, topGroupH,
             TRUE);
  MoveWindow(g.lblInputFile, margin + 16, margin + 24, labelW, 22, TRUE);
  MoveWindow(g.hInputFile, margin + 16 + labelW, margin + 20,
             clientW - 2 * margin - 16 - labelW - browseW - 36, 30, TRUE);
  MoveWindow(g.btnBrowse, clientW - margin - browseW - 16, margin + 20, browseW,
             30, TRUE);
  MoveWindow(g.txtSourceHint, margin + 16 + labelW, margin + 56,
             clientW - 2 * margin - labelW - 48, 28, TRUE);

  MoveWindow(g.grpParameters, xLeft, midTop, leftW, leftGroupH, TRUE);
  const int paramX = xLeft + 18;
  const int editX = paramX + labelW;
  const int hintX = editX + editW + 14;
  const int hintW = std::max(110, leftW - (hintX - xLeft) - 20);
  int y = midTop + 28;
  const int rowH = 30;
  const int hintH = 38;

  MoveWindow(g.lblMinDist, paramX, y + 4, labelW, 20, TRUE);
  MoveWindow(g.hMinDist, editX, y, editW, rowH, TRUE);
  MoveWindow(g.txtMinDist, hintX, y + 1, hintW, hintH, TRUE);
  y += 48;

  MoveWindow(g.lblPrecision, paramX, y + 4, labelW, 20, TRUE);
  MoveWindow(g.hPrecision, editX, y, editW, rowH, TRUE);
  MoveWindow(g.txtPrecision, hintX, y + 1, hintW, hintH, TRUE);
  y += 48;

  MoveWindow(g.lblPDMODE, paramX, y + 4, labelW, 20, TRUE);
  MoveWindow(g.hPDMODE, editX, y, editW, rowH, TRUE);
  MoveWindow(g.txtPDMODE, hintX, y + 1, hintW, hintH, TRUE);
  y += 48;

  MoveWindow(g.lblGridSpacing, paramX, y + 4, labelW, 20, TRUE);
  MoveWindow(g.hGridSpacing, editX, y, editW, rowH, TRUE);
  MoveWindow(g.txtGridSpacing, hintX, y + 1, hintW, hintH, TRUE);
  y += 50;

  MoveWindow(g.hCheckOutlier, paramX, y, leftW - 36, 24, TRUE);
  MoveWindow(g.txtCheckOutlier, paramX + 20, y + 24, leftW - 56, 38, TRUE);
  y += 72;

  MoveWindow(g.hCheckBoundary, paramX, y, leftW - 36, 24, TRUE);
  MoveWindow(g.txtCheckBoundary, paramX + 20, y + 24, leftW - 56, 42, TRUE);

  MoveWindow(g.grpMethod, xRight, midTop, rightW, rightGroupH, TRUE);
  const int methodX = xRight + 18;
  const int methodW = rightW - 36;
  int my = midTop + 24;
  MoveWindow(g.txtMethodIntro, methodX, my, methodW, 46, TRUE);
  my += 54;
  MoveWindow(g.hRadioBicubic, methodX, my, methodW, 22, TRUE);
  my += 24;
  MoveWindow(g.txtBicubic, methodX + 20, my, methodW - 20, 68, TRUE);
  my += 78;
  MoveWindow(g.hRadioTPS, methodX, my, methodW, 22, TRUE);
  my += 24;
  MoveWindow(g.txtTPS, methodX + 20, my, methodW - 20, 68, TRUE);
  my += 78;
  MoveWindow(g.hRadioClamped, methodX, my, methodW, 22, TRUE);
  my += 24;
  MoveWindow(g.txtClamped, methodX + 20, my, methodW - 20, 68, TRUE);

  MoveWindow(g.grpExecution, margin, yBottom, clientW - 2 * margin,
             bottomGroupH, TRUE);
  MoveWindow(g.txtExecutionHint, margin + 18, yBottom + 24,
             clientW - 2 * margin - 232, 24, TRUE);
  MoveWindow(g.btnRun, clientW - margin - 166, yBottom + 18, 148, 36, TRUE);
  MoveWindow(g.hwndStatus, margin + 18, yBottom + 58, clientW - 2 * margin - 36,
             52, TRUE);
}

/****************************************************************************
 * Window Procedure (GUI code)
 ****************************************************************************/

LRESULT CALLBACK WindowProc(HWND hwnd, UINT uMsg, WPARAM wParam,
                            LPARAM lParam) {
  if (uMsg == WM_NCCREATE) {
    auto *state = new (std::nothrow) GuiWindowState();
    if (state == nullptr) {
      return FALSE;
    }

    SetLastError(0);
    static_cast<void>(SetWindowLongPtrA(hwnd, GWLP_USERDATA,
                                        reinterpret_cast<LONG_PTR>(state)));
    if (GetLastError() != 0) {
      delete state;
      return FALSE;
    }
    return TRUE;
  }

  GuiWindowState *const state = GetGuiWindowState(hwnd);
  if (uMsg == WM_NCDESTROY) {
    if (state != nullptr) {
      SetWindowLongPtrA(hwnd, GWLP_USERDATA, 0);
      delete state;
    }
    return DefWindowProcA(hwnd, uMsg, wParam, lParam);
  }

  if (state == nullptr) {
    return DefWindowProcA(hwnd, uMsg, wParam, lParam);
  }

  GuiLayoutHandles &gui = state->gui;
  HFONT &hGuiFont = state->hGuiFont;
  HFONT &hGuiFontBold = state->hGuiFontBold;
  GuiWorkerState &workerState = state->workerState;

  switch (uMsg) {
  case WM_CREATE: {
    hGuiFont =
        CreateFontA(-15, 0, 0, 0, FW_NORMAL, FALSE, FALSE, FALSE,
                    DEFAULT_CHARSET, OUT_DEFAULT_PRECIS, CLIP_DEFAULT_PRECIS,
                    CLEARTYPE_QUALITY, DEFAULT_PITCH | FF_DONTCARE, "Segoe UI");
    hGuiFontBold =
        CreateFontA(-15, 0, 0, 0, FW_SEMIBOLD, FALSE, FALSE, FALSE,
                    DEFAULT_CHARSET, OUT_DEFAULT_PRECIS, CLIP_DEFAULT_PRECIS,
                    CLEARTYPE_QUALITY, DEFAULT_PITCH | FF_DONTCARE, "Segoe UI");
    if (hGuiFontBold == nullptr) {
      hGuiFontBold = hGuiFont;
    }

    gui.grpSource = CreateWindowA("BUTTON", " Source File ",
                                  WS_VISIBLE | WS_CHILD | BS_GROUPBOX, 0, 0, 10,
                                  10, hwnd, nullptr, nullptr, nullptr);
    gui.lblInputFile =
        CreateWindowA("STATIC", "Input XYZ file:", WS_VISIBLE | WS_CHILD, 0, 0,
                      10, 10, hwnd, nullptr, nullptr, nullptr);
    gui.hInputFile = CreateWindowExA(
        WS_EX_CLIENTEDGE, "EDIT", "",
        WS_CHILD | WS_VISIBLE | WS_BORDER | ES_AUTOHSCROLL | WS_TABSTOP, 0, 0,
        10, 10, hwnd, makeControlMenuHandle(IDC_INPUT_FILE), nullptr, nullptr);
    gui.btnBrowse = CreateWindowA(
        "BUTTON", "Browse...",
        WS_VISIBLE | WS_CHILD | WS_TABSTOP | BS_PUSHBUTTON, 0, 0, 10, 10, hwnd,
        makeControlMenuHandle(IDC_BROWSE_BUTTON), nullptr, nullptr);
    gui.txtSourceHint = CreateWindowA(
        "STATIC",
        "Select the source XYZ point cloud to be interpolated to DXF.",
        WS_VISIBLE | WS_CHILD, 0, 0, 10, 10, hwnd, nullptr, nullptr, nullptr);

    gui.grpParameters = CreateWindowA("BUTTON", " Processing Parameters ",
                                      WS_VISIBLE | WS_CHILD | BS_GROUPBOX, 0, 0,
                                      10, 10, hwnd, nullptr, nullptr, nullptr);
    gui.lblMinDist =
        CreateWindowA("STATIC", "Min Dist:", WS_VISIBLE | WS_CHILD, 0, 0, 10,
                      10, hwnd, nullptr, nullptr, nullptr);
    gui.hMinDist = CreateWindowExA(
        WS_EX_CLIENTEDGE, "EDIT", "5",
        WS_CHILD | WS_VISIBLE | WS_BORDER | ES_AUTOHSCROLL | WS_TABSTOP, 0, 0,
        10, 10, hwnd, makeControlMenuHandle(IDC_MIN_DIST), nullptr, nullptr);
    gui.txtMinDist = CreateWindowA(
        "STATIC", "Minimum XY spacing for point filtering.",
        WS_VISIBLE | WS_CHILD, 0, 0, 10, 10, hwnd, nullptr, nullptr, nullptr);

    gui.lblPrecision =
        CreateWindowA("STATIC", "Precision:", WS_VISIBLE | WS_CHILD, 0, 0, 10,
                      10, hwnd, nullptr, nullptr, nullptr);
    gui.hPrecision = CreateWindowExA(
        WS_EX_CLIENTEDGE, "EDIT", "2",
        WS_CHILD | WS_VISIBLE | WS_BORDER | ES_AUTOHSCROLL | WS_TABSTOP, 0, 0,
        10, 10, hwnd, makeControlMenuHandle(IDC_PRECISION), nullptr, nullptr);
    gui.txtPrecision = CreateWindowA(
        "STATIC", "Decimal places in XYZ and DXF outputs.",
        WS_VISIBLE | WS_CHILD, 0, 0, 10, 10, hwnd, nullptr, nullptr, nullptr);

    gui.lblPDMODE = CreateWindowA("STATIC", "PDMODE:", WS_VISIBLE | WS_CHILD, 0,
                                  0, 10, 10, hwnd, nullptr, nullptr, nullptr);
    gui.hPDMODE = CreateWindowExA(
        WS_EX_CLIENTEDGE, "EDIT", "3",
        WS_CHILD | WS_VISIBLE | WS_BORDER | ES_AUTOHSCROLL | WS_TABSTOP, 0, 0,
        10, 10, hwnd, makeControlMenuHandle(IDC_PDMODE), nullptr, nullptr);
    gui.txtPDMODE = CreateWindowA(
        "STATIC", "DXF point-style written into the DXF header.",
        WS_VISIBLE | WS_CHILD, 0, 0, 10, 10, hwnd, nullptr, nullptr, nullptr);

    gui.lblGridSpacing =
        CreateWindowA("STATIC", "Grid Spacing:", WS_VISIBLE | WS_CHILD, 0, 0,
                      10, 10, hwnd, nullptr, nullptr, nullptr);
    gui.hGridSpacing = CreateWindowExA(
        WS_EX_CLIENTEDGE, "EDIT", "20",
        WS_CHILD | WS_VISIBLE | WS_BORDER | ES_AUTOHSCROLL | WS_TABSTOP, 0, 0,
        10, 10, hwnd, makeControlMenuHandle(IDC_GRID_SPACING), nullptr,
        nullptr);
    gui.txtGridSpacing = CreateWindowA(
        "STATIC", "Regular interpolation spacing of surface grid.",
        WS_VISIBLE | WS_CHILD, 0, 0, 10, 10, hwnd, nullptr, nullptr, nullptr);

    gui.hCheckOutlier = CreateWindowA(
        "BUTTON", "Use robust outlier removal",
        WS_VISIBLE | WS_CHILD | BS_AUTOCHECKBOX | WS_TABSTOP, 0, 0, 10, 10,
        hwnd, makeControlMenuHandle(IDC_CHECK_OUTLIERS), nullptr, nullptr);
    SendMessageA(gui.hCheckOutlier, BM_SETCHECK, BST_CHECKED, 0);
    gui.txtCheckOutlier = CreateWindowA(
        "STATIC",
        "When enabled, robust local residual screening is applied to the "
        "original XYZ cloud before min-distance thinning.",
        WS_VISIBLE | WS_CHILD, 0, 0, 10, 10, hwnd, nullptr, nullptr, nullptr);

    gui.hCheckBoundary = CreateWindowA(
        "BUTTON", "Use outer boundary",
        WS_VISIBLE | WS_CHILD | BS_AUTOCHECKBOX | WS_TABSTOP, 0, 0, 10, 10,
        hwnd, makeControlMenuHandle(IDC_CHECK_BOUNDARY), nullptr, nullptr);
    SendMessageA(gui.hCheckBoundary, BM_SETCHECK, BST_CHECKED, 0);
    gui.txtCheckBoundary = CreateWindowA(
        "STATIC",
        "When enabled, a tight outer boundary from the filtered XYZ cloud "
        "clips outputs and is also written to DXF.",
        WS_VISIBLE | WS_CHILD, 0, 0, 10, 10, hwnd, nullptr, nullptr, nullptr);

    gui.grpMethod = CreateWindowA("BUTTON", " Interpolation Method ",
                                  WS_VISIBLE | WS_CHILD | BS_GROUPBOX, 0, 0, 10,
                                  10, hwnd, nullptr, nullptr, nullptr);
    gui.txtMethodIntro = CreateWindowA(
        "STATIC",
        "Choose the interpolation strategy that best matches the spatial "
        "character of the input cloud.",
        WS_VISIBLE | WS_CHILD, 0, 0, 10, 10, hwnd, nullptr, nullptr, nullptr);
    gui.hRadioBicubic = CreateWindowA(
        "BUTTON", "Bicubic Hermite with robust MLS (default)",
        WS_VISIBLE | WS_CHILD | BS_AUTORADIOBUTTON | WS_TABSTOP, 0, 0, 10, 10,
        hwnd, makeControlMenuHandle(IDC_RADIO_BICUBIC), nullptr, nullptr);
    SendMessageA(gui.hRadioBicubic, BM_SETCHECK, BST_CHECKED, 0);
    gui.txtBicubic = CreateWindowA(
        "STATIC",
        "Smooth general-purpose choice for regular grid production, and "
        "stable high-quality interpolation method. With sparse point clouds "
        "may overshoot.",
        WS_VISIBLE | WS_CHILD, 0, 0, 10, 10, hwnd, nullptr, nullptr, nullptr);
    gui.hRadioTPS = CreateWindowA(
        "BUTTON", "Thin Plate Spline (scattered data)",
        WS_VISIBLE | WS_CHILD | BS_AUTORADIOBUTTON | WS_TABSTOP, 0, 0, 10, 10,
        hwnd, makeControlMenuHandle(IDC_RADIO_TPS), nullptr, nullptr);
    gui.txtTPS = CreateWindowA(
        "STATIC",
        "Best for irregular scattered clouds when stronger local flexibility "
        "is needed, preserves smooth minimum-curvature and "
        "minimum-bending-energy style surface behaviour.",
        WS_VISIBLE | WS_CHILD, 0, 0, 10, 10, hwnd, nullptr, nullptr, nullptr);
    gui.hRadioClamped = CreateWindowA(
        "BUTTON", "Clamped local MLS (bounded)",
        WS_VISIBLE | WS_CHILD | BS_AUTORADIOBUTTON | WS_TABSTOP, 0, 0, 10, 10,
        hwnd, makeControlMenuHandle(IDC_RADIO_CLAMPED_MLS), nullptr, nullptr);
    gui.txtClamped = CreateWindowA(
        "STATIC",
        "Conservative bounded local method that reduces overshoot by "
        "constraining predictions to the nearby data envelope.",
        WS_VISIBLE | WS_CHILD, 0, 0, 10, 10, hwnd, nullptr, nullptr, nullptr);

    gui.grpExecution = CreateWindowA("BUTTON", " Execution ",
                                     WS_VISIBLE | WS_CHILD | BS_GROUPBOX, 0, 0,
                                     10, 10, hwnd, nullptr, nullptr, nullptr);
    gui.txtExecutionHint = CreateWindowA(
        "STATIC",
        "Review the parameter values, choose the interpolation method, and "
        "then run the conversion.",
        WS_VISIBLE | WS_CHILD, 0, 0, 10, 10, hwnd, nullptr, nullptr, nullptr);
    gui.btnRun = CreateWindowA(
        "BUTTON", "Run Conversion",
        WS_VISIBLE | WS_CHILD | BS_DEFPUSHBUTTON | WS_TABSTOP, 0, 0, 10, 10,
        hwnd, makeControlMenuHandle(IDC_RUN_BUTTON), nullptr, nullptr);
    gui.hwndStatus = CreateWindowExA(
        WS_EX_CLIENTEDGE, "STATIC", "Status: Idle",
        WS_VISIBLE | WS_CHILD | WS_BORDER, 0, 0, 10, 10, hwnd,
        makeControlMenuHandle(IDC_STATUS_STATIC), nullptr, nullptr);

    ApplyFontToWindowTree(gui, hGuiFont, hGuiFontBold);
    ApplyProfessionalGuiLayout(hwnd, gui);
    return 0;
  }

  case WM_ERASEBKGND: {
    RECT rc{};
    GetClientRect(hwnd, &rc);
    FillRect(reinterpret_cast<HDC>(wParam), &rc,
             GetSysColorBrush(COLOR_BTNFACE));
    return 1;
  }

  case WM_PAINT: {
    PAINTSTRUCT ps{};
    HDC hdc = BeginPaint(hwnd, &ps);
    FillRect(hdc, &ps.rcPaint, GetSysColorBrush(COLOR_BTNFACE));
    EndPaint(hwnd, &ps);
    return 0;
  }

  case WM_CTLCOLORSTATIC: {
    HDC hdc = reinterpret_cast<HDC>(wParam);
    SetBkMode(hdc, TRANSPARENT);
    SetTextColor(hdc, GetSysColor(COLOR_WINDOWTEXT));
    return reinterpret_cast<INT_PTR>(GetSysColorBrush(COLOR_BTNFACE));
  }

  case WM_CTLCOLORBTN: {
    HDC hdc = reinterpret_cast<HDC>(wParam);
    SetBkMode(hdc, TRANSPARENT);
    SetTextColor(hdc, GetSysColor(COLOR_WINDOWTEXT));
    return reinterpret_cast<INT_PTR>(GetSysColorBrush(COLOR_BTNFACE));
  }

  case WM_COMMAND: {
    if (LOWORD(wParam) == IDC_BROWSE_BUTTON) {
      const std::string filePath = openFileDialog(hwnd);
      if (!filePath.empty()) {
        SetWindowTextA(gui.hInputFile, filePath.c_str());
      }
    } else if (LOWORD(wParam) == IDC_RUN_BUTTON) {
      if (workerState.running.load()) {
        MessageBoxA(hwnd, "A conversion is already running.", "Busy",
                    MB_ICONINFORMATION | MB_OK);
        return 0;
      }
      EnableWindow(GetDlgItem(hwnd, IDC_RUN_BUTTON), FALSE);

      std::vector<char> inputFile(32768, '\0');
      char minDistStr[64], precisionStr[64], pdModeStr[64];
      char gridSpacingStr[64];
      GetWindowTextA(gui.hInputFile, inputFile.data(),
                     static_cast<int>(inputFile.size()));
      GetWindowTextA(gui.hMinDist, minDistStr, sizeof(minDistStr));
      GetWindowTextA(gui.hPrecision, precisionStr, sizeof(precisionStr));
      GetWindowTextA(gui.hPDMODE, pdModeStr, sizeof(pdModeStr));
      GetWindowTextA(gui.hGridSpacing, gridSpacingStr, sizeof(gridSpacingStr));

      if (std::strlen(inputFile.data()) == 0) {
        MessageBoxA(hwnd, "Please select an input file.", "Error",
                    MB_ICONERROR | MB_OK);
        EnableWindow(GetDlgItem(hwnd, IDC_RUN_BUTTON), TRUE);
        break;
      }

      double d_minDist = 0.0;
      int i_precision = 2;
      int i_pdmode = 3;
      double d_gridSpacing = 10.0;
      const bool useOutlierRemoval =
          (SendMessageA(gui.hCheckOutlier, BM_GETCHECK, 0, 0) == BST_CHECKED);
      const bool useBoundary =
          (SendMessageA(gui.hCheckBoundary, BM_GETCHECK, 0, 0) == BST_CHECKED);

      try {
        d_minDist = std::stod(minDistStr);
        i_precision = std::stoi(precisionStr);
        i_pdmode = std::stoi(pdModeStr);
        d_gridSpacing = std::stod(gridSpacingStr);
      } catch (...) {
        MessageBoxA(hwnd, "One or more numeric inputs are invalid.", "Error",
                    MB_ICONERROR | MB_OK);
        EnableWindow(GetDlgItem(hwnd, IDC_RUN_BUTTON), TRUE);
        break;
      }

      if (d_minDist < 0.0) {
        MessageBoxA(hwnd, "Min Dist must be >= 0.", "Error",
                    MB_ICONERROR | MB_OK);
        EnableWindow(GetDlgItem(hwnd, IDC_RUN_BUTTON), TRUE);
        break;
      }
      if (i_precision < 0 || i_precision > 15) {
        MessageBoxA(hwnd, "Precision must be between 0 and 15.", "Error",
                    MB_ICONERROR | MB_OK);
        EnableWindow(GetDlgItem(hwnd, IDC_RUN_BUTTON), TRUE);
        break;
      }
      if (!(d_gridSpacing > 0.0) || !std::isfinite(d_gridSpacing)) {
        MessageBoxA(hwnd, "Grid Spacing must be a positive finite value.",
                    "Error", MB_ICONERROR | MB_OK);
        EnableWindow(GetDlgItem(hwnd, IDC_RUN_BUTTON), TRUE);
        break;
      }

      InterpolationMethod methodUsed = METHOD_BICUBIC;
      if (SendMessageA(gui.hRadioTPS, BM_GETCHECK, 0, 0) == BST_CHECKED) {
        methodUsed = METHOD_TPS;
      } else if (SendMessageA(gui.hRadioClamped, BM_GETCHECK, 0, 0) ==
                 BST_CHECKED) {
        methodUsed = METHOD_CLAMPED_MLS;
      }

      updateStatus(gui.hwndStatus, "Starting conversion...");

      const std::string inputFileStr(inputFile.data());
      const auto aliveFlag = state->windowAlive;
      GuiWorkerState *const workerStatePtr = &workerState;
      workerState.running.store(true);
      workerState.worker =
          std::thread([hwnd, aliveFlag, workerStatePtr, inputFileStr, d_minDist,
                       i_precision, i_pdmode, d_gridSpacing, useOutlierRemoval,
                       useBoundary, methodUsed]() -> void {
            auto postStatus = [hwnd,
                               aliveFlag](const std::string &msg) -> void {
              if (!aliveFlag || !aliveFlag->load() || !IsWindow(hwnd)) {
                return;
              }
              auto sm = std::make_unique<GuiStatusMessage>();
              sm->message = msg;
              if (!PostMessageA(hwnd, WM_APP_STATUS_UPDATE, 0,
                                reinterpret_cast<LPARAM>(sm.get()))) {
                return;
              }
              static_cast<void>(sm.release());
            };

            bool success = false;
            std::string finalMessage;
            try {
              success = processXYZtoDXF(inputFileStr, d_minDist, i_precision,
                                        i_pdmode, d_gridSpacing,
                                        useOutlierRemoval, useBoundary,
                                        methodUsed, postStatus, finalMessage);
            } catch (const std::bad_alloc &) {
              success = false;
              finalMessage = "Error: insufficient memory during interpolation.";
            } catch (const std::exception &e) {
              success = false;
              finalMessage = std::string("Error: ") + e.what();
            } catch (...) {
              success = false;
              finalMessage = "Error: unexpected worker-thread failure.";
            }

            workerStatePtr->running.store(false);

            if (!aliveFlag || !aliveFlag->load() || !IsWindow(hwnd)) {
              return;
            }

            auto cm = std::make_unique<GuiCompletionMessage>();
            cm->success = success;
            cm->message = std::move(finalMessage);
            if (!PostMessageA(hwnd, WM_APP_PROCESSING_DONE, 0,
                              reinterpret_cast<LPARAM>(cm.get()))) {
              return;
            }
            static_cast<void>(cm.release());
          });
    }
    return 0;
  }

  case WM_APP_STATUS_UPDATE: {
    std::unique_ptr<GuiStatusMessage> sm(
        reinterpret_cast<GuiStatusMessage *>(lParam));
    if (sm && gui.hwndStatus && IsWindow(gui.hwndStatus)) {
      updateStatus(gui.hwndStatus, sm->message);
    }
    return 0;
  }

  case WM_APP_PROCESSING_DONE: {
    std::unique_ptr<GuiCompletionMessage> cm(
        reinterpret_cast<GuiCompletionMessage *>(lParam));
    EnableWindow(GetDlgItem(hwnd, IDC_RUN_BUTTON), TRUE);
    if (cm && gui.hwndStatus && IsWindow(gui.hwndStatus)) {
      const std::string fallback =
          cm->success ? "Conversion completed." : "Conversion failed.";
      updateStatus(gui.hwndStatus,
                   cm->message.empty() ? fallback : cm->message);
    }
    if (workerState.worker.joinable()) {
      workerState.worker.join();
    }
    return 0;
  }

  case WM_CLOSE:
    if (workerState.running.load()) {
      const int answer = MessageBoxA(hwnd,
                                     R"(A conversion is still running.

Closing now will terminate the full process immediately.
Any output files currently being written may be incomplete.

Do you want to abort the conversion and exit?)",
                                     "Abort conversion and exit",
                                     MB_ICONWARNING | MB_YESNO | MB_DEFBUTTON2);
      if (answer != IDYES) {
        return 0;
      }

      if (state->windowAlive) {
        state->windowAlive->store(false);
      }

      ExitProcess(0);
      return 0;
    }
    DestroyWindow(hwnd);
    return 0;

  case WM_DESTROY:
    if (state->windowAlive) {
      state->windowAlive->store(false);
    }
    if (workerState.worker.joinable()) {
      workerState.worker.join();
    }

    MSG pendingMessage{};
    while (PeekMessageA(&pendingMessage, hwnd, WM_APP_STATUS_UPDATE,
                        WM_APP_PROCESSING_DONE, PM_REMOVE)) {
      if (pendingMessage.message == WM_APP_STATUS_UPDATE) {
        delete reinterpret_cast<GuiStatusMessage *>(pendingMessage.lParam);
      } else if (pendingMessage.message == WM_APP_PROCESSING_DONE) {
        delete reinterpret_cast<GuiCompletionMessage *>(pendingMessage.lParam);
      }
    }

    PostQuitMessage(0);
    return 0;
  }

  return DefWindowProcA(hwnd, uMsg, wParam, lParam);
}

/*****************************************************************************
 * WinMain
 ****************************************************************************/
int WINAPI WinMain(HINSTANCE hInstance, HINSTANCE, LPSTR, int nCmdShow) {
  constexpr char CLASS_NAME[] = "xyz2dxfWindowClass";

  WNDCLASSA wc{};
  wc.style = CS_HREDRAW | CS_VREDRAW;
  wc.lpfnWndProc = WindowProc;
  wc.hInstance = hInstance;
  wc.lpszClassName = CLASS_NAME;
  wc.hbrBackground = reinterpret_cast<HBRUSH>(COLOR_BTNFACE + 1);
  wc.hCursor = LoadCursorA(nullptr, IDC_ARROW);

  if (RegisterClassA(&wc) == 0) {
    MessageBoxA(nullptr, "Failed to register the main window class.",
                "Startup Error", MB_ICONERROR | MB_OK);
    return 0;
  }

  RECT wr{0, 0, 950, 670};
  const DWORD windowStyle =
      WS_OVERLAPPED | WS_CAPTION | WS_SYSMENU | WS_MINIMIZEBOX;
  AdjustWindowRect(&wr, windowStyle, FALSE);

  HWND hwnd = CreateWindowExA(
      0, CLASS_NAME, "XYZ to DXF Converter - Professional Edition", windowStyle,
      CW_USEDEFAULT, CW_USEDEFAULT, wr.right - wr.left, wr.bottom - wr.top,
      nullptr, nullptr, hInstance, nullptr);
  if (!hwnd)
    return 0;

  ShowWindow(hwnd, nCmdShow);
  UpdateWindow(hwnd);

  MSG msg{};
  while (GetMessageA(&msg, nullptr, 0, 0)) {
    TranslateMessage(&msg);
    DispatchMessageA(&msg);
  }
  return 0;
}