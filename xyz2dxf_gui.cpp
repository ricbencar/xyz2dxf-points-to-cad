/******************************************************************************
 * XYZ to DXF Converter GUI
 * ------------------------
 * Windows GUI for reading large XYZ point clouds, cleaning them, interpolating
 * a regular surface grid, and exporting engineering-friendly XYZ and DXF
 * outputs.
 *
 * This is not just a format converter. The active pipeline performs
 * deterministic spatial thinning, robust local outlier rejection, interpolation,
 * confidence/diagnostic evaluation, and streamed output generation.
 *
 * Interpolation modes implemented in this build:
 *   1. Bicubic Hermite + robust anisotropic MLS
 *      - Primary smooth surface mode.
 *      - Uses robust local modeling to estimate nodal values/derivatives, then
 *        reconstructs bicubic patches on the output grid.
 *
 *   2. Local Thin Plate Spline (TPS)
 *      - Intended for scattered / irregular point sets.
 *      - Uses adaptive local neighborhoods, regularization, and blended local
 *        TPS patches.
 *
 *   3. Clamped local MLS (bounded)
 *      - Conservative local mode.
 *      - Evaluates a local MLS surface and clamps the result to the local data
 *        envelope to reduce overshoot.
 *
 * GUI inputs:
 *   - Input File      : source .xyz file selected with the standard Windows dialog
 *   - Min Dist        : minimum horizontal spacing for deterministic thinning
 *   - Precision       : decimal places written to output files
 *   - PDMODE          : DXF point display mode written to the DXF header
 *   - Grid Spacing    : spacing of the regular output grid
 *   - Max TPS Points  : optional control-point cap for TPS (0 = use all)
 *
 * Processing pipeline:
 *   1. Read XYZ input
 *   2. Deterministic minimum-distance thinning
 *   3. Robust local residual-based outlier rejection
 *   4. Build interpolation context / tuning data
 *   5. Interpolate onto a regular grid
 *   6. Compute local diagnostics / confidence
 *   7. Stream outputs to disk
 *
 * Output files:
 *   - input.xyz.filtered.xyz    cleaned point cloud after thinning/outlier removal
 *   - input.xyz.grid.xyz        interpolated regular grid
 *   - input.xyz.confidence.xyz  grid plus confidence / diagnostic values
 *   - input.xyz.dxf             DXF export for filtered points and optional grid
 *   - input.xyz.rpt.txt         run report, method summary, metrics, and status log
 *
 * Notes:
 *   - Precision affects file formatting only; internal calculations use double.
 *   - If Min Dist <= 0, the thinning stage falls back to exact duplicate removal.
 *   - If Max TPS Points > 0, the TPS branch may deterministically decimate the
 *     control set before model construction.
 *   - The code is designed to fail safely: singular/unstable local solves,
 *     insufficient neighborhoods, and I/O failures are checked and reported.
 *
 * Build (MinGW / g++):
 *   g++ -O3 -fopenmp -march=native -std=c++17 -Wall -Wextra -pedantic -Wconversion -Wsign-conversion -static -static-libgcc -static-libstdc++ -mwindows -o xyz2dxf_gui.exe xyz2dxf_gui.cpp -lkernel32 -lopengl32 -luuid -lcomdlg32 -lm
 *
 ******************************************************************************/

#include <windows.h>
#include <commdlg.h>
#include <string>
#include <sstream>
#include <iomanip>
#include <thread>
#include <mutex>
#include <future>
#include <chrono>
#include <functional>
#include <cmath>
#include <unordered_set>
#include <unordered_map>
#include <numeric>
#include <fstream>
#include <iostream>
#include <algorithm>
#include <random>
#include <array>
#include <vector>
#include <limits>
#include <memory>
#include <atomic>
#include <stdexcept>
#include <cstdint>
#include <omp.h>


/*****************************************************************************
 * GUI Component IDs, etc.
 ****************************************************************************/
#define IDC_INPUT_FILE     1001
#define IDC_MIN_DIST       1002
#define IDC_PRECISION      1003
#define IDC_PDMODE         1004
#define IDC_GRID_SPACING   1005
#define IDC_MAX_TPS_POINTS 1006
#define IDC_BROWSE_BUTTON  1007
#define IDC_RUN_BUTTON     1008
#define IDC_STATUS_STATIC  1009
#define IDC_RADIO_BICUBIC  1010
#define IDC_RADIO_TPS      1011
#define IDC_RADIO_CLAMPED_MLS 1012

enum InterpolationMethod
{
    METHOD_BICUBIC     = 0,
    METHOD_TPS         = 1,
    METHOD_CLAMPED_MLS = 2
};

/*****************************************************************************
 * Data Structures
 ****************************************************************************/
struct Point3D
{
    double x, y, z; // 3D coordinates

    bool operator==(const Point3D &other) const
    {
        return (x == other.x && y == other.y && z == other.z);
    }
};

// Hash functor for Point3D to enable usage in std::unordered_set
struct Point3DHash
{
    size_t operator()(const Point3D &p) const
    {
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

struct Cell2D
{
    std::int64_t ix = 0;
    std::int64_t iy = 0;

    bool operator==(const Cell2D &other) const
    {
        return ix == other.ix && iy == other.iy;
    }
};

struct Cell2DHash
{
    size_t operator()(const Cell2D &c) const
    {
        const auto h1 = std::hash<std::int64_t>{}(c.ix);
        const auto h2 = std::hash<std::int64_t>{}(c.iy);
        size_t seed = 0;
        seed ^= h1 + 0x9e3779b97f4a7c15ULL + (seed << 6) + (seed >> 2);
        seed ^= h2 + 0x9e3779b97f4a7c15ULL + (seed << 6) + (seed >> 2);
        return seed;
    }
};

static std::int64_t cellCoord(double value, double cellSize)
{
    return static_cast<std::int64_t>(std::floor(value / cellSize));
}

struct Grid2D
{
    size_t nx = 0;
    size_t ny = 0;
    std::vector<double> values;

    Grid2D() = default;
    Grid2D(size_t nx_, size_t ny_, double init = 0.0)
        : nx(nx_), ny(ny_), values(nx_ * ny_, init) {}

    double &operator()(size_t i, size_t j)
    {
        return values[i * ny + j];
    }

    const double &operator()(size_t i, size_t j) const
    {
        return values[i * ny + j];
    }
};

struct SpatialHash2D
{
    double cellSize = 1.0;
    double xMin = 0.0;
    double yMin = 0.0;
    std::unordered_map<long long, std::vector<size_t>> buckets;

    static long long cellKey(int ix, int iy)
    {
        return (static_cast<long long>(ix) << 32) ^ static_cast<unsigned int>(iy);
    }

    void build(const std::vector<Point3D> &points, double requestedCellSize)
    {
        buckets.clear();
        cellSize = std::max(requestedCellSize, 1e-9);
        xMin = std::numeric_limits<double>::max();
        yMin = std::numeric_limits<double>::max();
        for (const auto &p : points)
        {
            if (p.x < xMin) xMin = p.x;
            if (p.y < yMin) yMin = p.y;
        }
        if (!std::isfinite(xMin)) xMin = 0.0;
        if (!std::isfinite(yMin)) yMin = 0.0;
        for (size_t i = 0; i < points.size(); ++i)
        {
            const int ix = static_cast<int>(std::floor((points[i].x - xMin) / cellSize));
            const int iy = static_cast<int>(std::floor((points[i].y - yMin) / cellSize));
            buckets[cellKey(ix, iy)].push_back(i);
        }
    }

    std::vector<size_t> gatherNearest(const std::vector<Point3D> &points,
                                      double x,
                                      double y,
                                      size_t desiredCount,
                                      size_t minCount) const
    {
        std::vector<size_t> candidates;
        if (points.empty())
        {
            return candidates;
        }

        desiredCount = std::min(desiredCount, points.size());
        minCount = std::min(minCount, points.size());

        const int ix0 = static_cast<int>(std::floor((x - xMin) / cellSize));
        const int iy0 = static_cast<int>(std::floor((y - yMin) / cellSize));

        const int maxRing = 64;
        for (int ring = 0; ring <= maxRing; ++ring)
        {
            for (int dx = -ring; dx <= ring; ++dx)
            {
                for (int dy = -ring; dy <= ring; ++dy)
                {
                    if (std::max(std::abs(dx), std::abs(dy)) != ring)
                    {
                        continue;
                    }
                    const auto it = buckets.find(cellKey(ix0 + dx, iy0 + dy));
                    if (it != buckets.end())
                    {
                        candidates.insert(candidates.end(), it->second.begin(), it->second.end());
                    }
                }
            }
            if (candidates.size() >= minCount && (candidates.size() >= desiredCount * 2 || ring >= 4))
            {
                break;
            }
        }

        if (candidates.size() < minCount)
        {
            candidates.resize(points.size());
            std::iota(candidates.begin(), candidates.end(), static_cast<size_t>(0));
        }

        std::sort(candidates.begin(), candidates.end(), [&](size_t a, size_t b)
        {
            const double da = (points[a].x - x) * (points[a].x - x) + (points[a].y - y) * (points[a].y - y);
            const double db = (points[b].x - x) * (points[b].x - x) + (points[b].y - y) * (points[b].y - y);
            if (da != db)
            {
                return da < db;
            }
            return a < b;
        });

        if (candidates.size() > desiredCount)
        {
            candidates.resize(desiredCount);
        }
        return candidates;
    }
};

struct ValidationMetrics
{
    double rmse = std::numeric_limits<double>::infinity();
    double mae = std::numeric_limits<double>::infinity();
    double p95 = std::numeric_limits<double>::infinity();
    size_t count = 0;
    bool valid = false;
};

struct InterpolationDiagnostics
{
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

struct BicubicNodeData
{
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

struct ConvexHull2D
{
    std::vector<Point3D> vertices;
    double cellSize = 1.0;
    double xMin = 0.0;
    double yMin = 0.0;
    std::unordered_set<long long> occupied;
    bool valid = false;
};

struct TPSModel
{
    std::vector<Point3D> controlPoints;
    SpatialHash2D spatialIndex;
    ConvexHull2D hull;
    size_t neighborhoodSize = 48;
    double suggestedLambda = 1e-10;
    ValidationMetrics validation;
    bool valid = false;
};

struct GridDefinition
{
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

struct TuningChoice
{
    size_t neighborhoodTarget = 0;
    double supportMultiplier = 1.75;
    double ridgeFactor = 1.0;
    int basisOrder = 3;
    ValidationMetrics metrics;
};

struct ConfidenceSummary
{
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
struct BicubicSpline
{
    double a00, a01, a02, a03;
    double a10, a11, a12, a13;
    double a20, a21, a22, a23;
    double a30, a31, a32, a33;
};

/*****************************************************************************
 * Forward Declarations
 ****************************************************************************/
struct LocalFrame2D;
static bool solveDenseLinearSystem(std::vector<double> A,
                                   std::vector<double> b,
                                   size_t n,
                                   std::vector<double> &x);
static double averagePlanarSpacing(const std::vector<Point3D> &points);
static bool fitRobustLocalMLS(const std::vector<Point3D> &points,
                              const std::vector<size_t> &fitIndices,
                              const LocalFrame2D &frame,
                              double supportRadius,
                              int basisOrder,
                              double ridgeFactor,
                              std::vector<double> &coeffs,
                              ValidationMetrics *fitMetrics);
static bool solveLocalTPSAtPoint(double x,
                                 double y,
                                 const TPSModel &model,
                                 double &z,
                                 std::string *diagnostics,
                                 InterpolationDiagnostics *interpDiagnostics);
static void updateStatus(HWND hwndStatus, const std::string &message);
static std::string openFileDialog(HWND hwnd);
static void computeBoundingBox(const std::vector<Point3D> &points,
                               double &xMin, double &xMax,
                               double &yMin, double &yMax);

[[maybe_unused]]
static std::vector<Point3D> filterPointsGrid(const std::vector<Point3D> &points,
                                             double minDist);
static bool readAndFilterXYZFile(const std::string &inputFileName,
                                 double minDist,
                                 std::vector<Point3D> &filteredPoints,
                                 size_t &totalPointsRead,
                                 const std::function<void(const std::string &)> &progressCallback,
                                 std::string *errorMessage);
[[maybe_unused]]
static std::vector<Point3D> removeZOutliers(const std::vector<Point3D> &points,
                                            double neighborDist,
                                            double zThresholdFactor);

/*****************************************************************************
 * Basic GUI / geometry helpers
 ****************************************************************************/
static void updateStatus(HWND hwndStatus, const std::string &message)
{
    SetWindowTextA(hwndStatus, message.c_str());
}

static std::string openFileDialog(HWND hwnd)
{
    OPENFILENAMEA ofn;
    char szFile[260];
    ZeroMemory(&ofn, sizeof(ofn));
    ZeroMemory(szFile, sizeof(szFile));

    ofn.lStructSize = sizeof(ofn);
    ofn.hwndOwner = hwnd;
    ofn.lpstrFile = szFile;
    ofn.nMaxFile = static_cast<DWORD>(sizeof(szFile));
    ofn.lpstrFilter = "XYZ Files (*.xyz)\0*.xyz\0All Files (*.*)\0*.*\0";
    ofn.nFilterIndex = 1;
    ofn.Flags = OFN_PATHMUSTEXIST | OFN_FILEMUSTEXIST;

    if (GetOpenFileNameA(&ofn))
    {
        return std::string(ofn.lpstrFile);
    }
    return "";
}

static void computeBoundingBox(const std::vector<Point3D> &points,
                               double &xMin, double &xMax,
                               double &yMin, double &yMax)
{
    if (points.empty())
    {
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

    for (const auto &p : points)
    {
        xMin = std::min(xMin, p.x);
        xMax = std::max(xMax, p.x);
        yMin = std::min(yMin, p.y);
        yMax = std::max(yMax, p.y);
    }
}

[[maybe_unused]]
static std::vector<Point3D> filterPointsGrid(const std::vector<Point3D> &points,
                                             double minDist)
{
    if (points.empty())
    {
        return {};
    }

    if (!(minDist > 0.0))
    {
        std::unordered_set<Point3D, Point3DHash> seen;
        std::vector<Point3D> uniquePoints;
        uniquePoints.reserve(points.size());
        for (const auto &p : points)
        {
            if (seen.emplace(p).second)
            {
                uniquePoints.push_back(p);
            }
        }
        std::sort(uniquePoints.begin(), uniquePoints.end(), [](const Point3D &a, const Point3D &b)
        {
            if (a.x != b.x) return a.x < b.x;
            if (a.y != b.y) return a.y < b.y;
            return a.z < b.z;
        });
        return uniquePoints;
    }

    struct CandidatePoint
    {
        Point3D p;
        Cell2D cell;
        double score = 0.0;
    };

    std::vector<CandidatePoint> candidates;
    candidates.reserve(points.size());
    for (const auto &p : points)
    {
        const Cell2D cell{cellCoord(p.x, minDist), cellCoord(p.y, minDist)};
        const double cx = (static_cast<double>(cell.ix) + 0.5) * minDist;
        const double cy = (static_cast<double>(cell.iy) + 0.5) * minDist;
        const double dx = p.x - cx;
        const double dy = p.y - cy;
        candidates.push_back({p, cell, dx * dx + dy * dy});
    }
    std::stable_sort(candidates.begin(), candidates.end(), [](const CandidatePoint &a, const CandidatePoint &b)
    {
        if (a.cell.ix != b.cell.ix) return a.cell.ix < b.cell.ix;
        if (a.cell.iy != b.cell.iy) return a.cell.iy < b.cell.iy;
        if (a.score != b.score) return a.score < b.score;
        if (a.p.x != b.p.x) return a.p.x < b.p.x;
        if (a.p.y != b.p.y) return a.p.y < b.p.y;
        return a.p.z < b.p.z;
    });

    const double minDistSq = minDist * minDist;
    std::unordered_map<Cell2D, std::vector<Point3D>, Cell2DHash> sparseGrid;
    sparseGrid.reserve(std::min<size_t>(candidates.size(), 1U << 20));

    std::vector<Point3D> accepted;
    accepted.reserve(points.size());

    for (const auto &cand : candidates)
    {
        bool tooClose = false;
        for (std::int64_t dx = -1; dx <= 1 && !tooClose; ++dx)
        {
            for (std::int64_t dy = -1; dy <= 1 && !tooClose; ++dy)
            {
                const auto it = sparseGrid.find(Cell2D{cand.cell.ix + dx, cand.cell.iy + dy});
                if (it == sparseGrid.end())
                {
                    continue;
                }
                for (const auto &q : it->second)
                {
                    const double ddx = cand.p.x - q.x;
                    const double ddy = cand.p.y - q.y;
                    if (ddx * ddx + ddy * ddy < minDistSq)
                    {
                        tooClose = true;
                        break;
                    }
                }
            }
        }
        if (!tooClose)
        {
            sparseGrid[cand.cell].push_back(cand.p);
            accepted.push_back(cand.p);
        }
    }

    return accepted;
}

static bool readAndFilterXYZFile(const std::string &inputFileName,
                                 double minDist,
                                 std::vector<Point3D> &filteredPoints,
                                 size_t &totalPointsRead,
                                 const std::function<void(const std::string &)> &progressCallback,
                                 std::string *errorMessage)
{
    filteredPoints.clear();
    totalPointsRead = 0U;

    std::ifstream inFile(inputFileName);
    if (!inFile.is_open())
    {
        if (errorMessage)
        {
            *errorMessage = "Error: Unable to open input file.";
        }
        return false;
    }

    constexpr size_t progressStep = 1000000U;
    std::string line;
    std::vector<Point3D> rawPoints;
    rawPoints.reserve(1U << 20);

    while (std::getline(inFile, line))
    {
        if (line.empty() || line[0] == '#')
        {
            continue;
        }
        std::replace(line.begin(), line.end(), ',', ' ');
        std::istringstream ss(line);
        Point3D p{};
        if (!(ss >> p.x >> p.y >> p.z))
        {
            continue;
        }
        ++totalPointsRead;
        rawPoints.push_back(p);
        if (progressCallback && (totalPointsRead % progressStep) == 0U)
        {
            std::ostringstream msg;
            msg << "Read " << totalPointsRead << " points into deterministic pre-filter buffer ...";
            progressCallback(msg.str());
        }
    }

    if (rawPoints.empty())
    {
        if (errorMessage)
        {
            *errorMessage = "Error: No valid points found in file.";
        }
        return false;
    }

    filteredPoints = filterPointsGrid(rawPoints, minDist);
    if (filteredPoints.empty())
    {
        if (errorMessage)
        {
            *errorMessage = "Error: No valid points found in file after deterministic filtering.";
        }
        return false;
    }

    if (progressCallback)
    {
        std::ostringstream msg;
        msg << "Deterministic minDist filtering complete: retained " << filteredPoints.size() << " of " << rawPoints.size() << " points.";
        progressCallback(msg.str());
    }
    return true;
}

[[maybe_unused]]
static std::vector<Point3D> removeZOutliers(const std::vector<Point3D> &points,
                                            double neighborDist,
                                            double zThresholdFactor)
{
    if (points.empty())
    {
        return points;
    }

    if (!(zThresholdFactor > 0.0))
    {
        return points;
    }

    if (!(neighborDist > 0.0))
    {
        double sumZ = 0.0;
        double sumZ2 = 0.0;
        for (const auto &p : points)
        {
            sumZ += p.z;
            sumZ2 += p.z * p.z;
        }

        const double n = static_cast<double>(points.size());
        const double meanZ = sumZ / n;
        double varZ = sumZ2 / n - meanZ * meanZ;
        if (varZ < 0.0)
        {
            varZ = 0.0;
        }
        const double stdevZ = std::sqrt(varZ);
        if (stdevZ < 1e-15)
        {
            return points;
        }

        std::vector<Point3D> result;
        result.reserve(points.size());
        for (const auto &p : points)
        {
            if (std::fabs(p.z - meanZ) <= zThresholdFactor * stdevZ)
            {
                result.push_back(p);
            }
        }
        return result;
    }

    const double neighborDistSq = neighborDist * neighborDist;
    std::unordered_map<Cell2D, std::vector<size_t>, Cell2DHash> sparseGrid;
    sparseGrid.reserve(std::min<size_t>(points.size(), 1U << 20));

    for (size_t i = 0; i < points.size(); ++i)
    {
        sparseGrid[Cell2D{cellCoord(points[i].x, neighborDist), cellCoord(points[i].y, neighborDist)}].push_back(i);
    }

    std::vector<Point3D> filtered;
    filtered.reserve(points.size());

    for (size_t i = 0; i < points.size(); ++i)
    {
        const Point3D &pi = points[i];
        const Cell2D cell{cellCoord(pi.x, neighborDist), cellCoord(pi.y, neighborDist)};

        double sumZ = 0.0;
        double sumZ2 = 0.0;
        size_t count = 0U;

        for (std::int64_t dx = -1; dx <= 1; ++dx)
        {
            for (std::int64_t dy = -1; dy <= 1; ++dy)
            {
                const auto it = sparseGrid.find(Cell2D{cell.ix + dx, cell.iy + dy});
                if (it == sparseGrid.end())
                {
                    continue;
                }
                for (size_t j : it->second)
                {
                    const double ddx = pi.x - points[j].x;
                    const double ddy = pi.y - points[j].y;
                    if (ddx * ddx + ddy * ddy <= neighborDistSq)
                    {
                        sumZ += points[j].z;
                        sumZ2 += points[j].z * points[j].z;
                        ++count;
                    }
                }
            }
        }

        if (count < 2U)
        {
            filtered.push_back(pi);
            continue;
        }

        const double n = static_cast<double>(count);
        const double meanZ = sumZ / n;
        double varZ = sumZ2 / n - meanZ * meanZ;
        if (varZ < 0.0)
        {
            varZ = 0.0;
        }
        const double stdevZ = std::sqrt(varZ);
        if (stdevZ < 1e-15 || std::fabs(pi.z - meanZ) <= zThresholdFactor * stdevZ)
        {
            filtered.push_back(pi);
        }
    }

    return filtered;
}

static double computeMedianUnchecked(std::vector<double> values)
{
    if (values.empty())
    {
        return 0.0;
    }
    const size_t mid = values.size() / 2U;
    std::nth_element(values.begin(), values.begin() + static_cast<std::ptrdiff_t>(mid), values.end());
    double med = values[mid];
    if ((values.size() % 2U) == 0U)
    {
        std::nth_element(values.begin(), values.begin() + static_cast<std::ptrdiff_t>(mid - 1U), values.end());
        med = 0.5 * (med + values[mid - 1U]);
    }
    return med;
}

static ValidationMetrics makeValidationMetricsFromErrors(std::vector<double> errors)
{
    ValidationMetrics out;
    if (errors.empty())
    {
        return out;
    }
    double sumSq = 0.0;
    double sumAbs = 0.0;
    for (double e : errors)
    {
        const double ae = std::fabs(e);
        sumSq += ae * ae;
        sumAbs += ae;
    }
    std::sort(errors.begin(), errors.end(), [](double a, double b)
    {
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

static double validationScore(const ValidationMetrics &m)
{
    if (!m.valid)
    {
        return std::numeric_limits<double>::infinity();
    }
    return m.rmse + 0.35 * m.p95 + 0.15 * m.mae;
}

static double calibratedConfidenceFromObservedError(double rawConfidence,
                                                    double localResidual,
                                                    const ValidationMetrics &observedMetrics)
{
    if (!std::isfinite(rawConfidence) || !(rawConfidence > 0.0))
    {
        rawConfidence = 0.05;
    }
    if (!observedMetrics.valid)
    {
        return std::max(0.05, std::min(1.0, rawConfidence));
    }
    const double scale = std::max({observedMetrics.rmse, 0.5 * observedMetrics.p95, 1e-8});
    const double penalty = 1.0 / (1.0 + std::max(localResidual, 0.0) / scale);
    return std::max(0.02, std::min(1.0, rawConfidence * penalty));
}

static std::vector<size_t> makeDeterministicHoldoutIndices(size_t n, size_t maxCount)
{
    std::vector<size_t> out;
    if (n == 0U || maxCount == 0U)
    {
        return out;
    }
    const size_t target = std::min(maxCount, n);
    const size_t step = std::max<size_t>(1U, n / target);
    for (size_t i = step / 2U; i < n && out.size() < target; i += step)
    {
        out.push_back(i);
    }
    if (out.empty())
    {
        out.push_back(0U);
    }
    return out;
}

static long long supportMaskKey(int ix, int iy)
{
    return (static_cast<long long>(ix) << 32) ^ static_cast<unsigned int>(iy);
}

static double cross2D(const Point3D &a, const Point3D &b, const Point3D &c)
{
    return (b.x - a.x) * (c.y - a.y) - (b.y - a.y) * (c.x - a.x);
}

static ConvexHull2D buildConvexHull2D(std::vector<Point3D> points)
{
    ConvexHull2D hull;
    if (points.empty())
    {
        return hull;
    }

    double xMin = 0.0;
    double xMax = 0.0;
    double yMin = 0.0;
    double yMax = 0.0;
    computeBoundingBox(points, xMin, xMax, yMin, yMax);
    hull.cellSize = std::max(2.5 * averagePlanarSpacing(points), 1e-6);
    hull.xMin = xMin;
    hull.yMin = yMin;

    if (points.size() >= 3U)
    {
        std::sort(points.begin(), points.end(), [](const Point3D &a, const Point3D &b)
        {
            if (a.x != b.x) return a.x < b.x;
            if (a.y != b.y) return a.y < b.y;
            return a.z < b.z;
        });
        points.erase(std::unique(points.begin(), points.end(), [](const Point3D &a, const Point3D &b)
        {
            return a.x == b.x && a.y == b.y;
        }), points.end());

        if (points.size() >= 3U)
        {
            std::vector<Point3D> lower;
            std::vector<Point3D> upper;
            for (const auto &p : points)
            {
                while (lower.size() >= 2U && cross2D(lower[lower.size() - 2U], lower.back(), p) <= 0.0)
                {
                    lower.pop_back();
                }
                lower.push_back({p.x, p.y, 0.0});
            }
            for (size_t ii = points.size(); ii-- > 0U; )
            {
                const auto &p = points[ii];
                while (upper.size() >= 2U && cross2D(upper[upper.size() - 2U], upper.back(), p) <= 0.0)
                {
                    upper.pop_back();
                }
                upper.push_back({p.x, p.y, 0.0});
            }
            if (!lower.empty()) lower.pop_back();
            if (!upper.empty()) upper.pop_back();
            hull.vertices = lower;
            hull.vertices.insert(hull.vertices.end(), upper.begin(), upper.end());
        }
    }

    const int dilation = 1;
    for (const auto &p : points)
    {
        const int ix = static_cast<int>(std::floor((p.x - hull.xMin) / hull.cellSize));
        const int iy = static_cast<int>(std::floor((p.y - hull.yMin) / hull.cellSize));
        for (int dx = -dilation; dx <= dilation; ++dx)
        {
            for (int dy = -dilation; dy <= dilation; ++dy)
            {
                hull.occupied.insert(supportMaskKey(ix + dx, iy + dy));
            }
        }
    }
    hull.valid = !hull.occupied.empty();
    return hull;
}

static double distancePointToSegment2D(double x, double y,
                                       double x0, double y0,
                                       double x1, double y1)
{
    const double dx = x1 - x0;
    const double dy = y1 - y0;
    const double l2 = dx * dx + dy * dy;
    if (!(l2 > 0.0))
    {
        const double dx0 = x - x0;
        const double dy0 = y - y0;
        return std::sqrt(dx0 * dx0 + dy0 * dy0);
    }
    const double t = std::min(1.0, std::max(0.0, ((x - x0) * dx + (y - y0) * dy) / l2));
    const double px = x0 + t * dx;
    const double py = y0 + t * dy;
    const double ddx = x - px;
    const double ddy = y - py;
    return std::sqrt(ddx * ddx + ddy * ddy);
}

static bool pointInConvexHull2D(const ConvexHull2D &hull, double x, double y)
{
    if (!hull.valid)
    {
        return false;
    }
    const int ix = static_cast<int>(std::floor((x - hull.xMin) / std::max(hull.cellSize, 1e-6)));
    const int iy = static_cast<int>(std::floor((y - hull.yMin) / std::max(hull.cellSize, 1e-6)));
    if (hull.occupied.find(supportMaskKey(ix, iy)) != hull.occupied.end())
    {
        return true;
    }
    for (int dx = -1; dx <= 1; ++dx)
    {
        for (int dy = -1; dy <= 1; ++dy)
        {
            if (hull.occupied.find(supportMaskKey(ix + dx, iy + dy)) != hull.occupied.end())
            {
                return true;
            }
        }
    }
    if (hull.vertices.size() < 3U)
    {
        return false;
    }
    bool hasPos = false;
    bool hasNeg = false;
    for (size_t i = 0; i < hull.vertices.size(); ++i)
    {
        const Point3D &a = hull.vertices[i];
        const Point3D &b = hull.vertices[(i + 1U) % hull.vertices.size()];
        const double cp = (b.x - a.x) * (y - a.y) - (b.y - a.y) * (x - a.x);
        hasPos = hasPos || (cp > 1e-12);
        hasNeg = hasNeg || (cp < -1e-12);
        if (hasPos && hasNeg)
        {
            return false;
        }
    }
    return true;
}

static double distanceToConvexHull2D(const ConvexHull2D &hull, double x, double y)
{
    if (!hull.valid)
    {
        return 0.0;
    }
    if (pointInConvexHull2D(hull, x, y))
    {
        return 0.0;
    }

    double best = std::numeric_limits<double>::infinity();
    if (!hull.occupied.empty())
    {
        const double cs = std::max(hull.cellSize, 1e-6);
        const int ix = static_cast<int>(std::floor((x - hull.xMin) / cs));
        const int iy = static_cast<int>(std::floor((y - hull.yMin) / cs));
        for (int rx = -3; rx <= 3; ++rx)
        {
            for (int ry = -3; ry <= 3; ++ry)
            {
                const long long key = supportMaskKey(ix + rx, iy + ry);
                if (hull.occupied.find(key) == hull.occupied.end())
                {
                    continue;
                }
                const double cx = hull.xMin + (static_cast<double>(ix + rx) + 0.5) * cs;
                const double cy = hull.yMin + (static_cast<double>(iy + ry) + 0.5) * cs;
                const double dx = x - cx;
                const double dy = y - cy;
                best = std::min(best, std::sqrt(dx * dx + dy * dy));
            }
        }
    }
    if (hull.vertices.size() >= 2U)
    {
        for (size_t i = 0; i < hull.vertices.size(); ++i)
        {
            const Point3D &a = hull.vertices[i];
            const Point3D &b = hull.vertices[(i + 1U) % hull.vertices.size()];
            best = std::min(best, distancePointToSegment2D(x, y, a.x, a.y, b.x, b.y));
        }
    }
    return std::isfinite(best) ? best : 0.0;
}

static bool solveWeightedPlaneAtPoint(const std::vector<Point3D> &points,
                                      const std::vector<size_t> &neighbors,
                                      double x0,
                                      double y0,
                                      double &a,
                                      double &bx,
                                      double &by,
                                      ValidationMetrics *residualMetrics,
                                      size_t excludeIndex = std::numeric_limits<size_t>::max())
{
    a = 0.0;
    bx = 0.0;
    by = 0.0;
    if (neighbors.size() < 3U)
    {
        return false;
    }

    double meanR2 = 0.0;
    size_t usedCount = 0U;
    for (size_t idx : neighbors)
    {
        const double dx = points[idx].x - x0;
        const double dy = points[idx].y - y0;
        meanR2 += dx * dx + dy * dy;
        ++usedCount;
    }
    if (usedCount < 3U)
    {
        return false;
    }
    meanR2 /= static_cast<double>(usedCount);
    const double h = std::max(std::sqrt(std::max(meanR2, 1e-12)), 1e-6);

    std::vector<double> ATA(9U, 0.0);
    std::vector<double> ATb(3U, 0.0);
    for (size_t idx : neighbors)
    {
        const double dx = points[idx].x - x0;
        const double dy = points[idx].y - y0;
        const double r2 = dx * dx + dy * dy;
        const double w = 1.0 / (1.0 + r2 / (h * h));
        const double phi[3] = {1.0, dx, dy};
        for (size_t r = 0; r < 3U; ++r)
        {
            ATb[r] += w * phi[r] * points[idx].z;
            for (size_t c = 0; c < 3U; ++c)
            {
                ATA[r * 3U + c] += w * phi[r] * phi[c];
            }
        }
    }
    for (size_t d = 0; d < 3U; ++d)
    {
        ATA[d * 3U + d] += 1e-12;
    }

    std::vector<double> coeffs;
    if (!solveDenseLinearSystem(ATA, ATb, 3U, coeffs) || coeffs.size() < 3U)
    {
        return false;
    }
    a = coeffs[0];
    bx = coeffs[1];
    by = coeffs[2];

    if (residualMetrics)
    {
        std::vector<double> residuals;
        residuals.reserve(neighbors.size());
        for (size_t idx : neighbors)
        {
            if (idx == excludeIndex)
            {
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

static bool removeZOutliersInPlace(std::vector<Point3D> &points,
                                   double neighborDist,
                                   double zThresholdFactor,
                                   const std::function<void(const std::string &)> &progressCallback,
                                   std::string *errorMessage)
{
    if (points.empty())
    {
        return true;
    }
    if (!(zThresholdFactor > 0.0) || !(neighborDist > 0.0))
    {
        return true;
    }

    const double neighborDistSq = neighborDist * neighborDist;
    SpatialHash2D spatialIndex;
    spatialIndex.build(points, std::max(neighborDist, 1e-6));
    std::vector<unsigned char> keep(points.size(), 1U);

    const size_t reportStep = std::max<size_t>(250000U, points.size() / 100U + 1U);
    size_t kept = 0U;
    size_t rejected = 0U;

    for (size_t i = 0; i < points.size(); ++i)
    {
        const Point3D &pi = points[i];
        const auto candidates = spatialIndex.gatherNearest(points, pi.x, pi.y, std::min<size_t>(points.size(), 96U), std::min<size_t>(points.size(), 12U));
        std::vector<size_t> neighbors;
        neighbors.reserve(candidates.size());
        for (size_t idx : candidates)
        {
            const double dx = points[idx].x - pi.x;
            const double dy = points[idx].y - pi.y;
            if (idx != i && dx * dx + dy * dy <= neighborDistSq)
            {
                neighbors.push_back(idx);
            }
        }
        if (neighbors.size() < 6U)
        {
            keep[i] = 1U;
            ++kept;
            continue;
        }

        double a = 0.0;
        double bx = 0.0;
        double by = 0.0;
        ValidationMetrics residualMetrics;
        if (!solveWeightedPlaneAtPoint(points, neighbors, pi.x, pi.y, a, bx, by, &residualMetrics, i))
        {
            keep[i] = 1U;
            ++kept;
            continue;
        }

        std::vector<double> residuals;
        residuals.reserve(neighbors.size());
        for (size_t idx : neighbors)
        {
            const double dx = points[idx].x - pi.x;
            const double dy = points[idx].y - pi.y;
            residuals.push_back(points[idx].z - (a + bx * dx + by * dy));
        }
        const double medianResidual = computeMedianUnchecked(residuals);
        for (double &r : residuals)
        {
            r = std::fabs(r - medianResidual);
        }
        const double mad = computeMedianUnchecked(residuals);
        const double robustSigma = std::max(1.4826 * mad, 1e-8);
        const double selfResidual = std::fabs(pi.z - a - medianResidual);
        const bool keepPoint = (selfResidual <= zThresholdFactor * robustSigma);
        keep[i] = keepPoint ? 1U : 0U;
        if (keepPoint)
        {
            ++kept;
        }
        else
        {
            ++rejected;
        }

        if (progressCallback && ((i + 1U) % reportStep) == 0U)
        {
            std::ostringstream msg;
            const double pct = 100.0 * static_cast<double>(i + 1U) / static_cast<double>(points.size());
            msg << "Robust residual outlier scan " << (i + 1U) << "/" << points.size()
                << ", kept=" << kept << ", rejected=" << rejected
                << " (" << std::fixed << std::setprecision(1) << pct << "%)";
            progressCallback(msg.str());
        }
    }

    size_t writeIdx = 0U;
    for (size_t i = 0; i < points.size(); ++i)
    {
        if (keep[i] != 0U)
        {
            if (writeIdx != i)
            {
                points[writeIdx] = std::move(points[i]);
            }
            ++writeIdx;
        }
    }
    points.resize(writeIdx);
    points.shrink_to_fit();

    if (points.empty())
    {
        if (errorMessage)
        {
            *errorMessage = "Error: All points were removed during robust residual outlier filtering.";
        }
        return false;
    }
    if (progressCallback)
    {
        std::ostringstream msg;
        msg << "Robust residual outlier filtering complete: retained " << points.size()
            << " points, removed " << rejected << ".";
        progressCallback(msg.str());
    }
    return true;
}

static std::vector<Point3D> subsamplePointsUniformly(const std::vector<Point3D> &points,
                                                     size_t maxTPSPoints)
{
    if (maxTPSPoints == 0U || points.size() <= maxTPSPoints)
    {
        return points;
    }

    std::vector<Point3D> candidatePool = points;
    if (points.size() > maxTPSPoints * 12U)
    {
        double xMin, xMax, yMin, yMax;
        computeBoundingBox(points, xMin, xMax, yMin, yMax);
        const size_t gridCount = static_cast<size_t>(std::ceil(std::sqrt(static_cast<double>(maxTPSPoints * 8U))));
        const double cellWidth = std::max((xMax - xMin) / static_cast<double>(std::max<size_t>(gridCount, 1U)), 1e-9);
        const double cellHeight = std::max((yMax - yMin) / static_cast<double>(std::max<size_t>(gridCount, 1U)), 1e-9);
        std::unordered_map<long long, size_t> chosen;
        candidatePool.clear();
        candidatePool.reserve(maxTPSPoints * 8U);
        for (size_t i = 0; i < points.size(); ++i)
        {
            const int ix = static_cast<int>(std::floor((points[i].x - xMin) / cellWidth));
            const int iy = static_cast<int>(std::floor((points[i].y - yMin) / cellHeight));
            const long long key = SpatialHash2D::cellKey(ix, iy);
            const double cx = xMin + (static_cast<double>(ix) + 0.5) * cellWidth;
            const double cy = yMin + (static_cast<double>(iy) + 0.5) * cellHeight;
            const double score = (points[i].x - cx) * (points[i].x - cx) + (points[i].y - cy) * (points[i].y - cy);
            auto it = chosen.find(key);
            if (it == chosen.end())
            {
                chosen.emplace(key, candidatePool.size());
                candidatePool.push_back(points[i]);
            }
            else
            {
                const Point3D &current = candidatePool[it->second];
                const double currentScore = (current.x - cx) * (current.x - cx) + (current.y - cy) * (current.y - cy);
                if (score < currentScore)
                {
                    candidatePool[it->second] = points[i];
                }
            }
        }
        if (candidatePool.size() < maxTPSPoints)
        {
            candidatePool = points;
        }
    }

    if (candidatePool.size() <= maxTPSPoints)
    {
        return candidatePool;
    }

    std::vector<Point3D> selected;
    selected.reserve(maxTPSPoints);
    std::vector<double> minDist2(candidatePool.size(), std::numeric_limits<double>::infinity());

    double cx = 0.0;
    double cy = 0.0;
    for (const auto &p : candidatePool)
    {
        cx += p.x;
        cy += p.y;
    }
    cx /= static_cast<double>(candidatePool.size());
    cy /= static_cast<double>(candidatePool.size());

    size_t firstIdx = 0U;
    double bestFirst = std::numeric_limits<double>::infinity();
    for (size_t i = 0; i < candidatePool.size(); ++i)
    {
        const double d2 = (candidatePool[i].x - cx) * (candidatePool[i].x - cx) + (candidatePool[i].y - cy) * (candidatePool[i].y - cy);
        if (d2 < bestFirst)
        {
            bestFirst = d2;
            firstIdx = i;
        }
    }
    selected.push_back(candidatePool[firstIdx]);

    for (size_t i = 0; i < candidatePool.size(); ++i)
    {
        const double dx = candidatePool[i].x - candidatePool[firstIdx].x;
        const double dy = candidatePool[i].y - candidatePool[firstIdx].y;
        minDist2[i] = dx * dx + dy * dy;
    }
    minDist2[firstIdx] = -1.0;

    while (selected.size() < maxTPSPoints)
    {
        size_t bestIdx = 0U;
        double bestDist2 = -1.0;
        for (size_t i = 0; i < candidatePool.size(); ++i)
        {
            if (minDist2[i] > bestDist2)
            {
                bestDist2 = minDist2[i];
                bestIdx = i;
            }
        }
        if (!(bestDist2 > 0.0) || !std::isfinite(bestDist2))
        {
            break;
        }
        selected.push_back(candidatePool[bestIdx]);
        for (size_t i = 0; i < candidatePool.size(); ++i)
        {
            if (minDist2[i] < 0.0)
            {
                continue;
            }
            const double dx = candidatePool[i].x - candidatePool[bestIdx].x;
            const double dy = candidatePool[i].y - candidatePool[bestIdx].y;
            const double d2 = dx * dx + dy * dy;
            if (d2 < minDist2[i])
            {
                minDist2[i] = d2;
            }
        }
        minDist2[bestIdx] = -1.0;
    }

    if (selected.size() > maxTPSPoints)
    {
        selected.resize(maxTPSPoints);
    }
    return selected;
}

/****************************************************************************
 * 7) High-accuracy interpolation helpers
 ****************************************************************************/
static double averagePlanarSpacing(const std::vector<Point3D> &points)
{
    if (points.size() < 2)
    {
        return 1.0;
    }

    double xMin, xMax, yMin, yMax;
    computeBoundingBox(points, xMin, xMax, yMin, yMax);
    const double area = std::max((xMax - xMin) * (yMax - yMin), 1e-12);
    return std::sqrt(area / static_cast<double>(points.size()));
}

static bool solveDenseLinearSystem(std::vector<double> A,
                                   std::vector<double> b,
                                   size_t n,
                                   std::vector<double> &x)
{
    auto idx = [n](size_t r, size_t c) -> size_t
    {
        return r * n + c;
    };

    x.assign(n, 0.0);
    for (size_t k = 0; k < n; ++k)
    {
        size_t pivot = k;
        double pivotAbs = std::fabs(A[idx(k, k)]);
        for (size_t r = k + 1; r < n; ++r)
        {
            const double cand = std::fabs(A[idx(r, k)]);
            if (cand > pivotAbs)
            {
                pivotAbs = cand;
                pivot = r;
            }
        }
        if (!(pivotAbs > 1e-18) || !std::isfinite(pivotAbs))
        {
            return false;
        }
        if (pivot != k)
        {
            for (size_t c = k; c < n; ++c)
            {
                std::swap(A[idx(k, c)], A[idx(pivot, c)]);
            }
            std::swap(b[k], b[pivot]);
        }

        const double diag = A[idx(k, k)];
        for (size_t r = k + 1; r < n; ++r)
        {
            const double factor = A[idx(r, k)] / diag;
            A[idx(r, k)] = 0.0;
            if (factor == 0.0)
            {
                continue;
            }
            for (size_t c = k + 1; c < n; ++c)
            {
                A[idx(r, c)] -= factor * A[idx(k, c)];
            }
            b[r] -= factor * b[k];
        }
    }

    for (size_t ii = n; ii-- > 0; )
    {
        double sum = b[ii];
        for (size_t c = ii + 1; c < n; ++c)
        {
            sum -= A[idx(ii, c)] * x[c];
        }
        const double diag = A[idx(ii, ii)];
        if (!(std::fabs(diag) > 1e-18) || !std::isfinite(diag))
        {
            return false;
        }
        x[ii] = sum / diag;
    }
    return true;
}

static double compactWeightWendland(double distance, double supportRadius)
{
    if (!(supportRadius > 0.0))
    {
        return 1.0;
    }
    const double q = distance / supportRadius;
    if (q >= 1.0)
    {
        return 0.0;
    }
    const double oneMinus = 1.0 - q;
    return oneMinus * oneMinus * oneMinus * oneMinus * (4.0 * q + 1.0);
}

static double inverseDistanceFallback(const std::vector<Point3D> &points,
                                      const SpatialHash2D &spatialIndex,
                                      double x,
                                      double y)
{
    const auto neighbors = spatialIndex.gatherNearest(points, x, y, std::min<size_t>(16, points.size()), 1);
    double sumW = 0.0;
    double sumZ = 0.0;
    for (size_t idx : neighbors)
    {
        const double dx = points[idx].x - x;
        const double dy = points[idx].y - y;
        const double w = 1.0 / (1e-12 + dx * dx + dy * dy);
        sumW += w;
        sumZ += w * points[idx].z;
    }
    if (sumW <= 0.0)
    {
        return points.empty() ? 0.0 : points.front().z;
    }
    return sumZ / sumW;
}


struct LocalFrame2D
{
    double x0 = 0.0;
    double y0 = 0.0;
    double e1x = 1.0;
    double e1y = 0.0;
    double e2x = 0.0;
    double e2y = 1.0;
    double s1 = 1.0;
    double s2 = 1.0;
};

static void computeLocalPrincipalFrame(const std::vector<Point3D> &points,
                                       const std::vector<size_t> &neighbors,
                                       double x0,
                                       double y0,
                                       double referenceScale,
                                       LocalFrame2D &frame)
{
    frame = LocalFrame2D{};
    frame.x0 = x0;
    frame.y0 = y0;

    if (neighbors.empty())
    {
        frame.s1 = std::max(referenceScale, 1e-6);
        frame.s2 = std::max(referenceScale, 1e-6);
        return;
    }

    double meanR2 = 0.0;
    for (size_t idx : neighbors)
    {
        const double dx = points[idx].x - x0;
        const double dy = points[idx].y - y0;
        meanR2 += dx * dx + dy * dy;
    }
    meanR2 /= static_cast<double>(neighbors.size());
    const double h = std::max(std::sqrt(std::max(meanR2, 1e-12)), std::max(referenceScale, 1e-6));

    double sxx = 0.0;
    double syy = 0.0;
    double sxy = 0.0;
    double wsum = 0.0;
    for (size_t idx : neighbors)
    {
        const double dx = points[idx].x - x0;
        const double dy = points[idx].y - y0;
        const double r2 = dx * dx + dy * dy;
        const double w = 1.0 / (1.0 + r2 / (h * h));
        sxx += w * dx * dx;
        syy += w * dy * dy;
        sxy += w * dx * dy;
        wsum += w;
    }

    if (!(wsum > 0.0))
    {
        frame.s1 = std::max(referenceScale, 1e-6);
        frame.s2 = std::max(referenceScale, 1e-6);
        return;
    }

    sxx /= wsum;
    syy /= wsum;
    sxy /= wsum;

    const double tr = sxx + syy;
    const double det_term = std::sqrt(std::max((sxx - syy) * (sxx - syy) + 4.0 * sxy * sxy, 0.0));
    const double lambda1 = std::max(0.5 * (tr + det_term), 0.0);
    const double lambda2 = std::max(0.5 * (tr - det_term), 0.0);
    const double theta = 0.5 * std::atan2(2.0 * sxy, sxx - syy);

    frame.e1x = std::cos(theta);
    frame.e1y = std::sin(theta);
    frame.e2x = -frame.e1y;
    frame.e2y = frame.e1x;
    frame.s1 = std::max(std::sqrt(lambda1), 0.75 * std::max(referenceScale, 1e-6));
    frame.s2 = std::max(std::sqrt(lambda2), 0.35 * std::max(referenceScale, 1e-6));
}


static void toLocalUV(const LocalFrame2D &frame,
                      double x,
                      double y,
                      double &u,
                      double &v)
{
    const double dx = x - frame.x0;
    const double dy = y - frame.y0;
    u = (dx * frame.e1x + dy * frame.e1y) / frame.s1;
    v = (dx * frame.e2x + dy * frame.e2y) / frame.s2;
}

static double robustSupportRadiusFromLocalUV(const std::vector<Point3D> &points,
                                             const std::vector<size_t> &neighbors,
                                             const LocalFrame2D &frame,
                                             double minimumRadius)
{
    std::vector<double> radii;
    radii.reserve(neighbors.size());
    for (size_t idx : neighbors)
    {
        double u = 0.0;
        double v = 0.0;
        toLocalUV(frame, points[idx].x, points[idx].y, u, v);
        radii.push_back(std::sqrt(u * u + v * v));
    }
    if (radii.empty())
    {
        return std::max(1.5, minimumRadius);
    }
    std::sort(radii.begin(), radii.end());
    const size_t q80 = std::min(radii.size() - 1U, (radii.size() * 4U) / 5U);
    const size_t q90 = std::min(radii.size() - 1U, (radii.size() * 9U) / 10U);
    return std::max({minimumRadius, 1.5, radii[q80], 0.75 * radii[q90]});
}

static std::vector<size_t> selectAdaptiveNeighborhood(const std::vector<Point3D> &points,
                                                      const std::vector<size_t> &candidates,
                                                      const LocalFrame2D &frame,
                                                      size_t minKeep,
                                                      size_t maxKeep,
                                                      double supportRadius)
{
    struct Entry
    {
        size_t idx = 0;
        double rho = 0.0;
    };

    std::vector<Entry> ranked;
    ranked.reserve(candidates.size());
    for (size_t idx : candidates)
    {
        double u = 0.0;
        double v = 0.0;
        toLocalUV(frame, points[idx].x, points[idx].y, u, v);
        ranked.push_back({idx, std::sqrt(u * u + v * v)});
    }
    std::sort(ranked.begin(), ranked.end(), [](const Entry &a, const Entry &b)
    {
        if (a.rho != b.rho)
        {
            return a.rho < b.rho;
        }
        return a.idx < b.idx;
    });

    std::vector<size_t> selected;
    selected.reserve(std::min(maxKeep, ranked.size()));
    for (const auto &entry : ranked)
    {
        if (entry.rho <= supportRadius || selected.size() < minKeep)
        {
            selected.push_back(entry.idx);
            if (selected.size() >= maxKeep)
            {
                break;
            }
        }
        else
        {
            break;
        }
    }
    if (selected.size() < minKeep)
    {
        selected.clear();
        for (size_t k = 0; k < ranked.size() && k < maxKeep; ++k)
        {
            selected.push_back(ranked[k].idx);
        }
    }
    return selected;
}

static double pointDistanceXY(double x0, double y0, double x1, double y1)
{
    const double dx = x0 - x1;
    const double dy = y0 - y1;
    return std::sqrt(dx * dx + dy * dy);
}

static size_t polynomialBasisCount(int basisOrder)
{
    return (basisOrder >= 3) ? 10U : 6U;
}

static void fillPolynomialBasis(double u,
                                double v,
                                int basisOrder,
                                double *phi)
{
    phi[0] = 1.0;
    phi[1] = u;
    phi[2] = v;
    phi[3] = 0.5 * u * u;
    phi[4] = u * v;
    phi[5] = 0.5 * v * v;
    if (basisOrder >= 3)
    {
        phi[6] = (u * u * u) / 6.0;
        phi[7] = 0.5 * u * u * v;
        phi[8] = 0.5 * u * v * v;
        phi[9] = (v * v * v) / 6.0;
    }
}

static double evaluatePolynomialValueFromCoeffs(const std::vector<double> &coeffs,
                                                int basisOrder,
                                                double u,
                                                double v)
{
    double phi[10] = {};
    fillPolynomialBasis(u, v, basisOrder, phi);
    const size_t basisCount = polynomialBasisCount(basisOrder);
    double value = 0.0;
    for (size_t i = 0; i < basisCount; ++i)
    {
        value += coeffs[i] * phi[i];
    }
    return value;
}

static ValidationMetrics evaluateMLSValidation(const std::vector<Point3D> &points,
                                               const std::vector<size_t> &indices,
                                               const LocalFrame2D &frame,
                                               const std::vector<double> &coeffs,
                                               int basisOrder)
{
    std::vector<double> errors;
    errors.reserve(indices.size());
    for (size_t idx : indices)
    {
        double u = 0.0;
        double v = 0.0;
        toLocalUV(frame, points[idx].x, points[idx].y, u, v);
        const double pred = evaluatePolynomialValueFromCoeffs(coeffs, basisOrder, u, v);
        errors.push_back(pred - points[idx].z);
    }
    return makeValidationMetricsFromErrors(std::move(errors));
}

static ValidationMetrics evaluateLocalMLSPredictiveValidation(const std::vector<Point3D> &points,
                                                              const std::vector<size_t> &neighbors,
                                                              const LocalFrame2D &frame,
                                                              double supportRadius,
                                                              int basisOrder,
                                                              double ridgeFactor,
                                                              size_t maxChecks = 8U)
{
    std::vector<double> errors;
    if (neighbors.size() < polynomialBasisCount(basisOrder) + 2U)
    {
        return ValidationMetrics{};
    }

    const size_t step = std::max<size_t>(1U, neighbors.size() / std::max<size_t>(1U, maxChecks));
    for (size_t k = step / 2U; k < neighbors.size() && errors.size() < maxChecks; k += step)
    {
        const size_t holdoutIdx = neighbors[k];
        std::vector<size_t> fitIndices;
        fitIndices.reserve(neighbors.size() - 1U);
        for (size_t idx : neighbors)
        {
            if (idx != holdoutIdx)
            {
                fitIndices.push_back(idx);
            }
        }
        if (fitIndices.size() < polynomialBasisCount(basisOrder))
        {
            continue;
        }
        std::vector<double> coeffs;
        if (!fitRobustLocalMLS(points, fitIndices, frame, supportRadius, basisOrder, ridgeFactor, coeffs, nullptr))
        {
            continue;
        }
        double u = 0.0;
        double v = 0.0;
        toLocalUV(frame, points[holdoutIdx].x, points[holdoutIdx].y, u, v);
        const double pred = evaluatePolynomialValueFromCoeffs(coeffs, basisOrder, u, v);
        if (std::isfinite(pred))
        {
            errors.push_back(pred - points[holdoutIdx].z);
        }
    }
    return makeValidationMetricsFromErrors(std::move(errors));
}

static bool fitRobustLocalMLS(const std::vector<Point3D> &points,
                              const std::vector<size_t> &fitIndices,
                              const LocalFrame2D &frame,
                              double supportRadius,
                              int basisOrder,
                              double ridgeFactor,
                              std::vector<double> &coeffs,
                              ValidationMetrics *fitMetrics)
{
    const size_t basisCount = polynomialBasisCount(basisOrder);
    coeffs.assign(basisCount, 0.0);
    if (fitIndices.size() < basisCount)
    {
        return false;
    }

    std::vector<double> robustWeights(fitIndices.size(), 1.0);
    ValidationMetrics localFitMetrics;

    for (int iter = 0; iter < 3; ++iter)
    {
        std::vector<double> ATA(basisCount * basisCount, 0.0);
        std::vector<double> ATb(basisCount, 0.0);
        double baseTrace = 0.0;

        for (size_t k = 0; k < fitIndices.size(); ++k)
        {
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
            for (size_t r = 0; r < basisCount; ++r)
            {
                ATb[r] += w * phi[r] * points[idx].z;
                for (size_t c = 0; c < basisCount; ++c)
                {
                    ATA[r * basisCount + c] += w * phi[r] * phi[c];
                }
            }
        }
        for (size_t d = 0; d < basisCount; ++d)
        {
            baseTrace += ATA[d * basisCount + d];
        }
        const double ridge = std::max(1e-12, ridgeFactor * 1e-11 * baseTrace / static_cast<double>(basisCount));
        for (size_t d = 0; d < basisCount; ++d)
        {
            ATA[d * basisCount + d] += ridge;
        }
        if (!solveDenseLinearSystem(ATA, ATb, basisCount, coeffs))
        {
            return false;
        }

        localFitMetrics = evaluateMLSValidation(points, fitIndices, frame, coeffs, basisOrder);
        std::vector<double> residuals;
        residuals.reserve(fitIndices.size());
        for (size_t idx : fitIndices)
        {
            double u = 0.0;
            double v = 0.0;
            toLocalUV(frame, points[idx].x, points[idx].y, u, v);
            const double pred = evaluatePolynomialValueFromCoeffs(coeffs, basisOrder, u, v);
            residuals.push_back(points[idx].z - pred);
        }
        const double med = computeMedianUnchecked(residuals);
        for (double &r : residuals)
        {
            r = std::fabs(r - med);
        }
        const double mad = computeMedianUnchecked(residuals);
        const double sigma = std::max(1.4826 * mad, 1e-8);
        for (size_t k = 0; k < fitIndices.size(); ++k)
        {
            const size_t idx = fitIndices[k];
            double u = 0.0;
            double v = 0.0;
            toLocalUV(frame, points[idx].x, points[idx].y, u, v);
            const double pred = evaluatePolynomialValueFromCoeffs(coeffs, basisOrder, u, v);
            const double t = std::fabs(points[idx].z - pred) / (4.685 * sigma);
            if (t >= 1.0)
            {
                robustWeights[k] = 1e-4;
            }
            else
            {
                const double oneMinus = 1.0 - t * t;
                robustWeights[k] = std::max(oneMinus * oneMinus, 1e-4);
            }
        }
    }

    if (fitMetrics)
    {
        *fitMetrics = localFitMetrics;
    }
    return std::all_of(coeffs.begin(), coeffs.end(), [](double v) { return std::isfinite(v); });
}

static TuningChoice defaultTuningChoice()
{
    TuningChoice t;
    t.neighborhoodTarget = 64U;
    t.supportMultiplier = 1.75;
    t.ridgeFactor = 1.0;
    t.basisOrder = 3;
    return t;
}

static void limitHermiteNodeDerivatives(BicubicNodeData &node,
                                        double gridSpacing,
                                        double localZMin,
                                        double localZMax)
{
    const double zRange = std::max(localZMax - localZMin, 1e-8);
    const double slopeLimit = 3.0 * zRange / std::max(gridSpacing, 1e-6);
    const double twistLimit = 4.0 * zRange / std::max(gridSpacing * gridSpacing, 1e-6);
    node.fx = std::max(-slopeLimit, std::min(slopeLimit, node.fx));
    node.fy = std::max(-slopeLimit, std::min(slopeLimit, node.fy));
    node.fxy = std::max(-twistLimit, std::min(twistLimit, node.fxy));
}

static double minmod2(double a, double b)
{
    if (a * b <= 0.0)
    {
        return 0.0;
    }
    return (std::fabs(a) < std::fabs(b)) ? a : b;
}

static void limitHermiteCellDerivatives(BicubicNodeData &n00,
                                        BicubicNodeData &n10,
                                        BicubicNodeData &n01,
                                        BicubicNodeData &n11,
                                        double hx,
                                        double hy)
{
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

    n00.fx = limFx00; n10.fx = limFx10; n01.fx = limFx01; n11.fx = limFx11;
    n00.fy = limFy00; n01.fy = limFy01; n10.fy = limFy10; n11.fy = limFy11;

    const double crossBase = 0.5 * ((sy1 - sy0) / std::max(hx, 1e-6) + (sx1 - sx0) / std::max(hy, 1e-6));
    const double crossLimit = 2.0 * std::fabs(crossBase) + 1e-8;
    n00.fxy = std::max(-crossLimit, std::min(crossLimit, minmod2(n00.fxy, crossBase)));
    n10.fxy = std::max(-crossLimit, std::min(crossLimit, minmod2(n10.fxy, crossBase)));
    n01.fxy = std::max(-crossLimit, std::min(crossLimit, minmod2(n01.fxy, crossBase)));
    n11.fxy = std::max(-crossLimit, std::min(crossLimit, minmod2(n11.fxy, crossBase)));
}

static bool fitLocalPolynomialNode(const std::vector<Point3D> &points,
                                   const SpatialHash2D &spatialIndex,
                                   double x0,
                                   double y0,
                                   double gridSpacing,
                                   BicubicNodeData &node,
                                   const TuningChoice *tuningChoice = nullptr,
                                   InterpolationDiagnostics *diagnostics = nullptr,
                                   size_t excludeIndex = std::numeric_limits<size_t>::max())
{
    node = BicubicNodeData{};
    if (diagnostics)
    {
        *diagnostics = InterpolationDiagnostics{};
    }
    if (points.empty())
    {
        return false;
    }

    const TuningChoice tuning = tuningChoice ? *tuningChoice : defaultTuningChoice();
    const size_t probeCount = std::min<size_t>(points.size(), std::max<size_t>(tuning.neighborhoodTarget * 3U, 192U));
    const size_t minCount = std::min<size_t>(points.size(), 24U);
    auto candidates = spatialIndex.gatherNearest(points, x0, y0, probeCount, minCount);
    if (excludeIndex != std::numeric_limits<size_t>::max())
    {
        candidates.erase(std::remove(candidates.begin(), candidates.end(), excludeIndex), candidates.end());
    }
    if (candidates.empty())
    {
        return false;
    }

    LocalFrame2D probeFrame;
    computeLocalPrincipalFrame(points, candidates, x0, y0, gridSpacing, probeFrame);
    const double baseSupport = robustSupportRadiusFromLocalUV(points, candidates, probeFrame, 1.6);
    const size_t maxKeep = std::min<size_t>(points.size(), std::max<size_t>(tuning.neighborhoodTarget + 24U, 96U));
    auto neighbors = selectAdaptiveNeighborhood(points,
                                                candidates,
                                                probeFrame,
                                                minCount,
                                                maxKeep,
                                                std::max(1.5, baseSupport * tuning.supportMultiplier));
    if (neighbors.size() < std::max<size_t>(12U, polynomialBasisCount(tuning.basisOrder) + 4U))
    {
        node.z = inverseDistanceFallback(points, spatialIndex, x0, y0);
        node.localConfidence = 0.1;
        node.fallbackUsed = true;
        node.valid = true;
        if (diagnostics)
        {
            diagnostics->fallbackUsed = true;
            diagnostics->confidence = 0.1;
            diagnostics->neighborCount = neighbors.size();
            diagnostics->summary = "MLS fallback to IDW due to insufficient neighborhood.";
        }
        return true;
    }

    LocalFrame2D frame;
    computeLocalPrincipalFrame(points, neighbors, x0, y0, gridSpacing, frame);
    const double supportRadius = std::max(1.5, robustSupportRadiusFromLocalUV(points, neighbors, frame, 1.5) * tuning.supportMultiplier);

    std::vector<size_t> fitIndices = neighbors;

    std::vector<double> coeffs;
    ValidationMetrics fitMetrics;
    if (!fitRobustLocalMLS(points, fitIndices, frame, supportRadius, tuning.basisOrder, tuning.ridgeFactor, coeffs, &fitMetrics))
    {
        node.z = inverseDistanceFallback(points, spatialIndex, x0, y0);
        node.localConfidence = 0.1;
        node.fallbackUsed = true;
        node.valid = true;
        if (diagnostics)
        {
            diagnostics->fallbackUsed = true;
            diagnostics->confidence = 0.1;
            diagnostics->neighborCount = neighbors.size();
            diagnostics->summary = "MLS fallback to IDW because the robust solve failed.";
        }
        return true;
    }

    ValidationMetrics validationMetrics = evaluateLocalMLSPredictiveValidation(points, neighbors, frame, supportRadius, tuning.basisOrder, tuning.ridgeFactor);
    if (!validationMetrics.valid)
    {
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
    node.fxy = fuu * du_dx * du_dy
             + fuv * (du_dx * dv_dy + dv_dx * du_dy)
             + fvv * dv_dx * dv_dy;
    node.localResidual = validationMetrics.valid ? validationMetrics.rmse : fitMetrics.rmse;
    node.localConfidence = validationMetrics.valid ? calibratedConfidenceFromObservedError(1.0 / (1.0 + 4.0 * validationScore(validationMetrics)), validationMetrics.rmse, validationMetrics) : 0.25;
    node.supportRadius = supportRadius;
    node.neighborCount = neighbors.size();
    node.valid = std::isfinite(node.z) && std::isfinite(node.fx) && std::isfinite(node.fy) && std::isfinite(node.fxy);

    double localZMin = std::numeric_limits<double>::infinity();
    double localZMax = -std::numeric_limits<double>::infinity();
    for (size_t idx : neighbors)
    {
        localZMin = std::min(localZMin, points[idx].z);
        localZMax = std::max(localZMax, points[idx].z);
    }
    if (std::isfinite(localZMin) && std::isfinite(localZMax))
    {
        limitHermiteNodeDerivatives(node, gridSpacing, localZMin, localZMax);
    }

    if (!node.valid)
    {
        node.z = inverseDistanceFallback(points, spatialIndex, x0, y0);
        node.fx = 0.0;
        node.fy = 0.0;
        node.fxy = 0.0;
        node.localConfidence = 0.1;
        node.fallbackUsed = true;
        node.valid = true;
    }

    if (diagnostics)
    {
        diagnostics->confidence = node.localConfidence;
        diagnostics->residual = node.localResidual;
        diagnostics->neighborCount = node.neighborCount;
        diagnostics->patchCount = 1U;
        diagnostics->fallbackUsed = node.fallbackUsed;
        diagnostics->summary = validationMetrics.valid ? "robust anisotropic MLS" : "robust anisotropic MLS (fit-metric only)";
    }
    return true;
}

static double evaluateBicubicHermite(const BicubicNodeData &n00,
                                     const BicubicNodeData &n10,
                                     const BicubicNodeData &n01,
                                     const BicubicNodeData &n11,
                                     double u,
                                     double v,
                                     double hx,
                                     double hy)
{
    const double h00u = 2.0 * u * u * u - 3.0 * u * u + 1.0;
    const double h10u = u * u * u - 2.0 * u * u + u;
    const double h01u = -2.0 * u * u * u + 3.0 * u * u;
    const double h11u = u * u * u - u * u;

    const double h00v = 2.0 * v * v * v - 3.0 * v * v + 1.0;
    const double h10v = v * v * v - 2.0 * v * v + v;
    const double h01v = -2.0 * v * v * v + 3.0 * v * v;
    const double h11v = v * v * v - v * v;

    return
        n00.z * h00u * h00v +
        n10.z * h01u * h00v +
        n01.z * h00u * h01v +
        n11.z * h01u * h01v +
        hx * (n00.fx * h10u * h00v + n10.fx * h11u * h00v + n01.fx * h10u * h01v + n11.fx * h11u * h01v) +
        hy * (n00.fy * h00u * h10v + n10.fy * h01u * h10v + n01.fy * h00u * h11v + n11.fy * h01u * h11v) +
        hx * hy * (n00.fxy * h10u * h10v + n10.fxy * h11u * h10v + n01.fxy * h10u * h11v + n11.fxy * h11u * h11v);
}

static bool buildTPSModel(std::vector<Point3D> pts,
                          TPSModel &model,
                          std::string *diagnostics)
{
    model = TPSModel{};
    if (pts.empty())
    {
        if (diagnostics)
        {
            *diagnostics = "TPS solver received no control points.";
        }
        return false;
    }

    const size_t pointCount = pts.size();
    const double spacing = averagePlanarSpacing(pts);
    model.controlPoints = std::move(pts);
    model.spatialIndex.build(model.controlPoints, std::max(3.0 * spacing, 1e-6));
    model.neighborhoodSize = std::min<size_t>(std::max<size_t>(64U, static_cast<size_t>(3.0 * std::sqrt(static_cast<double>(pointCount)))), std::min<size_t>(160U, pointCount));
    model.suggestedLambda = 1e-9;
    model.hull = buildConvexHull2D(model.controlPoints);

    std::vector<double> learnedLambdas;
    std::vector<double> errors;
    const auto holdout = makeDeterministicHoldoutIndices(model.controlPoints.size(), 32U);
    if (holdout.size() >= 8U)
    {
        std::vector<Point3D> training;
        training.reserve(model.controlPoints.size() - holdout.size());
        std::vector<unsigned char> isHoldout(model.controlPoints.size(), 0U);
        for (size_t idx : holdout) isHoldout[idx] = 1U;
        for (size_t i = 0; i < model.controlPoints.size(); ++i)
        {
            if (isHoldout[i] == 0U)
            {
                training.push_back(model.controlPoints[i]);
            }
        }
        TPSModel trainingModel;
        trainingModel.controlPoints = training;
        trainingModel.spatialIndex.build(trainingModel.controlPoints, std::max(3.0 * averagePlanarSpacing(trainingModel.controlPoints), 1e-6));
        trainingModel.neighborhoodSize = std::min<size_t>(std::max<size_t>(64U, static_cast<size_t>(3.0 * std::sqrt(static_cast<double>(trainingModel.controlPoints.size())))), std::min<size_t>(160U, trainingModel.controlPoints.size()));
        trainingModel.suggestedLambda = model.suggestedLambda;
        trainingModel.hull = buildConvexHull2D(trainingModel.controlPoints);
        trainingModel.valid = true;
        for (size_t idx : holdout)
        {
            double pred = 0.0;
            InterpolationDiagnostics d;
            if (solveLocalTPSAtPoint(model.controlPoints[idx].x, model.controlPoints[idx].y, trainingModel, pred, nullptr, &d) && std::isfinite(pred))
            {
                errors.push_back(pred - model.controlPoints[idx].z);
                if (d.lambda > 0.0 && std::isfinite(d.lambda))
                {
                    learnedLambdas.push_back(d.lambda);
                }
            }
        }
    }
    model.validation = makeValidationMetricsFromErrors(std::move(errors));
    if (!learnedLambdas.empty())
    {
        std::sort(learnedLambdas.begin(), learnedLambdas.end());
        model.suggestedLambda = learnedLambdas[learnedLambdas.size() / 2U];
    }
    model.valid = true;

    if (diagnostics)
    {
        std::ostringstream oss;
        oss << "High-accuracy local TPS controls=" << model.controlPoints.size()
            << ", target neighborhood=" << model.neighborhoodSize
            << ", support-mask cells=" << model.hull.occupied.size()
            << ", learned lambda=" << model.suggestedLambda;
        if (model.validation.valid)
        {
            oss << ", holdout RMSE=" << model.validation.rmse;
        }
        *diagnostics = oss.str();
    }
    return true;
}

static bool solveLocalTPSAtPoint(double x,
                                 double y,
                                 const TPSModel &model,
                                 double &z,
                                 std::string *diagnostics,
                                 InterpolationDiagnostics *interpDiagnostics = nullptr)
{
    z = 0.0;
    if (interpDiagnostics)
    {
        *interpDiagnostics = InterpolationDiagnostics{};
    }
    if (!model.valid || model.controlPoints.empty())
    {
        return false;
    }

    auto phi = [](double r2) -> double
    {
        if (r2 <= 1e-24)
        {
            return 0.0;
        }
        return r2 * std::log(r2);
    };

    auto solvePatchSystem = [&](const std::vector<size_t> &neighbors,
                                const LocalFrame2D &frame,
                                const std::vector<size_t> &fitSelection,
                                double lambda,
                                std::vector<double> &solution,
                                std::vector<double> &un,
                                std::vector<double> &vn) -> bool
    {
        const size_t m = fitSelection.size();
        if (m < 8U)
        {
            return false;
        }
        const size_t N = m + 3U;
        un.assign(m, 0.0);
        vn.assign(m, 0.0);
        std::vector<double> rhs(N, 0.0);
        for (size_t i = 0; i < m; ++i)
        {
            const Point3D &p = model.controlPoints[neighbors[fitSelection[i]]];
            toLocalUV(frame, p.x, p.y, un[i], vn[i]);
            rhs[i] = p.z;
        }

        std::vector<double> A(N * N, 0.0);
        for (size_t i = 0; i < m; ++i)
        {
            A[i * N + i] = lambda;
            A[i * N + m] = 1.0;
            A[i * N + (m + 1U)] = un[i];
            A[i * N + (m + 2U)] = vn[i];
            A[m * N + i] = 1.0;
            A[(m + 1U) * N + i] = un[i];
            A[(m + 2U) * N + i] = vn[i];
        }
        for (size_t i = 0; i < m; ++i)
        {
            for (size_t j = i + 1U; j < m; ++j)
            {
                const double du = un[i] - un[j];
                const double dv = vn[i] - vn[j];
                const double value = phi(du * du + dv * dv);
                A[i * N + j] = value;
                A[j * N + i] = value;
            }
        }
        return solveDenseLinearSystem(A, rhs, N, solution)
            && std::all_of(solution.begin(), solution.end(), [](double v){ return std::isfinite(v); });
    };

    auto evaluateTPSPatch = [&](const std::vector<double> &solution,
                                const std::vector<double> &un,
                                const std::vector<double> &vn,
                                const LocalFrame2D &frame,
                                double qx,
                                double qy) -> double
    {
        const size_t m = un.size();
        const size_t N = m + 3U;
        if (solution.size() < N)
        {
            return std::numeric_limits<double>::quiet_NaN();
        }
        double uq = 0.0;
        double vq = 0.0;
        toLocalUV(frame, qx, qy, uq, vq);
        double value = solution[m] + solution[m + 1U] * uq + solution[m + 2U] * vq;
        for (size_t i = 0; i < m; ++i)
        {
            const double du = un[i] - uq;
            const double dv = vn[i] - vq;
            value += solution[i] * phi(du * du + dv * dv);
        }
        return value;
    };

    auto solveLocalPatch = [&](double anchorX,
                               double anchorY,
                               double &zPatch,
                               double &patchResidual,
                               double &patchSupport,
                               double &chosenLambda,
                               size_t &usedNeighbors) -> bool
    {
        const size_t probeCount = std::min<size_t>(model.controlPoints.size(), std::max<size_t>(model.neighborhoodSize * 2U, 96U));
        const size_t minCount = std::min<size_t>(model.controlPoints.size(), 24U);
        auto candidates = model.spatialIndex.gatherNearest(model.controlPoints, anchorX, anchorY, probeCount, minCount);
        if (candidates.size() < 8U)
        {
            return false;
        }

        LocalFrame2D probeFrame;
        computeLocalPrincipalFrame(model.controlPoints, candidates, anchorX, anchorY, model.spatialIndex.cellSize, probeFrame);
        const double probeSupport = robustSupportRadiusFromLocalUV(model.controlPoints, candidates, probeFrame, 1.75);
        const size_t maxKeep = std::min<size_t>(model.controlPoints.size(), model.neighborhoodSize);
        auto neighbors = selectAdaptiveNeighborhood(model.controlPoints, candidates, probeFrame, minCount, maxKeep, probeSupport);
        if (neighbors.size() < 8U)
        {
            return false;
        }

        LocalFrame2D frame;
        computeLocalPrincipalFrame(model.controlPoints, neighbors, anchorX, anchorY, model.spatialIndex.cellSize, frame);
        patchSupport = robustSupportRadiusFromLocalUV(model.controlPoints, neighbors, frame, 1.5);
        usedNeighbors = neighbors.size();

        std::vector<size_t> fitSelection;
        std::vector<size_t> validationSelection;
        for (size_t k = 0; k < neighbors.size(); ++k)
        {
            if (neighbors.size() >= 20U && (k % 5U) == 0U)
            {
                validationSelection.push_back(k);
            }
            else
            {
                fitSelection.push_back(k);
            }
        }
        if (fitSelection.size() < 8U)
        {
            fitSelection.clear();
            for (size_t k = 0; k < neighbors.size(); ++k)
            {
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
        for (double lambda : {1e-12, 1e-11, 1e-10, 1e-9, 1e-8, 1e-6})
        {
            std::vector<double> xvec;
            std::vector<double> un;
            std::vector<double> vn;
            if (!solvePatchSystem(neighbors, frame, fitSelection, lambda, xvec, un, vn))
            {
                continue;
            }

            ValidationMetrics metrics;
            if (!validationSelection.empty())
            {
                std::vector<double> errors;
                errors.reserve(validationSelection.size());
                for (size_t localIdx : validationSelection)
                {
                    const Point3D &vp = model.controlPoints[neighbors[localIdx]];
                    const double pred = evaluateTPSPatch(xvec, un, vn, frame, vp.x, vp.y);
                    if (std::isfinite(pred))
                    {
                        errors.push_back(pred - vp.z);
                    }
                }
                metrics = makeValidationMetricsFromErrors(std::move(errors));
            }
            if (!metrics.valid)
            {
                std::vector<double> errors;
                errors.reserve(fitSelection.size());
                for (size_t localIdx : fitSelection)
                {
                    const Point3D &fp = model.controlPoints[neighbors[localIdx]];
                    const double pred = evaluateTPSPatch(xvec, un, vn, frame, fp.x, fp.y);
                    if (std::isfinite(pred))
                    {
                        errors.push_back(pred - fp.z);
                    }
                }
                metrics = makeValidationMetricsFromErrors(std::move(errors));
            }
            const double score = validationScore(metrics);
            if (score < bestScore)
            {
                bestScore = score;
                bestMetrics = metrics;
                bestX = std::move(xvec);
                bestUn = std::move(un);
                bestVn = std::move(vn);
                chosenLambda = lambda;
                solved = true;
            }
        }

        if (!solved)
        {
            return false;
        }

        zPatch = evaluateTPSPatch(bestX, bestUn, bestVn, frame, x, y);
        patchResidual = bestMetrics.valid ? bestMetrics.rmse : 1.0;
        return std::isfinite(zPatch);
    };

    const size_t anchorProbe = std::min<size_t>(model.controlPoints.size(), 8U);
    const auto anchorCandidates = model.spatialIndex.gatherNearest(model.controlPoints, x, y, anchorProbe, 1U);
    if (anchorCandidates.empty())
    {
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
    for (size_t a = 0; a < blendCount; ++a)
    {
        const Point3D &anchor = model.controlPoints[anchorCandidates[a]];
        double zPatch = 0.0;
        double residual = std::numeric_limits<double>::infinity();
        double supportRadius = 1.5;
        double lambda = model.suggestedLambda;
        size_t usedNeighbors = 0U;
        if (!solveLocalPatch(anchor.x, anchor.y, zPatch, residual, supportRadius, lambda, usedNeighbors))
        {
            continue;
        }

        const double d = pointDistanceXY(x, y, anchor.x, anchor.y);
        const double supportXY = std::max(1.0, supportRadius * std::max(model.spatialIndex.cellSize, 1e-6));
        const double w = std::max(compactWeightWendland(d, supportXY * 1.001) / (1e-8 + residual), 1e-12);
        weightedZ += w * zPatch;
        weightSum += w;
        ++successes;
        meanNeighbors += usedNeighbors;

        if (residual < bestSingleResidual)
        {
            bestSingleResidual = residual;
            bestSingleZ = zPatch;
            bestLambda = lambda;
        }
    }

    bool fallbackUsed = false;
    if (weightSum > 0.0)
    {
        z = weightedZ / weightSum;
    }
    else if (std::isfinite(bestSingleZ))
    {
        z = bestSingleZ;
    }
    else
    {
        z = inverseDistanceFallback(model.controlPoints, model.spatialIndex, x, y);
        fallbackUsed = true;
    }

    const bool outsideHull = model.hull.valid ? !pointInConvexHull2D(model.hull, x, y) : false;
    const double outsideDistance = outsideHull ? distanceToConvexHull2D(model.hull, x, y) : 0.0;
    const double boundaryPenalty = outsideHull ? (1.0 / (1.0 + outsideDistance / std::max(model.spatialIndex.cellSize, 1e-6))) : 1.0;
    z *= 1.0;

    if (interpDiagnostics)
    {
        interpDiagnostics->confidence = calibratedConfidenceFromObservedError(boundaryPenalty / (1.0 + 4.0 * std::max(bestSingleResidual, 0.0)), std::max(bestSingleResidual, 0.0), model.validation);
        interpDiagnostics->residual = std::isfinite(bestSingleResidual) ? bestSingleResidual : 1.0;
        interpDiagnostics->lambda = bestLambda;
        interpDiagnostics->boundaryPenalty = boundaryPenalty;
        interpDiagnostics->outsideDistance = outsideDistance;
        interpDiagnostics->neighborCount = (successes > 0U) ? (meanNeighbors / successes) : 0U;
        interpDiagnostics->patchCount = successes;
        interpDiagnostics->fallbackUsed = fallbackUsed;
        interpDiagnostics->outsideHull = outsideHull;
        std::ostringstream oss;
        oss << "TPS patches=" << successes << ", lambda=" << bestLambda
            << ", residual=" << bestSingleResidual
            << (outsideHull ? ", outside-hull" : "");
        interpDiagnostics->summary = oss.str();
    }
    if (diagnostics)
    {
        if (interpDiagnostics)
        {
            *diagnostics = interpDiagnostics->summary;
        }
        else
        {
            std::ostringstream oss;
            oss << "TPS patches=" << successes << ", residual=" << bestSingleResidual;
            *diagnostics = oss.str();
        }
    }
    return std::isfinite(z);
}

/****************************************************************************
 * 8) Bicubic interpolation via local polynomial node estimation
 ****************************************************************************//****************************************************************************
 * 8) Bicubic interpolation via local polynomial node estimation
 ****************************************************************************/
[[maybe_unused]]
static std::vector<Point3D> generateGridPointsBicubic(const std::vector<Point3D> &points,
                                                      double gridSpacing)
{
    if (points.empty() || !(gridSpacing > 0.0))
    {
        return {};
    }

    double dataXMin, dataXMax, dataYMin, dataYMax;
    computeBoundingBox(points, dataXMin, dataXMax, dataYMin, dataYMax);

    const double margin = 1.5 * gridSpacing;
    const double xMin = dataXMin - margin;
    const double xMax = dataXMax + margin;
    const double yMin = dataYMin - margin;
    const double yMax = dataYMax + margin;

    const size_t nx = static_cast<size_t>(std::ceil(std::max(xMax - xMin, 1.0) / gridSpacing)) + 1;
    const size_t ny = static_cast<size_t>(std::ceil(std::max(yMax - yMin, 1.0) / gridSpacing)) + 1;
    if (nx < 2 || ny < 2)
    {
        return {};
    }

    SpatialHash2D spatialIndex;
    const double spacing = averagePlanarSpacing(points);
    spatialIndex.build(points, std::max(2.5 * spacing, 1.5 * gridSpacing));

    std::vector<BicubicNodeData> nodes(nx * ny);
#ifdef _OPENMP
#pragma omp parallel for collapse(2)
#endif
    for (size_t i = 0; i < nx; ++i)
    {
        for (size_t j = 0; j < ny; ++j)
        {
            const double gx = xMin + static_cast<double>(i) * gridSpacing;
            const double gy = yMin + static_cast<double>(j) * gridSpacing;
            BicubicNodeData node;
            fitLocalPolynomialNode(points, spatialIndex, gx, gy, gridSpacing, node);
            nodes[i * ny + j] = node;
        }
    }

    std::vector<Point3D> gridPoints(nx * ny);
#ifdef _OPENMP
#pragma omp parallel for collapse(2)
#endif
    for (size_t i = 0; i < nx; ++i)
    {
        for (size_t j = 0; j < ny; ++j)
        {
            const double gx = xMin + static_cast<double>(i) * gridSpacing;
            const double gy = yMin + static_cast<double>(j) * gridSpacing;

            size_t ci = (i == nx - 1) ? nx - 2 : i;
            size_t cj = (j == ny - 1) ? ny - 2 : j;
            const double u = (i == nx - 1) ? 1.0 : 0.0;
            const double v = (j == ny - 1) ? 1.0 : 0.0;

            const BicubicNodeData &n00 = nodes[ci * ny + cj];
            const BicubicNodeData &n10 = nodes[(ci + 1) * ny + cj];
            const BicubicNodeData &n01 = nodes[ci * ny + (cj + 1)];
            const BicubicNodeData &n11 = nodes[(ci + 1) * ny + (cj + 1)];

            double gz = 0.0;
            if (n00.valid && n10.valid && n01.valid && n11.valid)
            {
                gz = evaluateBicubicHermite(n00, n10, n01, n11, u, v, gridSpacing, gridSpacing);
            }
            else
            {
                gz = inverseDistanceFallback(points, spatialIndex, gx, gy);
            }
            gridPoints[i * ny + j] = {gx, gy, gz};
        }
    }

    return gridPoints;
}

/****************************************************************************
 * 9) Local TPS interpolation with adaptive neighborhoods
 ****************************************************************************/
[[maybe_unused]]
static std::vector<Point3D> generateGridPointsTPS(const std::vector<Point3D> &tpsPoints,
                                                  double gridSpacing,
                                                  std::string *errorMessage)
{
    if (tpsPoints.empty() || !(gridSpacing > 0.0))
    {
        if (errorMessage)
        {
            *errorMessage = "TPS received invalid input.";
        }
        return {};
    }

    TPSModel model;
    std::string diagnostics;
    if (!buildTPSModel(tpsPoints, model, &diagnostics))
    {
        if (errorMessage)
        {
            *errorMessage = diagnostics;
        }
        return {};
    }

    double xMin, xMax, yMin, yMax;
    computeBoundingBox(tpsPoints, xMin, xMax, yMin, yMax);
    const double margin = 1.5 * gridSpacing;
    xMin -= margin;
    xMax += margin;
    yMin -= margin;
    yMax += margin;

    const size_t nx = static_cast<size_t>(std::ceil(std::max(xMax - xMin, 1.0) / gridSpacing)) + 1;
    const size_t ny = static_cast<size_t>(std::ceil(std::max(yMax - yMin, 1.0) / gridSpacing)) + 1;
    std::vector<Point3D> gridPoints(nx * ny);

#ifdef _OPENMP
#pragma omp parallel for collapse(2)
#endif
    for (size_t i = 0; i < nx; ++i)
    {
        for (size_t j = 0; j < ny; ++j)
        {
            const double gx = xMin + static_cast<double>(i) * gridSpacing;
            const double gy = yMin + static_cast<double>(j) * gridSpacing;
            double gz = 0.0;
            if (!solveLocalTPSAtPoint(gx, gy, model, gz, nullptr))
            {
                gz = inverseDistanceFallback(model.controlPoints, model.spatialIndex, gx, gy);
            }
            gridPoints[i * ny + j] = {gx, gy, gz};
        }
    }

    if (errorMessage)
    {
        *errorMessage = diagnostics;
    }
    return gridPoints;
}

static std::string formatPointXYZLine(const Point3D &p, int precision)
{
    std::ostringstream oss;
    oss << std::fixed << std::setprecision(precision)
        << p.x << ' ' << p.y << ' ' << p.z << "\n";
    return oss.str();
}


static bool checkedMultiplySizeT(size_t a, size_t b, size_t &out)
{
    if (a == 0 || b == 0)
    {
        out = 0;
        return true;
    }
    if (a > (std::numeric_limits<size_t>::max() / b))
    {
        out = 0;
        return false;
    }
    out = a * b;
    return true;
}

static bool buildGridDefinition(const std::vector<Point3D> &points,
                                double gridSpacing,
                                GridDefinition &grid,
                                std::string *errorMessage)
{
    grid = GridDefinition{};
    if (points.empty())
    {
        if (errorMessage)
        {
            *errorMessage = "Error: no points available to define interpolation grid.";
        }
        return false;
    }
    if (!(gridSpacing > 0.0) || !std::isfinite(gridSpacing))
    {
        if (errorMessage)
        {
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
    if (grid.nx < 2U || grid.ny < 2U)
    {
        if (errorMessage)
        {
            *errorMessage = "Error: interpolation grid is degenerate.";
        }
        return false;
    }
    if (!checkedMultiplySizeT(grid.nx, grid.ny, grid.totalPoints))
    {
        if (errorMessage)
        {
            *errorMessage = "Error: interpolation grid dimensions overflow size_t.";
        }
        return false;
    }

    grid.valid = true;
    return true;
}

struct BicubicStreamingContext
{
    const std::vector<Point3D> *points = nullptr;
    SpatialHash2D spatialIndex;
    ConvexHull2D hull;
    GridDefinition grid;
    TuningChoice tuning;
    ValidationMetrics validation;
};

struct ClampedMLSStreamingContext
{
    const std::vector<Point3D> *points = nullptr;
    SpatialHash2D spatialIndex;
    ConvexHull2D hull;
    GridDefinition grid;
    TuningChoice tuning;
    ValidationMetrics validation;
};

struct TPSStreamingContext
{
    TPSModel model;
    GridDefinition grid;
};

static ValidationMetrics evaluateGlobalMLSHoldout(const std::vector<Point3D> &points,
                                                  const SpatialHash2D &spatialIndex,
                                                  double gridSpacing,
                                                  const TuningChoice &tuning,
                                                  bool clampedMode)
{
    std::vector<double> errors;
    const auto holdout = makeDeterministicHoldoutIndices(points.size(), 48U);
    errors.reserve(holdout.size());
    for (size_t idx : holdout)
    {
        BicubicNodeData node;
        InterpolationDiagnostics diag;
        if (!fitLocalPolynomialNode(points, spatialIndex, points[idx].x, points[idx].y, gridSpacing, node, &tuning, &diag, idx))
        {
            continue;
        }
        double pred = node.z;
        if (clampedMode)
        {
            auto neighbors = spatialIndex.gatherNearest(points, points[idx].x, points[idx].y, std::min<size_t>(points.size(), 32U), 1U);
            double zMin = std::numeric_limits<double>::infinity();
            double zMax = -std::numeric_limits<double>::infinity();
            for (size_t j : neighbors)
            {
                if (j == idx) continue;
                zMin = std::min(zMin, points[j].z);
                zMax = std::max(zMax, points[j].z);
            }
            if (std::isfinite(zMin) && std::isfinite(zMax))
            {
                pred = std::max(zMin, std::min(zMax, pred));
            }
        }
        if (std::isfinite(pred))
        {
            errors.push_back(pred - points[idx].z);
        }
    }
    return makeValidationMetricsFromErrors(std::move(errors));
}

static TuningChoice tuneBicubicChoice(const std::vector<Point3D> &points,
                                      const SpatialHash2D &spatialIndex,
                                      double gridSpacing)
{
    TuningChoice best = defaultTuningChoice();
    if (points.empty())
    {
        return best;
    }

    std::vector<size_t> neighborhoodOptions = {32U, 48U, 64U, 96U};
    std::vector<double> supportOptions = {1.5, 1.75, 2.25};
    std::vector<double> ridgeOptions = {1.0, 10.0};
    std::vector<int> orderOptions = {2, 3};

    double bestScore = std::numeric_limits<double>::infinity();
    for (size_t neigh : neighborhoodOptions)
    {
        neigh = std::min(neigh, points.size());
        for (double support : supportOptions)
        {
            for (double ridge : ridgeOptions)
            {
                for (int order : orderOptions)
                {
                    if (order >= 3 && points.size() < 18U)
                    {
                        continue;
                    }
                    TuningChoice candidate;
                    candidate.neighborhoodTarget = std::max<size_t>(neigh, 24U);
                    candidate.supportMultiplier = support;
                    candidate.ridgeFactor = ridge;
                    candidate.basisOrder = order;
                    candidate.metrics = evaluateGlobalMLSHoldout(points, spatialIndex, gridSpacing, candidate, false);
                    if (!candidate.metrics.valid)
                    {
                        continue;
                    }
                    const double score = validationScore(candidate.metrics);
                    if (score < bestScore)
                    {
                        bestScore = score;
                        best = candidate;
                    }
                }
            }
        }
    }
    return best;
}

static bool buildBicubicStreamingContext(const std::vector<Point3D> &points,
                                         double gridSpacing,
                                         BicubicStreamingContext &ctx,
                                         std::string *errorMessage)
{
    ctx = BicubicStreamingContext{};
    if (!buildGridDefinition(points, gridSpacing, ctx.grid, errorMessage))
    {
        return false;
    }
    const double spacing = averagePlanarSpacing(points);
    ctx.points = &points;
    ctx.spatialIndex.build(points, std::max(2.5 * spacing, 1.5 * gridSpacing));
    ctx.hull = buildConvexHull2D(points);
    ctx.tuning = tuneBicubicChoice(points, ctx.spatialIndex, gridSpacing);
    ctx.validation = evaluateGlobalMLSHoldout(points, ctx.spatialIndex, gridSpacing, ctx.tuning, false);
    return true;
}

static bool buildClampedMLSStreamingContext(const std::vector<Point3D> &points,
                                          double gridSpacing,
                                          ClampedMLSStreamingContext &ctx,
                                          std::string *errorMessage)
{
    ctx = ClampedMLSStreamingContext{};
    if (!buildGridDefinition(points, gridSpacing, ctx.grid, errorMessage))
    {
        return false;
    }
    const double spacing = averagePlanarSpacing(points);
    ctx.points = &points;
    ctx.spatialIndex.build(points, std::max(2.5 * spacing, 1.5 * gridSpacing));
    ctx.hull = buildConvexHull2D(points);
    ctx.tuning = tuneBicubicChoice(points, ctx.spatialIndex, gridSpacing);
    ctx.validation = evaluateGlobalMLSHoldout(points, ctx.spatialIndex, gridSpacing, ctx.tuning, true);
    return true;
}

static bool buildTPSStreamingContext(std::vector<Point3D> controlPoints,
                                     const GridDefinition &extentGrid,
                                     TPSStreamingContext &ctx,
                                     std::string *errorMessage)
{
    ctx = TPSStreamingContext{};
    ctx.grid = extentGrid;
    if (!ctx.grid.valid)
    {
        if (errorMessage)
        {
            *errorMessage = "Error: invalid TPS grid definition.";
        }
        return false;
    }
    std::string diagnostics;
    if (!buildTPSModel(std::move(controlPoints), ctx.model, &diagnostics))
    {
        if (errorMessage)
        {
            *errorMessage = diagnostics.empty() ? "Error: failed to build TPS model." : diagnostics;
        }
        return false;
    }
    return true;
}

static bool evaluateBicubicStreamPointDetailed(const BicubicStreamingContext &ctx,
                                               size_t i,
                                               size_t j,
                                               Point3D &outPoint,
                                               InterpolationDiagnostics *diagnostics,
                                               std::string *errorMessage)
{
    if (diagnostics)
    {
        *diagnostics = InterpolationDiagnostics{};
    }
    if (!ctx.points || !ctx.grid.valid || ctx.grid.nx < 2U || ctx.grid.ny < 2U)
    {
        if (errorMessage)
        {
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
    const bool ok00 = fitLocalPolynomialNode(*ctx.points, ctx.spatialIndex, xCell0,       yCell0,       gs, n00, &ctx.tuning, &d00);
    const bool ok10 = fitLocalPolynomialNode(*ctx.points, ctx.spatialIndex, xCell0 + gs,  yCell0,       gs, n10, &ctx.tuning, &d10);
    const bool ok01 = fitLocalPolynomialNode(*ctx.points, ctx.spatialIndex, xCell0,       yCell0 + gs,  gs, n01, &ctx.tuning, &d01);
    const bool ok11 = fitLocalPolynomialNode(*ctx.points, ctx.spatialIndex, xCell0 + gs,  yCell0 + gs,  gs, n11, &ctx.tuning, &d11);

    double gz = 0.0;
    bool fallback = false;
    if (ok00 && ok10 && ok01 && ok11 && n00.valid && n10.valid && n01.valid && n11.valid)
    {
        limitHermiteCellDerivatives(n00, n10, n01, n11, gs, gs);
        gz = evaluateBicubicHermite(n00, n10, n01, n11, u, v, gs, gs);
        const double zMin = std::min(std::min(n00.z, n10.z), std::min(n01.z, n11.z));
        const double zMax = std::max(std::max(n00.z, n10.z), std::max(n01.z, n11.z));
        const double pad = 0.05 * std::max(zMax - zMin, 1e-8);
        gz = std::max(zMin - pad, std::min(zMax + pad, gz));
    }
    else
    {
        gz = inverseDistanceFallback(*ctx.points, ctx.spatialIndex, gx, gy);
        fallback = true;
    }

    if (!std::isfinite(gz))
    {
        if (errorMessage)
        {
            *errorMessage = "Error: bicubic Hermite evaluation produced a non-finite value.";
        }
        return false;
    }

    const bool outsideHull = ctx.hull.valid ? !pointInConvexHull2D(ctx.hull, gx, gy) : false;
    const double outsideDistance = outsideHull ? distanceToConvexHull2D(ctx.hull, gx, gy) : 0.0;
    const double boundaryPenalty = outsideHull ? (1.0 / (1.0 + outsideDistance / std::max(gs, 1e-6))) : 1.0;

    outPoint = {gx, gy, gz};
    if (diagnostics)
    {
        diagnostics->residual = 0.25 * (d00.residual + d10.residual + d01.residual + d11.residual);
        diagnostics->confidence = calibratedConfidenceFromObservedError(0.25 * (d00.confidence + d10.confidence + d01.confidence + d11.confidence) * boundaryPenalty, diagnostics->residual, ctx.validation);
        diagnostics->neighborCount = (d00.neighborCount + d10.neighborCount + d01.neighborCount + d11.neighborCount) / 4U;
        diagnostics->patchCount = 4U;
        diagnostics->fallbackUsed = fallback;
        diagnostics->outsideHull = outsideHull;
        diagnostics->outsideDistance = outsideDistance;
        diagnostics->boundaryPenalty = boundaryPenalty;
        diagnostics->summary = fallback ? "bicubic fallback" : "bicubic Hermite with robust MLS";
    }
    return std::isfinite(outPoint.x) && std::isfinite(outPoint.y) && std::isfinite(outPoint.z);
}

static bool evaluateBicubicStreamPoint(const BicubicStreamingContext &ctx,
                                       size_t i,
                                       size_t j,
                                       Point3D &outPoint,
                                       std::string *errorMessage)
{
    return evaluateBicubicStreamPointDetailed(ctx, i, j, outPoint, nullptr, errorMessage);
}

static bool evaluateClampedMLSStreamPointDetailed(const ClampedMLSStreamingContext &ctx,
                                                size_t i,
                                                size_t j,
                                                Point3D &outPoint,
                                                InterpolationDiagnostics *diagnostics,
                                                std::string *errorMessage)
{
    if (diagnostics)
    {
        *diagnostics = InterpolationDiagnostics{};
    }
    if (!ctx.points || !ctx.grid.valid)
    {
        if (errorMessage)
        {
            *errorMessage = "Error: clamped local MLS context is invalid.";
        }
        return false;
    }
    const double gx = ctx.grid.xMin + static_cast<double>(i) * ctx.grid.gridSpacing;
    const double gy = ctx.grid.yMin + static_cast<double>(j) * ctx.grid.gridSpacing;

    BicubicNodeData node;
    InterpolationDiagnostics localDiag;
    if (!fitLocalPolynomialNode(*ctx.points, ctx.spatialIndex, gx, gy, ctx.grid.gridSpacing, node, &ctx.tuning, &localDiag))
    {
        if (errorMessage)
        {
            *errorMessage = "Error: local monotone MLS evaluation failed.";
        }
        return false;
    }

    auto neighbors = ctx.spatialIndex.gatherNearest(*ctx.points, gx, gy, std::min<size_t>(ctx.points->size(), 32U), 1U);
    double zMin = std::numeric_limits<double>::infinity();
    double zMax = -std::numeric_limits<double>::infinity();
    for (size_t idx : neighbors)
    {
        zMin = std::min(zMin, (*ctx.points)[idx].z);
        zMax = std::max(zMax, (*ctx.points)[idx].z);
    }
    double gz = node.z;
    bool clamped = false;
    if (std::isfinite(zMin) && std::isfinite(zMax))
    {
        const double clampedZ = std::max(zMin, std::min(zMax, gz));
        clamped = (clampedZ != gz);
        gz = clampedZ;
    }

    const bool outsideHull = ctx.hull.valid ? !pointInConvexHull2D(ctx.hull, gx, gy) : false;
    const double outsideDistance = outsideHull ? distanceToConvexHull2D(ctx.hull, gx, gy) : 0.0;
    const double boundaryPenalty = outsideHull ? (1.0 / (1.0 + outsideDistance / std::max(ctx.grid.gridSpacing, 1e-6))) : 1.0;

    outPoint = {gx, gy, gz};
    if (diagnostics)
    {
        *diagnostics = localDiag;
        diagnostics->outsideHull = outsideHull;
        diagnostics->outsideDistance = outsideDistance;
        diagnostics->boundaryPenalty = boundaryPenalty;
        diagnostics->confidence = calibratedConfidenceFromObservedError(diagnostics->confidence * boundaryPenalty, diagnostics->residual, ctx.validation);
        diagnostics->clampedLocal = clamped;
        diagnostics->summary = clamped ? "clamped local MLS (clamped)" : "clamped local MLS";
    }
    return std::isfinite(outPoint.z);
}

static bool evaluateClampedMLSStreamPoint(const ClampedMLSStreamingContext &ctx,
                                        size_t i,
                                        size_t j,
                                        Point3D &outPoint,
                                        std::string *errorMessage)
{
    return evaluateClampedMLSStreamPointDetailed(ctx, i, j, outPoint, nullptr, errorMessage);
}

static bool evaluateTPSStreamPointDetailed(const TPSStreamingContext &ctx,
                                           size_t i,
                                           size_t j,
                                           Point3D &outPoint,
                                           InterpolationDiagnostics *diagnostics,
                                           std::string *errorMessage)
{
    if (diagnostics)
    {
        *diagnostics = InterpolationDiagnostics{};
    }
    if (!ctx.grid.valid || !ctx.model.valid)
    {
        if (errorMessage)
        {
            *errorMessage = "Error: TPS context is invalid.";
        }
        return false;
    }

    const double gx = ctx.grid.xMin + static_cast<double>(i) * ctx.grid.gridSpacing;
    const double gy = ctx.grid.yMin + static_cast<double>(j) * ctx.grid.gridSpacing;
    double gz = 0.0;
    InterpolationDiagnostics localDiag;
    if (!solveLocalTPSAtPoint(gx, gy, ctx.model, gz, nullptr, &localDiag) || !std::isfinite(gz))
    {
        gz = inverseDistanceFallback(ctx.model.controlPoints, ctx.model.spatialIndex, gx, gy);
        localDiag.fallbackUsed = true;
        localDiag.confidence = 0.05;
    }
    outPoint = {gx, gy, gz};
    if (diagnostics)
    {
        *diagnostics = localDiag;
    }
    return std::isfinite(outPoint.x) && std::isfinite(outPoint.y) && std::isfinite(outPoint.z);
}

static bool evaluateTPSStreamPoint(const TPSStreamingContext &ctx,
                                   size_t i,
                                   size_t j,
                                   Point3D &outPoint,
                                   std::string *errorMessage)
{
    return evaluateTPSStreamPointDetailed(ctx, i, j, outPoint, nullptr, errorMessage);
}

static void maybeReportRowProgress(const std::function<void(const std::string &)> &progressCallback,
                                   const char *phaseLabel,
                                   size_t rowIndex,
                                   size_t rowCount)
{
    if (!progressCallback || rowCount == 0U)
    {
        return;
    }
    const size_t step = std::max<size_t>(1U, rowCount / 50U);
    const bool shouldReport = (rowIndex == 0U) || ((rowIndex + 1U) == rowCount) || (((rowIndex + 1U) % step) == 0U);
    if (!shouldReport)
    {
        return;
    }

    const double pct = (100.0 * static_cast<double>(rowIndex + 1U)) / static_cast<double>(rowCount);
    std::ostringstream oss;
    oss << phaseLabel << " row " << (rowIndex + 1U) << "/" << rowCount
        << " (" << std::fixed << std::setprecision(1) << pct << "%)";
    progressCallback(oss.str());
}

static bool writeGridXYZStream(const std::string &outputFileName,
                               const GridDefinition &grid,
                               int precision,
                               const std::function<bool(size_t, size_t, Point3D &, std::string *)> &evaluator,
                               const std::function<void(const std::string &)> &progressCallback,
                               std::string *errorMessage)
{
    std::ofstream outFile(outputFileName, std::ios::out | std::ios::binary);
    if (!outFile.is_open())
    {
        if (errorMessage)
        {
            *errorMessage = "Error creating grid XYZ file: " + outputFileName;
        }
        return false;
    }
    std::vector<char> fileBuffer(8U * 1024U * 1024U);
    outFile.rdbuf()->pubsetbuf(fileBuffer.data(), static_cast<std::streamsize>(fileBuffer.size()));
    outFile << std::fixed << std::setprecision(precision);

    if (progressCallback)
    {
        std::ostringstream msg;
        msg << "grid.xyz writer opened. Rows=" << grid.nx << ", cols=" << grid.ny << ", total nodes=" << grid.totalPoints;
        progressCallback(msg.str());
    }

    std::string buffer;
    buffer.reserve(4U * 1024U * 1024U);

    for (size_t i = 0; i < grid.nx; ++i)
    {
        for (size_t j = 0; j < grid.ny; ++j)
        {
            Point3D p{};
            std::string localError;
            if (!evaluator(i, j, p, &localError))
            {
                if (errorMessage)
                {
                    *errorMessage = localError.empty() ? ("Error evaluating interpolation point for: " + outputFileName) : localError;
                }
                return false;
            }
            buffer += formatPointXYZLine(p, precision);
            if (buffer.size() >= (4U * 1024U * 1024U))
            {
                outFile << buffer;
                if (!outFile.good())
                {
                    if (errorMessage)
                    {
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
    if (!ok && errorMessage)
    {
        *errorMessage = "Error writing grid XYZ file: " + outputFileName;
    }
    return ok;
}

static bool writeConfidenceXYZStream(const std::string &outputFileName,
                                     const GridDefinition &grid,
                                     int precision,
                                     const std::function<bool(size_t, size_t, Point3D &, InterpolationDiagnostics &, std::string *)> &evaluator,
                                     const std::function<void(const std::string &)> &progressCallback,
                                     ConfidenceSummary &summary,
                                     std::string *errorMessage)
{
    summary = ConfidenceSummary{};
    std::ofstream outFile(outputFileName, std::ios::out | std::ios::binary);
    if (!outFile.is_open())
    {
        if (errorMessage)
        {
            *errorMessage = "Error creating confidence XYZ file: " + outputFileName;
        }
        return false;
    }
    outFile << std::fixed << std::setprecision(precision);
    outFile << "# x y z confidence residual lambda boundaryPenalty fallback outsideHull clampedLocal neighborCount patchCount\n";

    for (size_t i = 0; i < grid.nx; ++i)
    {
        for (size_t j = 0; j < grid.ny; ++j)
        {
            Point3D p{};
            InterpolationDiagnostics d;
            std::string localError;
            if (!evaluator(i, j, p, d, &localError))
            {
                if (errorMessage)
                {
                    *errorMessage = localError.empty() ? ("Error evaluating confidence point for: " + outputFileName) : localError;
                }
                return false;
            }
            outFile << p.x << ' ' << p.y << ' ' << p.z << ' '
                    << d.confidence << ' ' << d.residual << ' ' << d.lambda << ' '
                    << d.boundaryPenalty << ' ' << (d.fallbackUsed ? 1 : 0) << ' '
                    << (d.outsideHull ? 1 : 0) << ' ' << (d.clampedLocal ? 1 : 0) << ' '
                    << d.neighborCount << ' ' << d.patchCount << "\n";
            ++summary.count;
            summary.meanConfidence += d.confidence;
            summary.meanResidual += d.residual;
            summary.fallbackCount += d.fallbackUsed ? 1U : 0U;
            summary.outsideHullCount += d.outsideHull ? 1U : 0U;
            summary.clampedLocalCount += d.clampedLocal ? 1U : 0U;
        }
        maybeReportRowProgress(progressCallback, "Writing confidence.xyz:", i, grid.nx);
    }
    if (summary.count > 0U)
    {
        summary.meanConfidence /= static_cast<double>(summary.count);
        summary.meanResidual /= static_cast<double>(summary.count);
    }
    const bool ok = outFile.good();
    outFile.close();
    if (!ok && errorMessage)
    {
        *errorMessage = "Error writing confidence XYZ file: " + outputFileName;
    }
    return ok;
}

static bool computeBoundingBoxFromXYZFile(const std::string &inputFileName,
                                          double &xmin,
                                          double &xmax,
                                          double &ymin,
                                          double &ymax)
{
    std::ifstream inFile(inputFileName);
    if (!inFile.is_open())
    {
        return false;
    }

    xmin = std::numeric_limits<double>::max();
    xmax = -std::numeric_limits<double>::max();
    ymin = std::numeric_limits<double>::max();
    ymax = -std::numeric_limits<double>::max();

    std::string line;
    bool found = false;
    while (std::getline(inFile, line))
    {
        if (line.empty() || line[0] == '#')
        {
            continue;
        }
        std::replace(line.begin(), line.end(), ',', ' ');
        std::istringstream ss(line);
        Point3D p{};
        if (!(ss >> p.x >> p.y >> p.z))
        {
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

static bool writeDXFStreamFromXYZFile(const std::string &outputFileName,
                                      const std::string &xyzInputFileName,
                                      const GridDefinition &grid,
                                      int precision,
                                      int pdmode,
                                      bool hasGrid,
                                      const std::function<bool(size_t, size_t, Point3D &, std::string *)> &evaluator,
                                      const std::function<void(const std::string &)> &progressCallback,
                                      std::string *errorMessage)
{
    double xmin = 0.0;
    double xmax = 0.0;
    double ymin = 0.0;
    double ymax = 0.0;
    if (!computeBoundingBoxFromXYZFile(xyzInputFileName, xmin, xmax, ymin, ymax))
    {
        if (errorMessage)
        {
            *errorMessage = "Error reading filtered XYZ file for DXF bounding box: " + xyzInputFileName;
        }
        return false;
    }

    if (hasGrid && grid.valid)
    {
        const double gridXMax = grid.xMin + static_cast<double>(grid.nx - 1U) * grid.gridSpacing;
        const double gridYMax = grid.yMin + static_cast<double>(grid.ny - 1U) * grid.gridSpacing;
        xmin = std::min(xmin, grid.xMin);
        xmax = std::max(xmax, gridXMax);
        ymin = std::min(ymin, grid.yMin);
        ymax = std::max(ymax, gridYMax);
    }

    std::ofstream outFile(outputFileName, std::ios::out | std::ios::binary);
    if (!outFile.is_open())
    {
        if (errorMessage)
        {
            *errorMessage = "Error creating DXF file: " + outputFileName;
        }
        return false;
    }
    std::vector<char> fileBuffer(8U * 1024U * 1024U);
    outFile.rdbuf()->pubsetbuf(fileBuffer.data(), static_cast<std::streamsize>(fileBuffer.size()));
    outFile << std::fixed << std::setprecision(precision);

    if (progressCallback)
    {
        std::ostringstream msg;
        msg << "DXF writer opened. Filtered source=" << xyzInputFileName << ", grid nodes=" << grid.totalPoints;
        progressCallback(msg.str());
    }

    const double centerX = 0.5 * (xmin + xmax);
    const double centerY = 0.5 * (ymin + ymax);
    double viewSize = std::max(xmax - xmin, ymax - ymin) * 1.1;
    if (!(viewSize > 0.0) || !std::isfinite(viewSize))
    {
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
            << "0\nLAYER\n2\nxyz_labels\n70\n0\n62\n3\n6\nCONTINUOUS\n";
    if (hasGrid)
    {
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

    auto appendPointAndLabel = [&](std::string &buffer,
                                   const Point3D &p,
                                   const std::string &layerPoints,
                                   const std::string &layerLabels)
    {
        std::ostringstream ss;
        ss << std::fixed << std::setprecision(precision);
        ss << "0\nPOINT\n8\n"
           << layerPoints
           << "\n10\n"
           << p.x
           << "\n20\n"
           << p.y
           << "\n30\n"
           << p.z
           << "\n0\nTEXT\n8\n"
           << layerLabels
           << "\n10\n"
           << (p.x + 0.2)
           << "\n20\n"
           << (p.y + 0.2)
           << "\n30\n0.0\n40\n1.0\n1\n"
           << p.z << "\n";
        buffer += ss.str();
    };

    std::ifstream inFile(xyzInputFileName);
    if (!inFile.is_open())
    {
        if (errorMessage)
        {
            *errorMessage = "Error reading filtered XYZ file for DXF: " + xyzInputFileName;
        }
        return false;
    }

    std::string buffer;
    buffer.reserve(64U * 1024U);
    std::string line;
    while (std::getline(inFile, line))
    {
        if (line.empty() || line[0] == '#')
        {
            continue;
        }
        std::replace(line.begin(), line.end(), ',', ' ');
        std::istringstream ss(line);
        Point3D p{};
        if (!(ss >> p.x >> p.y >> p.z))
        {
            continue;
        }
        appendPointAndLabel(buffer, p, "xyz_points", "xyz_labels");
        if (buffer.size() >= (64U * 1024U))
        {
            outFile << buffer;
            if (!outFile.good())
            {
                if (errorMessage)
                {
                    *errorMessage = "Error writing DXF file: " + outputFileName;
                }
                return false;
            }
            buffer.clear();
        }
    }

    if (hasGrid)
    {
        for (size_t i = 0; i < grid.nx; ++i)
        {
            for (size_t j = 0; j < grid.ny; ++j)
            {
                Point3D p{};
                std::string localError;
                if (!evaluator(i, j, p, &localError))
                {
                    if (errorMessage)
                    {
                        *errorMessage = localError.empty() ? ("Error evaluating interpolation point for DXF: " + outputFileName) : localError;
                    }
                    return false;
                }
                appendPointAndLabel(buffer, p, "grid_points", "grid_labels");
                if (buffer.size() >= (4U * 1024U * 1024U))
                {
                    outFile << buffer;
                    if (!outFile.good())
                    {
                        if (errorMessage)
                        {
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
    if (!ok && errorMessage)
    {
        *errorMessage = "Error writing DXF file: " + outputFileName;
    }
    return ok;
}

[[maybe_unused]]
static bool writeDXFStream(const std::string &outputFileName,
                           const std::vector<Point3D> &xyzPoints,
                           const GridDefinition &grid,
                           int precision,
                           int pdmode,
                           bool hasGrid,
                           const std::function<bool(size_t, size_t, Point3D &, std::string *)> &evaluator,
                           const std::function<void(const std::string &)> &progressCallback,
                           std::string *errorMessage)
{
    std::ofstream outFile(outputFileName);
    if (!outFile.is_open())
    {
        if (errorMessage)
        {
            *errorMessage = "Error creating DXF file: " + outputFileName;
        }
        return false;
    }
    outFile << std::fixed << std::setprecision(precision);

    double xmin = std::numeric_limits<double>::max();
    double xmax = -std::numeric_limits<double>::max();
    double ymin = std::numeric_limits<double>::max();
    double ymax = -std::numeric_limits<double>::max();

    for (const auto &p : xyzPoints)
    {
        xmin = std::min(xmin, p.x);
        xmax = std::max(xmax, p.x);
        ymin = std::min(ymin, p.y);
        ymax = std::max(ymax, p.y);
    }
    if (hasGrid && grid.valid)
    {
        const double gridXMax = grid.xMin + static_cast<double>(grid.nx - 1U) * grid.gridSpacing;
        const double gridYMax = grid.yMin + static_cast<double>(grid.ny - 1U) * grid.gridSpacing;
        xmin = std::min(xmin, grid.xMin);
        xmax = std::max(xmax, gridXMax);
        ymin = std::min(ymin, grid.yMin);
        ymax = std::max(ymax, gridYMax);
    }

    const double centerX = 0.5 * (xmin + xmax);
    const double centerY = 0.5 * (ymin + ymax);
    double viewSize = std::max(xmax - xmin, ymax - ymin) * 1.1;
    if (!(viewSize > 0.0) || !std::isfinite(viewSize))
    {
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
            << "0\nLAYER\n2\nxyz_labels\n70\n0\n62\n3\n6\nCONTINUOUS\n";
    if (hasGrid)
    {
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

    auto appendPointAndLabel = [&](std::string &buffer,
                                   const Point3D &p,
                                   const std::string &layerPoints,
                                   const std::string &layerLabels)
    {
        std::ostringstream ss;
        ss << std::fixed << std::setprecision(precision);
        ss << "0\nPOINT\n8\n"
           << layerPoints
           << "\n10\n"
           << p.x
           << "\n20\n"
           << p.y
           << "\n30\n"
           << p.z
           << "\n0\nTEXT\n8\n"
           << layerLabels
           << "\n10\n"
           << (p.x + 0.2)
           << "\n20\n"
           << (p.y + 0.2)
           << "\n30\n0.0\n40\n1.0\n1\n"
           << p.z << "\n";
        buffer += ss.str();
    };

    std::string buffer;
    buffer.reserve(64U * 1024U);
    for (const auto &p : xyzPoints)
    {
        appendPointAndLabel(buffer, p, "xyz_points", "xyz_labels");
        if (buffer.size() >= (64U * 1024U))
        {
            outFile << buffer;
            if (!outFile.good())
            {
                if (errorMessage)
                {
                    *errorMessage = "Error writing DXF file: " + outputFileName;
                }
                return false;
            }
            buffer.clear();
        }
    }

    if (hasGrid)
    {
        for (size_t i = 0; i < grid.nx; ++i)
        {
            for (size_t j = 0; j < grid.ny; ++j)
            {
                Point3D p{};
                std::string localError;
                if (!evaluator(i, j, p, &localError))
                {
                    if (errorMessage)
                    {
                        *errorMessage = localError.empty() ? ("Error evaluating interpolation point for DXF: " + outputFileName) : localError;
                    }
                    return false;
                }
                appendPointAndLabel(buffer, p, "grid_points", "grid_labels");
                if (buffer.size() >= (4U * 1024U * 1024U))
                {
                    outFile << buffer;
                    if (!outFile.good())
                    {
                        if (errorMessage)
                        {
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
    if (!ok && errorMessage)
    {
        *errorMessage = "Error writing DXF file: " + outputFileName;
    }
    return ok;
}

/****************************************************************************
 * 11) Writers
 ****************************************************************************/
static bool writeXYZWithProgress(const std::string &outputFileName,
                                 const std::vector<Point3D> &points,
                                 int precision,
                                 const char *stageName,
                                 std::string *errorMessage,
                                 const std::function<void(const std::string &)> *statusUpdate = nullptr)
{
    std::ofstream outFile(outputFileName, std::ios::out | std::ios::binary);
    if (!outFile.is_open())
    {
        if (errorMessage)
        {
            *errorMessage = std::string("Error creating ") + stageName + " file: " + outputFileName;
        }
        return false;
    }

    std::vector<char> fileBuffer(8U * 1024U * 1024U);
    outFile.rdbuf()->pubsetbuf(fileBuffer.data(), static_cast<std::streamsize>(fileBuffer.size()));
    outFile << std::fixed << std::setprecision(precision);

    std::string chunk;
    chunk.reserve(4U * 1024U * 1024U);

    const size_t total = points.size();
    const size_t reportStep = std::max<size_t>(50000U, total / 100U + 1U);
    size_t nextReport = reportStep;

    for (size_t i = 0; i < total; ++i)
    {
        chunk += formatPointXYZLine(points[i], precision);
        if (chunk.size() >= 4U * 1024U * 1024U)
        {
            outFile.write(chunk.data(), static_cast<std::streamsize>(chunk.size()));
            chunk.clear();
            if (!outFile.good())
            {
                if (errorMessage)
                {
                    *errorMessage = std::string("Error writing ") + stageName + " file: " + outputFileName;
                }
                return false;
            }
        }

        if (statusUpdate && i + 1U >= nextReport)
        {
            std::ostringstream msg;
            const double pct = total > 0U ? (100.0 * static_cast<double>(i + 1U) / static_cast<double>(total)) : 100.0;
            msg << "Writing " << stageName << " ... "
                << (i + 1U) << "/" << total
                << " (" << std::fixed << std::setprecision(1) << pct << "%)";
            (*statusUpdate)(msg.str());
            nextReport += reportStep;
        }
    }

    if (!chunk.empty())
    {
        outFile.write(chunk.data(), static_cast<std::streamsize>(chunk.size()));
    }

    const bool ok = outFile.good();
    outFile.close();
    if (!ok)
    {
        if (errorMessage)
        {
            *errorMessage = std::string("Error writing ") + stageName + " file: " + outputFileName;
        }
        return false;
    }

    if (statusUpdate)
    {
        std::ostringstream msg;
        msg << "Writing " << stageName << " ... " << total << "/" << total << " (100.0%)";
        (*statusUpdate)(msg.str());
    }
    return true;
}

static bool writeFilteredXYZ(const std::string &outputFileName,
                             const std::vector<Point3D> &filteredPoints,
                             int precision,
                             std::string *errorMessage,
                             const std::function<void(const std::string &)> *statusUpdate = nullptr)
{
    return writeXYZWithProgress(outputFileName,
                                filteredPoints,
                                precision,
                                "filtered.xyz",
                                errorMessage,
                                statusUpdate);
}

[[maybe_unused]]
static bool writeGridXYZ(const std::string &outputFileName,
                         const std::vector<Point3D> &gridPoints,
                         int precision,
                         std::string *errorMessage,
                         const std::function<void(const std::string &)> *statusUpdate = nullptr)
{
    return writeXYZWithProgress(outputFileName,
                                gridPoints,
                                precision,
                                "grid.xyz",
                                errorMessage,
                                statusUpdate);
}

[[maybe_unused]]
static bool writeDXF(const std::string &outputFileName,
                     const std::vector<Point3D> &xyzPoints,
                     const std::vector<Point3D> &gridPoints,
                     int precision,
                     int pdmode,
                     bool hasGrid,
                     std::string *errorMessage)
{
    std::ofstream outFile(outputFileName);
    if (!outFile.is_open())
    {
        if (errorMessage)
        {
            *errorMessage = "Error creating DXF file: " + outputFileName;
        }
        return false;
    }
    outFile << std::fixed << std::setprecision(precision);

    outFile << "0\nSECTION\n2\nHEADER\n"
            << "9\n$PDMODE\n70\n"
            << pdmode << "\n"
            << "9\n$PDSIZE\n40\n0.5\n"
            << "0\nENDSEC\n";

    double xmin = std::numeric_limits<double>::max();
    double xmax = -std::numeric_limits<double>::max();
    double ymin = std::numeric_limits<double>::max();
    double ymax = -std::numeric_limits<double>::max();

    for (const auto &p : xyzPoints)
    {
        if (p.x < xmin) xmin = p.x;
        if (p.x > xmax) xmax = p.x;
        if (p.y < ymin) ymin = p.y;
        if (p.y > ymax) ymax = p.y;
    }
    if (hasGrid)
    {
        for (const auto &p : gridPoints)
        {
            if (p.x < xmin) xmin = p.x;
            if (p.x > xmax) xmax = p.x;
            if (p.y < ymin) ymin = p.y;
            if (p.y > ymax) ymax = p.y;
        }
    }

    const double centerX = 0.5 * (xmin + xmax);
    const double centerY = 0.5 * (ymin + ymax);
    double viewSize = std::max(xmax - xmin, ymax - ymin) * 1.1;
    if (!(viewSize > 0.0) || !std::isfinite(viewSize))
    {
        viewSize = 1.0;
    }

    outFile << "0\nSECTION\n2\nTABLES\n"
            << "0\nTABLE\n2\nLAYER\n"
            << "0\nLAYER\n2\nxyz_points\n70\n0\n62\n7\n6\nCONTINUOUS\n"
            << "0\nLAYER\n2\nxyz_labels\n70\n0\n62\n3\n6\nCONTINUOUS\n";
    if (hasGrid)
    {
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

    auto writePointAndLabel = [&](const Point3D &p,
                                  const std::string &layerPoints,
                                  const std::string &layerLabels)
    {
        outFile << "0\nPOINT\n8\n"
                << layerPoints
                << "\n10\n"
                << p.x
                << "\n20\n"
                << p.y
                << "\n30\n"
                << p.z
                << "\n0\nTEXT\n8\n"
                << layerLabels
                << "\n10\n"
                << (p.x + 0.2)
                << "\n20\n"
                << (p.y + 0.2)
                << "\n30\n0.0\n40\n1.0\n1\n"
                << p.z << "\n";
    };

    for (const auto &p : xyzPoints)
    {
        writePointAndLabel(p, "xyz_points", "xyz_labels");
    }
    if (hasGrid)
    {
        for (const auto &p : gridPoints)
        {
            writePointAndLabel(p, "grid_points", "grid_labels");
        }
    }

    outFile << "0\nENDSEC\n0\nEOF\n";
    const bool ok = outFile.good();
    outFile.close();
    if (!ok && errorMessage)
    {
        *errorMessage = "Error writing DXF file: " + outputFileName;
    }
    return ok;
}

/*****************************************************************************
 * 12) processXYZtoDXF
 ****************************************************************************/
static bool processXYZtoDXF(const std::string &inputFileName,
                            double minDist,
                            int precision,
                            int pdmode,
                            double gridSpacing,
                            size_t maxTPSPoints,
                            InterpolationMethod methodUsed,
                            std::function<void(const std::string &)> statusUpdate,
                            std::string &finalMessage)
{
    std::vector<std::string> reportLines;

    std::function<void(const std::string &)> reportStatus = [&](const std::string &msg)
    {
        reportLines.push_back(msg);
        statusUpdate(msg);
    };

    auto writeReportFile = [&]()
    {
        const std::string reportFileName = inputFileName + ".rpt.txt";
        std::ofstream rptFile(reportFileName);
        if (rptFile.is_open())
        {
            for (const auto &line : reportLines)
            {
                rptFile << line << "\n";
            }
        }
    };

    auto finish = [&](bool success, const std::string &msg) -> bool
    {
        finalMessage = msg;
        if (!msg.empty())
        {
            reportStatus(msg);
        }
        writeReportFile();
        return success;
    };

    try
    {
        const auto startTime = std::chrono::high_resolution_clock::now();

        if (minDist < 0.0)
        {
            return finish(false, "Error: minDist must be >= 0.");
        }
        if (precision < 0 || precision > 15)
        {
            return finish(false, "Error: precision must be between 0 and 15.");
        }
        if (!(gridSpacing > 0.0) || !std::isfinite(gridSpacing))
        {
            return finish(false, "Error: gridSpacing must be a positive finite value.");
        }

        reportStatus("Phase 1/6: reading input file and applying minDist filter ...");
        size_t totalPointsRead = 0U;
        std::vector<Point3D> filteredPoints;
        std::string readError;
        if (!readAndFilterXYZFile(inputFileName,
                                  minDist,
                                  filteredPoints,
                                  totalPointsRead,
                                  [&](const std::string &msg) { reportStatus(msg); },
                                  &readError))
        {
            return finish(false, readError.empty() ? "Error: No valid points found in file." : readError);
        }
        {
            std::ostringstream msg;
            msg << "Total valid points read: " << totalPointsRead;
            reportStatus(msg.str());
        }
        {
            std::ostringstream msg;
            msg << "Points after minDist filter: " << filteredPoints.size();
            reportStatus(msg.str());
        }
        if (filteredPoints.empty())
        {
            return finish(false, "Error: All points were removed by the minDist filter.");
        }

        reportStatus("Phase 2/6: robust local residual outlier filtering ...");
        const double neighborDist = std::max(5.0 * minDist, 0.01);
        const double zThresholdFactor = 3.5;
        {
            std::ostringstream msg;
            msg << "Robust outlier parameters: neighborDist=" << neighborDist
                << ", zThresholdFactor=" << zThresholdFactor;
            reportStatus(msg.str());
        }
        {
            std::string outlierError;
            if (!removeZOutliersInPlace(filteredPoints,
                                        neighborDist,
                                        zThresholdFactor,
                                        [&](const std::string &msg) { reportStatus(msg); },
                                        &outlierError))
            {
                return finish(false, outlierError.empty() ? "Error: All points were removed during robust outlier filtering." : outlierError);
            }
        }
        {
            std::ostringstream msg;
            msg << "Points after robust outlier removal: " << filteredPoints.size();
            reportStatus(msg.str());
        }

        reportStatus("Phase 3/6: writing filtered.xyz to disk before interpolation ...");
        {
            std::string ioError;
            if (!writeFilteredXYZ(inputFileName + ".filtered.xyz", filteredPoints, precision, &ioError, &reportStatus))
            {
                return finish(false, ioError);
            }
        }
        reportStatus("filtered.xyz successfully written.");

        GridDefinition grid;
        std::function<bool(size_t, size_t, Point3D &, std::string *)> evaluator;
        std::function<bool(size_t, size_t, Point3D &, InterpolationDiagnostics &, std::string *)> diagEvaluator;
        std::string methodLabel;
        ValidationMetrics methodValidation;

        if (methodUsed == METHOD_TPS)
        {
            methodLabel = "TPS";
            GridDefinition fullExtentGrid;
            std::string gridError;
            if (!buildGridDefinition(filteredPoints, gridSpacing, fullExtentGrid, &gridError))
            {
                return finish(false, gridError);
            }

            std::vector<Point3D> tpsSubset;
            const bool useAllFilteredPointsForTPS = (maxTPSPoints == 0U) || (filteredPoints.size() <= maxTPSPoints);
            if (useAllFilteredPointsForTPS)
            {
                reportStatus("TPS control set: using all filtered points without making an extra copy.");
                tpsSubset.swap(filteredPoints);
            }
            else
            {
                reportStatus("TPS control set: Poisson/FPS-style decimation in progress ...");
                tpsSubset = subsamplePointsUniformly(filteredPoints, maxTPSPoints);
            }
            {
                std::ostringstream msg;
                msg << "TPS subset size: " << tpsSubset.size();
                reportStatus(msg.str());
            }
            if (tpsSubset.empty())
            {
                return finish(false, "Error: TPS received an empty control-point subset.");
            }

            TPSStreamingContext tpsCtx;
            std::string ctxError;
            if (!buildTPSStreamingContext(std::move(tpsSubset), fullExtentGrid, tpsCtx, &ctxError))
            {
                return finish(false, ctxError);
            }
            grid = tpsCtx.grid;
            methodValidation = tpsCtx.model.validation;
            evaluator = [tpsCtx](size_t i, size_t j, Point3D &p, std::string *err) -> bool
            {
                return evaluateTPSStreamPoint(tpsCtx, i, j, p, err);
            };
            diagEvaluator = [tpsCtx](size_t i, size_t j, Point3D &p, InterpolationDiagnostics &d, std::string *err) -> bool
            {
                return evaluateTPSStreamPointDetailed(tpsCtx, i, j, p, &d, err);
            };
            reportStatus("Phase 4/6: TPS context ready with cross-validated local lambda selection.");
        }
        else if (methodUsed == METHOD_CLAMPED_MLS)
        {
            methodLabel = "Clamped local MLS (bounded)";
            ClampedMLSStreamingContext monoCtx;
            std::string ctxError;
            if (!buildClampedMLSStreamingContext(filteredPoints, gridSpacing, monoCtx, &ctxError))
            {
                return finish(false, ctxError);
            }
            grid = monoCtx.grid;
            methodValidation = monoCtx.validation;
            evaluator = [monoCtx](size_t i, size_t j, Point3D &p, std::string *err) -> bool
            {
                return evaluateClampedMLSStreamPoint(monoCtx, i, j, p, err);
            };
            diagEvaluator = [monoCtx](size_t i, size_t j, Point3D &p, InterpolationDiagnostics &d, std::string *err) -> bool
            {
                return evaluateClampedMLSStreamPointDetailed(monoCtx, i, j, p, &d, err);
            };
            reportStatus("Phase 4/6: clamped local MLS context ready with predictive holdout validation.");
        }
        else
        {
            methodLabel = "Bicubic Hermite + robust MLS";
            BicubicStreamingContext bicCtx;
            std::string ctxError;
            if (!buildBicubicStreamingContext(filteredPoints, gridSpacing, bicCtx, &ctxError))
            {
                return finish(false, ctxError);
            }
            grid = bicCtx.grid;
            methodValidation = bicCtx.validation;
            {
                std::ostringstream msg;
                msg << "Selected bicubic tuning: neighbors=" << bicCtx.tuning.neighborhoodTarget
                    << ", supportMultiplier=" << bicCtx.tuning.supportMultiplier
                    << ", ridgeFactor=" << bicCtx.tuning.ridgeFactor
                    << ", basisOrder=" << bicCtx.tuning.basisOrder;
                reportStatus(msg.str());
            }
            evaluator = [bicCtx](size_t i, size_t j, Point3D &p, std::string *err) -> bool
            {
                return evaluateBicubicStreamPoint(bicCtx, i, j, p, err);
            };
            diagEvaluator = [bicCtx](size_t i, size_t j, Point3D &p, InterpolationDiagnostics &d, std::string *err) -> bool
            {
                return evaluateBicubicStreamPointDetailed(bicCtx, i, j, p, &d, err);
            };
            reportStatus("Phase 4/6: bicubic context ready with predictive robust anisotropic MLS and 2D derivative limiting.");
        }

        if (methodValidation.valid)
        {
            std::ostringstream msg;
            msg << "Validation summary for " << methodLabel << ": RMSE=" << methodValidation.rmse
                << ", MAE=" << methodValidation.mae
                << ", P95=" << methodValidation.p95
                << ", n=" << methodValidation.count;
            reportStatus(msg.str());
        }

        {
            std::ostringstream msg;
            msg << "Grid nodes planned: " << grid.totalPoints << " (" << grid.nx << " x " << grid.ny << ")"
                << ", spacing=" << grid.gridSpacing;
            reportStatus(msg.str());
        }

        reportStatus("Phase 5/6: starting streamed grid.xyz writing ...");
        {
            std::string ioError;
            if (!writeGridXYZStream(inputFileName + ".grid.xyz",
                                    grid,
                                    precision,
                                    evaluator,
                                    [&](const std::string &msg) { reportStatus(msg); },
                                    &ioError))
            {
                return finish(false, ioError);
            }
        }

        reportStatus("Phase 5b/6: writing confidence diagnostics grid ...");
        ConfidenceSummary confidenceSummary;
        {
            std::string ioError;
            if (!writeConfidenceXYZStream(inputFileName + ".confidence.xyz",
                                          grid,
                                          precision,
                                          diagEvaluator,
                                          [&](const std::string &msg) { reportStatus(msg); },
                                          confidenceSummary,
                                          &ioError))
            {
                return finish(false, ioError);
            }
        }
        {
            std::ostringstream msg;
            msg << "Confidence summary: meanConfidence=" << confidenceSummary.meanConfidence
                << ", meanResidual=" << confidenceSummary.meanResidual
                << ", fallbacks=" << confidenceSummary.fallbackCount
                << ", outsideHull=" << confidenceSummary.outsideHullCount
                << ", clampedLocal=" << confidenceSummary.clampedLocalCount;
            reportStatus(msg.str());
        }

        reportStatus("Phase 6/6: starting streamed DXF generation ...");
        {
            std::string ioError;
            if (!writeDXFStreamFromXYZFile(inputFileName + ".dxf",
                                           inputFileName + ".filtered.xyz",
                                           grid,
                                           precision,
                                           pdmode,
                                           grid.valid && grid.totalPoints > 0U,
                                           evaluator,
                                           [&](const std::string &msg) { reportStatus(msg); },
                                           &ioError))
            {
                return finish(false, ioError);
            }
        }

        const auto endTime = std::chrono::high_resolution_clock::now();
        const double elapsedSec = std::chrono::duration<double>(endTime - startTime).count();
        const int totalSeconds = static_cast<int>(std::round(elapsedSec));

        std::ostringstream finishMsg;
        finishMsg << "Done. Total time: " << totalSeconds << " sec.";
        return finish(true, finishMsg.str());
    }
    catch (const std::bad_alloc &)
    {
        return finish(false, "Error: insufficient memory during interpolation. The process now aborts safely instead of crashing.");
    }
    catch (const std::exception &e)
    {
        return finish(false, std::string("Error: ") + e.what());
    }
    catch (...)
    {
        return finish(false, "Error: unexpected failure during interpolation.");
    }
}

struct GuiStatusMessage
{
    std::string message;
};

struct GuiCompletionMessage
{
    bool success = false;
    std::string message;
};


/****************************************************************************
 * Professional GUI layout helpers
 ****************************************************************************/
struct GuiLayoutHandles
{
    HWND grpSource = NULL;
    HWND lblInputFile = NULL;
    HWND hInputFile = NULL;
    HWND btnBrowse = NULL;
    HWND txtSourceHint = NULL;

    HWND grpParameters = NULL;
    HWND lblMinDist = NULL;
    HWND hMinDist = NULL;
    HWND txtMinDist = NULL;
    HWND lblPrecision = NULL;
    HWND hPrecision = NULL;
    HWND txtPrecision = NULL;
    HWND lblPDMODE = NULL;
    HWND hPDMODE = NULL;
    HWND txtPDMODE = NULL;
    HWND lblGridSpacing = NULL;
    HWND hGridSpacing = NULL;
    HWND txtGridSpacing = NULL;
    HWND lblMaxTPSPoints = NULL;
    HWND hMaxTPSPoints = NULL;
    HWND txtMaxTPSPoints = NULL;

    HWND grpMethod = NULL;
    HWND txtMethodIntro = NULL;
    HWND hRadioBicubic = NULL;
    HWND txtBicubic = NULL;
    HWND hRadioTPS = NULL;
    HWND txtTPS = NULL;
    HWND hRadioClamped = NULL;
    HWND txtClamped = NULL;

    HWND grpExecution = NULL;
    HWND btnRun = NULL;
    HWND txtExecutionHint = NULL;
    HWND hwndStatus = NULL;
};

static void ApplyFontToWindowTree(const GuiLayoutHandles &g, HFONT font, HFONT fontBold)
{
    const HWND handles[] = {
        g.grpSource, g.lblInputFile, g.hInputFile, g.btnBrowse, g.txtSourceHint,
        g.grpParameters, g.lblMinDist, g.hMinDist, g.txtMinDist,
        g.lblPrecision, g.hPrecision, g.txtPrecision,
        g.lblPDMODE, g.hPDMODE, g.txtPDMODE,
        g.lblGridSpacing, g.hGridSpacing, g.txtGridSpacing,
        g.lblMaxTPSPoints, g.hMaxTPSPoints, g.txtMaxTPSPoints,
        g.grpMethod, g.txtMethodIntro, g.hRadioBicubic, g.txtBicubic,
        g.hRadioTPS, g.txtTPS, g.hRadioClamped, g.txtClamped,
        g.grpExecution, g.btnRun, g.txtExecutionHint, g.hwndStatus
    };
    for (HWND h : handles)
    {
        if (h != NULL)
        {
            SendMessageA(h, WM_SETFONT, reinterpret_cast<WPARAM>(font), TRUE);
        }
    }
    const HWND groupBoxes[] = {g.grpSource, g.grpParameters, g.grpMethod, g.grpExecution};
    for (HWND h : groupBoxes)
    {
        if (h != NULL)
        {
            SendMessageA(h, WM_SETFONT, reinterpret_cast<WPARAM>(fontBold != NULL ? fontBold : font), TRUE);
        }
    }
    if (g.btnRun != NULL && fontBold != NULL)
    {
        SendMessageA(g.btnRun, WM_SETFONT, reinterpret_cast<WPARAM>(fontBold), TRUE);
    }
}

static void ApplyProfessionalGuiLayout(HWND hwnd, const GuiLayoutHandles &g)
{
    RECT rc{};
    GetClientRect(hwnd, &rc);
    const int clientW = rc.right - rc.left;

    const int margin = 16;
    const int sectionGap = 12;
    const int labelW = 118;
    const int editW = 126;
    const int browseW = 102;
    const int topGroupH = 96;
    const int leftGroupH = 282;
    const int rightGroupH = 282;
    const int bottomGroupH = 126;
    const int midTop = margin + topGroupH + sectionGap;
    const int leftW = std::max(455, (clientW - 3 * margin) / 2);
    const int rightW = clientW - 3 * margin - leftW;

    const int xLeft = margin;
    const int xRight = margin + leftW + margin;
    const int yBottom = midTop + std::max(leftGroupH, rightGroupH) + sectionGap;

    MoveWindow(g.grpSource, margin, margin, clientW - 2 * margin, topGroupH, TRUE);
    MoveWindow(g.lblInputFile, margin + 16, margin + 24, labelW, 22, TRUE);
    MoveWindow(g.hInputFile, margin + 16 + labelW, margin + 20,
               clientW - 2 * margin - 16 - labelW - browseW - 36, 30, TRUE);
    MoveWindow(g.btnBrowse, clientW - margin - browseW - 16, margin + 20, browseW, 30, TRUE);
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
    y += 48;

    MoveWindow(g.lblMaxTPSPoints, paramX, y + 4, labelW, 20, TRUE);
    MoveWindow(g.hMaxTPSPoints, editX, y, editW, rowH, TRUE);
    MoveWindow(g.txtMaxTPSPoints, hintX, y + 1, hintW, hintH + 12, TRUE);

    MoveWindow(g.grpMethod, xRight, midTop, rightW, rightGroupH, TRUE);
    const int methodX = xRight + 18;
    const int methodW = rightW - 36;
    int my = midTop + 26;
    MoveWindow(g.txtMethodIntro, methodX, my, methodW, 40, TRUE);
    my += 44;
    MoveWindow(g.hRadioBicubic, methodX, my, methodW, 22, TRUE);
    my += 22;
    MoveWindow(g.txtBicubic, methodX + 20, my, methodW - 20, 40, TRUE);
    my += 52;
    MoveWindow(g.hRadioTPS, methodX, my, methodW, 22, TRUE);
    my += 22;
    MoveWindow(g.txtTPS, methodX + 20, my, methodW - 20, 40, TRUE);
    my += 52;
    MoveWindow(g.hRadioClamped, methodX, my, methodW, 22, TRUE);
    my += 22;
    MoveWindow(g.txtClamped, methodX + 20, my, methodW - 20, 46, TRUE);

    MoveWindow(g.grpExecution, margin, yBottom, clientW - 2 * margin, bottomGroupH, TRUE);
    MoveWindow(g.txtExecutionHint, margin + 18, yBottom + 24, clientW - 2 * margin - 232, 24, TRUE);
    MoveWindow(g.btnRun, clientW - margin - 166, yBottom + 18, 148, 36, TRUE);
    MoveWindow(g.hwndStatus, margin + 18, yBottom + 58, clientW - 2 * margin - 36, 52, TRUE);
}

/****************************************************************************
 * Window Procedure (GUI code)
 ****************************************************************************/

LRESULT CALLBACK WindowProc(HWND hwnd, UINT uMsg, WPARAM wParam, LPARAM lParam)
{
    static GuiLayoutHandles gui;
    static HFONT hGuiFont = NULL;
    static HFONT hGuiFontBold = NULL;
    static std::shared_ptr<std::atomic_bool> windowAlive;

    switch (uMsg)
    {
    case WM_CREATE:
    {
        windowAlive = std::make_shared<std::atomic_bool>(true);

        hGuiFont = CreateFontA(-15, 0, 0, 0, FW_NORMAL, FALSE, FALSE, FALSE,
                               DEFAULT_CHARSET, OUT_DEFAULT_PRECIS, CLIP_DEFAULT_PRECIS,
                               CLEARTYPE_QUALITY, DEFAULT_PITCH | FF_DONTCARE, "Segoe UI");
        hGuiFontBold = CreateFontA(-15, 0, 0, 0, FW_SEMIBOLD, FALSE, FALSE, FALSE,
                                   DEFAULT_CHARSET, OUT_DEFAULT_PRECIS, CLIP_DEFAULT_PRECIS,
                                   CLEARTYPE_QUALITY, DEFAULT_PITCH | FF_DONTCARE, "Segoe UI");

        gui.grpSource = CreateWindowA("BUTTON", " Source File ", WS_VISIBLE | WS_CHILD | BS_GROUPBOX,
                                      0, 0, 10, 10, hwnd, NULL, NULL, NULL);
        gui.lblInputFile = CreateWindowA("STATIC", "Input XYZ file:", WS_VISIBLE | WS_CHILD,
                                         0, 0, 10, 10, hwnd, NULL, NULL, NULL);
        gui.hInputFile = CreateWindowExA(WS_EX_CLIENTEDGE, "EDIT", "",
                                         WS_CHILD | WS_VISIBLE | WS_BORDER | ES_AUTOHSCROLL | WS_TABSTOP,
                                         0, 0, 10, 10, hwnd, (HMENU)IDC_INPUT_FILE, NULL, NULL);
        gui.btnBrowse = CreateWindowA("BUTTON", "Browse...",
                                      WS_VISIBLE | WS_CHILD | WS_TABSTOP | BS_PUSHBUTTON,
                                      0, 0, 10, 10, hwnd, (HMENU)IDC_BROWSE_BUTTON, NULL, NULL);
        gui.txtSourceHint = CreateWindowA("STATIC",
                                          "Select the source XYZ point cloud. The full path remains visible while editing.",
                                          WS_VISIBLE | WS_CHILD,
                                          0, 0, 10, 10, hwnd, NULL, NULL, NULL);

        gui.grpParameters = CreateWindowA("BUTTON", " Processing Parameters ", WS_VISIBLE | WS_CHILD | BS_GROUPBOX,
                                          0, 0, 10, 10, hwnd, NULL, NULL, NULL);
        gui.lblMinDist = CreateWindowA("STATIC", "Min Dist:", WS_VISIBLE | WS_CHILD,
                                       0, 0, 10, 10, hwnd, NULL, NULL, NULL);
        gui.hMinDist = CreateWindowExA(WS_EX_CLIENTEDGE, "EDIT", "5",
                                       WS_CHILD | WS_VISIBLE | WS_BORDER | ES_AUTOHSCROLL | WS_TABSTOP,
                                       0, 0, 10, 10, hwnd, (HMENU)IDC_MIN_DIST, NULL, NULL);
        gui.txtMinDist = CreateWindowA("STATIC", "Minimum XY spacing for point filtering.", WS_VISIBLE | WS_CHILD,
                                       0, 0, 10, 10, hwnd, NULL, NULL, NULL);

        gui.lblPrecision = CreateWindowA("STATIC", "Precision:", WS_VISIBLE | WS_CHILD,
                                         0, 0, 10, 10, hwnd, NULL, NULL, NULL);
        gui.hPrecision = CreateWindowExA(WS_EX_CLIENTEDGE, "EDIT", "2",
                                         WS_CHILD | WS_VISIBLE | WS_BORDER | ES_AUTOHSCROLL | WS_TABSTOP,
                                         0, 0, 10, 10, hwnd, (HMENU)IDC_PRECISION, NULL, NULL);
        gui.txtPrecision = CreateWindowA("STATIC", "Decimal places in XYZ and DXF outputs.", WS_VISIBLE | WS_CHILD,
                                         0, 0, 10, 10, hwnd, NULL, NULL, NULL);

        gui.lblPDMODE = CreateWindowA("STATIC", "PDMODE:", WS_VISIBLE | WS_CHILD,
                                      0, 0, 10, 10, hwnd, NULL, NULL, NULL);
        gui.hPDMODE = CreateWindowExA(WS_EX_CLIENTEDGE, "EDIT", "3",
                                      WS_CHILD | WS_VISIBLE | WS_BORDER | ES_AUTOHSCROLL | WS_TABSTOP,
                                      0, 0, 10, 10, hwnd, (HMENU)IDC_PDMODE, NULL, NULL);
        gui.txtPDMODE = CreateWindowA("STATIC", "DXF point-style written into the DXF header.", WS_VISIBLE | WS_CHILD,
                                      0, 0, 10, 10, hwnd, NULL, NULL, NULL);

        gui.lblGridSpacing = CreateWindowA("STATIC", "Grid Spacing:", WS_VISIBLE | WS_CHILD,
                                           0, 0, 10, 10, hwnd, NULL, NULL, NULL);
        gui.hGridSpacing = CreateWindowExA(WS_EX_CLIENTEDGE, "EDIT", "20",
                                           WS_CHILD | WS_VISIBLE | WS_BORDER | ES_AUTOHSCROLL | WS_TABSTOP,
                                           0, 0, 10, 10, hwnd, (HMENU)IDC_GRID_SPACING, NULL, NULL);
        gui.txtGridSpacing = CreateWindowA("STATIC", "Regular interpolation spacing of surface grid.", WS_VISIBLE | WS_CHILD,
                                           0, 0, 10, 10, hwnd, NULL, NULL, NULL);

        gui.lblMaxTPSPoints = CreateWindowA("STATIC", "Max TPS Points:", WS_VISIBLE | WS_CHILD,
                                            0, 0, 10, 10, hwnd, NULL, NULL, NULL);
        gui.hMaxTPSPoints = CreateWindowExA(WS_EX_CLIENTEDGE, "EDIT", "0",
                                            WS_CHILD | WS_VISIBLE | WS_BORDER | ES_AUTOHSCROLL | WS_TABSTOP,
                                            0, 0, 10, 10, hwnd, (HMENU)IDC_MAX_TPS_POINTS, NULL, NULL);
        gui.txtMaxTPSPoints = CreateWindowA("STATIC", "Maximum TPS points. 0 to keep all points.",
                                            WS_VISIBLE | WS_CHILD,
                                            0, 0, 10, 10, hwnd, NULL, NULL, NULL);

        gui.grpMethod = CreateWindowA("BUTTON", " Interpolation Method ", WS_VISIBLE | WS_CHILD | BS_GROUPBOX,
                                      0, 0, 10, 10, hwnd, NULL, NULL, NULL);
        gui.txtMethodIntro = CreateWindowA("STATIC", "Choose the interpolation strategy that best matches the spatial character of the input cloud.",
                                           WS_VISIBLE | WS_CHILD,
                                           0, 0, 10, 10, hwnd, NULL, NULL, NULL);
        gui.hRadioBicubic = CreateWindowA("BUTTON", "Bicubic Hermite + robust MLS (default)",
                                          WS_VISIBLE | WS_CHILD | BS_AUTORADIOBUTTON | WS_TABSTOP,
                                          0, 0, 10, 10, hwnd, (HMENU)IDC_RADIO_BICUBIC, NULL, NULL);
        SendMessageA(gui.hRadioBicubic, BM_SETCHECK, BST_CHECKED, 0);
        gui.txtBicubic = CreateWindowA("STATIC", "Best general-purpose option for smooth surfaces and regular grid production.",
                                       WS_VISIBLE | WS_CHILD,
                                       0, 0, 10, 10, hwnd, NULL, NULL, NULL);
        gui.hRadioTPS = CreateWindowA("BUTTON", "Thin Plate Spline (scattered data)",
                                      WS_VISIBLE | WS_CHILD | BS_AUTORADIOBUTTON | WS_TABSTOP,
                                      0, 0, 10, 10, hwnd, (HMENU)IDC_RADIO_TPS, NULL, NULL);
        gui.txtTPS = CreateWindowA("STATIC", "Preferred for irregularly distributed points and locally flexible scattered-data interpolation.",
                                   WS_VISIBLE | WS_CHILD,
                                   0, 0, 10, 10, hwnd, NULL, NULL, NULL);
        gui.hRadioClamped = CreateWindowA("BUTTON", "Clamped local MLS (bounded)",
                                          WS_VISIBLE | WS_CHILD | BS_AUTORADIOBUTTON | WS_TABSTOP,
                                          0, 0, 10, 10, hwnd, (HMENU)IDC_RADIO_CLAMPED_MLS, NULL, NULL);
        gui.txtClamped = CreateWindowA("STATIC", "Conservative local option that limits overshoot by constraining values to the local data envelope.",
                                       WS_VISIBLE | WS_CHILD,
                                       0, 0, 10, 10, hwnd, NULL, NULL, NULL);

        gui.grpExecution = CreateWindowA("BUTTON", " Execution ", WS_VISIBLE | WS_CHILD | BS_GROUPBOX,
                                         0, 0, 10, 10, hwnd, NULL, NULL, NULL);
        gui.txtExecutionHint = CreateWindowA("STATIC", "Review the parameter values, choose the interpolation method, and then run the conversion.",
                                             WS_VISIBLE | WS_CHILD,
                                             0, 0, 10, 10, hwnd, NULL, NULL, NULL);
        gui.btnRun = CreateWindowA("BUTTON", "Run Conversion",
                                   WS_VISIBLE | WS_CHILD | BS_DEFPUSHBUTTON | WS_TABSTOP,
                                   0, 0, 10, 10, hwnd, (HMENU)IDC_RUN_BUTTON, NULL, NULL);
        gui.hwndStatus = CreateWindowExA(WS_EX_CLIENTEDGE, "STATIC", "Status: Idle",
                                         WS_VISIBLE | WS_CHILD | WS_BORDER,
                                         0, 0, 10, 10, hwnd, (HMENU)IDC_STATUS_STATIC, NULL, NULL);

        ApplyFontToWindowTree(gui, hGuiFont, hGuiFontBold);
        ApplyProfessionalGuiLayout(hwnd, gui);
        return 0;
    }

    case WM_COMMAND:
    {
        if (LOWORD(wParam) == IDC_BROWSE_BUTTON)
        {
            const std::string filePath = openFileDialog(hwnd);
            if (!filePath.empty())
            {
                SetWindowTextA(gui.hInputFile, filePath.c_str());
            }
        }
        else if (LOWORD(wParam) == IDC_RUN_BUTTON)
        {
            EnableWindow(GetDlgItem(hwnd, IDC_RUN_BUTTON), FALSE);

            std::vector<char> inputFile(32768, '\0');
            char minDistStr[64], precisionStr[64], pdModeStr[64];
            char gridSpacingStr[64], maxTPSStr[64];
            GetWindowTextA(gui.hInputFile, inputFile.data(), static_cast<int>(inputFile.size()));
            GetWindowTextA(gui.hMinDist, minDistStr, sizeof(minDistStr));
            GetWindowTextA(gui.hPrecision, precisionStr, sizeof(precisionStr));
            GetWindowTextA(gui.hPDMODE, pdModeStr, sizeof(pdModeStr));
            GetWindowTextA(gui.hGridSpacing, gridSpacingStr, sizeof(gridSpacingStr));
            GetWindowTextA(gui.hMaxTPSPoints, maxTPSStr, sizeof(maxTPSStr));

            if (strlen(inputFile.data()) == 0)
            {
                MessageBoxA(hwnd, "Please select an input file.", "Error", MB_ICONERROR | MB_OK);
                EnableWindow(GetDlgItem(hwnd, IDC_RUN_BUTTON), TRUE);
                break;
            }

            double d_minDist = 0.0;
            int i_precision = 2;
            int i_pdmode = 3;
            double d_gridSpacing = 10.0;
            size_t s_maxTPSPoints = 0;

            try
            {
                d_minDist = std::stod(minDistStr);
                i_precision = std::stoi(precisionStr);
                i_pdmode = std::stoi(pdModeStr);
                d_gridSpacing = std::stod(gridSpacingStr);
                s_maxTPSPoints = static_cast<size_t>(std::stoull(maxTPSStr));
            }
            catch (...)
            {
                MessageBoxA(hwnd, "One or more numeric inputs are invalid.", "Error", MB_ICONERROR | MB_OK);
                EnableWindow(GetDlgItem(hwnd, IDC_RUN_BUTTON), TRUE);
                break;
            }

            if (d_minDist < 0.0)
            {
                MessageBoxA(hwnd, "Min Dist must be >= 0.", "Error", MB_ICONERROR | MB_OK);
                EnableWindow(GetDlgItem(hwnd, IDC_RUN_BUTTON), TRUE);
                break;
            }
            if (i_precision < 0 || i_precision > 15)
            {
                MessageBoxA(hwnd, "Precision must be between 0 and 15.", "Error", MB_ICONERROR | MB_OK);
                EnableWindow(GetDlgItem(hwnd, IDC_RUN_BUTTON), TRUE);
                break;
            }
            if (!(d_gridSpacing > 0.0) || !std::isfinite(d_gridSpacing))
            {
                MessageBoxA(hwnd, "Grid Spacing must be a positive finite value.", "Error", MB_ICONERROR | MB_OK);
                EnableWindow(GetDlgItem(hwnd, IDC_RUN_BUTTON), TRUE);
                break;
            }

            InterpolationMethod methodUsed = METHOD_BICUBIC;
            if (SendMessageA(gui.hRadioTPS, BM_GETCHECK, 0, 0) == BST_CHECKED)
            {
                methodUsed = METHOD_TPS;
            }
            else if (SendMessageA(gui.hRadioClamped, BM_GETCHECK, 0, 0) == BST_CHECKED)
            {
                methodUsed = METHOD_CLAMPED_MLS;
            }

            updateStatus(gui.hwndStatus, "Starting conversion...");

            const std::string inputFileStr(inputFile.data());
            const auto aliveFlag = windowAlive;
            std::thread processingThread([hwnd, aliveFlag, inputFileStr, d_minDist, i_precision, i_pdmode,
                                          d_gridSpacing, s_maxTPSPoints, methodUsed]()
            {
                auto postStatus = [hwnd, aliveFlag](const std::string &msg)
                {
                    if (!aliveFlag || !aliveFlag->load() || !IsWindow(hwnd))
                    {
                        return;
                    }
                    GuiStatusMessage *sm = new GuiStatusMessage{msg};
                    if (!PostMessageA(hwnd, WM_APP_STATUS_UPDATE, 0, reinterpret_cast<LPARAM>(sm)))
                    {
                        delete sm;
                    }
                };

                bool success = false;
                std::string finalMessage;
                try
                {
                    success = processXYZtoDXF(inputFileStr,
                                              d_minDist,
                                              i_precision,
                                              i_pdmode,
                                              d_gridSpacing,
                                              s_maxTPSPoints,
                                              methodUsed,
                                              postStatus,
                                              finalMessage);
                }
                catch (const std::bad_alloc &)
                {
                    success = false;
                    finalMessage = "Error: insufficient memory during interpolation.";
                }
                catch (const std::exception &e)
                {
                    success = false;
                    finalMessage = std::string("Error: ") + e.what();
                }
                catch (...)
                {
                    success = false;
                    finalMessage = "Error: unexpected worker-thread failure.";
                }

                if (!aliveFlag || !aliveFlag->load() || !IsWindow(hwnd))
                {
                    return;
                }

                GuiCompletionMessage *cm = new GuiCompletionMessage{success, finalMessage};
                if (!PostMessageA(hwnd, WM_APP_PROCESSING_DONE, 0, reinterpret_cast<LPARAM>(cm)))
                {
                    delete cm;
                }
            });
            processingThread.detach();
        }
        return 0;
    }

    case WM_APP_STATUS_UPDATE:
    {
        std::unique_ptr<GuiStatusMessage> sm(reinterpret_cast<GuiStatusMessage *>(lParam));
        if (sm && gui.hwndStatus && IsWindow(gui.hwndStatus))
        {
            updateStatus(gui.hwndStatus, sm->message);
        }
        return 0;
    }

    case WM_APP_PROCESSING_DONE:
    {
        std::unique_ptr<GuiCompletionMessage> cm(reinterpret_cast<GuiCompletionMessage *>(lParam));
        EnableWindow(GetDlgItem(hwnd, IDC_RUN_BUTTON), TRUE);
        if (cm && gui.hwndStatus && IsWindow(gui.hwndStatus))
        {
            const std::string fallback = cm->success ? "Conversion completed." : "Conversion failed.";
            updateStatus(gui.hwndStatus, cm->message.empty() ? fallback : cm->message);
        }
        return 0;
    }

    case WM_DESTROY:
        if (windowAlive)
        {
            windowAlive->store(false);
        }
        if (hGuiFont != NULL)
        {
            DeleteObject(hGuiFont);
            hGuiFont = NULL;
        }
        if (hGuiFontBold != NULL)
        {
            DeleteObject(hGuiFontBold);
            hGuiFontBold = NULL;
        }
        PostQuitMessage(0);
        return 0;
    }

    return DefWindowProcA(hwnd, uMsg, wParam, lParam);
}

/*****************************************************************************
 * WinMain
 ****************************************************************************/
int WINAPI WinMain(HINSTANCE hInstance, HINSTANCE, LPSTR, int nCmdShow)
{
    const char CLASS_NAME[] = "xyz2dxfWindowClass";

    WNDCLASSA wc = {};
    wc.lpfnWndProc   = WindowProc;
    wc.hInstance     = hInstance;
    wc.lpszClassName = CLASS_NAME;

    RegisterClassA(&wc);

    RECT wr{0, 0, 950, 580};
    const DWORD windowStyle = WS_OVERLAPPED | WS_CAPTION | WS_SYSMENU | WS_MINIMIZEBOX;
    AdjustWindowRect(&wr, windowStyle, FALSE);

    HWND hwnd = CreateWindowExA(0, CLASS_NAME, "XYZ to DXF Converter - Professional Edition",
                                windowStyle,
                                CW_USEDEFAULT, CW_USEDEFAULT, wr.right - wr.left, wr.bottom - wr.top,
                                NULL, NULL, hInstance, NULL);
    if (!hwnd)
        return 0;

    ShowWindow(hwnd, nCmdShow);

    MSG msg;
    while (GetMessageA(&msg, NULL, 0, 0))
    {
        TranslateMessage(&msg);
        DispatchMessageA(&msg);
    }
    return 0;
}
