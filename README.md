# XYZ to DXF Converter GUI

## Abstract

The **XYZ to DXF Converter GUI** is a Windows-based scientific and engineering application for the ingestion, filtering, regularization, interpolation, and CAD export of large irregular three-dimensional point clouds in XYZ format. The program is not a mere file converter. Its numerical purpose is the construction of a **stable, quality-controlled, interpolated surface representation** suitable for engineering inspection, grid export, and DXF-based visualization.

The active code implements three interpolation modes:

1. **Bicubic Hermite + robust anisotropic MLS**, intended as the primary smooth high-accuracy method.
2. **Local Thin Plate Spline (TPS)** with adaptive neighborhoods, local regularization, and blended patches, intended for highly scattered and irregular point distributions.
3. **Clamped local MLS (bounded)**, a conservative local alternative that constrains predictions to the neighborhood data envelope in order to reduce overshoot.

The software also includes:

- deterministic minimum-distance thinning independent of file order;
- robust local outlier rejection using local plane residuals and MAD-based scale estimation;
- predictive validation for parameter selection and method quality reporting;
- support-mask / hull awareness to reduce overconfidence near extrapolative regions;
- streamed writing of filtered XYZ, interpolated grid XYZ, confidence diagnostics, and DXF;
- GUI-based parameter entry with Windows-native file selection and runtime status reporting.

This document is written as a **technical scientific manual**. It explains the mathematics, software architecture, numerical workflow, data structures, diagnostics, outputs, and implementation philosophy of the active code base.

---

## 1. Scope and technical objective

Let the raw input point cloud be

```math
\mathcal{P} = \{(x_i, y_i, z_i)\}_{i=1}^{N}.
```

The computational workflow transforms this cloud into a hierarchy of engineering products:

```math
\mathcal{P}
\to
\mathcal{P}_{d_{\min}}
\to
\mathcal{P}_{\mathrm{clean}}
\to
\mathcal{G}.
```

with the following stage meanings:

- `\mathcal{P} \to \mathcal{P}_{d_{\min}}`: deterministic minimum-distance filtering;
- `\mathcal{P}_{d_{\min}} \to \mathcal{P}_{\mathrm{clean}}`: robust local residual outlier rejection;
- `\mathcal{P}_{\mathrm{clean}} \to \mathcal{G}`: regular-grid interpolation;
- final exports: `filtered.xyz`, `grid.xyz`, `confidence.xyz`, `DXF`, and `report`.

The numerical objectives are:

- to preserve geometric fidelity while removing redundant sampling;
- to suppress local vertical anomalies without over-smoothing the cloud;
- to interpolate on a regular Cartesian grid suitable for downstream engineering use;
- to quantify local interpolation confidence and failure or fallback conditions;
- to export the processed surface to standard XYZ and DXF products.

The code is therefore a coupled **preprocessing + interpolation + diagnostics + export** system.

---

## 2. High-level software workflow

The active program logic is organized in the following stages:

1. **Input reading and deterministic filtering**
2. **Robust local residual outlier rejection**
3. **Interpolation context construction**
4. **Predictive tuning / validation**
5. **Streamed writing of `filtered.xyz`**
6. **Streamed writing of `grid.xyz`**
7. **Streamed writing of `confidence.xyz`**
8. **Streamed writing of `dxf`**
9. **Textual execution report and GUI status messages**

The code is implemented in a single C++ source file and uses:

- Win32 GUI primitives (`windows.h`, `commdlg.h`);
- OpenMP for selected parallel loops;
- explicit dense linear algebra without external matrix libraries;
- streamed file I/O for large grid and DXF products.

The current build command is:

```bash
g++ -O3 -fopenmp -march=native -std=c++17 -Wall -Wextra -pedantic \
    -Wconversion -Wsign-conversion -static -static-libgcc -static-libstdc++ \
    -mwindows -o xyz2dxf_gui.exe \
    xyz2dxf_gui.cpp -lkernel32 -lopengl32 -luuid -lcomdlg32 -lm
```

---

## 3. Graphical interface and exposed parameters

The GUI exposes the following primary numerical inputs:

- `Min Dist`
- `Precision`
- `PDMODE`
- `Grid Spacing`
- `Max TPS Points`

It also exposes the interpolation method selector:

- **Bicubic Hermite + robust MLS**
- **Local TPS**
- **Clamped local MLS (bounded)**

The default values in the code are:

- `d_min = 5`
- `Precision = 2`
- `PDMODE = 3`
- `Grid Spacing = 20`
- `Max TPS Points = 0`

### 3.1 Minimum distance

The quantity `d_min` is the minimum admissible horizontal separation between retained points after thinning. If `d_min <= 0`, the code switches to exact duplicate removal only.

### 3.2 Precision

The parameter `Precision` controls only the number of decimal places written to output files. It does **not** change internal floating-point precision, which remains `double`.

### 3.3 PDMODE

`PDMODE` is written into the DXF header and controls the CAD point-display style:

`PDMODE` is a user-specified integer.

The DXF point size is written as the constant

```math
\mathrm{PDSIZE} = 0.5.
```

### 3.4 Grid spacing

The interpolation grid spacing is

```math
h = \mathrm{GridSpacing}.
```

This is the Cartesian node spacing of the regular output grid.

### 3.5 Maximum TPS points

`Max TPS Points` limits the number of control points used in the TPS branch. If the value is zero, the full filtered cloud is used. If the value is positive and smaller than the filtered point count, a deterministic decimation stage is applied before TPS model construction.

---

## 4. Deterministic minimum-distance filtering

### 4.1 Motivation

A simple “first-come-first-kept” filter is order-dependent: the result changes if the lines in the XYZ file are reordered. The active code avoids that problem by constructing a **deterministic pre-filter sequence** before acceptance testing.

### 4.2 Exact duplicate mode

If

```math
d_{\min} \le 0,
```

the code performs only exact duplicate removal. Two points are considered duplicates if and only if

```math
(x_i, y_i, z_i) = (x_j, y_j, z_j)
```

in exact `double` comparison.

### 4.3 Cell assignment

If `d_min > 0`, each point is assigned to a spatial cell

```math
(i_x, i_y) =
\left(
\left\lfloor \frac{x}{d_{\min}} \right\rfloor,
\left\lfloor \frac{y}{d_{\min}} \right\rfloor
\right).
```

For the cell center,

```math
x_c = \left(i_x + \frac{1}{2}\right)d_{\min},
\qquad
y_c = \left(i_y + \frac{1}{2}\right)d_{\min},
```

the point receives the deterministic score

```math
s = (x-x_c)^2 + (y-y_c)^2.
```

The candidate list is then sorted primarily by cell index and secondarily by this score, followed by deterministic coordinate tie-breaks. This makes the thinning stage independent of the original file order.

### 4.4 Acceptance criterion

For a candidate point

```math
p_i = (x_i, y_i, z_i),
```

the code searches only the `3 x 3` neighborhood of adjacent cells already accepted in the sparse hash table. A point is rejected if any previously retained point `p_j` satisfies

```math
(x_i-x_j)^2 + (y_i-y_j)^2 < d_{\min}^2.
```

Otherwise the point is accepted.

### 4.5 Properties

This strategy gives the filter the following properties:

- deterministic under line reordering;
- sparse in memory;
- local in spatial access;
- consistent with minimum-distance thinning in the horizontal plane.

---

## 5. Robust local outlier rejection

### 5.1 Philosophy

The active code does **not** rely on a simple global `z` threshold. Instead it uses a **local plane residual test** with robust scale estimation. This is substantially better for sloping topography, bathymetry, embankments, breaklines, and irregular local trends.

### 5.2 Neighborhood radius

The local radius is defined as

```math
r_n = \max(5d_{\min},\,0.01).
```

Only neighboring points within that radius are considered candidates for the local residual model.

### 5.3 Local weighted plane fit

For each tested point `(x_0, y_0)`, the local surface is approximated as

```math
z(x,y) \approx a + b_x(x-x_0) + b_y(y-y_0).
```

The weighted least-squares system uses the basis

```math
\phi(x,y) = \left(1,\ x-x_0,\ y-y_0\right)^{\mathsf T}.
```

Let

```math
\Delta x_k = x_k - x_0,
\qquad
\Delta y_k = y_k - y_0.
```

The local reference scale is

```math
h =
\max\!\left(
\sqrt{\frac{1}{n}\sum_{k=1}^{n}(\Delta x_k^2+\Delta y_k^2)},
\,10^{-6}
\right).
```

The geometric weight is

```math
w_k =
\frac{1}{
1 + \dfrac{\Delta x_k^2+\Delta y_k^2}{h^2}
}.
```

The normal equations are

```math
(\mathbf{A}^{\mathsf T}\mathbf{W}\mathbf{A})\,\mathbf{c}
=
\mathbf{A}^{\mathsf T}\mathbf{W}\mathbf{z},
```

with

```math
\mathbf{c} = \left(a,\ b_x,\ b_y\right)^{\mathsf T}.
```

A very small diagonal stabilization is added:

```math
\mathbf{A}^{\mathsf T}\mathbf{W}\mathbf{A}
\leftarrow
\mathbf{A}^{\mathsf T}\mathbf{W}\mathbf{A} + 10^{-12}\mathbf{I}.
```

### 5.4 Leave-one-out logic

When testing a point `p_i`, the point itself is excluded from the fitting neighborhood. This prevents the tested point from masking its own outlier status.

### 5.5 Residual and robust scale

Once the local plane is fitted, local residuals are formed as

```math
r_k = z_k - \left(a + b_x\Delta x_k + b_y\Delta y_k\right).
```

The residual center is the sample median:

```math
\widetilde{r} = \mathrm{median}\{r_k\}.
```

The absolute deviations are

```math
d_k = |r_k - \widetilde{r}|.
```

The MAD scale is

```math
\mathrm{MAD} = \mathrm{median}\{d_k\}.
```

The robust scale is then

```math
\sigma_{\mathrm{rob}} =
\max\!\left(1.4826\,\mathrm{MAD},\,10^{-8}\right).
```

### 5.6 Rejection criterion

The tested point residual is evaluated against the fitted local plane and median-centered residual field. The point is retained if

```math
|r_i - \widetilde{r}| \le \alpha \sigma_{\mathrm{rob}},
```

where in the active pipeline

```math
\alpha = 3.5.
```

If the neighborhood is too small for a stable local estimate, the point is conservatively retained.

---

## 6. Grid construction

After filtering and outlier removal, the code constructs a regular Cartesian grid over an expanded bounding box.

Let the clean cloud bounding box be

```math
x_{\min}^{d},\ x_{\max}^{d},\ y_{\min}^{d},\ y_{\max}^{d}.
```

The grid margin is

```math
m = 1.5h.
```

Hence the interpolation domain is

```math
x_{\min} = x_{\min}^{d} - m,
\qquad
x_{\max} = x_{\max}^{d} + m,
```

```math
y_{\min} = y_{\min}^{d} - m,
\qquad
y_{\max} = y_{\max}^{d} + m.
```

The grid sizes are

```math
n_x = \left\lceil \frac{\max(x_{\max}-x_{\min},\,1)}{h} \right\rceil + 1,
```

```math
n_y = \left\lceil \frac{\max(y_{\max}-y_{\min},\,1)}{h} \right\rceil + 1.
```

The grid nodes are

```math
x_i = x_{\min} + ih,
\qquad
i = 0,\dots,n_x-1,
```

```math
y_j = y_{\min} + jh,
\qquad
j = 0,\dots,n_y-1.
```

The total node count is

```math
N_g = n_x n_y.
```

Overflow checks are applied when forming this product.

---

## 7. Auxiliary geometric quantities

### 7.1 Average planar spacing

The average spacing used repeatedly in the code is

```math
\bar{\Delta}_{xy} = \sqrt{\frac{A}{N}},
```

with

```math
A = \max\!\left((x_{\max}-x_{\min})(y_{\max}-y_{\min}),\,10^{-12}\right).
```

This quantity is used to define spatial hash cell sizes and local neighborhood scales.

### 7.2 Inverse-distance fallback

Whenever a preferred interpolant fails, the code falls back to inverse-distance weighting:

```math
z(x,y) =
\frac{\sum_{k=1}^{m} w_k z_k}{\sum_{k=1}^{m} w_k},
```

with

```math
w_k = \frac{1}{10^{-12} + (x-x_k)^2 + (y-y_k)^2}.
```

This fallback is a safeguard, not the preferred method.

---

## 8. Local anisotropic coordinate system

Both the robust MLS and TPS branches rely on a local anisotropic coordinate system obtained from a weighted covariance analysis.

### 8.1 Preliminary scale

For a query point `(x_0, y_0)` and neighborhood `\mathcal{N}`,

```math
\bar{r}^{\,2}
=
\frac{1}{n}\sum_{k=1}^{n}\left[(x_k-x_0)^2 + (y_k-y_0)^2\right].
```

The preliminary scale is

```math
h_0 =
\max\!\left(
\sqrt{\max(\bar{r}^{\,2},\,10^{-12})},
\ \max(h,\,10^{-6})
\right).
```

### 8.2 Preliminary weights

The covariance weights are

```math
w_k^{\mathrm{frame}} =
\frac{1}{
1 + \dfrac{(x_k-x_0)^2 + (y_k-y_0)^2}{h_0^2}
}.
```

### 8.3 Weighted covariance matrix

The code forms the weighted second moments

```math
S_{xx} = \frac{1}{W}\sum_k w_k^{\mathrm{frame}}(x_k-x_0)^2,
```

```math
S_{yy} = \frac{1}{W}\sum_k w_k^{\mathrm{frame}}(y_k-y_0)^2,
```

```math
S_{xy} = \frac{1}{W}\sum_k w_k^{\mathrm{frame}}(x_k-x_0)(y_k-y_0),
```

with

```math
W = \sum_k w_k^{\mathrm{frame}}.
```

### 8.4 Principal directions

The trace is

```math
\mathrm{tr}(\mathbf{S}) = S_{xx}+S_{yy}.
```

The spectral gap term is

```math
\Delta = \sqrt{(S_{xx}-S_{yy})^2 + 4S_{xy}^2}.
```

The principal eigenvalues are

```math
\lambda_1 =
\max\!\left(
\frac{\mathrm{tr}(\mathbf{S})+\Delta}{2},
\,0
\right),
```

```math
\lambda_2 =
\max\!\left(
\frac{\mathrm{tr}(\mathbf{S})-\Delta}{2},
\,0
\right).
```

The principal-frame angle is obtained from the two-argument arctangent:

```math
\theta = \frac{1}{2}\,\mathrm{atan2}(2S_{xy},\,S_{xx}-S_{yy}).
```

Therefore,

```math
\mathbf{e}_1 = (\cos\theta,\ \sin\theta),
\qquad
\mathbf{e}_2 = (-\sin\theta,\ \cos\theta).
```

### 8.5 Local scales

The anisotropic normalization scales are

```math
s_1 = \max\!\left(\sqrt{\lambda_1},\ 0.75\max(h,10^{-6})\right),
```

```math
s_2 = \max\!\left(\sqrt{\lambda_2},\ 0.35\max(h,10^{-6})\right).
```

### 8.6 Local coordinates

For any point `(x,y)`,

```math
u = \frac{(x-x_0)e_{1x} + (y-y_0)e_{1y}}{s_1},
```

```math
v = \frac{(x-x_0)e_{2x} + (y-y_0)e_{2y}}{s_2}.
```

The normalized radial coordinate is

```math
\rho = \sqrt{u^2 + v^2}.
```

---

## 9. Adaptive neighborhood support

The active code uses quantile-based support selection. The ordered local radii `\rho_k` are computed, and quantiles are extracted, including approximately the 80th and 90th percentiles.

The final support radius is of the form

```math
\rho_s = \max\!\left(\rho_{\min},\ 1.5,\ \rho_{80},\ 0.75\rho_{90}\right),
```

possibly multiplied by a tuning factor in the MLS branches.

The neighborhood is then selected adaptively from the spatial-hash candidate pool under minimum and maximum cardinality requirements.

---

## 10. Robust anisotropic MLS

### 10.1 Basis order

The code supports two polynomial basis orders:

- **quadratic**, with 6 coefficients;
- **cubic**, with 10 coefficients.

Thus,

```math
n_b = 6 \quad (\mathrm{quadratic\ MLS}),
```

```math
n_b = 10 \quad (\mathrm{cubic\ MLS}).
```

### 10.2 Polynomial basis

The local polynomial basis used in the code is

```math
\phi_0 = 1,
\qquad
\phi_1 = u,
\qquad
\phi_2 = v,
```

```math
\phi_3 = \frac{1}{2}u^2,
\qquad
\phi_4 = uv,
\qquad
\phi_5 = \frac{1}{2}v^2.
```

If cubic mode is enabled, the basis is extended by

```math
\phi_6 = \frac{1}{6}u^3,
\qquad
\phi_7 = \frac{1}{2}u^2v,
```

```math
\phi_8 = \frac{1}{2}uv^2,
\qquad
\phi_9 = \frac{1}{6}v^3.
```

The local approximation is

```math
p(u,v) = \sum_{k=0}^{n_b-1} c_k\phi_k(u,v).
```

### 10.3 Compact support kernel

For `q < 1`, the compact kernel is

```math
w_c(q) = (1-q)^4(4q+1).
```

For `q \ge 1`, the compact kernel is zero.

The normalized support coordinate is

```math
q = \frac{\rho}{\rho_s}.
```

The geometric MLS weight is

```math
w_{\mathrm{geom}}(\rho) =
\max\!\left(
\frac{w_c(q)}{1+\rho^2},
\,10^{-14}
\right).
```

### 10.4 Robust reweighting

The MLS solve is repeated for three robust IRLS-like passes. Let `r_k` be the residual at a fitted point. The code computes the median residual center and MAD scale:

```math
\widetilde{r} = \mathrm{median}\{r_k\},
```

```math
\sigma = \max\!\left(1.4826\,\mathrm{MAD},\,10^{-8}\right).
```

The normalized robust coordinate is

```math
t_k = \frac{|r_k|}{4.685\,\sigma}.
```

For `t_k < 1`, the robust weight is

```math
w_k^{\mathrm{rob}} = \max\!\left((1-t_k^2)^2,\,10^{-4}\right).
```

For `t_k \ge 1`, the robust weight is

```math
w_k^{\mathrm{rob}} = 10^{-4}.
```

The total fitting weight becomes

```math
w_k =
\max\!\left(
w_{\mathrm{geom},k}\,w_k^{\mathrm{rob}},
\,10^{-14}
\right).
```

### 10.5 Regularized normal equations

The MLS coefficients are obtained from

```math
\left(\mathbf{A}^{\mathsf T}\mathbf{W}\mathbf{A} + \lambda\mathbf{I}\right)\mathbf{c}
=
\mathbf{A}^{\mathsf T}\mathbf{W}\mathbf{z}.
```

The code computes a trace-based ridge parameter

```math
\lambda =
\max\!\left(
10^{-12},
\ \mathrm{ridgeFactor}\times 10^{-11}
\frac{\mathrm{tr}(\mathbf{A}^{\mathsf T}\mathbf{W}\mathbf{A})}{n_b}
\right).
```

### 10.6 Validation metrics

For a set of errors `e_k`, the code reports

```math
\mathrm{RMSE} =
\sqrt{\frac{1}{n}\sum_{k=1}^{n} e_k^2},
```

```math
\mathrm{MAE} =
\frac{1}{n}\sum_{k=1}^{n}|e_k|,
```

and the 95th percentile absolute error `P95`.

The combined selection score is

```math
S = \mathrm{RMSE} + 0.35\,P95 + 0.15\,\mathrm{MAE}.
```

### 10.7 Predictive local validation

The code also performs limited local predictive validation by removing selected neighborhood points and predicting them from the remaining set. This provides a **local holdout-like quality estimate** for node construction.

---

## 11. Bicubic Hermite + robust MLS

### 11.1 Nodal quantities

At each grid node, the MLS fit supplies a local value and local derivatives:

```math
z = c_0,
\qquad
f_u = c_1,
\qquad
f_v = c_2.
```

The second-order local quantities are

```math
f_{uu} = c_3,
\qquad
f_{uv} = c_4,
\qquad
f_{vv} = c_5.
```

### 11.2 Chain-rule mapping to global coordinates

The transformation derivatives are

```math
\frac{\partial u}{\partial x} = \frac{e_{1x}}{s_1},
\qquad
\frac{\partial u}{\partial y} = \frac{e_{1y}}{s_1},
```

```math
\frac{\partial v}{\partial x} = \frac{e_{2x}}{s_2},
\qquad
\frac{\partial v}{\partial y} = \frac{e_{2y}}{s_2}.
```

Hence,

```math
f_x = f_u\frac{\partial u}{\partial x} + f_v\frac{\partial v}{\partial x},
```

```math
f_y = f_u\frac{\partial u}{\partial y} + f_v\frac{\partial v}{\partial y},
```

```math
f_{xy}
=
f_{uu}\frac{\partial u}{\partial x}\frac{\partial u}{\partial y}
+
f_{uv}
\left(
\frac{\partial u}{\partial x}\frac{\partial v}{\partial y}
+
\frac{\partial v}{\partial x}\frac{\partial u}{\partial y}
\right)
+
f_{vv}\frac{\partial v}{\partial x}\frac{\partial v}{\partial y}.
```

### 11.3 Nodewise derivative limiting

Before nodal data are used in Hermite patches, the code limits nodal derivatives using a local vertical range

```math
\Delta z_{\mathrm{loc}} =
\max\!\left(z_{\max}^{\mathrm{loc}} - z_{\min}^{\mathrm{loc}},\,10^{-8}\right).
```

The slope and twist bounds are

```math
L_s = \frac{3\Delta z_{\mathrm{loc}}}{\max(h,10^{-6})},
```

```math
L_t = \frac{4\Delta z_{\mathrm{loc}}}{\max(h^2,10^{-6})}.
```

The clamped node derivatives satisfy

```math
f_x \in [-L_s,L_s],
\qquad
f_y \in [-L_s,L_s],
\qquad
f_{xy} \in [-L_t,L_t].
```

### 11.4 Cellwise minmod derivative limiting

At the cell level, the code computes directional slopes

```math
s_{x0} = \frac{z_{10}-z_{00}}{h_x},
\qquad
s_{x1} = \frac{z_{11}-z_{01}}{h_x},
```

```math
s_{y0} = \frac{z_{01}-z_{00}}{h_y},
\qquad
s_{y1} = \frac{z_{11}-z_{10}}{h_y}.
```

The minmod operator is defined as follows:

- `minmod(a,b) = 0` if `ab <= 0`;
- `minmod(a,b) = a` if `ab > 0` and `|a| < |b|`;
- `minmod(a,b) = b` if `ab > 0` and `|b| <= |a|`.

The cellwise derivative limits are applied consistently to `f_x`, `f_y`, and to the cross derivative via a cross-base estimate. This is intended to reduce overshoot and oscillatory behavior.

### 11.5 Hermite basis functions

The bicubic Hermite evaluation uses the classical cubic basis functions

```math
h_{00}(u) = 2u^3 - 3u^2 + 1,
\qquad
h_{10}(u) = u^3 - 2u^2 + u,
```

```math
h_{01}(u) = -2u^3 + 3u^2,
\qquad
h_{11}(u) = u^3 - u^2,
```

and analogously in `v`.

### 11.6 Bicubic cell surface

For node data `(z, f_x, f_y, f_{xy})` at corners `00`, `10`, `01`, and `11`, the interpolated value is

```math
f(u,v) = T_0 + T_x + T_y + T_{xy}.
```

The value term is

```math
T_0 =
z_{00}h_{00}(u)h_{00}(v)
+ z_{10}h_{01}(u)h_{00}(v)
+ z_{01}h_{00}(u)h_{01}(v)
+ z_{11}h_{01}(u)h_{01}(v).
```

The `x`-derivative contribution is

```math
T_x =
h_x\!\left[
f_{x,00}h_{10}(u)h_{00}(v)
+ f_{x,10}h_{11}(u)h_{00}(v)
+ f_{x,01}h_{10}(u)h_{01}(v)
+ f_{x,11}h_{11}(u)h_{01}(v)
\right].
```

The `y`-derivative contribution is

```math
T_y =
h_y\!\left[
f_{y,00}h_{00}(u)h_{10}(v)
+ f_{y,10}h_{01}(u)h_{10}(v)
+ f_{y,01}h_{00}(u)h_{11}(v)
+ f_{y,11}h_{01}(u)h_{11}(v)
\right].
```

The mixed-derivative contribution is

```math
T_{xy} =
h_x h_y\!\left[
f_{xy,00}h_{10}(u)h_{10}(v)
+ f_{xy,10}h_{11}(u)h_{10}(v)
+ f_{xy,01}h_{10}(u)h_{11}(v)
+ f_{xy,11}h_{11}(u)h_{11}(v)
\right].
```

### 11.7 Final value envelope clamp

After Hermite evaluation, the code constrains the result to the corner envelope enlarged by a small pad:

```math
z_{\min} = \min(z_{00},z_{10},z_{01},z_{11}),
```

```math
z_{\max} = \max(z_{00},z_{10},z_{01},z_{11}),
```

```math
p = 0.05\max(z_{\max}-z_{\min},\,10^{-8}),
```

```math
f(u,v) \leftarrow
\max\!\left(
z_{\min}-p,\,
\min(z_{\max}+p,\,f(u,v))
\right).
```

This is an additional safety clamp against excessive oscillation.

### 11.8 Parameter tuning

The bicubic or MLS branch tunes:

- target neighborhood size;
- support multiplier;
- ridge factor;
- polynomial order.

The tested values are:

- neighborhood target in `{32, 48, 64, 96}`;
- support multiplier in `{1.5, 1.75, 2.25}`;
- ridge factor in `{1, 10}`;
- basis order in `{2, 3}`.

Selection is based on predictive global holdout metrics.

---

## 12. Clamped local MLS (bounded)

This method reuses the robust anisotropic MLS node evaluation, but instead of constructing a bicubic Hermite patch, it evaluates a single local MLS prediction at the query location and then constrains the value to the local data envelope.

For the neighborhood `\mathcal{N}`,

```math
z_{\min}^{\mathcal{N}} = \min_{k\in\mathcal{N}} z_k,
```

```math
z_{\max}^{\mathcal{N}} = \max_{k\in\mathcal{N}} z_k.
```

The final value is

```math
z^{*} =
\max\!\left(
z_{\min}^{\mathcal{N}},
\min\!\left(z_{\max}^{\mathcal{N}},\,z_{\mathrm{MLS}}\right)
\right).
```

This mode is intended for conservative local evaluation with reduced overshoot risk. It is not a full monotonicity-preserving finite-volume scheme; rather, it is a bounded local regression method.

---

## 13. Local Thin Plate Spline (TPS)

### 13.1 Control set preparation

If the filtered point count exceeds `Max TPS Points`, the control set is reduced in two stages:

1. a coarse spatial preselection by grid cell;
2. a farthest-point-like deterministic diversification stage.

This reduces the TPS control count while preserving coverage.

### 13.2 Model construction

The TPS model stores:

- the control points;
- a spatial hash index;
- a support-mask or hull proxy;
- a target neighborhood size;
- a learned suggested regularization parameter;
- holdout validation metrics.

The target neighborhood size is chosen as

```math
N_{\mathrm{TPS}}
=
\min\!\left(
\max\!\left(64,\ 3\sqrt{N_c}\right),
\ \min(160, N_c)
\right),
```

where `N_c` is the control-point count.

### 13.3 Holdout learning of `\lambda`

The code performs deterministic holdout validation on a subset of control points, builds a training model without those points, predicts the holdouts, and records local `\lambda` values returned by the TPS solver. The median of the learned `\lambda` values becomes the model-level suggested regularization.

### 13.4 TPS kernel

For `r^2 <= 10^{-24}`, the thin-plate kernel is zero.

For `r^2 > 10^{-24}`, the thin-plate kernel is

```math
\phi(r^2) = r^2 \log(r^2).
```

### 13.5 Local patch system

For a local patch with `m` selected control points, the TPS system has dimension

```math
N = m + 3.
```

The affine tail is

```math
a_0 + a_1u + a_2v.
```

Instead of writing the full block matrix, the local TPS system can be expressed equivalently as

```math
(\mathbf{K} + \lambda\mathbf{I})\mathbf{w} + \mathbf{P}\mathbf{a} = \mathbf{z},
```

```math
\mathbf{P}^{\mathsf T}\mathbf{w} = \mathbf{0}.
```

Here,

```math
K_{ij} = \phi\!\left((u_i-u_j)^2 + (v_i-v_j)^2\right),
```

and

`\mathbf{P}` is the `m x 3` matrix whose `i`-th row is `(1, u_i, v_i)`.

### 13.6 Patch evaluation

For a query location `(u_q, v_q)`, the local TPS patch evaluates

```math
z_q =
a_0 + a_1u_q + a_2v_q
+ \sum_{i=1}^{m} w_i
\phi\!\left((u_i-u_q)^2 + (v_i-v_q)^2\right).
```

### 13.7 Local cross-validation for `\lambda`

Within each candidate local patch, the code tries the regularization set

```math
\lambda \in \{10^{-12},10^{-11},10^{-10},10^{-9},10^{-8},10^{-6}\}.
```

A subset of local points is reserved for validation whenever the neighborhood is sufficiently large. The selected `\lambda` is the one minimizing the validation score

```math
S = \mathrm{RMSE} + 0.35\,P95 + 0.15\,\mathrm{MAE}.
```

### 13.8 Blended multi-anchor evaluation

The query point uses up to four nearby anchor patches. For each patch, let `d` be the anchor-query distance and `r_s` the support scale. The blending weight is

```math
w_{\mathrm{blend}} =
\max\!\left(
\frac{w_c(d/r_s)}{10^{-8} + \mathrm{residual}},
\,10^{-12}
\right).
```

The final value is

```math
z = \frac{\sum_a w_a z_a}{\sum_a w_a},
```

if at least one patch succeeds. Otherwise the code falls back to IDW.

---

## 14. Support mask, hull awareness, and extrapolation control

The code uses a hybrid support-awareness structure:

1. a convex boundary polygon built from the clean point set;
2. a dilated occupancy mask on a regular support grid.

This is used as an inexpensive practical approximation of a concave support mask.

### 14.1 Occupancy grid spacing

The support-mask cell size is

```math
c_s = \max(2.5\bar{\Delta}_{xy},\,10^{-6}).
```

### 14.2 Support occupancy

For a point `(x,y)`, the support mask stores dilated occupied cells around

```math
i_x = \left\lfloor \frac{x-x_{\min}}{c_s} \right\rfloor,
\qquad
i_y = \left\lfloor \frac{y-y_{\min}}{c_s} \right\rfloor.
```

### 14.3 Outside-support distance

If a query point lies outside the support mask or hull proxy, a support distance is computed and converted into a boundary penalty

```math
\beta =
\frac{1}{1 + d_{\mathrm{out}}/L},
```

where `L` is a characteristic local scale such as grid spacing or spatial-hash cell size.

This penalty enters the local confidence score.

---

## 15. Confidence diagnostics and calibration

The code writes a dedicated `confidence.xyz` product. For each grid node it records:

- `x`
- `y`
- `z`
- confidence
- residual
- selected `\lambda` (TPS only; otherwise often zero or nominal)
- boundary penalty
- fallback flag
- outside-hull flag
- clamped-local flag
- neighborhood count
- patch count

The file header is:

```text
# x y z confidence residual lambda boundaryPenalty fallback outsideHull clampedLocal neighborCount patchCount
```

### 15.1 Confidence calibration

The raw method confidence is not used directly. It is calibrated against observed validation error using

```math
C^{*}
=
C_{\mathrm{raw}}
\cdot
\frac{1}{1 + r_{\mathrm{loc}}/s_{\mathrm{obs}}}.
```

The reported confidence is then clipped to the interval `[0.02, 1]`.

The observation scale is

```math
s_{\mathrm{obs}} =
\max\!\left(\mathrm{RMSE}_{\mathrm{obs}},\ 0.5P95_{\mathrm{obs}},\ 10^{-8}\right).
```

This ties the reported confidence to predictive error evidence rather than only internal residual surrogates.

---

## 16. Dense linear solver

The code uses a custom dense linear solver based on Gaussian elimination with partial pivoting.

At elimination stage `k`, the pivot row is chosen from

```math
\max_{r \ge k}|a_{rk}|.
```

If the pivot magnitude is too small or non-finite, the system is rejected.

The elimination multiplier is

```math
\gamma_{rk} = \frac{a_{rk}}{a_{kk}}.
```

The row update is

```math
a_{rc} \leftarrow a_{rc} - \gamma_{rk}a_{kc},
```

```math
b_r \leftarrow b_r - \gamma_{rk}b_k.
```

Back substitution recovers

```math
x_i =
\frac{1}{a_{ii}}
\left(
b_i - \sum_{c=i+1}^{n-1} a_{ic}x_c
\right).
```

This solver is reused by:

- the weighted local plane fit;
- the robust MLS fit;
- the local TPS patch solves.

---

## 17. Streamed grid and DXF generation

### 17.1 Streamed grid writing

The grid writer traverses the regular grid row by row and writes `x y z` records directly to disk. This avoids holding the entire grid product in memory.

### 17.2 Confidence stream

The confidence writer traverses the same grid and writes, per node,

```math
(x,y,z,C,r,\lambda,\beta,\mathrm{flags},N_{\mathrm{nbr}},N_{\mathrm{patch}}).
```

It also accumulates a global summary:

- mean confidence;
- mean residual;
- fallback count;
- outside-support count;
- clamped-local count.

### 17.3 DXF viewport geometry

The DXF writer constructs the drawing center as

```math
x_c = \frac{x_{\min}+x_{\max}}{2},
\qquad
y_c = \frac{y_{\min}+y_{\max}}{2},
```

and the view size as

```math
V = 1.1\max(x_{\max}-x_{\min},\ y_{\max}-y_{\min}).
```

If `V` is not finite and positive, the fallback is

```math
V = 1.
```

### 17.4 DXF layers

The writer creates layers:

- `xyz_points`
- `xyz_labels`
- `grid_points`
- `grid_labels`

Each retained input point and each grid point can be exported as both `POINT` and `TEXT` entities.

---

## 18. Main execution pipeline

The active runtime pipeline in `processXYZtoDXF(...)` consists of six user-visible phases:

1. deterministic input reading and `minDist` filtering;
2. robust local residual outlier rejection;
3. writing `filtered.xyz`;
4. building the interpolation context;
5. writing `grid.xyz` and `confidence.xyz`;
6. writing `dxf`.

The active outlier settings are

```math
r_n = \max(5d_{\min},\,0.01),
\qquad
\alpha = 3.5.
```

The method labels used in the report are:

- **Bicubic Hermite + robust MLS**
- **TPS**
- **Clamped local MLS (bounded)**

The code also reports predictive validation metrics for the selected interpolation branch whenever available.

---

## 19. GUI architecture

The Win32 GUI is organized into four grouped panels:

1. **Source File**
2. **Processing Parameters**
3. **Interpolation Method**
4. **Execution**

The window is fixed-size and non-resizable. The GUI uses `Segoe UI` fonts, dedicated group boxes, a file-browse button, parameter edit boxes, method radio buttons, a prominent `Run Conversion` button, and a status box that is updated asynchronously by the worker thread.

The method descriptions are deliberately textual because the interpolation branches have materially different numerical philosophies:

- smooth Hermite surface reconstruction from robust local nodal derivatives;
- locally adaptive radial-basis patching for scattered data;
- bounded local MLS for conservative overshoot control.

---

## 20. Output products

The program writes the following files:

1. `input.xyz.filtered.xyz`
2. `input.xyz.grid.xyz`
3. `input.xyz.confidence.xyz`
4. `input.xyz.dxf`
5. `input.xyz.rpt.txt`

### 20.1 `filtered.xyz`

Contains the deterministically thinned and robustly de-outliered point cloud.

### 20.2 `grid.xyz`

Contains the final interpolated regular grid.

### 20.3 `confidence.xyz`

Contains the interpolated grid plus local diagnostics.

### 20.4 `dxf`

Contains DXF entities for filtered points and optionally grid points.

### 20.5 `rpt.txt`

Contains a textual report or log of the run, including method label, selected tuning, validation metrics, confidence summary, and phase progress messages.

---

## 21. Practical method selection guidance

### 21.1 Bicubic Hermite + robust MLS

Recommended when:

- a smooth regular surface is desired;
- local derivative quality matters;
- the surface is reasonably continuous;
- a structured grid product is the primary target.

### 21.2 Local TPS

Recommended when:

- the point cloud is strongly scattered;
- point spacing is irregular;
- local radial-basis flexibility is more important than direct derivative control.

### 21.3 Clamped local MLS (bounded)

Recommended when:

- overshoot control is more important than maximal smoothness;
- local values must remain inside the neighborhood data envelope;
- a conservative interpolated surface is preferred.

---

## 22. Numerical safeguards and failure philosophy

The code is designed to **fail safely** rather than silently return unstable algebraic results. The following protections are implemented:

- dense solves reject nearly singular pivots;
- MLS and TPS local fits check finiteness;
- insufficient neighborhoods trigger conservative fallback;
- boundary penalties reduce confidence outside support;
- derivative limiters suppress excessive Hermite oscillation;
- output streams are written incrementally with explicit error checks;
- allocation failures are caught and reported.

The overall philosophy is therefore **robustness first, smoothness second, and silent numerical optimism never**.

---

## 23. Limitations and intended use

This software is a strong engineering interpolator and exporter, but it remains a deterministic local-surface tool rather than a universal geostatistical model. In particular:

- it does not estimate kriging variance;
- it does not explicitly detect semantic breaklines;
- its support mask is a practical hull or occupancy surrogate, not a full computational-geometry alpha-shape library;
- the clamped MLS mode is bounded, but it is not a formal TVD finite-volume scheme.

Within these limits, the code is nevertheless highly suitable for technical surface reconstruction from engineering XYZ clouds.

---

## 24. Summary

The active code base implements a complete and scientifically defensible workflow for transforming irregular XYZ point clouds into quality-controlled regular grids and DXF visualization products. Its principal scientific ingredients are:

- deterministic minimum-distance thinning;
- robust local residual-based outlier rejection;
- anisotropic local-coordinate construction;
- robust MLS with predictive validation;
- bicubic Hermite cell reconstruction with derivative limiting;
- local TPS patches with cross-validated regularization;
- bounded local MLS as a conservative alternative;
- support-aware confidence diagnostics;
- streamed engineering outputs.

The software is therefore best understood as a **numerically robust surface reconstruction and export framework**, rather than a simple format converter.
