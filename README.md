# XYZ to DXF Converter GUI

## Abstract

The **XYZ to DXF Converter GUI** is a Win32 desktop application for the ingestion, conditioning, interpolation, diagnostic assessment, and CAD-oriented export of large irregular XYZ point clouds. The code is not a simple file-format converter. Its active numerical role is the construction of a **cleaned, regularized, diagnostically traceable surface representation** from scattered planar samples

$$
\mathcal{P} = \{(x_i,y_i,z_i)\}_{i=1}^{N},
$$

where $z$ is typically elevation, depth, or another scalar field defined over $(x,y)$.

The current C++ implementation exposes three interpolation branches:

1. **Bicubic Hermite + robust anisotropic MLS**;
2. **Local TPS** with adaptive local patches and local cross-validated regularization;
3. **Clamped local MLS (bounded)**.

The active code also implements:

- robust local residual outlier rejection on the **original** XYZ cloud;
- deterministic minimum-distance thinning applied **after** optional outlier rejection;
- raster-occupancy outer-boundary construction and clipping;
- predictive holdout validation for MLS and TPS branches;
- streamed writing of filtered XYZ, interpolated grid XYZ, confidence diagnostics, report, and DXF;
- a fixed-size Win32 GUI with checkboxes for outlier removal and outer-boundary usage.

**Keywords:** scattered-data interpolation, moving least squares, bicubic Hermite interpolation, thin plate spline, robust regression, point-cloud conditioning, DXF export, engineering surface reconstruction.

---

## 1. Scope and technical objective

The active processing chain transforms the raw point set into a hierarchy of working products:

$$
\mathcal{P}
\rightarrow
\mathcal{P}_{\mathrm{outlier\ filtered}}
\rightarrow
\mathcal{P}_{d_{\min}}
\rightarrow
\mathcal{G}
\rightarrow
\{\mathrm{filtered.xyz},\ \mathrm{grid.xyz},\ \mathrm{confidence.xyz},\ \mathrm{dxf},\ \mathrm{rpt.txt}\}.
$$

The stage meanings are:

- $\mathcal{P} \rightarrow \mathcal{P}_{\mathrm{outlier\ filtered}}$: optional robust local residual screening on the original XYZ cloud;
- $\mathcal{P}_{\mathrm{outlier\ filtered}} \rightarrow \mathcal{P}_{d_{\min}}$: deterministic horizontal thinning or exact duplicate suppression;
- $\mathcal{P}_{d_{\min}} \rightarrow \mathcal{G}$: regular-grid interpolation by the selected method;
- $\mathcal{G} \rightarrow$ exports: streamed writing of the surface, diagnostics, and CAD representation.

The technical objectives of the code are:

- to reduce redundant sampling without introducing file-order dependence;
- to reject locally inconsistent vertical anomalies without using a crude global $z$ cutoff;
- to interpolate a regular Cartesian grid suitable for engineering inspection and downstream numerical use;
- to attach practical confidence and fallback diagnostics to each predicted node;
- to produce traceable exported products rather than an opaque black-box result.

---

## 2. High-level software workflow

The runtime pipeline implemented in `processXYZtoDXF(...)` is:

1. **Read original XYZ input**;
2. **Optional robust outlier removal on the original cloud**;
3. **Deterministic `minDist` thinning**;
4. **Optional outer-boundary construction**;
5. **Write `filtered.xyz`**;
6. **Build interpolation context and stream `grid.xyz` + `confidence.xyz`**;
7. **Stream `dxf` and write the textual report**.

For a successful run, the program therefore behaves as a coupled **conditioning + interpolation + diagnostics + export** system rather than a single interpolation call.

---

## 3. Graphical interface and exposed parameters

The GUI exposes the following active controls:

- `Input XYZ file`
- `Min Dist`
- `Precision`
- `PDMODE`
- `Grid Spacing`
- `Use robust outlier removal`
- `Use outer boundary`
- interpolation method radio buttons

The default GUI values created by the active code are:

$$
d_{\min}=5,
\qquad
\mathrm{Precision}=2,
\qquad
\mathrm{PDMODE}=3,
\qquad
h=20.
$$

There is **no active `Max TPS Points` GUI parameter** in the current code. The TPS branch uses the full filtered point cloud as its control set.

### 3.1 Minimum distance

If $d_{\min} > 0$, the program performs deterministic horizontal thinning. If

$$
d_{\min} \le 0,
$$

it falls back to exact duplicate removal only.

### 3.2 Precision

`Precision` controls only text formatting in exported XYZ and DXF files. Internal arithmetic is performed in `double` precision throughout.

### 3.3 PDMODE

`PDMODE` is written into the DXF header variable `$PDMODE`. The DXF writer also emits the fixed point-size constant

$$
\mathrm{PDSIZE}=0.5.
$$

### 3.4 Grid spacing

The regular interpolation spacing is the Cartesian grid increment

$$
h > 0.
$$

This parameter governs the grid extent, node count, Hermite cell size, derivative limiting scales, and boundary-penalty length scale in the bicubic and clamped MLS branches.

### 3.5 Method options and toggles

The method selector chooses one of:

- **Bicubic Hermite + robust MLS**;
- **TPS**;
- **Clamped local MLS (bounded)**.

The checkboxes enable or disable:

- robust outlier removal on the original XYZ cloud;
- outer-boundary construction and clipping of grid and DXF outputs.

---

## 4. Deterministic minimum-distance filtering

### 4.1 Motivation

The code deliberately avoids a first-come-first-kept thinning rule. Instead, it constructs a deterministic candidate ordering so that the retained cloud does not depend on the row order of the source file.

### 4.2 Exact duplicate mode

If

$$
d_{\min} \le 0,
$$

the code removes only exact duplicates under exact `double` comparison:

$$
(x_i,y_i,z_i)=(x_j,y_j,z_j).
$$

### 4.3 Cell assignment

For $d_{\min}>0$, each point is assigned to a cell

$$
(i_x,i_y)=
\left(
\left\lfloor \frac{x}{d_{\min}} \right\rfloor,
\left\lfloor \frac{y}{d_{\min}} \right\rfloor
\right).
$$

The corresponding cell center is

$$
x_c = \left(i_x+\frac12\right)d_{\min},
\qquad
y_c = \left(i_y+\frac12\right)d_{\min}.
$$

Each candidate receives the deterministic score

$$
s=(x-x_c)^2+(y-y_c)^2.
$$

The candidate list is stably sorted by cell index, then by $s$, then by coordinate tie-breaks.

### 4.4 Acceptance criterion

Let $p_i=(x_i,y_i,z_i)$ be a candidate point. The code inspects accepted points in the surrounding $3\times 3$ neighborhood of sparse cells. The candidate is rejected if any previously accepted point $p_j$ satisfies

$$
(x_i-x_j)^2+(y_i-y_j)^2 < d_{\min}^2.
$$

Otherwise it is accepted.

### 4.5 Properties

The implemented thinning rule is therefore:

- deterministic;
- local in memory access;
- horizontal rather than full 3D distance based;
- suitable as a stable preconditioner for the subsequent interpolation stages.

---

## 5. Robust local outlier rejection

### 5.1 Position in the active pipeline

In the active code, robust outlier removal is applied **before** deterministic `minDist` thinning. This is materially important, because the residual test is evaluated on the original cloud rather than on an already-thinned subset.

### 5.2 Neighborhood radius used at runtime

The runtime pipeline computes the outlier search distance as

$$
r_n = \max\left(5d_{\min},\ 3\bar{\Delta}_{xy},\ 0.01\right),
$$

where $\bar{\Delta}_{xy}$ is the average planar spacing defined later.

The robust threshold factor used by the code is

$$
\alpha = 3.5.
$$

### 5.3 Local weighted plane fit

For a tested point $(x_0,y_0,z_0)$, the program collects neighboring points within radius $r_n$ and fits the local plane

$$
z(x,y) \approx a + b_x(x-x_0) + b_y(y-y_0).
$$

The basis is

$$
\phi(x,y)=\begin{bmatrix}1 & x-x_0 & y-y_0\end{bmatrix}^{T}.
$$

The local scale is

$$
h = \max\left(\sqrt{\frac1n\sum_{k=1}^{n}\big[(x_k-x_0)^2+(y_k-y_0)^2\big]},\ 10^{-6}\right).
$$

The geometric weights are

$$
w_k = \frac{1}{1 + \dfrac{(x_k-x_0)^2+(y_k-y_0)^2}{h^2}}.
$$

The weighted normal equations are

$$
(\mathbf{A}^{T}\mathbf{W}\mathbf{A})\,\mathbf{c} = \mathbf{A}^{T}\mathbf{W}\mathbf{z},
$$

with

$$
\mathbf{c}=\begin{bmatrix}a & b_x & b_y\end{bmatrix}^{T}.
$$

The code adds a tiny diagonal stabilization:

$$
\mathbf{A}^{T}\mathbf{W}\mathbf{A}
\leftarrow
\mathbf{A}^{T}\mathbf{W}\mathbf{A}+10^{-12}\mathbf{I}.
$$

### 5.4 Neighborhood sufficiency logic

If fewer than 6 usable neighbors remain after range filtering, the point is conservatively retained rather than aggressively removed.

### 5.5 Residual center and robust scale

After fitting the plane, the neighborhood residuals are

$$
r_k = z_k - \big(a + b_x(x_k-x_0) + b_y(y_k-y_0)\big).
$$

The residual center is the median

$$
\widetilde{r}=\mathrm{median}(r_k).
$$

The median absolute deviation is

$$
\mathrm{MAD} = \mathrm{median}\big(|r_k-\widetilde{r}|\big).
$$

The robust scale is

$$
\sigma_{\mathrm{rob}} = \max\big(1.4826\,\mathrm{MAD},\ 10^{-8}\big).
$$

### 5.6 Pointwise decision rule

The code evaluates the tested point using

$$
r_{\mathrm{self}} = \left|z_0-a-\widetilde{r}\right|.
$$

The point is kept when

$$
r_{\mathrm{self}} \le \alpha\,\sigma_{\mathrm{rob}}.
$$

Otherwise it is rejected.

---

## 6. Grid construction

After conditioning, the code builds a regular Cartesian grid over an expanded bounding box of the filtered point cloud.

Let the filtered-data bounds be

$$
x_{\min}^{d},\ x_{\max}^{d},\ y_{\min}^{d},\ y_{\max}^{d}.
$$

The margin is

$$
m = 1.5h.
$$

Hence the grid extent is

$$
x_{\min}=x_{\min}^{d}-m,
\qquad
x_{\max}=x_{\max}^{d}+m,
$$

$$
y_{\min}=y_{\min}^{d}-m,
\qquad
y_{\max}=y_{\max}^{d}+m.
$$

The code then forms

$$
W = \max(x_{\max}-x_{\min}, 1),
\qquad
H = \max(y_{\max}-y_{\min}, 1),
$$

$$
n_x = \left\lceil \frac{W}{h} \right\rceil + 1,
\qquad
n_y = \left\lceil \frac{H}{h} \right\rceil + 1.
$$

The node coordinates are

$$
x_i = x_{\min} + ih,
\qquad
y_j = y_{\min} + jh.
$$

The planned node count is

$$
N_g = n_x n_y,
$$

with explicit overflow checks before multiplication.

---

## 7. Auxiliary geometric quantities

### 7.1 Average planar spacing

The code repeatedly uses the characteristic spacing

$$
\bar{\Delta}_{xy} = \sqrt{\frac{A}{N}},
$$

where

$$
A = \max\big((x_{\max}-x_{\min})(y_{\max}-y_{\min}),\ 10^{-12}\big).
$$

This is used to size spatial-hash cells, boundary raster cells, and local support scales.

### 7.2 Inverse-distance fallback

Whenever a preferred local interpolant fails, the code falls back to inverse-distance weighting over up to 16 nearest points:

$$
z(x,y) = \frac{\sum_{k=1}^{m} w_k z_k}{\sum_{k=1}^{m} w_k},
$$

with

$$
w_k = \frac{1}{10^{-12} + (x-x_k)^2 + (y-y_k)^2}.
$$

This fallback is intentionally conservative and is explicitly recorded in the diagnostics.

---

## 8. Local anisotropic coordinate system

The MLS and TPS branches both build a local anisotropic principal frame.

### 8.1 Preliminary scale

For query point $(x_0,y_0)$ and neighborhood $\mathcal{N}$,

$$
\bar{r}^2 = \frac1n\sum_{k=1}^{n}\big[(x_k-x_0)^2+(y_k-y_0)^2\big],
$$

$$
h_0 = \max\big(\sqrt{\max(\bar{r}^2,10^{-12})},\ \max(h,10^{-6})\big).
$$

### 8.2 Preliminary weights

The covariance weights are

$$
w_k^{\mathrm{frame}} = \frac{1}{1 + \dfrac{(x_k-x_0)^2+(y_k-y_0)^2}{h_0^2}}.
$$

### 8.3 Weighted covariance matrix

The weighted second moments are

$$
S_{xx} = \frac{1}{W}\sum_k w_k^{\mathrm{frame}}(x_k-x_0)^2,
$$

$$
S_{yy} = \frac{1}{W}\sum_k w_k^{\mathrm{frame}}(y_k-y_0)^2,
$$

$$
S_{xy} = \frac{1}{W}\sum_k w_k^{\mathrm{frame}}(x_k-x_0)(y_k-y_0),
$$

with

$$
W = \sum_k w_k^{\mathrm{frame}}.
$$

### 8.4 Principal directions

The code computes

$$
\mathrm{tr}(\mathbf{S}) = S_{xx}+S_{yy},
$$

$$
\Delta = \sqrt{(S_{xx}-S_{yy})^2 + 4S_{xy}^2},
$$

$$
\lambda_1 = \max\left(\frac{\mathrm{tr}(\mathbf{S})+\Delta}{2},0\right),
\qquad
\lambda_2 = \max\left(\frac{\mathrm{tr}(\mathbf{S})-\Delta}{2},0\right).
$$

The frame angle is

$$
\theta = \frac12\mathrm{atan2}(2S_{xy}, S_{xx}-S_{yy}).
$$

The orthonormal axes are

$$
\mathbf{e}_1=(\cos\theta,\sin\theta),
\qquad
\mathbf{e}_2=(-\sin\theta,\cos\theta).
$$

### 8.5 Local scales

The anisotropic scaling lengths are

$$
s_1 = \max\big(\sqrt{\lambda_1},\ 0.75\max(h,10^{-6})\big),
$$

$$
s_2 = \max\big(\sqrt{\lambda_2},\ 0.35\max(h,10^{-6})\big).
$$

### 8.6 Local coordinates

For any point $(x,y)$,

$$
u = \frac{(x-x_0)e_{1x}+(y-y_0)e_{1y}}{s_1},
\qquad
v = \frac{(x-x_0)e_{2x}+(y-y_0)e_{2y}}{s_2}.
$$

The normalized local radius is

$$
\rho = \sqrt{u^2+v^2}.
$$

---

## 9. Adaptive neighborhood support

For a candidate neighborhood, the code computes the transformed radii $\rho_k$ and sorts them. The support radius is then chosen as

$$
r_s = \max\big(r_{\min},\ 1.5,\ \rho_{80},\ 0.75\rho_{90}\big),
$$

where $\rho_{80}$ and $\rho_{90}$ are the 80th and 90th percentile radii.

The final adaptive neighborhood is selected by retaining all candidates with

$$
\rho_k \le r_s,
$$

while enforcing minimum and maximum cardinality constraints. If too few points satisfy the support-radius test, the nearest-ranked points are retained until the minimum count is met.

---

## 10. Robust anisotropic MLS

### 10.1 Basis order

The MLS branch supports two polynomial orders:

- quadratic: 6 coefficients;
- cubic: 10 coefficients.

The active helper function is therefore

$n_b = 6$ for $\mathrm{basisOrder}=2$, and $n_b = 10$ for $\mathrm{basisOrder}\ge 3$.

### 10.2 Polynomial basis

The implemented basis is

$$
\phi_0 = 1,
\quad
\phi_1 = u,
\quad
\phi_2 = v,
$$

$$
\phi_3 = \frac12u^2,
\quad
\phi_4 = uv,
\quad
\phi_5 = \frac12v^2.
$$

If cubic mode is enabled, the code appends

$$
\phi_6 = \frac16u^3,
\quad
\phi_7 = \frac12u^2v,
\quad
\phi_8 = \frac12uv^2,
\quad
\phi_9 = \frac16v^3.
$$

The fitted polynomial value is

$$
f(u,v) = \sum_{k=0}^{n_b-1} c_k\phi_k(u,v).
$$

### 10.3 Compact support kernel

The geometric support kernel is the compact Wendland weight

$w_c(q) = (1-q)^4(4q+1)$ for $q<1$, and $w_c(q)=0$ for $q\ge 1$.

with

$$
q = \frac{\rho}{r_s}.
$$

The final geometric factor used in the MLS fit is

$$
w_{\mathrm{geom}} = \max\left(\frac{w_c(\rho/r_s)}{1+\rho^2},\ 10^{-14}\right).
$$

### 10.4 Robust reweighting

The fit uses three iteratively reweighted passes. After each pass, residuals are formed, their median-centered MAD scale is computed,

$$
\sigma = \max(1.4826\,\mathrm{MAD}, 10^{-8}),
$$

and a Tukey-type reweighting is applied with

$$
t = \frac{|z_i-f_i|}{4.685\sigma}.
$$

The robust weight update is

$w_i^{\mathrm{rob}} = 10^{-4}$ for $t \ge 1$, and $w_i^{\mathrm{rob}} = \max\big((1-t^2)^2,\ 10^{-4}\big)$ for $t < 1$.

### 10.5 Regularized normal equations

The local normal equations are assembled as

$$
\mathbf{A}^{T}\mathbf{W}\mathbf{A}\,\mathbf{c} = \mathbf{A}^{T}\mathbf{W}\mathbf{z},
$$

with combined weight

$$
w_i = \max\big(w_{\mathrm{geom},i}\,w_i^{\mathrm{rob}},\ 10^{-14}\big).
$$

The diagonal ridge is scaled from the matrix trace:

$$
\lambda_{\mathrm{ridge}} = \max\left(10^{-12},\ \mathrm{ridgeFactor}\times 10^{-11}\frac{\mathrm{tr}(\mathbf{A}^{T}\mathbf{W}\mathbf{A})}{n_b}\right).
$$

Hence the solved system is

$$
\big(\mathbf{A}^{T}\mathbf{W}\mathbf{A}+\lambda_{\mathrm{ridge}}\mathbf{I}\big)\mathbf{c} = \mathbf{A}^{T}\mathbf{W}\mathbf{z}.
$$

### 10.6 Validation metrics

For MLS and TPS, the code uses the metrics

$$
\mathrm{RMSE} = \sqrt{\frac1n\sum_{i=1}^{n} e_i^2}.
$$

$$
\mathrm{MAE} = \frac1n\sum_{i=1}^{n}|e_i|.
$$

$P95$ denotes the 95th percentile of $|e_i|$.

The scalar tuning score is

$$
S = \mathrm{RMSE} + 0.35\,P95 + 0.15\,\mathrm{MAE}.
$$

### 10.7 Predictive local validation

For local MLS node estimation, the code performs deterministic local holdout checks. If those checks fail to produce a valid metric set, it falls back to the fit residual metrics.

---

## 11. Bicubic Hermite + robust MLS

### 11.1 Nodal quantities

At each grid node, the MLS fit provides

$$
z = c_0,
\qquad
f_u = c_1,
\qquad
f_v = c_2,
\qquad
f_{uu} = c_3,
\qquad
f_{uv} = c_4,
\qquad
f_{vv} = c_5.
$$

### 11.2 Chain-rule mapping to global coordinates

The local-to-global derivatives used by the Hermite patch are

$$
\frac{\partial u}{\partial x}=\frac{e_{1x}}{s_1},
\qquad
\frac{\partial u}{\partial y}=\frac{e_{1y}}{s_1},
$$

$$
\frac{\partial v}{\partial x}=\frac{e_{2x}}{s_2},
\qquad
\frac{\partial v}{\partial y}=\frac{e_{2y}}{s_2}.
$$

The code then computes

$$
f_x = f_u\frac{\partial u}{\partial x} + f_v\frac{\partial v}{\partial x},
$$

$$
f_y = f_u\frac{\partial u}{\partial y} + f_v\frac{\partial v}{\partial y},
$$

$$
\begin{aligned}
f_{xy}
&= f_{uu}\,\frac{\partial u}{\partial x}\,\frac{\partial u}{\partial y}
\\[4pt]
&\quad + f_{uv}\!\left(
\frac{\partial u}{\partial x}\,\frac{\partial v}{\partial y}
+
\frac{\partial v}{\partial x}\,\frac{\partial u}{\partial y}
\right)
\\[4pt]
&\quad + f_{vv}\,\frac{\partial v}{\partial x}\,\frac{\partial v}{\partial y}.
\end{aligned}
$$

### 11.3 Nodewise derivative limiting

Before Hermite reconstruction, the nodal derivatives are bounded using the local $z$ range

$$
\Delta z = \max(z_{\max}^{\mathrm{loc}}-z_{\min}^{\mathrm{loc}}, 10^{-8}).
$$

The nodewise slope and twist caps are

$$
L_s = \frac{3\Delta z}{\max(h,10^{-6})},
\qquad
L_t = \frac{4\Delta z}{\max(h^2,10^{-6})}.
$$

The code clamps

$$
f_x \leftarrow \mathrm{clip}(f_x,-L_s,L_s),
\qquad
f_y \leftarrow \mathrm{clip}(f_y,-L_s,L_s),
$$

$$
f_{xy} \leftarrow \mathrm{clip}(f_{xy},-L_t,L_t).
$$

### 11.4 Cellwise minmod derivative limiting

For a cell with corner values $z_{00},z_{10},z_{01},z_{11}$ and spacings $h_x,h_y$, the finite-difference slopes are

$$
s_{x0} = \frac{z_{10}-z_{00}}{\max(h_x,10^{-6})},
\qquad
s_{x1} = \frac{z_{11}-z_{01}}{\max(h_x,10^{-6})},
$$

$$
s_{y0} = \frac{z_{01}-z_{00}}{\max(h_y,10^{-6})},
\qquad
s_{y1} = \frac{z_{11}-z_{10}}{\max(h_y,10^{-6})}.
$$

The limiter is

$\mathrm{minmod}(a,b)=0$ when $ab \le 0$, and $\mathrm{minmod}(a,b)=\mathrm{sign}(a)\min(|a|,|b|)$ when $ab>0$.

The mixed-derivative baseline is

$$
c_{\times} = \frac12\left(
\frac{s_{y1}-s_{y0}}{\max(h_x,10^{-6})}
+
\frac{s_{x1}-s_{x0}}{\max(h_y,10^{-6})}
\right),
$$

with cap

$$
L_{\times} = 2|c_{\times}| + 10^{-8}.
$$

### 11.5 Hermite basis functions

The one-dimensional cubic Hermite functions are

$$
h_{00}(t)=2t^3-3t^2+1,
\qquad
h_{10}(t)=t^3-2t^2+t,
$$

$$
h_{01}(t)=-2t^3+3t^2,
\qquad
h_{11}(t)=t^3-t^2.
$$

### 11.6 Bicubic cell surface

For normalized local coordinates $(u,v)\in[0,1]^2$, the implemented patch is

$$
\begin{aligned}
f(u,v) ={}&
 z_{00}h_{00}(u)h_{00}(v)
+z_{10}h_{01}(u)h_{00}(v)
+z_{01}h_{00}(u)h_{01}(v)
+z_{11}h_{01}(u)h_{01}(v) \\
&+ h_x\Big(f_{x,00}h_{10}(u)h_{00}(v)+f_{x,10}h_{11}(u)h_{00}(v)
+f_{x,01}h_{10}(u)h_{01}(v)+f_{x,11}h_{11}(u)h_{01}(v)\Big) \\
&+ h_y\Big(f_{y,00}h_{00}(u)h_{10}(v)+f_{y,10}h_{01}(u)h_{10}(v)
+f_{y,01}h_{00}(u)h_{11}(v)+f_{y,11}h_{01}(u)h_{11}(v)\Big) \\
&+ h_xh_y\Big(f_{xy,00}h_{10}(u)h_{10}(v)+f_{xy,10}h_{11}(u)h_{10}(v)
+f_{xy,01}h_{10}(u)h_{11}(v)+f_{xy,11}h_{11}(u)h_{11}(v)\Big).
\end{aligned}
$$

### 11.7 Final value envelope clamp

After evaluating the bicubic patch, the result is softly clamped to the four-corner value envelope. Define

$$
z_{\min}=\min(z_{00},z_{10},z_{01},z_{11}),
\qquad
z_{\max}=\max(z_{00},z_{10},z_{01},z_{11}),
$$

$$
p = 0.05\max(z_{\max}-z_{\min},10^{-8}).
$$

The final bicubic value is limited by

$$
f \leftarrow \max(z_{\min}-p,\ \min(z_{\max}+p, f)).
$$

### 11.8 Parameter tuning

The bicubic/MLS tuning search uses the exact active option sets

$$
\mathrm{neighborhoodTarget}\in\{48,64,96,128\},
$$

$$
\mathrm{supportMultiplier}\in\{1.75,2.25,2.75\},
$$

$$
\mathrm{ridgeFactor}\in\{1,10,50\},
$$

$$
\mathrm{basisOrder}\in\{2,3\}.
$$

Selection is based on the smallest tuning score $S$ from deterministic global holdout validation.

---

## 12. Clamped local MLS (bounded)

This branch reuses the robust anisotropic MLS node evaluation at the query point, but replaces bicubic cell reconstruction by a direct local prediction followed by an envelope clamp.

For the local envelope extracted from up to 32 nearest neighbors,

$$
z_{\min}^{\mathcal{N}} = \min_{k\in\mathcal{N}} z_k,
\qquad
z_{\max}^{\mathcal{N}} = \max_{k\in\mathcal{N}} z_k.
$$

If $z_{\mathrm{MLS}}$ is the raw local MLS prediction, the final value is

$$
z^{*} = \max\left(z_{\min}^{\mathcal{N}},\ \min\left(z_{\max}^{\mathcal{N}}, z_{\mathrm{MLS}}\right)\right).
$$

This branch is therefore not a monotone finite-volume method. It is a **bounded local regression** method.

---

## 13. Local Thin Plate Spline (TPS)

### 13.1 Control set preparation

The active code uses **all filtered points** as TPS control points. There is no GUI-side control-count cap in the current implementation.

### 13.2 Model construction

Let $N_c$ be the number of TPS control points. The target TPS neighborhood size is

$$
N_{\mathrm{TPS}} = \min\left(\max\left(64,\ 3\sqrt{N_c}\right),\ \min(160,N_c)\right).
$$

The TPS spatial hash uses cell size

$$
c_{\mathrm{TPS}} = \max(3\bar{\Delta}_{xy},10^{-6}).
$$

The model-level initial regularization is

$$
\lambda_{0}=10^{-9}.
$$

### 13.3 Holdout learning of $\lambda$

The code builds a deterministic holdout subset of at most 32 control points. Predictions on the training-only model are used to compute holdout errors and collect locally selected $\lambda$ values. The model-level suggested regularization becomes the median learned $\lambda$.

### 13.4 TPS kernel

The radial basis kernel used by the code is

$\Phi(r^2)=0$ for $r^2 \le 10^{-24}$, and $\Phi(r^2)=r^2\log(r^2)$ for $r^2 > 10^{-24}$.

### 13.5 Local patch system

For a local patch with $m$ retained control points in local coordinates $(u_i,v_i)$, the code solves the standard TPS augmented system in the equivalent split form

$$
(\mathbf{K}+\lambda\mathbf{I})\mathbf{w}+\mathbf{P}\mathbf{a}=\mathbf{z},
$$

$$
\mathbf{P}^{T}\mathbf{w}=\mathbf{0},
$$

with

$$
K_{ij}=\Phi\big((u_i-u_j)^2+(v_i-v_j)^2\big).
$$

The polynomial design matrix $\mathbf{P}$ has rows

$$
\mathbf{P}_i = [1\ \ u_i\ \ v_i].
$$

### 13.6 Patch evaluation

At query coordinates $(u_q,v_q)$, the patch value is

$$
z_q = a_0 + a_1u_q + a_2v_q + \sum_{i=1}^{m} w_i\,\Phi\big((u_i-u_q)^2+(v_i-v_q)^2\big).
$$

### 13.7 Local cross-validation for $\lambda$

For each local patch, the code tests the exact set

$$
\lambda \in \{10^{-12},10^{-11},10^{-10},10^{-9},10^{-8},10^{-6}\}.
$$

Each candidate is scored by predictive residual metrics on a local validation split if available, otherwise by fit residual metrics. The smallest score $S$ wins.

### 13.8 Blended multi-anchor evaluation

The TPS evaluation blends up to four nearby anchor patches. If patch $a$ has support radius $r_{s,a}$, anchor distance $d_a$, and residual estimate $r_a$, the blend weight is

$$
w_a = \max\left(\frac{w_c\big(d_a/(1.001\,r_{s,a}^{xy})\big)}{\sqrt{10^{-8}+r_a}},\ 10^{-12}\right),
$$

where

$$
r_{s,a}^{xy}=\max\big(1,\ r_{s,a}\max(c_{\mathrm{TPS}},10^{-6})\big).
$$

The blended TPS value is

$$
z = \frac{\sum_a w_a z_a}{\sum_a w_a}.
$$

If no valid blend is available, the code falls back first to the best single patch and finally to inverse-distance weighting.

---

## 14. Support mask, hull awareness, and extrapolation control

The structure name `ConvexHull2D` is historical. The implemented algorithm is **not** a mathematical convex hull. It is a raster-occupancy outer-boundary model built from the filtered point cloud.

### 14.1 Occupancy grid spacing

The boundary raster cell size is

$c_s = \max(d_{\min},10^{-6})$ if a positive minimum-distance hint exists, and $c_s = \max(\bar{\Delta}_{xy},10^{-6})$ otherwise.

Inside `buildConvexHull2D(...)`, if no explicit request is supplied, the cell size defaults to

$$
c_s = \max(2.5\bar{\Delta}_{xy},10^{-6}).
$$

### 14.2 Occupancy support and loop extraction

The filtered points are rasterized into occupied cells. The code then applies morphological closing, fills interior holes, extracts boundary edges, traces loops, and keeps the largest loop by area as the active outer boundary.

When the boundary is built from the runtime pipeline, the cell size is iteratively relaxed by

$$
c_s \leftarrow 1.35\,c_s
$$

until all filtered points lie inside the resulting boundary or the iteration budget is exhausted.

### 14.3 Outside-support distance and penalty

For a query point $(x,y)$ outside the active boundary, the code computes the shortest distance to the polygon loop,

$$
d_{\mathrm{out}} = \mathrm{dist}\big((x,y),\partial\Omega\big).
$$

The boundary penalty used by bicubic, clamped MLS, and TPS diagnostics is

$$
\beta = \frac{1}{1 + d_{\mathrm{out}}/L},
$$

where the active length scale is:

- $L=h$ for bicubic Hermite and clamped MLS;
- $L=\max(c_{\mathrm{TPS}},10^{-6})$ for TPS.

---

## 15. Confidence diagnostics and calibration

The file `confidence.xyz` writes the fields

```text
x y z confidence residual lambda boundaryPenalty fallback outsideHull clampedLocal neighborCount patchCount
```

The interpretation is branch dependent, but the common calibration routine is implemented as follows.

### 15.1 Confidence calibration

If the raw confidence is invalid or non-positive, the code resets it to

$$
C_{\mathrm{raw}} = 0.05.
$$

If no observed validation metrics are available, the returned confidence is simply clamped:

$$
C = \mathrm{clip}(C_{\mathrm{raw}}, 0.05, 1).
$$

Otherwise the observed error scale is

$$
s_{\mathrm{obs}} = \max\big(\mathrm{RMSE},\ 0.5P95,\ 10^{-8}\big),
$$

and the calibrated confidence is

$$
C = \mathrm{clip}\left(C_{\mathrm{raw}}\frac{1}{1+r_{\mathrm{loc}}/s_{\mathrm{obs}}},\ 0.02, 1\right).
$$

The raw branch-specific formulas are:

- **Bicubic Hermite:** $C_{\mathrm{raw}} = \frac14(C_{00}+C_{10}+C_{01}+C_{11})\,\beta.$
- **Local MLS node estimate:** $C_{\mathrm{raw}} = \mathrm{clip}\left(\frac{1}{1+1.5S},\ 0.05,1\right).$
- **Clamped local MLS stream:** $C_{\mathrm{raw}} = C_{\mathrm{local}}\,\beta.$
- **TPS stream:** $C_{\mathrm{raw}} = \frac{\beta}{1+4r_{\mathrm{best}}}.$

The residual field written to `confidence.xyz` is the branch-local fit or validation residual estimate. The `lambda` field is active for TPS and zero for the MLS-based branches.

---

## 16. Dense linear solver

All local algebra is solved by an internal dense direct solver with partial pivoting. At elimination step $k$, the solver selects the pivot row with largest

$$
|a_{rk}|.
$$

If the pivot magnitude is not finite or does not satisfy

$$
|a_{kk}| > 10^{-18},
$$

the solve fails safely.

The elimination factor is

$$
\gamma_{rk} = \frac{a_{rk}}{a_{kk}},
$$

and the back-substitution stage computes

$$
x_i = \frac{1}{a_{ii}}\left(b_i - \sum_{c=i+1}^{n-1} a_{ic}x_c\right).
$$

This internal solver is used by the weighted plane fit, MLS normal equations, and TPS patch systems.

---

## 17. Streamed grid and DXF generation

### 17.1 Streamed grid writing

The code writes `grid.xyz` row by row. If outer-boundary usage is enabled, nodes outside the active boundary are skipped rather than written.

### 17.2 Confidence stream

The confidence writer streams the same grid node pattern and accumulates the global summary:

- node count;
- mean confidence;
- mean residual;
- fallback count;
- outside-hull count;
- clamped-local count.

### 17.3 DXF viewport geometry

The DXF writer computes the viewport center as

$$
x_c = \frac{x_{\min}+x_{\max}}{2},
\qquad
y_c = \frac{y_{\min}+y_{\max}}{2}.
$$

The view size is

$$
V = 1.1\max(x_{\max}-x_{\min},\ y_{\max}-y_{\min}).
$$

If $V$ is invalid or non-positive, the fallback is

$$
V=1.
$$

### 17.4 DXF layers

The writer creates the layers:

- `xyz_points`
- `xyz_labels`
- `boundary`
- `grid_points`
- `grid_labels` (when a grid is written)

The DXF therefore contains both point geometry and a practical boundary outline suitable for CAD inspection.

---

## 18. Main execution pipeline

The runtime report produced by the code records the following major stages:

1. read original XYZ input;
2. robust outlier removal on the original cloud, or explicit notice that it was disabled;
3. deterministic `minDist` thinning;
4. boundary construction, or explicit notice that it was disabled;
5. writing `filtered.xyz`;
6. interpolation-context construction and streamed writing of `grid.xyz` and `confidence.xyz`;
7. streamed DXF generation.

The method labels written by the active code are:

- **Bicubic Hermite + robust MLS**
- **TPS**
- **Clamped local MLS (bounded)**

For the bicubic branch, the report also records the selected tuning tuple

$(\mathrm{neighbors},\ \mathrm{supportMultiplier},\ \mathrm{ridgeFactor},\ \mathrm{basisOrder})$.

---

## 19. GUI architecture

The Win32 interface is organized into four group boxes:

1. **Source File**
2. **Processing Parameters**
3. **Interpolation Method**
4. **Execution**

The window is fixed-size and non-resizable. The active logical creation size is approximately

$$
950 \times 670
$$

window units before `AdjustWindowRect(...)` expands it to the full outer window size.

The GUI uses `Segoe UI` fonts, a file-browse button, edit boxes for the scalar parameters, method radio buttons, two checkboxes for preprocessing/domain options, a `Run Conversion` button, and a status box updated through a joinable worker thread.

---

## 20. Output products

The current code writes the following files:

1. `input.xyz.filtered.xyz`
2. `input.xyz.grid.xyz`
3. `input.xyz.confidence.xyz`
4. `input.xyz.dxf`
5. `input.xyz.rpt.txt`

### 20.1 `filtered.xyz`

Contains the optional outlier-filtered and deterministically thinned working cloud used as interpolation support.

### 20.2 `grid.xyz`

Contains the final interpolated regular-grid nodes written inside the active boundary, if boundary clipping is enabled.

### 20.3 `confidence.xyz`

Contains the interpolated nodes plus confidence, residual, lambda, boundary penalty, fallback flag, outside-hull flag, clamped-local flag, neighborhood count, and patch count.

### 20.4 `dxf`

Contains DXF point and label entities for the filtered XYZ cloud, optional grid entities, and the active boundary polyline.

### 20.5 `rpt.txt`

Contains the textual execution report with phase messages, validation summaries, tuning information, and confidence summary.

---

## 21. Practical method selection guidance

### 21.1 Bicubic Hermite + robust MLS

Recommended when:

- a smooth regular surface is desired;
- derivative continuity matters;
- the cloud is sufficiently sampled for stable local polynomial recovery;
- a general-purpose engineering grid is the primary target.

### 21.2 TPS

Recommended when:

- the point cloud is strongly scattered or geometrically irregular;
- local patch flexibility is more important than direct derivative control;
- a local radial-basis reconstruction is preferable to a cellwise Hermite surface.

### 21.3 Clamped local MLS (bounded)

Recommended when:

- overshoot suppression is more important than maximal smoothness;
- values should remain close to the nearby observed envelope;
- a conservative engineering surface is preferred.

---

## 22. Numerical safeguards and failure philosophy

The code is explicitly written to fail conservatively rather than silently produce unstable algebraic output. The principal safeguards are:

- pivot thresholding in all dense solves;
- finite-value checks after local solves;
- inverse-distance fallback when local polynomial or TPS evaluation fails;
- explicit fallback flags in the diagnostics stream;
- derivative limiting in the bicubic branch;
- boundary penalties and outside-hull diagnostics;
- allocation failure handling with safe abort and report output.

The implementation philosophy is therefore **robustness first, smoothness second, and silent numerical optimism never**.

---

## 23. Limitations and intended use

This program is a strong engineering interpolator and exporter, but it is not a full geostatistical or feature-aware terrain-modelling framework. In particular:

- it does not compute kriging variance;
- it does not infer semantic breaklines or discontinuities;
- its boundary object is a practical raster-derived support mask, not an exact alpha-shape or constrained triangulation engine;
- the bounded MLS mode is a conservative clamp, not a formal monotonicity-preserving PDE scheme.

Within those limits, the code is well suited to technical surface reconstruction from large engineering XYZ clouds.

---

## 24. Summary

The active `xyz2dxf_gui.cpp` implementation performs a numerically explicit workflow consisting of robust original-cloud screening, deterministic horizontal thinning, optional raster-derived outer-boundary construction, method-specific interpolation, diagnostic confidence estimation, and streamed CAD/text export. The README presented here is aligned with the current code path, current GUI controls, current tuning sets, and current formulas implemented in the MLS, bicubic Hermite, TPS, confidence, boundary, and DXF branches.
