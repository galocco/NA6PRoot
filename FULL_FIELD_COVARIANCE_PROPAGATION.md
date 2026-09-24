# Track and covariance propagation in a three-dimensional magnetic field

## Scope and conventions

The NA6P track state at a detector plane orthogonal to `Z` is

$$
\mathbf a = (x,\ y,\ t_x,\ t_y,\ q/p_{XZ})^T,
$$

with

$$
t_x = \frac{p_x}{p_{XZ}},\qquad
t_y = \frac{p_y}{p_{XZ}},\qquad
p_{XZ}=\sqrt{p_x^2+p_z^2},\qquad p_z>0.
$$

Consequently,

$$
t_x=\sin\psi,\qquad
c\equiv\cos\psi=\sqrt{1-t_x^2},
$$

and the geometrical slopes are

$$
\frac{dx}{dz}=\frac{t_x}{c},\qquad
\frac{dy}{dz}=\frac{t_y}{c}.
$$

The covariance matrix is stored as the packed lower triangle `mC[15]`, in
the parameter order above. Propagation applies

$$
C_1 = F C_0 F^T,
\qquad
F_{ij}=\frac{\partial a_{1,i}}{\partial a_{0,j}},
$$

before process noise from material is added by the material-correction code.

## Previous propagation: full-field state and pure-`By` Jacobian

The previous full-field function transported the reference state with all
three field components, but transported state deviations and covariance with
a Jacobian derived for

$$
\mathbf B=(0,B_y,0).
$$

For a pure `By` field, define

$$
K=k_{B2C}B_y,
\qquad
\kappa=K\frac{q}{p_{XZ}},
\qquad
\Delta z=z_1-z_0.
$$

The exact bending relation in the `XZ` plane is

$$
s_0=t_{x,0},\qquad
s_1=s_0+\kappa\Delta z,
$$

with

$$
c_0=\sqrt{1-s_0^2},\qquad
c_1=\sqrt{1-s_1^2}.
$$

The position change in the bending plane is

$$
x_1=x_0+\Delta z\frac{s_0+s_1}{c_0+c_1}.
$$

To linearize it, denote

$$
r=\frac q{p_{XZ}},\qquad
u=\frac{s_0+s_1}{c_0+c_1}.
$$

The elementary variations are

$$
\delta s_1=\delta s_0+K\Delta z\,\delta r,
$$

$$
\delta c_0=-\frac{s_0}{c_0}\delta s_0,
\qquad
\delta c_1=-\frac{s_1}{c_1}\delta s_1,
$$

and

$$
\delta u=
\frac{(\delta s_0+\delta s_1)(c_0+c_1)
-(s_0+s_1)(\delta c_0+\delta c_1)}{(c_0+c_1)^2}.
$$

Therefore,

$$
\delta x_1=\delta x_0+\Delta z\,\delta u.
$$

Collecting separately the coefficient multiplying `delta tx0` and the
coefficient multiplying `delta(q/pXZ)0` gives `f02` and `f04` below. The
direction variation follows immediately from

$$
\delta t_{x,1}=\delta t_{x,0}+K\Delta z\,\delta r,
$$

which gives `f24`.

For Y, the motion is the XZ arc length multiplied by `ty`:

$$
y_1=y_0+t_{y,0}L_{XZ},
$$

where, in a constant pure-`By` field,

$$
L_{XZ}=\frac{\psi_1-\psi_0}{\kappa},
\qquad
\psi_i=\arcsin s_i,
$$

with its continuous straight-line limit used as `kappa` approaches zero.
Thus

$$
\delta y_1=\delta y_0+L_{XZ}\delta t_{y,0}
            +t_{y,0}\delta L_{XZ}.
$$

Expanding `delta LXZ` through `delta s0` and `delta r` produces `f12`,
`f13`, and `f14` below.

In a pure `By` field, `ty` and `q/pXZ` are invariant. The old transport
therefore used a sparse Jacobian of the form

$$
F_{B_y}=
\begin{pmatrix}
1&0&f_{02}&0&f_{04}\\
0&1&f_{12}&f_{13}&f_{14}\\
0&0&1&0&f_{24}\\
0&0&0&1&0\\
0&0&0&0&1
\end{pmatrix}.
$$

The coefficients implemented in `NA6PTrackParCov` are obtained by defining

$$
d=c_0+c_1,\qquad
u=\frac{s_0+s_1}{d},
$$

$$
h=\frac{\Delta z}{d\,c_1}
  \left(1+c_0c_1+s_0s_1\right),
\qquad
j=\Delta z\left(u-\frac{s_1}{c_1}\right),
$$

and then

$$
f_{02}=\frac{h}{c_0},
\qquad
f_{24}=\Delta z K,
\qquad
f_{04}=h\frac{\Delta z}{d}K,
$$

$$
f_{12}=t_{y,0}\left(f_{02}s_1+j\right),
$$

$$
f_{13}=\Delta z\left(c_1+s_1u\right),
$$

$$
f_{14}=t_{y,0}\left(f_{04}s_1+jf_{24}\right).
$$

This Jacobian is analytic, fast, and correct for a locally constant pure
`By` field. It remains available through the `byOnly` propagation option.

## Selecting the covariance transport

The reconstruction configurable

```ini
useFullFieldJacobian=true
```

selects the new dense 3D-field Jacobian. Setting

```ini
useFullFieldJacobian=false
```

restores the legacy mixed behavior for direct comparisons: the reference
state is propagated in the full field, while deviations and covariance use
the sparse `By`-only Jacobian derived above.

The additional configurable

```ini
useFieldGradientJacobian=true
```

adds the response of the field map to transverse displacement. It is used only
when `useFullFieldJacobian=true` and `byOnly=false`. It defaults to `false`, so
the constant-field dense Jacobian remains the default full-field mode.

These switches are independent of `PropOpt::byOnly`. The available modes are:

| `byOnly` | `useFullFieldJacobian` | `useFieldGradientJacobian` | State transport | Covariance transport |
|---:|---:|---:|---|---|
| `true` | either | either | pure `By` | analytic pure `By` |
| `false` | `false` | either | full `(Bx,By,Bz)` | legacy analytic `By` approximation |
| `false` | `true` | `false` | full `(Bx,By,Bz)` | dense, locally constant full-field Jacobian |
| `false` | `true` | `true` | full `(Bx,By,Bz)` | dense full-field Jacobian with field gradients |

Both options are applied to VT, MuonSpec, and matching
`NA6PFastTrackFitter` instances when they are configured from
`NA6PRecoParam`.

### Why the old Jacobian is incomplete for a general field

For nonzero `Bx` or `Bz`, the assumptions represented by the last two rows of
the matrix are no longer valid:

- `Bx` bends the trajectory in a plane involving the Y momentum;
- `Bz` rotates the transverse momentum components and couples X and Y;
- `ty` is not generally invariant;
- `pXZ` is not generally invariant, even though the total momentum is
  invariant in a magnetic field;
- correlations such as `Cov(ty,q/pXZ)`, `Cov(tx,ty)`, and additional
  position-slope correlations are generated.

The mean trajectory could therefore follow the full field while its
uncertainty ellipse followed pure-`By` dynamics. Since the Kalman gain depends
on that uncertainty, the mismatch can change the fitted parameters and not
only their reported errors.

## New propagation: full-field state and full-field Jacobian

The reference state is still propagated through a locally constant field

$$
\mathbf B=(B_x,B_y,B_z)
$$

using the helix transport already used by `NA6PTrackPar`:

1. Convert `(tx,ty,q/pXZ)` to the normalized momentum direction and total
   momentum.
2. Rotate position and direction into a frame whose third axis is parallel
   to the local magnetic field.
3. Apply the constant-field GEANT3 helix step in that frame.
4. Rotate back to the laboratory frame.
5. Apply the final linear correction to land exactly on the requested
   constant-Z plane.
6. Convert the resulting direction back to `(tx,ty,q/pXZ)`.

The full transport can be written abstractly as

$$
\mathbf a_1=f(\mathbf a_0,\mathbf B,\Delta z).
$$

### Full-field equations of motion

Let

$$
\mathbf u=\frac{\mathbf p}{p}
          =(u_x,u_y,u_z)
$$

be the normalized momentum direction and let `s` be signed path length. In a
magnetic field, the momentum magnitude is constant and the direction obeys

$$
\frac{d\mathbf u}{ds}
=\frac{qk_{B2C}}{p}\,\mathbf u\times\mathbf B.
$$

Written component by component,

$$
\frac{du_x}{ds}=\frac{qk_{B2C}}p(u_yB_z-u_zB_y),
$$

$$
\frac{du_y}{ds}=\frac{qk_{B2C}}p(u_zB_x-u_xB_z),
$$

$$
\frac{du_z}{ds}=\frac{qk_{B2C}}p(u_xB_y-u_yB_x).
$$

These equations show the missing couplings directly. A pure `By` field
couples `ux` and `uz`; `Bx` additionally couples `uy` and `uz`; and `Bz`
couples `ux` and `uy`.

From the NA6P parameters at the beginning of a step,

$$
c=\sqrt{1-t_x^2},\qquad
n=\sqrt{1+t_y^2},
$$

the direction is

$$
u_x=\frac{t_x}{n},\qquad
u_y=\frac{t_y}{n},\qquad
u_z=\frac{c}{n}.
$$

The momentum is

$$
p=p_{XZ}n,
\qquad
p_{XZ}=\frac{q}{q/p_{XZ}}
$$

for a charged track.

### Rotation to the field frame

For a locally constant field, define

$$
B=|\mathbf B|,
\qquad
\widehat{\mathbf b}=\frac{\mathbf B}{B}.
$$

The implementation constructs a rotation `R` whose third axis is
`b-hat`. Position and direction are transformed as

$$
\mathbf r'=R\mathbf r,
\qquad
\mathbf u'=R\mathbf u.
$$

In this frame the field is `(0,0,B)`, so the motion is a simple helix. Define

$$
\rho=\frac{qk_{B2C}B}{p},
\qquad
\theta=\rho\ell,
$$

where `ell` is the signed path length used for the substep. The exact direction
update in the field frame is

$$
u'_{x,1}=u'_{x,0}\cos\theta-u'_{y,0}\sin\theta,
$$

$$
u'_{y,1}=u'_{x,0}\sin\theta+u'_{y,0}\cos\theta,
\qquad
u'_{z,1}=u'_{z,0}.
$$

The corresponding position update is

$$
x'_1=x'_0+
\frac{\sin\theta}{\rho}u'_{x,0}
-\frac{1-\cos\theta}{\rho}u'_{y,0},
$$

$$
y'_1=y'_0+
\frac{\sin\theta}{\rho}u'_{y,0}
+\frac{1-\cos\theta}{\rho}u'_{x,0},
$$

$$
z'_1=z'_0+\ell u'_{z,0}.
$$

The code evaluates the ratios involving `theta` with small-angle expansions
when necessary, avoiding divisions by a tiny curvature.

After the helix step,

$$
\mathbf r_1=R^T\mathbf r'_1,
\qquad
\mathbf u_1=R^T\mathbf u'_1.
$$

Because the helix path-length estimate does not generally land exactly on the
requested plane, a final local correction is made:

$$
\delta z=z_{target}-z_1,
$$

$$
x_1\leftarrow x_1+\delta z\frac{u_{x,1}}{u_{z,1}},
\qquad
y_1\leftarrow y_1+\delta z\frac{u_{y,1}}{u_{z,1}}.
$$

Finally, with

$$
D=\sqrt{u_{x,1}^2+u_{z,1}^2}=\frac{p_{XZ,1}}p,
$$

the direction is converted back to the NA6P state:

$$
t_{x,1}=\frac{u_{x,1}}D,
\qquad
t_{y,1}=\frac{u_{y,1}}D,
\qquad
\left(\frac q{p_{XZ}}\right)_1=\frac q{pD}.
$$

This last conversion is important: even though `p` is constant in a magnetic
field, `pXZ` and therefore `q/pXZ` can change when `Bx` or `Bz` rotates
momentum into or out of the Y direction.

### Direction derivatives generated by the full field

The conversion above gives, for infinitesimal direction changes,

$$
dt_x=\frac{u_z^2}{D^3}du_x
     -\frac{u_xu_z}{D^3}du_z,
$$

$$
dt_y=\frac{du_y}{D}
     -\frac{u_y}{D^3}(u_xdu_x+u_zdu_z),
$$

and

$$
d\left(\frac q{p_{XZ}}\right)
=-\frac{q}{pD^3}(u_xdu_x+u_zdu_z)
$$

for fixed total momentum. Combining these expressions with the component
Lorentz equations explains why a general field produces nonzero derivatives
such as

$$
\frac{\partial t_{y,1}}{\partial t_{x,0}},\quad
\frac{\partial t_{y,1}}{\partial(q/p_{XZ})_0},\quad
\frac{\partial(q/p_{XZ})_1}{\partial t_{y,0}},
$$

which are zero in the old sparse pure-`By` Jacobian.

### Numerical evaluation of the full Jacobian

Instead of using the sparse pure-`By` derivatives, the new implementation
evaluates the derivatives of this full transport in double precision:

$$
F_{ij}\simeq
\frac{f_i(\mathbf a_0+h_j\mathbf e_j)-
      f_i(\mathbf a_0-h_j\mathbf e_j)}{2h_j}.
$$

Central differences are evaluated for `tx`, `ty`, and `q/pXZ`. When field
gradients are disabled, the transport is translationally invariant in X and
Y, so the first two Jacobian columns are known exactly:

$$
\frac{\partial\mathbf a_1}{\partial x_0}
=(1,0,0,0,0)^T,
\qquad
\frac{\partial\mathbf a_1}{\partial y_0}
=(0,1,0,0,0)^T.
$$

The nominal finite-difference scale is

$$
h_j=\sqrt[3]{\epsilon_{double}}
    \max(1,|a_j|),
$$

with additional bounds that prevent a perturbed `tx` from reaching
`|tx|=1` or a perturbed `q/pXZ` from crossing zero.

### Explicit field-gradient contribution

Let the constant-field substep be written as

$$
a_1 = f\left(a_0,B(r_0)\right),
\qquad
r_0=(x_0,y_0,z_0).
$$

The total derivative with respect to the input state is

$$
F = \frac{d a_1}{d a_0}
  = \left.\frac{\partial f}{\partial a_0}\right|_B
    + \frac{\partial f}{\partial B}
      \frac{\partial B}{\partial r_0}
      \frac{\partial r_0}{\partial a_0}.
$$

The first term is the constant-field Jacobian. The second term is the explicit
field-gradient correction. Since the state is defined on a fixed-Z plane,

$$
\frac{\partial r_0}{\partial a_0}
=
\begin{pmatrix}
1&0&0&0&0\\
0&1&0&0&0\\
0&0&0&0&0
\end{pmatrix}.
$$

Therefore the local explicit correction contributes to the `x` and `y`
columns:

$$
F_{i0}=\left.\frac{\partial f_i}{\partial x_0}\right|_B
       +\sum_{k=x,y,z}
        \frac{\partial f_i}{\partial B_k}
        \frac{\partial B_k}{\partial x},
$$

$$
F_{i1}=\left.\frac{\partial f_i}{\partial y_0}\right|_B
       +\sum_{k=x,y,z}
        \frac{\partial f_i}{\partial B_k}
        \frac{\partial B_k}{\partial y}.
$$

The implementation evaluates this total derivative directly. For each
perturbed input state, it first constructs the linearly displaced field

$$
B_k^{\pm}=B_k
 +\frac{\partial B_k}{\partial x}\left(x_0^{\pm}-x_0\right)
 +\frac{\partial B_k}{\partial y}\left(y_0^{\pm}-y_0\right),
$$

then propagates the perturbed state with `B+` or `B-`. Central differences of
the complete result automatically contain both chain-rule terms.

The map derivatives are calculated with a 0.1 cm central stencil:

$$
\frac{\partial B_k}{\partial x}
\simeq\frac{B_k(x+0.1\,\mathrm{cm},y,z)
             -B_k(x-0.1\,\mathrm{cm},y,z)}{0.2\,\mathrm{cm}},
$$

and analogously for Y. A `dB/dz` map derivative is not needed in this local
Jacobian because all states are compared on the same Z plane, so
`delta z=0`. Longitudinal field variation of the nominal trajectory is still
included by resampling the map at every propagation substep.

Once the dense Jacobian has been constructed, covariance propagation uses a
general fixed-size matrix multiplication:

$$
C_1=F C_0 F^T.
$$

The result is copied back into the packed lower-triangle layout. All temporary
state, Jacobian, and covariance matrices use fixed-size stack storage; no heap
allocation is performed.

## How the Jacobian changes the fitted track direction

It is important to distinguish the nominal propagated trajectory from the
fitted trajectory after measurement updates.

### Propagation without an external linearization reference

If the track itself is used as the reference, its deviation from the reference
is zero. The nominal prediction is simply

$$
\mathbf a_1^-=f(\mathbf a_0^+,\mathbf B,\Delta z).
$$

In that case, changing only the Jacobian does **not** directly rotate this
single predicted central trajectory. It changes its predicted covariance:

$$
C_1^-=F C_0^+ F^T+Q,
$$

where `Q` is the material process-noise contribution.

At a detector plane, the measurement is

$$
\mathbf m=(x_m,y_m)^T,
$$

with measurement matrix

$$
H=
\begin{pmatrix}
1&0&0&0&0\\
0&1&0&0&0
\end{pmatrix}.
$$

The residual covariance and Kalman gain are

$$
S=HC_1^-H^T+V,
$$

$$
K=C_1^-H^TS^{-1},
$$

and the updated state is

$$
\mathbf a_1^+=\mathbf a_1^-+
K(\mathbf m-H\mathbf a_1^-).
$$

The direction updates are therefore

$$
\Delta t_x=K_{t_xx}\,r_x+K_{t_xy}\,r_y,
$$

$$
\Delta t_y=K_{t_yx}\,r_x+K_{t_yy}\,r_y,
$$

where

$$
\mathbf r=(r_x,r_y)^T=\mathbf m-H\mathbf a_1^-.
$$

The relevant rows of `K` are built from propagated covariance elements such
as

$$
C_{t_xx},\ C_{t_xy},\ C_{t_yx},\ C_{t_yy}.
$$

Those elements are generated by `F C F^T`. A wrong Jacobian therefore gives
the wrong position-direction correlations, the wrong Kalman gain, and the
wrong correction to `tx` and `ty`. The next propagation starts from these
updated direction parameters, so the effect accumulates across detector
planes.

The same mechanism changes the momentum parameter through

$$
\Delta(q/p_{XZ})=
K_{(q/p_{XZ})x}r_x+K_{(q/p_{XZ})y}r_y.
$$

This is how a covariance Jacobian can affect the reconstructed momentum and
mass resolution even when the nominal propagation already uses the full
field.

### Propagation with an external linearization reference

When `useLinRef` is enabled, the fitted track is represented by its deviation
from a separately propagated reference trajectory. At the initial plane,

$$
\delta a_0 = a_0 - a_{\mathrm{ref},0}.
$$

The reference trajectory is propagated with the nonlinear transport function:

$$
a_{\mathrm{ref},1} = f(a_{\mathrm{ref},0}).
$$

The fitted state at the destination plane is then predicted by linearizing
around that reference:

$$
a_1^- = a_{\mathrm{ref},1} + F\,\delta a_0
      = f(a_{\mathrm{ref},0})
        + F\left(a_0-a_{\mathrm{ref},0}\right).
$$

Therefore, the Jacobian does not act only on the covariance. It also acts
directly on the difference between the fitted track and the reference track.
In particular, the propagated direction corrections are

$$
\delta t_{x,1} = \sum_{j=0}^{4} F_{2j}\,\delta a_{0,j},
\qquad
\delta t_{y,1} = \sum_{j=0}^{4} F_{3j}\,\delta a_{0,j}.
$$

The `tx` and `ty` rows of `F` therefore determine how a mismatch in any input
parameter changes the predicted track direction. The old sparse Jacobian
omitted direction couplings generated by `Bx` and `Bz`; the new dense Jacobian
includes them.

### Indirect effects through selection

The predicted cluster chi-square is

$$
\chi^2=\mathbf r^T S^{-1}\mathbf r.
$$

Changing `F` changes `S`, which can change cluster rejection, accepted track
candidates, and which candidate survives shared-cluster selection. Thus the
Jacobian can affect the final track sample as well as the parameters of an
individual accepted track.

The transport also now handles the case

$$
B_y\simeq0,\qquad B_x\ne0\ \text{or}\ B_z\ne0.
$$

Previously, the small-`By` test incorrectly selected straight-line transport
even when another field component was nonzero.

## Expected effect on MuonSpec tracking

The largest effects are expected where `Bx/By` or `Bz/By` is significant:

- magnet fringe regions;
- tracks far from the field-map symmetry axis;
- large-acceptance configurations;
- long extrapolations with few measurement constraints;
- layouts with only one station on one side of the magnet.

Expected improvements include:

- pull widths closer to one;
- more realistic `q/pXZ` uncertainty;
- correctly generated X-Y and momentum-slope correlations;
- better-calibrated predicted cluster chi-square;
- more appropriate Kalman gains and cluster weights;
- potentially improved momentum and dimuon-mass resolution when the previous
  covariance mismatch was significant.

The central trajectory was already using all three field components, so large
changes to the mean residual are not expected in a central, nearly pure-`By`
region. The primary expected change is in covariance consistency and in the
Kalman updates driven by that covariance.

This correction cannot recover information lost through detector geometry.
In particular, it does not compensate for removing the second upstream
MuonSpec station. The asymmetric station layout remains a separate and likely
larger source of momentum-resolution degradation.

## Limitations and computational cost

The state propagation still treats the field as constant during each
propagation substep. The outer propagator samples the field again after every
substep, normally at most 2 cm, so variation along the reference trajectory
is included. With `useFieldGradientJacobian=true`, the Jacobian also describes
how a transverse displacement changes the field sampled at the start of each
substep. It does not integrate continuous field variation inside one
substep. Reducing `maxPropagationStep` reduces that remaining approximation.

The full numerical Jacobian requires six double-precision perturbed
propagations per substep: positive and negative perturbations for each of the
three nontrivial input parameters. Enabling field gradients requires ten
perturbed propagations, covering all five state parameters, plus four field-map
queries for the X and Y central differences. It is therefore the most
expensive mode. The `byOnly` option retains the original fast path for studies
where the dipole approximation is sufficient.

The most useful physics validation is to compare the old and new modes using
the same generated events, clusters, geometry, and track candidates, checking:

1. residual widths for all five native parameters;
2. pull means and widths;
3. pull underflow and overflow counts;
4. cluster and track chi-square distributions;
5. reconstruction efficiency;
6. single-muon momentum resolution;
7. dimuon invariant-mass resolution.