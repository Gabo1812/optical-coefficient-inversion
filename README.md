# Optical Coefficient Inversion in Thin Films

## Overview

This project focuses on the numerical inversion of optical properties of thin films from transmittance data. The goal is to recover the refractive index $n(\lambda)$, extinction coefficient $k(\lambda)$, and film thickness $d$ by solving a nonlinear inverse problem.

The implementation is based on the **PUMA (Program for Unconstrained Minimization Applied to optical characterization)** framework, applied to experimental and simulated UV-Vis transmittance spectra of thin films such as:

- ITO (Indium Tin Oxide)
- WO₃ (Tungsten Oxide)

deposited on transparent substrates such as soda-lime glass (SLG) or amorphous quartz.

---

## Physical Model

The forward model describes the optical transmittance of a thin absorbing film deposited on a transparent substrate, accounting for:

- Multiple internal reflections (interference effects)
- Absorption within the film
- Wavelength-dependent refractive indices

The transmittance is modeled as:

$$
T(\lambda) = \frac{A x}{B - Cx\cos\phi + Dx^2}
$$

where:

$$
x = e^{-\alpha d}, \quad \alpha = \frac{4\pi k}{\lambda}, \quad \phi = \frac{4\pi n d}{\lambda}
$$

and $A$, $B$, $C$, $D$ depend on $n$, $k$, and the refractive index of the substrate. Here $x$ is the attenuation factor and $\phi$ is the optical phase difference introduced by the film.

This model is **nonlinear, highly coupled, and oscillatory**, making inversion nontrivial.

### Key assumptions

- Planar, homogeneous thin films
- Known substrate refractive index $s(\lambda)$
- Normal incidence
- Coherent multiple reflections
- No scattering (pure absorption via $k(\lambda)$)

---

## Inverse Problem Formulation

Given measured transmittance data:

$$
T_{\text{meas}}(\lambda_i)
$$

we solve:

$$
\min_{d,\,n(\lambda),\,k(\lambda)} \sum_i \left[T_{\text{model}}(\lambda_i) - T_{\text{meas}}(\lambda_i)\right]^2
$$

### Ill-posedness

This problem is fundamentally underdetermined.

- For each $\lambda$, there is 1 equation and 2 unknowns $(n, k)$
- Without additional constraints, infinitely many solutions fit the data equally well

To regularize the problem, **physical constraints are imposed**:

- $n(\lambda) \geq 1$, $k(\lambda) \geq 0$
- Monotonic behavior away from the absorption edge: $n'(\lambda) \leq 0$, $k'(\lambda) \leq 0$
- Convexity conditions on $n(\lambda)$ and $k(\lambda)$
- A single inflection point in $k(\lambda)$, marking where the absorption edge changes curvature

These constraints encode **physically expected dispersion behavior** near the absorption edge.

---

## Numerical Method

### PUMA Framework

Instead of directly solving a constrained optimization problem, PUMA reformulates it as an **unconstrained problem** via a change of variables:

- Positivity enforced via squares (e.g. $n = 1 + u^2$)
- Convexity enforced through second derivatives

This transforms the problem into a high-dimensional nonlinear minimization:

$$
\min_x f(x)
$$

where $x$ includes:

- Film thickness $d$
- Inflection point $\lambda_{\text{infl}}$
- Discretized representations of $n(\lambda)$, $k(\lambda)$

---

### Optimization Algorithm

The minimization is performed using the **Spectral Projected Gradient method (SPG)**:

- First-order method (gradient-based)
- Adaptive step size (Barzilai–Borwein type)

Key properties:

- Fast per iteration
- Sensitive to scaling and initialization
- Can converge to local minima

---

## Practical Implementation Notes

From actual simulations:

- Run time depends on how fine the initial parameter exploration is, ranging from a few minutes to several hours (cluster execution)
- Strong sensitivity to:
  - Initial parameter ranges
  - Inflection point $\lambda_{\text{infl}}$
- Noise in experimental data significantly affects stability

### Observed issues

- Non-physical solutions (e.g. $k < 0$) when constraints are poorly enforced
- Slow convergence due to ill-conditioning
- Multiple local minima

### Mitigation strategies

- Savitzky–Golay filtering of experimental data
- Careful tuning of parameter ranges
- Using the PUMA solution as an **initial seed** for a complementary local optimization (an SPG-by-spectral-windows scheme implemented in Fortran), to cross-check robustness

---

## Results

Preliminary results show:

- Good agreement between simulated and experimental transmittance
- Physically consistent $n(\lambda)$, $k(\lambda)$ after parameter tuning
- Final squared errors as low as ~$10^{-6}$–$10^{-5}$ in favorable cases

However:

- Full convergence is not always achieved on the first attempt
- There is a trade-off between fit quality and physical plausibility

---

## Installation

```bash
conda env create -f environment.yml
conda activate optical-inversion
```

## Future Work

- Extend to multilayer systems (ITO–WO₃ stacks)
- Incorporate reflectance data (better conditioning)
- Replace PUMA with hybrid methods (global + local optimization)
- Parallelize parameter sweeps

## Author

José Gabriel Álvarez Castrillo
Physics Student, Universidad de Costa Rica

Guided by Dr. Edgar Rojas González (CICIMA, Universidad de Costa Rica)