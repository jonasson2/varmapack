# varmasim

**VARMA time-series simulation library** — intended for CRAN submission but also usable as a standalone C program.

---

## Overview

This project provides exact, spin-up-free simulation of **VARMA(p, q)** and related models.  
It can be used in two ways:

- **As an R package** (for CRAN) providing efficient simulation and likelihood utilities.
- **As a standalone C library** or demo executable (`RunVarmaSim`) built via [Meson](https://mesonbuild.com).

---

## Model

A **VARMA(p, q)** process is defined as

$$
x_t - \mu = \sum_{i=1}^{p} A_i (x_{t-i} - \mu)
           + \sum_{j=0}^{q} B_j \varepsilon_{t-j},
\qquad
\varepsilon_t \sim \mathcal{N}(0, \Sigma).
$$

- \(x_t\) is an \(r\)-dimensional time series.  
- The innovations \(\varepsilon_t\) are uncorrelated in time and have covariance \(\Sigma\).  
- When \(r = 1\) the model reduces to an **ARMA(p, q)** series.  
- When \(q = 0\) and \(r > 1\) it becomes a **VAR(p)** model.  
- When \(q = 0\) and \(r = 1\) it is a simple **AR(p)** model.

The simulation avoids any “burn-in” period: the generated series has the correct covariance structure from the first observation.

---

## Usage (standalone C build)

```bash
meson setup build
cd build
meson compile
./RunVarmaSim smallMA     # named testcase
./RunVarmaSim 5           # numbered testcase
./RunVarmaSim 3,2,2       # dimensions p,q,r