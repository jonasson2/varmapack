<!-- -*- poly-markdown -*- -->
# Varmapack
<https://github.com/jonasson2/varmapack>

## Overview

Varmapack is a C library for simulation and analysis of Gaussian VAR, VMA,
VARMA, and VARMAX time-series models. Its primary purpose is to generate
simulated series from supplied model parameters, but it also provides model
testcases, theoretical and sample autocovariances, covariance-to-correlation
conversion, spectral radii, and impulse response functions. In contrast to some
other packages, the simulated series have the correct distribution from the
start; they are burn-in (or spin-up) free. The package uses Randompack
(https://github.com/jonasson2/randompack), [xx], to generate random numbers for
the simulation. C, Python, and MATLAB interfaces are described in the
respective README files in the GitHub repository.

Varmapack is based on a part of Algorithm 878 published in ACM TOMS in 2008
[xx]. That algorithm is a collection of MATLAB functions to evaluate the
likelihood of a VARMA model, which also includes functions to simulate such
models and provide testcase data.

The public header file `varmapack.h` serves as a compact reference for the C
API: all user-facing functions are declared there, with comments describing the
role of each parameter.

A fixed seed gives reproducible simulation paths for a fixed compiler, BLAS
library, and Randompack build. Bitwise agreement across platforms is not
guaranteed when a conditional startup covariance is rank deficient: floating
point factorization can select a different numerical rank and therefore consume
a different number of normal variates. The resulting paths remain draws from
the same conditional distribution.

A manuscript about Varmapack is being prepared and expected to be submitted to
the journal SoftwareX. [TODO: fill in details]

## Mathematical Description

For more details of the mathematics, see the report [1]. (TODO: put in the
right reference).

### VARMA and VARMAX models

The models considered are either VARMA $(p,q)$:

$$
\tag{1}
x_t = \sum_{i=1}^{p} A_i x_{t-i} + \varepsilon_t
      + \sum_{j=1}^{q} B_j \varepsilon_{t-j},
      \quad \varepsilon_t \sim N(0,\Sigma).
$$

or VARMAX $(p,q,s)$:

$$
\tag{2}
x_t = \varepsilon_t + \sum_{i=1}^{p} A_i x_{t-i}
      + \sum_{j=1}^{q} B_j \varepsilon_{t-j}
      + \sum_{k=1}^{s} C_k z_{t-k+1},
      \quad \varepsilon_t \sim N(0,\Sigma),
$$

where in both cases $x_t$ is $r$-dimensional, $\varepsilon_t$ are shocks or
innovations,
$A_i$ and $B_j$ are $r$ by $r$ autoregressive and moving-average coefficient
matrices, $C_k$ is an $r$ by $d$ exogenous coefficient matrix, $z_t$ is a
$d$-dimensional exogenous vector, and $\Sigma$ is the innovation covariance
matrix. Varmapack supports a time-dependent mean path $\mu_t$ for VARMA
simulation. Explicit mean paths are unsupported for VARMAX. For VARMA, the
recursion is applied to the centered series
$x_t - \mu_t$, which replaces $x_t$ in equation (1). The mean may also be
fixed, i.e. independent of the time step.

### VAR, VMA, VARX, and VMAX models

A model with $q=s=0$ is a pure vector-autoregressive (VAR) model, and when
$p=s=0$ the model is a pure moving-average (VMA) model. Similarly, when $q$ or
$p$ are zero and $s > 0$ the models are designated VARX or VMAX.

### Simulation start

There are two possibilities to start VARMA simulation with Varmapack: (a) by
drawing both shocks $\varepsilon_t$ and states $x_t$ from the exact joint
distribution of $(x,\varepsilon)$ for the initial segment $t=0,\ldots,h-1$,
where $h=\max(p,q)$, and (b) by specifying $h \geq\max(p,q)$ initial values of
the series and drawing the first $h$ shocks from the conditional distribution of
$(\varepsilon|x)$. In case (b), the covariance matrix of the startup states must
be positive definite.

For a stationary model, both (a) and (b) are possible, and the simulated series
will have the correct distribution from the first generated term, so that no
burn-in segment needs to be discarded.

For a nonstationary model, starting values must be supplied. For a pure VAR
model ($q=0$), the recurrence simply runs forward from those values. If MA terms
are present, startup innovations are drawn from their theoretical distribution,
conditional on constraints imposed by the model and supplied states. $\Sigma$
must be positive semidefinite and is allowed to be singular except for
nonstationary models with MA terms.

With VARMAX simulation, startup states must be provided: $x_t$ for
$t=0,\ldots,h-1$, where $h\geq\max(p,q,s-1)$, as well as the whole sequence of
exogenous values $z_t$, $t=0,\ldots,n-1$. As for nonstationary VARMA, $\Sigma$
must be positive definite, and the startup innovations are then drawn from their
theoretical distribution, conditional on constraints imposed by the model and
the supplied states and exogenous sequence. When more than $\max(p,q,s-1)$
values of $x_t$ and corresponding $z_t$ are available, they should be supplied,
to improve the information used to draw the innovations.

### Autocovariances

Varmapack can compute the theoretical autocovariance function of model (1),

$$
\Gamma_k = \operatorname{Cov}(x_t, x_{t-k}), \qquad k=0,1,2,\ldots,
$$

up to a user-specified maximum lag. It can also compute sample autocovariance
matrices from observed data. The default maximum-likelihood normalization is

$$
\widehat{\Gamma}_k =
\frac{1}{n}\sum_{t=k}^{n-1}(x_t-\bar{x})(x_{t-k}-\bar{x})^T,
$$

and the corrected normalization replaces the denominator $n$ by $n-k$.

### Spectral radii

The spectral radius $\rho$ of the model (1) is defined as the spectral radius
(maximum absolute eigenvalue) of the autoregressive companion matrix, which is a
square block matrix with $A_1,\ldots,A_p$ on the first block row and identity
matrices on the block-subdiagonal. Stationary models have $\rho < 1$. The
spectral radius of the moving average companion matrix, which has
$-B_1,\ldots,-B_q$ on the top block row, plays a different role: it provides
an invertibility diagnostic. Models with MA companion spectral radius $<1$ are
invertible, in the sense that the moving-average polynomial of the backshift
operator $L$, $B(L)$ has a convergent inverse,
$B(L)^{-1} = I + D_1L + D_2 L^2 + \ldots$. It follows that the model can be
written as an equivalent infinite pure VAR model

$$
x_t = E_1 x_{t-1} + E_2 x_{t-2} + \ldots + \varepsilon_t,
$$

where $I-E_1L-E_2L^2-\ldots = (I + D_1L + D_2 L^2 + \ldots) A(L)$.

### Impulse response functions

An important concept in time-series theory is that of impulse response
functions. The impulse response matrix $\Psi_j$ maps a change in the shock at
time $t$ to the resulting change in the process at time $t + j$. If the time
series is stationary, it can be expressed as an infinite VMA series:

$$
x_t = \sum_{j=0}^{\infty} \Psi_j \varepsilon_{t-j},
$$

which has the impulse response matrices as coefficients. With $\Psi_0=I$, the
remaining $\Psi$ matrices can be computed with the recursion

$$
\Psi_j = B_j + \sum_{i=1}^{\min(p,j)} A_i\Psi_{j-i},
$$

for $j\geq 1$, with $B_j=0$ for $j>q$. Varmapack provides functions to compute
these impulse responses as well as *orthogonalized impulse responses*, defined
as

$$
\Theta_j=\Psi_jL,
$$

where $L L^T=\Sigma$; for positive-definite $\Sigma$, $L$ is its lower
Cholesky factor. Only the lower triangle of $\Sigma$ is used.

## Installation

The supported installation method is to build Varmapack from source. Building
the library requires Meson, Ninja, `pkg-config`, C and Fortran compilers, and
an optimized BLAS library. macOS provides BLAS through the Accelerate framework.
On Linux and Windows, OpenBLAS or MKL must be installed separately. Varmapack
also requires the Randompack C library; see below.

Install or load these prerequisites using the normal method for your system.
For example, an HPC cluster may supply them through environment-module systems,
using `module load` commands.

### Install Randompack

Install Randompack 0.1.10 or later system-wide or under `~/.local` as described
in its [README](https://github.com/jonasson2/randompack#installation). Before
configuring Varmapack, verify that `pkg-config` can find a suitable version:

```sh
    pkg-config --atleast-version=0.1.10 randompack
```

For source installations, specifying `--libdir=lib` for both Randompack and
Varmapack is recommended to ensure consistent library paths across platforms.
The Randompack README and the Varmapack commands below use this setting. With
`--prefix=/usr/local`, the libraries are installed under `/usr/local/lib`;
without an explicit `libdir`, Meson may choose a platform-specific multiarch
directory.

For a user-local installation, this requires the setting described in the
Randompack README:

```sh
    export PKG_CONFIG_PATH="$HOME/.local/lib/pkgconfig${PKG_CONFIG_PATH:+:$PKG_CONFIG_PATH}"
```

This export sets `PKG_CONFIG_PATH` if it is undefined, or prepends the new
directory to it if it exists. Add the export to the appropriate shell startup
file if the user-local installation should be available in future shells. No
such setting is normally needed for a system-wide installation.

### Install Varmapack

Clone Varmapack and choose either a system-wide or user-local installation
before configuring the build:

```sh
    git clone https://github.com/jonasson2/varmapack.git
    cd varmapack
```

For a system-wide installation, configure, build, test, and install with:

```sh
    meson setup build --buildtype=release --prefix=/usr/local --libdir=lib
    ninja -C build
    meson test -C build
    sudo meson install -C build
```

For a user-local installation that does not require `sudo`, use `~/.local`:

```sh
    meson setup build --buildtype=release --prefix="$HOME/.local" --libdir=lib
    ninja -C build
    meson test -C build
    meson install -C build
```

The `PKG_CONFIG_PATH` setting above also makes the user-local Varmapack
installation visible. In either case, verify the installed version and inspect
the compiler and linker flags with:

```sh
    pkg-config --modversion varmapack
    pkg-config --cflags --libs varmapack
```

Programs using Varmapack can then be compiled with:

```sh
    cc -o myprog myprog.c $(pkg-config --cflags --libs varmapack)
```

Varmapack's `pkg-config` file declares Randompack as a dependency, so this
command supplies the flags for both libraries; Randompack need not be added
separately.

## C API overview by example

Varmapack stores matrices in column-major order. Coefficient matrices belonging
to successive lags are contiguous blocks, and simulated series are stored as an
`r` by `n` by `M` array.

The code blocks through the VARMAX simulation section are successive parts of one
`main` function in `QuickStart.c` in the `examples` folder. For
readability, error checking is omitted; the final section of the chapter shows
how to check errors.

### VAR simulation

The program begins by simulating 10 replicates of length 200 from a
two-dimensional VAR(1) model. Passing `B=0` and `q=0` selects a pure VAR model;
similarly, `A=0` and `p=0` would select a VMA model. Once Varmapack is
installed, compile the complete `QuickStart.c` example using the `cc`
command given above.

```c
    #include <stdio.h>
    #include <stdlib.h>
    #include <varmapack.h>

    int main(void) {
      int p = 1, q = 0, r = 2, n = 200, M = 10;
      double Avar[] = {0.6, 0, 0.1, 0.4}; // r by r by p
      double Sigvar[] = {2, 0, 0, 1};     // r by r
      double X[2*200*10];                 // r by n by M
      randompack_rng *rng = randompack_create(0); // select default Randompack engine
      randompack_seed(123, 0, 0, rng); // optional, but gives reproducible results
      varmapack_sim(Avar, 0, Sigvar, 0, 0, p, q, r, n, M, 0, 0, 1, X, 0, rng);
      for (int t=0; t<5; t++) printf("%8.4f %8.4f\n", X[r*t], X[1 + r*t]);
```

### Testcase construction

The next part selects the fixed `smallARMA1` testcase. First inquire about its
dimensions, then allocate and fill its coefficient arrays. This testcase has
`p=q=1` and `r=2`.

```c
      char name[VARMAPACK_TESTCASE_NAME_LEN] = "smallARMA1";
      int index = 0;

      varmapack_testcase(name, &index, 0, &p, &q, &r, 0, 0, 0, 0);
      double *A = malloc((p ? p : 1)*r*r*sizeof(*A));
      double *B = malloc((q ? q : 1)*r*r*sizeof(*B));
      double *Sig = malloc(r*r*sizeof(*Sig));
      varmapack_testcase(name, &index, 0, &p, &q, &r, A, B, Sig, 0);
```

### VARMA simulation

The testcase is simulated with a time-dependent mean and a supplied startup
path. The last `mu` column repeats to the end, so the mean at t=0 is (0,0) and
at every subsequent time step is (1,1). Here `nX0=2`, while `MX0=1` broadcasts
the single startup path. Set `MX0=M` and supply an `r` by `nX0` by `M` array to
give each replicate its own path.

```c
      double mu[] = {0, 0, 1, 1};        // r by nmu (time dependent means)
      double X0[] = {0, 0, 0.1, -0.1};   // r by nX0 (common start for all replicates)
      double E[2*200*10];                // r by n by M (returns shocks)

      varmapack_sim(A, B, Sig, mu, 2, p, q, r, n, M, X0, 2, 1, X, E, rng);
```

### Model analysis

Next, `QuickStart.c` simulates a zero-mean realization of `smallARMA1` and
computes theoretical and sample autocovariances, theoretical correlations,
ordinary and orthogonalized impulse responses, and AR and MA spectral radii.

```c
      varmapack_sim(A, B, Sig, 0, 0, p, q, r, n, M, 0, 0, 1, X, 0, rng);
      enum { maxlag = 3 };
      double Gamma[2*2*(maxlag+1)], Corr[2*2*(maxlag+1)];
      double Psi[2*2*(maxlag+1)], Theta[2*2*(maxlag+1)], GammaHat[2*2*(maxlag+1)];
      double rho = varmapack_specrad(A, r, p);
      double rhoMA = varmapack_ma_specrad(B, r, q);
      varmapack_acvf(A, B, Sig, p, q, r, Gamma, maxlag);
      varmapack_cov2corr(Gamma, r, maxlag, Corr);
      varmapack_autocov("N", "ML", r, n, X, maxlag, GammaHat);
      varmapack_psi(A, B, p, q, r, maxlag, Psi);
      varmapack_irf(A, B, Sig, p, q, r, maxlag, Theta);
      printf("AR spectral radius: %.4f\n", rho);
      printf("MA spectral radius: %.4f\n", rhoMA);
      for (int k=0; k<=maxlag; k++) printf("Gamma[%d]: %.4f %.4f %.4f %.4f\n",
        k, Gamma[4*k], Gamma[4*k+1], Gamma[4*k+2], Gamma[4*k+3]);
      for (int k=0; k<=maxlag; k++) printf("Corr[%d]: %.4f %.4f %.4f %.4f\n",
        k, Corr[4*k], Corr[4*k+1], Corr[4*k+2], Corr[4*k+3]);
      for (int k=0; k<=maxlag; k++) printf("GammaHat[%d]: %.4f %.4f %.4f %.4f\n",
        k, GammaHat[4*k], GammaHat[4*k+1], GammaHat[4*k+2], GammaHat[4*k+3]);
      for (int k=0; k<=maxlag; k++) printf("Psi[%d]: %.4f %.4f %.4f %.4f\n",
        k, Psi[4*k], Psi[4*k+1], Psi[4*k+2], Psi[4*k+3]);
      for (int k=0; k<=maxlag; k++) printf("Theta[%d]: %.4f %.4f %.4f %.4f\n",
        k, Theta[4*k], Theta[4*k+1], Theta[4*k+2], Theta[4*k+3]);
```

`Gamma`, `Corr`, `Psi`, `Theta`, and `GammaHat` each contain `maxlag + 1`
contiguous `r` by `r` matrices. The sample call above analyzes the first
replicate in `X`; use normalization `"C"` instead of `"ML"` for division by
`n - k`.

### Covariance-to-correlation conversion

`varmapack_cov2corr(cov, r, maxlag, corr)` converts a sequence of covariance
matrices such as `Gamma` or `GammaHat` to correlations. Both arrays have shape
`r` by `r` by `maxlag + 1` in column-major order. Each entry at lag `k` is
divided by the product of the corresponding standard deviations obtained from
the diagonal of the lag-zero covariance matrix.

The call above preserves the covariances in `Gamma` and writes the correlations
to `Corr`. Pass the same array for both arguments to convert it in place.
Distinct arrays must not overlap. The diagonal lag-zero correlations are set
exactly to one.

The function deliberately does not clip other values to `[-1,1]`. In
particular, the corrected sample autocovariances produced by
`varmapack_autocov` with normalization `"C"` use lag-dependent divisors and can
therefore yield values outside that interval after conversion.

Like the other Boolean-valued analysis functions, `varmapack_cov2corr` returns
`false` and sets `varmapack_last_error()` when an argument is invalid, an input
entry is not finite, or a lag-zero variance is not positive.

### VARMAX simulation

Finally, the example creates deterministic two-dimensional exogenous testcase
data. Its `d` by `n` input sequence is broadcast to every replicate. To use a
different sequence for each replicate, set `Mz=M` and store `z` as `d` by `n`
by `M`.

```c
      int s = 2, d = 2;
      double C[2*2*2];
      double z[2*200];
      varmapack_testcasex(s, d, r, n, C, z);
      varmapack_simx(A, B, C, Sig, z, 1, p, q, s, d, r, n, M, X0, 2, 1, X, E, rng);
      randompack_free(rng);
      free(Sig);
      free(B);
      free(A);
      return 0;
    }
```

### Error checking

Most Varmapack functions return `false` on failure and make a diagnostic
available through `varmapack_last_error()`. The exceptions are
`varmapack_specrad()` and `varmapack_ma_specrad()`, which return `NAN` on
failure. The diagnostic is thread-local, so a failure in one thread does not
replace the error reported in another. Randompack reports errors through its
RNG object. The standalone `ErrorCheck.c` in the `examples` folder checks
both libraries:

```c
    #include <stdio.h>
    #include <varmapack.h>

    int main(void) {
      int status = 0;
      int r = 2, n = 200, M = 10;
      double A[] = {0.6, 0, 0.1, 0.4};
      double Sig[] = {2, 0, 0, 1};
      double X[2*200*10];
      randompack_rng *rng = randompack_create(0);
      if (!rng) {
        fprintf(stderr, "could not create RNG\n");
        return 1;
      }
      if (!randompack_seed(123, 0, 0, rng)) {
        fprintf(stderr, "Randompack: %s\n", randompack_last_error(rng));
        status = 1;
      }
      else if (!varmapack_sim(A, 0, Sig, 0, 0, 1, 0, r, n, M, 0, 0, 1, X, 0,
                               rng)) {
        fprintf(stderr, "Varmapack: %s\n", varmapack_last_error());
        status = 1;
      }
      randompack_free(rng);
      return status;
    }
```

## Testing

Varmapack includes three Meson test programs:

- `RunTests` provides direct C unit and edge-case tests of the public functions.
- `AgainstMatlab` compares C results with MATLAB reference outputs.
- `TestLyapunov` tests the state-space covariance solver used for simulation
  startup.

Run the C test suite with:

```sh
    meson test -C build --print-errorlogs
```

`AgainstMatlab` compares the public numerical functions and simulated series,
using identical Randompack streams where applicable. The standard suite uses
the checked-in `tests/matlabcompare.txt` fixture, so MATLAB is not needed to
run it.

Most of the MATLAB reference functions are based on the thoroughly tested
MATLAB code from the 2008 ACM TOMS Algorithm 878 package. New reference
functions were developed mostly independently of the C implementations. MATLAB
functions are generally easier to write from mathematical algorithm
descriptions than C functions, making implementation errors less likely. This
makes the reference implementation a valuable independent cross-check.

When MATLAB is available, see `DEVELOPMENT.md` for regenerating the reference
fixture and running the full comparison.

## Timing

The `benchmark` directory contains four C timing programs. Build them with:

```sh
    ninja -C build benchmark/TimeSimulate benchmark/TimeScalability \
        benchmark/TimeSetup benchmark/TimeBreakEven
```

`TimeSimulate` measures direct C-library simulation for all named testcases and
an unnamed `(p,q,r,rho)=(3,3,10,.98)` model. It uses length 100, 1000
replicates, and a 0.1-second timing target per testcase by default. Its output
is directly comparable with `matlab/examples/TimeSimulate.sh` and
`python/examples/TimeSimulate.py`:

```sh
    build/benchmark/TimeSimulate
```

Use `-h` to list its timing and workload options.

`TimeScalability` varies one of the model order, dimension, series length,
replicate count, or spectral radius at a time around a configurable reference
model. It separates the amortized setup estimate from the remaining generation
time:

```sh
    build/benchmark/TimeScalability
```

`TimeSetup` and `TimeBreakEven` are development benchmarks for analyzing the
two startup solvers. See `DEVELOPMENT.md` for their use.
