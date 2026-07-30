# Varmapack for MATLAB

Varmapack is a MATLAB interface to the Varmapack C library for simulation and
analysis of Gaussian VAR, VMA, VARMA, and VARMAX time-series models. Its primary
purpose is to generate simulated series from supplied model parameters, but it
also provides model testcases, theoretical and sample autocovariances, spectral
radii, and impulse response functions. In contrast to some other packages, the
simulated series have the correct distribution from the start; they are burn-in
(or spin-up) free. The numerical work is performed by the Varmapack C library,
with random number generation provided by the companion library Randompack.

## Installation

The MATLAB interface is included in the Varmapack GitHub repository. The
following instructions build it from source as a set of MEX gateways; it is not
a separately distributed MATLAB toolbox.

Use Git to clone both source repositories. Then follow the complete [Varmapack C
installation guide](../README.md#installation): first perform its **Install
Randompack** section, then its **Install Varmapack** section. Those two sections
build the libraries and install Randompack where Varmapack's build can find it.
There is no need to follow a separate Randompack guide. Keep both cloned
repositories after this step. The MATLAB MEX build uses their headers and built
libraries.

Then start MATLAB, configure a supported MEX C compiler if necessary, and run:

```matlab
VP_root = "path/to/varmapack";         % top level folder of varmapack clone
RP_root = "path/to/randompack";        % top level folder of randompack clone
mexpath = fullfile(VP_root, "matlab"); % joins paths for the current platform
addpath(mexpath)                       % add to MATLAB search path
mex -setup C                           % only needed once
build_varmapack_mex(RP_root)           % build MEX gateways into mexpath
savepath                               % make path persistent across sessions
```

The C build creates the Varmapack library in `varmapack_root/build/src` and the
Randompack library in `randompack_root/build/src`. The MEX gateways link against
both libraries. On Linux and macOS, `build_varmapack_mex` records these two
paths in each gateway, so MATLAB can find the libraries without further
configuration. On Windows, add both `build/src` directories to the system `PATH`
before starting MATLAB, so it can find `varmapack.dll` and `randompack.dll` when
a gateway loads.

Check the installation with:

```matlab
test_varmapack_rng
test_varmapack_testcasex
test_varmapack_simx
```

## More information

[``README.md``](https://github.com/jonasson2/varmapack/blob/main/README.md) in
the repository root contains detailed information on the underlying C library as
well as mathematical description of the supported models and computations. For C
library builds, MATLAB reference tests, and the current MEX build helper, see
[``DEVELOPMENT.md``](https://github.com/jonasson2/varmapack/blob/main/DEVELOPMENT.md),
also in the repository root.

TODO: add reference to paper.

## Timing

On macOS and Linux, `matlab/examples/TimeSimulate.sh` measures Varmapack
simulation for all named testcases and an unnamed `(p,q,r,rho)=(3,3,10,.98)`
model. It uses length 100, 1000 replicates, and a 0.1-second timing target per
testcase by default. Its output can be compared with the C and Python
`TimeSimulate` benchmarks:

```sh
    matlab/examples/TimeSimulate.sh
```

Use `-h` to list its timing and workload options.

On Windows, run the same benchmark directly from MATLAB:

```matlab
root = "C:\path\to\varmapack";
addpath(fullfile(root, "matlab", "examples"))
addpath(fullfile(root, "matlab"))
addpath(fullfile(root, "matlab-reference"))
addpath(fullfile(root, "tests", "matlab"))
TimeSimulate
```

Pass options as name-value pairs, for example `TimeSimulate("M", 5000)`.

## Usage

### Quick start

A runnable version of this sequence is available as `QuickStart.m` in the
`matlab/examples` directory.

```matlab
% Create a reproducible generator using Randompack's default engine.
vrng = varmapack.Rng(123);
cleanup = onCleanup(@() delete(vrng));
n = 200; % all simulated series have length 200
M = 10;  % number of replicates for most simulations

% Simulate one two-dimensional, zero-mean VAR(1) series
Sig1 = [2, 0; 0, 1];
A = [0.6, 0.1; 0, 0.4];
X = varmapack.sim(A, [], Sig1, [], n, 1, [], vrng); % one replicate
disp(X(:, 1:5))

% For this multivariate model, the first dimension is the series dimension,
% the second is time, and the third selects the replicate. Request a second
% output from varmapack.sim to return the shocks as E.

% Display all named testcases; then construct smallARMA1; a 2d VARMA(1,1) model
varmapack.testcase('summary')
[A, B, Sig] = varmapack.testcase('smallARMA1');
X = varmapack.sim(A, B, Sig, [], n, M, [], vrng);

% Supply a time-dependent mean path; its last column repeats to the series end.
mu = [0, 1; 0, 1];
X0 = [0, 0.1; 0, -0.1];  % startup path shared by all replicates.
X0multi = cat(3, X0, X0 + [0.05, 0.05; -0.05, -0.05]); % different per replicate
Xcommon = varmapack.sim(A, B, Sig, mu, n, M, X0, vrng);
Xseparate = varmapack.sim(A, B, Sig, mu, n, 2, X0multi, vrng); % two replicates

% Simulate VAR(1), VMA(2), and VARMA(2,1) models.
A1 = [0.6, 0.1; 0, 0.4];
B1 = [0.4, 0.3; 0.1, 0.2];
A2 = 0.2*ones(2);
B2 = 0.5*eye(2);
X1 = varmapack.sim(A1, [], Sig, [], n, M, [], vrng);
X2 = varmapack.sim([], [B1, B2], Sig, [], n, M, [], vrng);
[X3, E] = varmapack.sim([A1, A2], B1, Sig, [], n, M, [], vrng); % E = shocks

% Add two-dimensional exogenous input z with coefficient matrix C.
C = [0.8, -0.3; 0.2, 0.4];
z = [cos((1:n)/10); sin((1:n)/10)]; % d by n, broadcast to all replicates
zmulti = cat(3, z, z + 0.01);       % d by n by 2, different paths
[X4,E4] = varmapack.simx(A1, B1, C, z, Sig, n, M, zeros(2), 2, vrng); %E4=shocks
X5 = varmapack.simx(A1, B1, C, zmulti, Sig, n, 2, zeros(2), 2, vrng);

% Compute AR and MA spectral radii, and theoretical and data autocovariances
maxlag = 3;
rho = varmapack.specrad([A1, A2]);
rhoMA = varmapack.ma_specrad(B1);
Gamma = varmapack.acvf([A1, A2], B1, Sig, maxlag);
GammaHat = varmapack.autocov(X3(:, :, 1), maxlag, "ML"); % Normalize by 1/n
GammaHatC = varmapack.autocov(X3(:, :, 1), maxlag, "C"); % Normalize by 1/(n-k)

% Compute orthogonalized input responses and input response coefficients
Theta = varmapack.irf([A1, A2], B1, Sig, maxlag);
Psi = varmapack.psi([A1, A2], B1, maxlag);
```

### Functions

The MATLAB interface is a package of functions in the `varmapack` namespace.
It does not use a model object: pass model parameters directly to the relevant
function.

| Function | Purpose |
| :--- | :--- |
| [`varmapack.Rng`](+varmapack/Rng.m) | Create random number generator. |
| [`varmapack.sim`](+varmapack/sim.m) | Simulate one or more VARMA series. |
| [`varmapack.simx`](+varmapack/simx.m) | Simulate one or more VARMAX series. |
| [`varmapack.specrad`](+varmapack/specrad.m) | Compute the AR spectral radius. |
| [`varmapack.ma_specrad`](+varmapack/ma_specrad.m) | Compute the MA spectral radius. |
| [`varmapack.acvf`](+varmapack/acvf.m) | Compute theoretical autocovariances. |
| [`varmapack.autocov`](+varmapack/autocov.m) | Compute sample autocovariances. |
| [`varmapack.psi`](+varmapack/psi.m) | Impulse-response matrices. |
| [`varmapack.irf`](+varmapack/irf.m) | Orthogonalized impulse responses. |
| [`varmapack.testcase`](+varmapack/testcase.m) | Construct testcase models. |
| [`varmapack.testcasex`](+varmapack/testcasex.m) | Construct exogenous testcase data. |

### Array shapes

MATLAB stores matrices in column-major order. Coefficient matrices at successive
lags are concatenated horizontally:

- `A` has shape `r` by `r*p`: `[A1, A2, ..., Ap]`.
- `B` has shape `r` by `r*q`: `[B1, B2, ..., Bq]`.
- `C` has shape `r` by `d*s`: `[C1, C2, ..., Cs]`, where each `Ck` is
  `r` by `d` and multiplies the `d`-vector `z(t-k+1)`.
- `Sig` has shape `r` by `r`.
- The mean `mu` has shape `r` by `nmu`, where `nmu <= n`; the last supplied
  column repeats to the end. Use `[]` or `0` for zero mean.
- Startup values `X0` have shape `r` by `nX0` or `r` by `nX0` by `M`.
- Exogenous inputs `z` have shape `d` by `n` or `d` by `n` by `M`.
  The scalar `d=1` forms `n` by `1` and `n` by `M` are also accepted.
- For `r > 1`, simulated series and returned shocks have shape `r` by `n` by
  `M`. For scalar models, `sim` and `simx` return `1` by `n` when `M=1` and
  `n` by `M` otherwise.
- `varmapack.autocov` accepts one observed series `X` with shape `r` by `n`.
- `acvf`, `autocov`, `psi`, and `irf` return an `r` by `r` by `(maxlag+1)`
  array, with the lag-zero matrix in its first page.

### Errors

Invalid inputs and C-library failures are reported as MATLAB errors by the MEX
gateways.
