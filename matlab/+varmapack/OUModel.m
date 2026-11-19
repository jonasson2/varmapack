classdef OUModel
  %VARMAPACK.OUMODEL  Stationary Gaussian Ornstein--Uhlenbeck model.
  %
  %   model = varmapack.OUModel(tau,cov,mu) has finite positive persistence
  %   times tau, stationary process covariance cov (not innovation covariance),
  %   and optional constant mean vector mu (default zero). Only the lower
  %   triangle of cov is used. cov and D(i,j)=cov(i,j)*(1/tau(i)+1/tau(j)) must
  %   be positive semidefinite; singular covariances are allowed.
  %
  %   model = varmapack.OUModel(K,cov,mu,'drift') uses the full drift term
  %   -K*(X-mu). Pass [] for a zero mean. K must have eigenvalues with
  %   strictly positive real parts, and K*cov+cov*K' must be positive
  %   semidefinite. The unused parameter property (tau or drift) is [].
  %
  %   Over a gap d, the autoregressive coefficient matrix is expm(-K*d).
  %   With tau, K is diag(1./tau). At a constant gap the sampled process
  %   is a VAR(1).
  %
  %   X = model.sim(t,M,rng) returns exact stationary paths without burn-in.
  %   t is a finite strictly increasing vector of observation times.
  %   M defaults to 1; rng is an optional varmapack.Rng. The result has shape
  %   r-by-numel(t)-by-M (MATLAB omits trailing singleton dimensions).
  %
  %   See also varmapack.gamma_times, varmapack.sim.
  properties (SetAccess = private)
    tau
    drift
    cov
    mu
    r
  end
  methods
    function obj = OUModel(parameter, cov, mu, mode)
      matrix_drift = nargin == 4;
      if matrix_drift && ~strcmp(mode, 'drift')
        error('varmapack:OUModel:mode', 'The fourth argument must be ''drift''');
      end
      if ~isnumeric(cov) || ~isreal(cov) || ~ismatrix(cov)
        error('varmapack:OUModel:cov', 'cov must be a real matrix');
      end
      r = size(cov, 1);
      if matrix_drift
        if ~isnumeric(parameter) || ~isreal(parameter) || ...
            ~isequal(size(parameter), [r r])
          error('varmapack:OUModel:drift', 'drift must be an r-by-r real matrix');
        end
      else
        if ~isnumeric(parameter) || ~isreal(parameter) || ...
            ~isvector(parameter) || isempty(parameter)
          error('varmapack:OUModel:tau', 'tau must be a nonempty real vector');
        end
        parameter = parameter(:);
      end
      if nargin < 3 || isempty(mu), mu = zeros(r, 1); end
      if ~isnumeric(mu) || ~isreal(mu) || ~isvector(mu)
        error('varmapack:OUModel:mu', 'mu must be a real vector');
      end
      parameter = full(double(parameter));
      mu = full(double(mu(:)));
      cov = full(double(cov));
      varmapack_ou_gateway(parameter, cov, mu, matrix_drift);
      if matrix_drift
        obj.drift = parameter;
      else
        obj.tau = parameter;
      end
      obj.cov = tril(cov) + tril(cov, -1)';
      obj.mu = mu;
      obj.r = r;
    end
    function rho = specrad(obj, gap)
      %SPECRAD  Sampled VAR(1) spectral radius at a positive finite gap.
      if ~isnumeric(gap) || ~isreal(gap) || ~isscalar(gap)
        error('varmapack:OUModel:gap', 'gap must be a real scalar');
      end
      matrix_drift = ~isempty(obj.drift);
      if matrix_drift, parameter = obj.drift; else, parameter = obj.tau; end
      rho = varmapack_ou_gateway(parameter, obj.cov, obj.mu, matrix_drift, ...
                                 full(double(gap)), 'specrad');
    end
    function G = acvf(obj, lags)
      %ACVF  Cov(X(t+h), X(t)) at real lags h, as r-by-r-by-numel(lags).
      if ~isnumeric(lags) || ~isreal(lags) || ~(isvector(lags) || isempty(lags))
        error('varmapack:OUModel:lags', 'lags must be a real vector');
      end
      matrix_drift = ~isempty(obj.drift);
      if matrix_drift, parameter = obj.drift; else, parameter = obj.tau; end
      G = varmapack_ou_gateway(parameter, obj.cov, obj.mu, matrix_drift, ...
                               full(double(lags(:))));
    end
    function X = sim(obj, t, M, rng)
      %SIM  Simulate M independent stationary paths at observation times t.
      if nargin < 3 || isempty(M), M = 1; end
      if ~isnumeric(t) || ~isreal(t) || ~isvector(t) || isempty(t)
        error('varmapack:OUModel:t', 't must be a nonempty real vector');
      end
      if nargin < 4 || isempty(rng)
        rng = varmapack.Rng();
        cleanup = onCleanup(@() delete(rng));
      end
      if ~isa(rng, 'varmapack.Rng') || ~isscalar(rng)
        error('varmapack:OUModel:rng', 'rng must be a varmapack.Rng object');
      end
      matrix_drift = ~isempty(obj.drift);
      if matrix_drift, parameter = obj.drift; else, parameter = obj.tau; end
      X = varmapack_ou_gateway(parameter, obj.cov, obj.mu, matrix_drift, ...
                               full(double(t(:))), M, rng.mex_handle());
    end
  end
end
