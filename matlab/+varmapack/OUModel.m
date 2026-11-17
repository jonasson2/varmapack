classdef OUModel
  %VARMAPACK.OUMODEL  Stationary Gaussian Ornstein--Uhlenbeck model.
  %
  %   model = varmapack.OUModel(tau,cov,mu) has finite positive persistence
  %   times tau, stationary process covariance cov (not innovation covariance),
  %   and optional constant mean vector mu (default zero). Only the lower
  %   triangle of cov is used. cov and D(i,j)=cov(i,j)*(1/tau(i)+1/tau(j)) must
  %   be positive semidefinite; singular covariances are allowed.
  %
  %   X = model.sim(t,M,rng) returns exact stationary paths without burn-in.
  %   t is a finite strictly increasing vector, in the same units as tau.
  %   M defaults to 1; rng is an optional varmapack.Rng. The result has shape
  %   r-by-numel(t)-by-M (MATLAB omits trailing singleton dimensions).
  %   Each component has autocorrelation exp(-gap/tau(i)). The drift is
  %   diagonal: dependence between components comes from correlated noise.
  %
  %   See also varmapack.gamma_times, varmapack.sim.
  properties (SetAccess = private)
    tau
    cov
    mu
    r
  end
  methods
    function obj = OUModel(tau, cov, mu)
      if nargin < 3 || isempty(mu), mu = zeros(numel(tau), 1); end
      if ~isnumeric(tau) || ~isreal(tau) || ~isvector(tau) || isempty(tau)
        error('varmapack:OUModel:tau', 'tau must be a nonempty real vector');
      end
      if ~isnumeric(mu) || ~isreal(mu) || ~isvector(mu)
        error('varmapack:OUModel:mu', 'mu must be a real vector');
      end
      if ~isnumeric(cov) || ~isreal(cov) || ~ismatrix(cov)
        error('varmapack:OUModel:cov', 'cov must be a real matrix');
      end
      tau = full(double(tau(:)));
      mu = full(double(mu(:)));
      cov = full(double(cov));
      varmapack_ou_gateway(tau, cov, mu);
      obj.tau = tau;
      obj.cov = tril(cov) + tril(cov, -1)';
      obj.mu = mu;
      obj.r = numel(tau);
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
      X = varmapack_ou_gateway(obj.tau, obj.cov, obj.mu, full(double(t(:))), ...
                               M, rng.mex_handle());
    end
  end
end
