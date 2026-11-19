classdef SlepianModel
  %VARMAPACK.SLEPIANMODEL  Stationary Slepian-based moving-average model.
  %
  %   model = varmapack.SlepianModel(theta,cov,mu) takes the real unit-step
  %   MA(1) coefficient matrix theta, stationary covariance cov, and an
  %   optional constant mean vector mu (default zero). Only the lower
  %   triangle of cov is used. The driving covariance Sig must be positive
  %   semidefinite and determined by cov = Sig + theta*Sig*theta'.
  %
  %   X = model.sim(t,M,rng) simulates exact stationary paths without burn-in
  %   using unit-width Brownian noise windows. t is a finite strictly
  %   increasing vector. M defaults to 1 and rng is an optional
  %   varmapack.Rng. X has shape r-by-numel(t)-by-M.
  %
  %   See also varmapack.OUModel, varmapack.gamma_times.
  properties (SetAccess = private)
    theta
    cov
    mu
    r
  end
  properties (Access = private)
    factor
  end
  methods
    function obj = SlepianModel(theta, cov, mu)
      r = size(theta, 1);
      if ~isnumeric(theta) || ~isreal(theta) || r < 1 || ...
          ~isequal(size(theta), [r r])
        error('varmapack:SlepianModel:theta', 'theta must be a nonempty square matrix');
      end
      if ~isnumeric(cov) || ~isreal(cov) || ~isequal(size(cov), [r r])
        error('varmapack:SlepianModel:cov', 'cov must be an r-by-r real matrix');
      end
      if nargin < 3 || isempty(mu), mu = zeros(r, 1); end
      if ~isnumeric(mu) || ~isreal(mu) || ~isvector(mu) || numel(mu) ~= r
        error('varmapack:SlepianModel:mu', 'mu must be a real vector of length r');
      end
      theta = full(double(theta));
      cov = full(double(cov));
      mu = full(double(mu(:)));
      obj.factor = varmapack_slepian_gateway(theta, cov, mu);
      obj.theta = theta;
      obj.cov = tril(cov) + tril(cov, -1)';
      obj.mu = mu;
      obj.r = r;
    end
    function G = acvf(obj, lags)
      %ACVF  Cov(X(t+h), X(t)) at real lags h, as r-by-r-by-numel(lags).
      if ~isnumeric(lags) || ~isreal(lags) || ~(isvector(lags) || isempty(lags))
        error('varmapack:SlepianModel:lags', 'lags must be a real vector');
      end
      G = varmapack_slepian_gateway(obj.theta, obj.factor, obj.mu, obj.cov, ...
                                    full(double(lags(:))));
    end
    function X = sim(obj, t, M, rng)
      %SIM  Simulate M independent stationary paths at observation times t.
      if nargin < 3 || isempty(M), M = 1; end
      if ~isnumeric(t) || ~isreal(t) || ~isvector(t) || isempty(t)
        error('varmapack:SlepianModel:t', 't must be a nonempty real vector');
      end
      if nargin < 4 || isempty(rng)
        rng = varmapack.Rng();
        cleanup = onCleanup(@() delete(rng));
      end
      if ~isa(rng, 'varmapack.Rng') || ~isscalar(rng)
        error('varmapack:SlepianModel:rng', 'rng must be a varmapack.Rng object');
      end
      X = varmapack_slepian_gateway(obj.theta, obj.factor, obj.mu, ...
                                    full(double(t(:))), M, rng.mex_handle());
    end
  end
end
