function t = gamma_times(n, shape, mean_gap, rng)
%VARMAPACK.GAMMA_TIMES  Observation times with independent Gamma gaps.
%
%   t = varmapack.gamma_times(n,shape,mean_gap,rng) returns an n-by-1 vector
%   starting at zero. The gaps have shape (default 16), scale mean_gap/shape,
%   and CV 1/sqrt(shape). mean_gap defaults to 1. The realized mean spacing
%   is not normalized. shape=Inf gives constant gaps. Unrepresentable positive
%   gaps or times produce an error. rng is an optional varmapack.Rng.
%
%   See also varmapack.OUModel.
  if nargin < 2 || isempty(shape), shape = 16; end
  if nargin < 3 || isempty(mean_gap), mean_gap = 1; end
  if nargin < 4 || isempty(rng)
    rng = varmapack.Rng();
    cleanup = onCleanup(@() delete(rng));
  end
  if ~isa(rng, 'varmapack.Rng') || ~isscalar(rng)
    error('varmapack:gamma_times:rng', 'rng must be a varmapack.Rng object');
  end
  t = varmapack_gamma_times_gateway(n, shape, mean_gap, rng.mex_handle());
end
