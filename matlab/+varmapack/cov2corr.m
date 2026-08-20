function corr = cov2corr(cov)
%VARMAPACK.COV2CORR  Convert covariance matrices to correlations.
%
%   corr = VARMAPACK.COV2CORR(cov) converts an r-by-r covariance matrix or
%   r-by-r-by-(maxlag+1) covariance sequence. Every entry is divided by the
%   product of the corresponding lag-zero marginal standard deviations.
%
%   corr has the same shape as cov. Its lag-zero diagonal entries are exactly
%   one. Other entries are not clipped to [-1,1].
%
%   See also VARMAPACK.ACVF, VARMAPACK.AUTOCOV.
  corr = varmapack_cov2corr_gateway(cov);
end
