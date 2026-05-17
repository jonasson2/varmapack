% Let E0 = [eps0...eps{h–1}] where each epst is N(0, Sigma). From the observed x’s and
% known z’s we compute:
%
%   rt = xt – sum{i=1:p} Ai*x{t–i} – sum{i=0:s} Ci*z{t–i}, for t = t0...h–1,
%
% where t0 = max(p,s). The model implies:
%
%   rt = epst + sum{i=1:q} Bi*eps{t–i}.
%
% Stacking these equations gives H*E0 = r, where r = [r{t0}...r{h–1}]. Thus
%
%   Y = H*E0
%
% is Gaussian with covariance W = H*D*H’. Equivalently, W is the block Toeplitz covariance
% matrix of the MA(q) process:
%
%   rt = epst + sum{i=1:q} Bi*eps{t–i}.
%
% Its lag–k covariance blocks are  Wk = sum{i=0:q–k} Bi*Sigma*B{i+k}’, where B{0} = I. Now
%
%   E0
%   Y
%
% is jointly Gaussian with covariance
%
%    D   D*H’
%    H*D   W
%
% where D is the block–diagonal matrix with Dii = Sigma. Conditioning on the observed
% value Y = r gives E0hat = D*H’*(W\r) and:
%
%   Sigma0 = D – D*H’*W\H*D.
%
% Thus:
%
%   E0 | H*E0=r ~ N(E0hat, Sigma0).
