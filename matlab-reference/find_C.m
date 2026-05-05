% FIND_CGW  Determine C matrices for VARMA log-likelihood
%
%  C = FIND_C(A,B,Sig,k) calculates the covariances Cj of xt and eps{t-k} in the
%  cell array C = {C0 C1...C{k-1}}. CC_build can subsequently be used to find CC 
%  = cov(x{1:k}, eps{1:k}):
%
%       CC = C0      0 ...      0
%            C1     C0          :
%            :         ..       :
%            :             ..   0
%            C{k-1} ...        C0

function C = find_C(A,B,Sig,k)
  [p,q,r] = get_dimensions(A,B,Sig);
  A = makecell(A);
  B = makecell(B);
  C = cell(1,k);
  C{1} = Sig;
  for j=1:k-1
    if j<= q
      C{j+1} = B{j}*Sig;
    else
      C{j+1} = zeros(r,r);
    end
    for i=1:min(j,p)
      C{j+1} = C{j+1} + A{i}*C{j-i+1};
    end
  end
end
