function [r,p,q,h] = getdims(A, B, Sig)
  if nargin >= 3
    r = size(Sig, 1);
  elseif ~isempty(A)
    r = size(A, 1);
  elseif ~isempty(B)
    r = size(B, 1);
  else
    r = 1;
  end
  p = size(A,2)/r;
  q = size(B,2)/r;
  h = max(p,q);
end
