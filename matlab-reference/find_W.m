%FIND_W  Determine covariance blocks for the MA(q) residual process
%
%  W = FIND_W(B,Sig) calculates W = [W0 W1...Wq], where
%
%    Wk = cov(y(t),y(t-k))
%
%  for y(t) = eps(t) + B1*eps(t-1) + ... + Bq*eps(t-q), with
%  eps(t) having covariance Sig.

function W = find_W(B, Sig)
  r = size(Sig, 1);
  if isempty(B), B = zeros(r,0); end
  q = size(B, 2)/r;
  Bc = [{eye(r)} makecell(B)];
  BSig = cell(1, q+1);
  Wc = cell(1, q+1);
  for i = 0:q
    BSig{i+1} = Bc{i+1}*Sig;
  end
  for k = 0:q
    Wc{k+1} = BSig{k+1};
    for i = 1:q-k
      Wc{k+1} = Wc{k+1} + BSig{i+k+1}*Bc{i+1}';
    end
  end
  W = cell2mat(Wc);
end
