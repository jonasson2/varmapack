%W_BUILD  Build block Toeplitz covariance matrix from MA covariance blocks
%
%  WW = W_BUILD(W,n) builds the covariance matrix of n consecutive values of
%  the MA(q) process whose lag blocks are W = [W0 W1...Wq].

function WW = W_build(W, n)
  r = size(W, 1);
  q = size(W, 2)/r - 1;
  Wc = makecell(W);
  WW = cell(n, n);
  for i = 1:n
    for j = 1:n
      k = abs(i - j);
      if k <= q
        if i >= j
          WW{i,j} = Wc{k+1};
        else
          WW{i,j} = Wc{k+1}';
        end
      else
        WW{i,j} = zeros(r, r);
      end
    end
  end
  WW = cell2mat(WW);
end
