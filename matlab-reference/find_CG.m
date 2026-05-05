%FIND_CG  Determine C and G matrices for VARMA simulation
%
%  [C,G] = FIND_CGW(A,B,Sig) calculates the matrices Ci and Gi for varma_sim. On
%  entry A=[A1 A2...Ap] and B=[B1 B2...Bq] are r × p·r and r × q·r arrays with
%  the r×r matrices Ai and Bi and Sig is also r × r. On exit
%                C = [C0 C1...Cq] with Cj = cov(x(t),eps(t-j)),
%                G = [G0 G1...Gq] with Gj = cov(y(t),x(t-j)), and 
%  These matrices are given by the formulae:
%                Cj = A1·C(j-1) + ... + A(j-1)·C1 + Aj·C0 + Bj·Sig 
%                Gj = Bj·C0' + ... + Bq·C(q-j)'
%  where C0 = Sig and B0 = I.

function [C, G] = find_CG(A,B,Sig)
  [p,q,r] = get_dimensions(A,B,Sig);
  C = zeros(r, r*(q+1));
  G = zeros(r, r*(q+1));
  J = 1:r;
  C(:, J) = Sig;
  for j=1:q
    C(:, J + r) = B(:, J)*Sig;
    I = 1:r;
    JmI = J;
    J = J + r;
    for i = 1:min(j,p)
      C(:, J) = C(:, J) + A(:,I)*C(:,JmI);
      I = I + r;
      JmI = JmI - r;
    end
  end
  J = 1:r;
  for j = 0:q
    I = J - r;
    ImJ = 1:r;
    G(:,J) = zeros(r,r);
    for i = j:q
      if i==0
        G(:, J) = G(:, J) + C(:,ImJ)';
      else
        G(:, J) = G(:, J) + B(:, I)*C(:, ImJ)'; 
      end
      I = I + r;
      ImJ = ImJ + r;
    end
    J = J + r;
  end
end