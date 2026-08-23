% REF_VARMA_TESTCASE  Create testcases for VARMA calculations
%
%   [A,B,Sig] = REF_VARMA_TESTCASE(NAME) returns a named testcase in matrices
%   A = [A1...Ap], B = [B1...Bq], and Sig, suitable for VARMA calculations.
%   REF_VARMA_TESTCASE('summary') prints a list of the named testcases with
%   possible values for NAME.
%
%   [A,B,Sig,p,q,r] = REF_VARMA_TESTCASE(NAME) also returns the dimensions of
%   the case.
%
%   [A,B,Sig,name] = REF_VARMA_TESTCASE(i) returns the one-based i-th named
%   testcase. To also get dimensions, use [A,B,Sig,p,q,r,name] =
%   REF_VARMA_TESTCASE(i).
%
%   n = REF_VARMA_TESTCASE('COUNT') returns the number of named testcases.
%
%   REF_VARMA_TESTCASE('summary') prints a list of the named testcases.
%
%   [A,B,Sig] = REF_VARMA_TESTCASE(p,q,r) returns a deterministic unnamed
%   testcase with dimensions p,q,r.
%
%   [A,B,Sig] = REF_VARMA_TESTCASE(p,q,r,rng) returns a random unnamed testcase
%   using the supplied Randompack RNG.

function [A,B,Sig,varargout] = ref_varma_testcase(varargin)
  p15 = 5; r15 = 7; seed15 = 42; c15 = 0.05;
  name = '';
  if nargin == 1
    type = varargin{1};
  else
    p = varargin{1};
    q = varargin{2}; 
    r = varargin{3};
    type = 'unnamed';
  end
  inputWasIndex = isnumeric(type);
  testcases = {
    %              p  q  r       
    'tinyAR'     % 1  0  1
    'tinyMA'
    'tinyARMA'
    'smallAR1'
    'smallAR2'
    'smallMA1'
    'smallMA2'
    'smallARMA1'
    'smallARMA2'
    'mediumAR'
    'mediumMA1'
    'mediumARMA1'
    'mediumARMA2'
    'mediumMA2'
    'largeAR'
    'largeARMA'
    %'pivotfailure'
  };
  N = length(testcases);
  if isstring(type), type = char(type); end
  if ischar(type) && strcmpi(type, 'count'), type = 'count'; end
  if ischar(type) && strcmpi(type, 'summary'), type = 'summary'; end
  if ~ischar(type), name = testcases{type}; 
  else, name=type; end
  A33 = [
    0.15 0.10 0.05  0.11 0.14 0.17  0.01 0.04 0.06;
    0.16 0.11 0.06  0.12 0.15 0.18  0.02 0.05 0.08;
    0.17 0.12 0.07  0.13 0.16 0.19  0.03 0.06 0.09];
  switch type
    case 'count'
      A = N; % count of named cases
      return
    case 'summary'
      fprintf('No. Name          p  q  r  spec.rad\n');
      for i = 1:N
        fmt = '%-2i  %-12s %2d %2d %2d   %6.4f\n';
        [A,B,Sig,name] = ref_varma_testcase(i);
        [p,q,r] = get_dimensions(A,B,Sig);
        if p == 0, rho = 0; else, rho = ref_varma_specrad(A); end
        fprintf(fmt, i, name, p, q, r, rho);
      end
      clear A B Sig name
      return
    case 'unnamed'
      Sig = hilb(r) + 0.2*eye(r);
      if nargin >= 4 && ~isempty(varargin{4})
        rng = varargin{4};
        A = zeros(r, r*p);
        B = zeros(r, r*q);
        if p > 0
          A = reshape(randompack_u01(rng, r*r*p), r, r*p)/(2*p*r);
          while ref_varma_specrad(A) >= 1, A=A/2; end
        end
        if q > 0
          B = reshape(randompack_u01(rng, r*r*q), r, r*q)/(q*r);
        end
      else
        A = zeros(r, r*p);
        B = zeros(r, r*q);
        if p > 0, A(:) = 0.5/(p*r); end
        if q > 0, B(:) = 1/(q*r); end
      end
    case {1,'tinyAR'}
      A = 0.5;
      B = [];
      Sig = 0.8;
    case {2,'tinyMA'}
      A = [];
      B = 0.5;
      Sig = 0.8;
    case {3,'tinyARMA'}
      A = 0.4;
      B = 0.4;
      Sig = 0.8;
    case {4,'smallAR1'}
      %A = [[0.3 0.1; 0.4 0.2],[0.1 0.1;0.1 0.2]];
      A = [0.1 0.1; 0.1 0.100000001];
      B = [];
      Sig = [2 1;1 3];
    case {5,'smallAR2'}
      A = [0.60 0.05 0.30 0.02; 0.04 0.60 0.02 0.30];
      B = [];
      Sig = [2 1;1 3];
    case {6,'smallMA1'}
      A = [];
      B = [0.3 0.1; 0.1 0.2];
      Sig = [2 1;1 2];
    case {7,'smallMA2'}
      A = [];
      B = [0.3 0.1 0.1 0.1; 0.3 0.1 0.1 0.1];
      Sig = [2 1;1 2];
    case {8,'smallARMA1'}
      A = [0.3 0.1; 0.4 0.2];
      B = [0.2 0.3; 0.2 0.3];
      Sig = [2 1;1 2];
    case {9,'smallARMA2'}
      B = [[0.4 0.1; 0.2 0.1], [0.2 0.1; 0.3 0.1]];
      A = [0.2 0.3; 0.2 0.3];
      Sig = [2 1;1 2];
    case {10,'mediumAR'}
%       A = [
%         0.35 0.15 0.15  0.11 0.14 0.07;
%         0.25 0.15 0.05  0.12 0.15 0.08;
%         0.15 0.05 0.01  0.13 0.16 0.09];
      A = [
        0.35 0.15 0.15;
        0.25 0.15 0.05;
        0.15 0.05 0.01];
      B = [];
      Sig=[2.0  0.5 0.0;
           0.5  2.0 0.5;
           0.0  0.5 1.0];
    case {11,'mediumMA1'}
      A = [];
      B = [
        0.35 0.25 0.15;
        0.25 0.15 0.05;
        0.15 0.05 0.01];
      Sig=[2.0  0.5 0.0;
           0.5  2.0 0.5;
           0.0  0.5 1.0];
    case {12,'mediumARMA1'}
      A = A33; 
      B = fliplr(A33);
      Sig=[2.0  0.5 0.0;
           0.5  2.0 0.5;
           0.0  0.5 1.0];
    case {13,'mediumARMA2'}
      A = fliplr(A33);
      A(:, 4:9) = A(:, 4:9) + 0.01;
      B = A33;
      Sig=[2.0  0.5 0.0;
           0.5  2.0 0.5;
           0.0  0.5 1.0];
    case {14,'mediumMA2'}
      A = [];
      B = [
        0.35 0.25 0.15  0.11 0.14 0.17;
        0.25 0.15 0.05  0.12 0.15 0.18;
        0.15 0.05 0.01  0.13 0.16 0.19];
      Sig=[2.0  0.5 0.0;
           0.5  2.0 0.5;
           0.0  0.5 1.0];
    case {15,'largeAR'}
      r = r15;
      p = p15;
      rng = randompack_create();
      cleanup = onCleanup(@() randompack_free(rng));
      randompack_seed(rng, seed15);
      A = reshape(randompack_u01(rng, r*r*p), r, r*p)*c15;
      Sig = hilb(r) + eye(r);
      B = [];
    case {16,'largeARMA'}
      p = 3;
      q = 3;
      r = 7;
      A = rho_coefficients(p, r, 0.98);
      B = ones(r, r*q)/(q*r);
      Sig = hilb(r) + 0.2*eye(r);
    case {'pivotfailure'} % create almost singular vyw equations
      A = [
        -0.6250  0.2400    0.2500    0.1250    0.6250    0.3750
        0             0    0.1250    0.5000    0.2500    0.1250
        0.2500        0    0.1430         0    0.2500    0.2500];
      A(1,1) = -0.62999813363195223; %% as singular as can be made
      B=[];
      Sig = eye(3);
    otherwise
      error('no such testcase')
  end
  if nargout > 4
    r = size(Sig,1);
    p = size(A,2)/r;
    q = size(B,2)/r;
    [varargout{1:3}] = deal(p, q, r);
  end
  if nargout == 4 || nargout >= 7
    k = find(strcmp(testcases, name));
    if inputWasIndex
      extra = name;
    else
      extra = k;
    end
    varargout{nargout-3} = extra;
  end
  ascertain(isequal(Sig,Sig'));
  ascertain(min(eig(Sig))>0);
end

function A = rho_coefficients(p, r, target)
  scale = 2/(r*p^(1/3));
  A = zeros(r, p*r);
  for j = 1:p*r
    for i = 1:r
      A(i,j) = scale*min(i,j)/max(i,j);
    end
  end
  lo = 0;
  hi = 1;
  while ref_varma_specrad(hi*A) <= target
    hi = 2*hi;
  end
  for iter = 1:100
    mid = (lo + hi)/2;
    rho = ref_varma_specrad(mid*A);
    if abs(rho - target) < 1e-6
      A = mid*A;
      return
    end
    if rho < target
      lo = mid;
    else
      hi = mid;
    end
  end
  A = ((lo + hi)/2)*A;
end
