%REF_VARMA_TESTCASEX  Create VARMAX testcases for reference simulation tests
%
%   [A,B,C,Sig,z] = REF_VARMA_TESTCASEX(NAME,n) returns a named VARMAX
%   testcase with an r by n exogenous sequence z.
%
%   [A,B,C,Sig,z,p,q,s,r] = REF_VARMA_TESTCASEX(NAME,n) also returns dimensions.
%
%   [A,B,C,Sig,z,name] = REF_VARMA_TESTCASEX(i,n) returns the i-th named
%   testcase name. To get dimensions use [A,B,C,Sig,z,p,q,s,r,name].
%
%   count = REF_VARMA_TESTCASEX('count') returns the number of named testcases.
%
%   REF_VARMA_TESTCASEX('summary',n) prints a list of the named testcases.

function [A, B, C, Sig, z, varargout] = ref_varma_testcasex(type, n)
  testcases = {
    'tinyARMAX'
    'tinyMAX'
    'tinyARMAX2'
    'smallARX1'
    'smallARX2'
    'smallMAX1'
    'smallMAX2'
    'smallARMAX1'
    'smallARMAX2'
    'mediumARX'
    'mediumMAX1'
    'mediumARMAX1'
    'mediumARMAX2'
    'mediumMAX2'
    'largeARX'
  };
  N = length(testcases);
  if nargin < 1, type = 1; end
  inputWasIndex = isnumeric(type);
  if isstring(type), type = char(type); end
  if ischar(type) && strcmpi(type, 'count')
    A = N;
    return
  end
  if ischar(type) && strcmpi(type, 'summary')
    if nargin < 2 || isempty(n), n = 10; end
    fprintf('No. Name          p  q  s  r  spec.rad\n');
    for i = 1:N
      [A, B, C, Sig, z, p, q, s, r, name] = ref_varma_testcasex(i, n);
      if p == 0, rho = 0; else, rho = ref_varma_specrad(A); end
      fprintf('%-2i  %-12s %2d %2d %2d %2d   %6.4f\n', i, name, p, q, s, r, rho);
    end
    clear A B C Sig z
    return
  end
  if isnumeric(type)
    if type < 1 || type > N, error('unknown testcase index'); end
    name = testcases{type};
  else
    name = type;
  end
  if nargin < 2 || isempty(n), error('n must be specified'); end
  icase = find(strcmpi(name, testcases));
  if isempty(icase), error('unknown testcase name'); end
  name = testcases{icase};
  if icase == 1
    A = 0.6;
    B = 0.2;
    C = 0.8;
    Sig = 1;
    p = 1; q = 1; s = 1; r = 1;
    z = (-1).^(0:n-1);
  else
    [A, B, Sig, p, q, r] = ref_varma_testcase(icase);
    s = max(p, q);
    C = xcoefficients(icase, s);
    z = zsequence(n);
  end
  if nargout > 5
    [varargout{1:4}] = deal(p, q, s, r);
  end
  if nargout == 6 || nargout >= 10
    if inputWasIndex
      extra = name;
    else
      extra = icase;
    end
    varargout{nargout-5} = extra;
  end
end

function C = xcoefficients(icase, s)
  switch icase
    case 2
      C = 0.70;
    case 3
      C = 0.80;
    case 4
      C = [0.70; -0.45];
    case 5
      C = [0.70 0.25; -0.45 0.18];
    case 6
      C = [0.65; 0.30];
    case 7
      C = [0.65 0.20; 0.30 -0.12];
    case 8
      C = [0.70; -0.45];
    case 9
      C = [0.70 0.25; -0.45 0.18];
    case 10
      C = [0.55; -0.40; 0.25];
    case 11
      C = [0.50; 0.25; -0.20];
    case 12
      C = [0.55 0.18 -0.08; -0.40 0.22 0.10; 0.25 -0.15 0.07];
    case 13
      C = [0.50 -0.20 0.09; 0.35 0.16 -0.08; -0.25 0.12 0.06];
    case 14
      C = [0.50 0.16; 0.25 -0.10; -0.20 0.08];
    case 15
      C = [
        0.45  0.18 -0.08  0.04;
       -0.35  0.14  0.09 -0.03;
        0.25 -0.12  0.06  0.02;
       -0.18  0.10 -0.05  0.03;
        0.12 -0.07  0.04 -0.02];
    otherwise
      error('missing x coefficients')
  end
  C = C(:, 1:s);
end

function z = zsequence(n)
  pattern = [3 -1 4 0 -5 2 1 -3 5 -2 0 4 -4 1 -1 2 -3]/5;
  z = pattern(mod(0:n-1, length(pattern)) + 1);
end
