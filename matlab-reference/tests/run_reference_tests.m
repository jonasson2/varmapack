% RUN_REFERENCE_TESTS  Test Varmapack reference m-functions
%
% RUN_REFERENCE_TESTS tests the Varmapack Matlab reference functions by calling all the
% defined individual tests, checking both the ref_varma_* functions (which have direct
% match in the varmapack_* C functions), as well as doing unit-test like checks of
% find_CG, S_build, CC_build, and the vector-Yule-Walker solver.
% 
% RUN_REFERENCE_TESTS(CASES) runs statistical tests of ref_varma_sim for the specified
% test cases. Default is 1:15 to test all the named test cases.
%
% RUN_REFERENCE_TESTS(CASES, M) computes M replicates for the statistical tests (default
% is M=10000)
%
% RUN_REFERENCE_TESTS(CASES, M, n) computes series of length n (default 20) for the
% statiistical tests.
%
% RUN_REFERENCE_TESTS("table",...) displays a table with the results of the
% statistical tests.

function run_reference_tests(varargin)
  showTable = false;
  if nargin >= 1 && (ischar(varargin{1}) || isstring(varargin{1})) && ...
      strcmpi(varargin{1}, 'table')
    showTable = true;
    varargin = varargin(2:end);
  end
  if length(varargin) >= 1, cases = varargin{1}; else, cases = []; end
  if length(varargin) >= 2, M = varargin{2}; else, M = []; end
  if length(varargin) >= 3, n = varargin{3}; else, n = []; end
  if nargin < 1 || isempty(cases), cases = 1:15; end
  if isempty(M), M = 10000; end
  if isempty(n), n = 20; end
  here = fileparts(mfilename('fullpath'));
  addpath(here);
  addpath(fileparts(here));
  addpath(fullfile(fileparts(fileparts(here)), 'tests', 'matlab'));
  test_ref_varma_sim_first;
  test_ref_varma_testcase;
  test_ref_varma_sim(cases, M, n, showTable);
  test_ref_acvf_autocov_specrad;
  test_ref_psi_irf;
  test_ref_varma_simx;
  test_ref_find_CG;
  test_ref_S_build;
  test_ref_CC_build;
  test_ref_vyw;
  fprintf('matlab-reference tests OK\n');
end
