function run_reference_tests
  here = fileparts(mfilename('fullpath'));
  addpath(here);
  addpath(fileparts(here));
  addpath(fullfile(fileparts(fileparts(here)), 'tests', 'matlab'));
  test_ref_varma_sim_first;
  test_ref_varma_testcase;
  test_ref_varma_sim;
  test_ref_find_CG;
  test_ref_S_build;
  test_ref_vyw;
  fprintf('matlab-reference tests OK\n');
end
