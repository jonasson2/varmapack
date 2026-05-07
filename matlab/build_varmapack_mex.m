function build_varmapack_mex
  root = fileparts(fileparts(mfilename('fullpath')));
  rp = "/Users/jonasson/randompack";
  inc = fullfile(root, "src");
  lib = fullfile(root, "build", "src");
  rpinc = fullfile(rp, "src");
  rplib = fullfile(rp, "release", "src");
  gateways = [
    "varmapack_sim_gateway"
    "varmapack_specrad_gateway"
    "varmapack_acvf_gateway"
    "varmapack_autocov_gateway"
    "varmapack_testcase_gateway"
  ];
  for k = 1:numel(gateways)
    mex("-I" + inc, "-I" + rpinc, "-L" + lib, "-lvarmapack", ...
        "-L" + rplib, "-lrandompack", ...
        "LDFLAGS=$LDFLAGS -Wl,-rpath," + lib + " -Wl,-rpath," + rplib, ...
        fullfile(root, "matlab", gateways(k) + ".c"), ...
        "-outdir", fullfile(root, "matlab"));
  end
end
