run(2)
function run(ncase)
  folder = "~/dropbox/varma/varmasim/dev";
  rand_init('ParkMillerPolar');
  fid = fopen(folder + "/matlabcompare.txt", "w");
  fprintf(fid, "  int ncase = %d;\n", ncase);
  for k=1:ncase
    if k < 12, ndec = 3; else, ndec = 15; end
    if k < 12, n = 4; else, n = 6; end
    rand_init(42)
    [A, B, Sig, p, q, r] = testcase(k);
    X = varma_sim(A, B, Sig, n, 0, 1);
    printdims(fid, k-1, p, q, r, n)
    fprintf(fid, "  double\n");
    printmat(fid, "A", A, k-1, ndec, ",")
    printmat(fid, "B", B, k-1, ndec, ",")
    printmat(fid, "Sig", Sig, k-1, ndec, ",")
    printmat(fid, "X", X, k-1, 15, ";")
  end
  printcomb(fid, "int", "p", ncase);
  printcomb(fid, "int", "q", ncase);
  printcomb(fid, "int", "r", ncase);
  printcomb(fid, "int", "n", ncase);
  printcomb(fid, "double*", "A", ncase);
  printcomb(fid, "double*", "B", ncase);
  printcomb(fid, "double*", "Sig", ncase);
  printcomb(fid, "double*", "X", ncase);
  fclose(fid);
end

function printdims(fid, icase, p, q, r, n)
  fmt = "  int p%d = %d, q%d = %d, r%d = %d, n%d = %d;\n";
  fprintf(fid, fmt, icase, p, icase, q, icase, r, icase, n);
end

function printcomb(fid, type, letter, ncase)
  list = strjoin(letter + string(0:ncase-1), ",");
  fprintf(fid, '  %s %s[] = {%s};\n', type, letter, list);
end

function printmat(fid, matrixName, M, icase, ndec, lineend)
  Mfmt = "%." + ndec + "g";
  Mstring = strjoin(string(compose(Mfmt, M)), ",");
  linefmt = "    %s%d[] = {%s}%s\n";
  fprintf(fid, linefmt, matrixName, icase, Mstring, lineend);
end