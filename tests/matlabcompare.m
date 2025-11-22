run(8,8)
function run(case1, caseN)
  fid = fopen("matlabcompare.inc", "w");
  fprintf(fid, "  int case1 = %d, caseN = %d;\n", case1, caseN);
  for k=case1:caseN
    if k < 12, ndec = 3; else, ndec = 15; end
    if k < 12, n = 2; else, n = 10; end
    [A, B, Sig, p, q, r] = testcase(k);
    rand_init('ParkMillerPolar', 42);
    X = new_varma_sim(A, B, Sig, n, 0, 1);
    printdims(fid, k-1, p, q, r, n)
    fprintf(fid, "  double\n");
    printmat(fid, "A", A, k-1, ndec, ",")
    printmat(fid, "B", B, k-1, ndec, ",")
    printmat(fid, "Sig", Sig, k-1, ndec, ",")
    printmat(fid, "X", X, k-1, 15, ";")
  end
  printcomb(fid, "int", "p", case1, caseN);
  printcomb(fid, "int", "q", case1, caseN);
  printcomb(fid, "int", "r", case1, caseN);
  printcomb(fid, "int", "n", case1, caseN);
  printcomb(fid, "double*", "A", case1, caseN);
  printcomb(fid, "double*", "B", case1, caseN);
  printcomb(fid, "double*", "Sig", case1, caseN);
  printcomb(fid, "double*", "X", case1, caseN);
  fclose(fid);
end

function printdims(fid, icase, p, q, r, n)
  fmt = "  int p%d = %d, q%d = %d, r%d = %d, n%d = %d;\n";
  fprintf(fid, fmt, icase, p, icase, q, icase, r, icase, n);
end

function printcomb(fid, type, letter, case1, caseN)
  list = strjoin(letter + string(case1-1:caseN-1), ",");
  fprintf(fid, '  %s %s[] = {%s};\n', type, letter, list);
end

function printmat(fid, matrixName, M, icase, ndec, lineend)
  if isempty(M)
    Mstring = "0";
  else
    Mfmt = "%." + ndec + "g";
    Mstring = strjoin(string(compose(Mfmt, M)), ",");
    if isempty(Mstring), Mstring = "0"; end
  end
  linefmt = "    %s%d[] = {%s}%s\n";
  fprintf(fid, linefmt, matrixName, icase, Mstring, lineend);
end
