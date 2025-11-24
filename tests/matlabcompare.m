function matlabcompare(cases)
  if nargin < 1, cases = [1:3,6:10,12]; end
  fid = fopen("matlabcompare.inc", "w");
  folder = """" + pwd() + filesep + """";
  cmd = "  strcpy(matmatfolder, " + folder + ");\n";
  fprintf(fid, cmd);
  char *matmatfolder[] = """ + pwd() + filesep + """");
  fprintf(fid, "  int cases[] = {%d", cases(1));
  for k = cases(2:end), fprintf(fid, ", %d", k); end
  fprintf(fid, "};\n");
  for k=cases
    fprintf("Case %d\n", k)
    if k < 12, ndec = 4; else, ndec = 15; end
    if k < 12, n = 3; else, n = 6; end
    [A, B, Sig, p, q, r] = testcase(k);
    rand_init('ParkMillerPolar', 42);
    X = new_varma_sim(A, B, Sig, n, 0, 1);
    printdims(fid, k, p, q, r, n)
    fprintf(fid, "  double\n");
    printmat(fid, "A", A, k, ndec, ",")
    printmat(fid, "B", B, k, ndec, ",")
    printmat(fid, "Sig", Sig, k, ndec, ",")
    printmat(fid, "X", X, k, 15, ";")
  end
  printcomb(fid, "int", "p", cases);
  printcomb(fid, "int", "q", cases);
  printcomb(fid, "int", "r", cases);
  printcomb(fid, "int", "n", cases);
  printcomb(fid, "double*", "A", cases);
  printcomb(fid, "double*", "B", cases);
  printcomb(fid, "double*", "Sig", cases);
  printcomb(fid, "double*", "X", cases);
  fclose(fid);
end

function printdims(fid, icase, p, q, r, n)
  fmt = "  int p%d = %d, q%d = %d, r%d = %d, n%d = %d;\n";
  fprintf(fid, fmt, icase, p, icase, q, icase, r, icase, n);
end

function printcomb(fid, type, letter, cases)
  list = strjoin(letter + string(cases), ",");
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
