function matlabcompare(cases, Ms, n)
  if nargin < 3, n = 5; end
  if nargin < 2, Ms = 1; end
  if nargin < 1 || isempty(cases), cases = [1:3,6:10,12]; end
  fid = fopen("../matlabcompare.txt", "w");
  mat2file(fid, length(Ms), "#Ms");
  mat2file(fid, length(cases), "#cases");
  mat2file(fid, n);
  mat2file(fid, Ms);
  mat2file(fid, cases);
  % folder = """" + pwd() + filesep + """";
  % cmd = "  strcpy(matmatfolder, " + folder + ");\n";
  % fprintf(fid, cmd);
  % char *matmatfolder[] = """ + pwd() + filesep + """");
  for k = cases
    icase = k + 1;
    [A, B, Sig, p, q, r] = testcase(k);    
    x0 = startmat(r, max(p,q));
    mat2file(fid, p, "p" + icase);
    mat2file(fid, q, "q" + icase);
    mat2file(fid, r, "r" + icase);
    mat2file(fid, A, "A" + icase);
    mat2file(fid, B, "B" + icase);
    mat2file(fid, Sig, "Sig" + icase);
    mat2file(fid, x0, "x0" + icase);
    for M = Ms
      rand_init('ParkMillerPolar', 42);
      [X, E, condR] = new_varma_sim(A, B, Sig, n, 0, M);
      [X0, E0]      = new_varma_sim(A, B, Sig, n, 0, M, x0);
      if M == Ms(1), mat2file(fid, condR, "condR-" + k); end
      mat2file(fid, X, "X" + M + "-" + k);
      mat2file(fid, E, "E" + M + "-" + k);
      mat2file(fid, X0, "X0" + M + "-" + k);
      mat2file(fid, E0, "E0" + M + "-" + k);
    end
  end
  fclose(fid);
end

function x0 = startmat(r, h)
  x0 = repmat([-2,-1,0,1,2], 1, ceil(r*h/5));
  x0 = reshape(x0(1:r*h), r, h);
end
