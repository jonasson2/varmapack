function test_ref_find_CG
  fprintf('TESTING FIND_CG... ');
  A1 = [1 1; 1 1];
  A2 = [2 2; 1 1];
  B1 = [1 1; 2 2];
  B2 = [2 2; 2 2];
  B3 = [3 3; 2 2];
  Sig = [3 1; 1 3];
  [C, G] = find_CG([A1, A2], [B1, B2, B3], Sig);
  Cb = blocks(C);
  Gb = blocks(G);
  ascertain(length(Gb) == 4);
  ascertain(almostequal(Gb{1}, Sig + B1*Cb{2}' + B2*Cb{3}' + B3*Cb{4}'));
  ascertain(almostequal(Gb{2}, B1*Sig + B2*Cb{2}' + B3*Cb{3}'));
  ascertain(almostequal(Gb{3}, B2*Sig + B3*Cb{2}'));
  ascertain(almostequal(Gb{4}, B3*Sig));
  [C, G] = find_CG([], [B1, B2, B3], Sig);
  Cb = blocks(C);
  Gb = blocks(G);
  ascertain(almostequal(Cb{1}, Sig));
  ascertain(length(Gb) == 4);
  ascertain(almostequal(Gb{1}, Sig + B1*Cb{2}' + B2*Cb{3}' + B3*Cb{4}'));
  [C, G] = find_CG([A1, A2], [], Sig);
  Cb = blocks(C);
  Gb = blocks(G);
  ascertain(length(Cb) == 1 && length(Gb) == 1);
  ascertain(almostequal(Cb{1}, Sig));
  ascertain(almostequal(Gb{1}, Sig));
  fprintf('OK\n');
end

function C = blocks(A)
  r = size(A, 1);
  C = makecell(A);
  if isempty(C), C = {zeros(r, 0)}; end
end
