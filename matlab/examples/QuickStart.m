% Basic VAR simulation
vrng = varmapack.Rng(123);
cleanup = onCleanup(@() delete(vrng));
n = 200;
M = 10;
Sig1 = [2, 0; 0, 1];
A = [0.6, 0.1; 0, 0.4];
X = varmapack.sim(A, [], Sig1, [], n, 1, [], vrng);
fprintf("VAR simulation, first five values:\n");
disp(X(:, 1:5))

% Named testcases
fprintf("\nNamed testcases:\n");
varmapack.testcase('summary')
[A, B, Sig] = varmapack.testcase('smallARMA1');
X = varmapack.sim(A, B, Sig, [], n, M, [], vrng);
fprintf("smallARMA1 simulation size: %s\n", mat2str(size(X)));

% Time-dependent means and fixed startup paths
mu = [0, 1; 0, 1];
X0 = [0, 0.1; 0, -0.1];
X0multi = cat(3, X0, X0 + [0.05, 0.05; -0.05, -0.05]);
Xcommon = varmapack.sim(A, B, Sig, mu, n, M, X0, vrng);
Xseparate = varmapack.sim(A, B, Sig, mu, n, 2, X0multi, vrng);
fprintf("fixed-start simulation sizes: %s and %s\n", ...
        mat2str(size(Xcommon)), mat2str(size(Xseparate)));

% VAR, VMA, and VARMA simulation
A1 = [0.6, 0.1; 0, 0.4];
B1 = [0.4, 0.3; 0.1, 0.2];
A2 = 0.2*ones(2);
B2 = 0.5*eye(2);
X1 = varmapack.sim(A1, [], Sig, [], n, M, [], vrng);
X2 = varmapack.sim([], [B1, B2], Sig, [], n, M, [], vrng);
[X3, E] = varmapack.sim([A1, A2], B1, Sig, [], n, M, [], vrng);
fprintf("VAR, VMA, VARMA, and shock sizes: %s, %s, %s, %s\n", ...
        mat2str(size(X1)), mat2str(size(X2)), mat2str(size(X3)), ...
        mat2str(size(E)));

% VARMAX simulation with common and replicate-specific vector exogenous inputs
C = [0.8, -0.3; 0.2, 0.4];
z = [cos((1:n)/10); sin((1:n)/10)];
zmulti = cat(3, z, z + 0.01);
[X4, E4] = varmapack.simx(A1, B1, C, z, Sig, n, M, zeros(2), 2, vrng);
X5 = varmapack.simx(A1, B1, C, zmulti, Sig, n, 2, zeros(2), 2, vrng);
fprintf("VARMAX and shock sizes: %s, %s, %s\n", ...
        mat2str(size(X4)), mat2str(size(E4)), mat2str(size(X5)));

% Model analysis
maxlag = 3;
rho = varmapack.specrad([A1, A2]);
rhoMA = varmapack.ma_specrad(B1);
Gamma = varmapack.acvf([A1, A2], B1, Sig, maxlag);
GammaHat = varmapack.autocov(X3(:, :, 1), maxlag, "ML");
GammaHatC = varmapack.autocov(X3(:, :, 1), maxlag, "C");
Theta = varmapack.irf([A1, A2], B1, Sig, maxlag);
Psi = varmapack.psi([A1, A2], B1, maxlag);
fprintf("\nVARMA model analysis:\n");
fprintf("AR spectral radius: %.6f\n", rho);
fprintf("MA spectral radius: %.6f\n", rhoMA);
fprintf("Psi, Theta, and covariance sizes: %s, %s, %s, %s, %s\n", ...
        mat2str(size(Psi)), mat2str(size(Theta)), mat2str(size(Gamma)), ...
        mat2str(size(GammaHat)), mat2str(size(GammaHatC)));
