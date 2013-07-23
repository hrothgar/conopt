% Example for computing pseudospectral abscissa.
%
%   >> ex_psabs
%   f =
%        1.000000000000063e-04
%   x =
%      0.463436192938173
%                      0
%   true_abs =
%      0.463436192938173

rng(1);
pars.A = rand(2);               % matrix to calc the psabs of
pars.ep = 1e-2;                 % epsilon-radius of the pseudospectrum
K = kron(-[2 2i; 2i 2],eye(2)); % 
gam = norm(K);                  % bound on the hessian
c  = [1; 0];                    % objective function, maximize real(lam)
x0 = [0.4; 0];                  % initial point
tol = 1e-15;                    % we can it real good

[f,x] = conopt(@f_psabs, pars, c, gam, x0, tol)
true_abs = 0.463436192938173    % true answer (from eigtool)