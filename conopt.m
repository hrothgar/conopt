% constrained optimization
% hrothgar, 20 july 2013
function [approx_abs, true_abs] = conopt
rng(1);
A = rand(2);
true_abs = 0.463436192938173; % true answer (from eigtool)

K   = kron([-2 -2i; -2i -2], eye(size(A)));
gam = norm(K);

ep = 0.01;          % epsilon-radius of the pseudospectrum
c  = [1; 0];        % objective function, maximize real(lam)
xx = [0.4; 0];      % initial point x0

% should be a while loop but I am just testing
for k = 2:16,
    x       = xx(:,k-1);
    [f,g]   = del(A, x);
    lambda  = norm(c)/sqrt(2*gam*(ep^2-f+1/2/gam*norm(g)^2));
    xx(:,k) = x + (c/lambda - g)/gam;
end

approx_abs = xx(1,end);
return

% eigenvalues and their derivatives
function [lam,dlam] = del(A, x),
z  = x(1) + 1i*x(2);
I  = eye(size(A));
A2 = (A - z*I)'*(A - z*I);

% just approximating the derivative for now
e   = 1e-9;
ew  = sortreal(eig((A-z*I)'*(A-z*I)));
ewr = sortreal(eig((A-(z+e)*I)'*(A-(z+e)*I)));
ewi = sortreal(eig((A-(z+1i*e)*I)'*(A-(z+1i*e)*I)));
ew  = ew(1);
ewr = ewr(1);
ewi = ewi(1);

lam  = real(ew);
dlam = [real(ewr-ew); imag(ewi-ew)]/e;
return

% sort by real part
function xx = sortreal(xx)
[~,indx] = sort(real(xx));
xx = xx(indx);
return