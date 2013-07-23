
% for computing the pseudospectral abscissa
% eigenvalues and their derivatives
function [f, g] = f_psabs(x, pars),

A  = pars.A;
z  = x(1) + 1i*x(2);
I  = eye(size(A));
A2 = (A - z*I)'*(A - z*I);

% FIXME: just approximating the derivative for now
%        which doesn't work when pars.ep shrinks much
e   = 1e-10;
ew  = sortreal(eig((A-z*I)'*(A-z*I)));
ewr = sortreal(eig((A-(z+e)*I)'*(A-(z+e)*I)));
ewi = sortreal(eig((A-(z+1i*e)*I)'*(A-(z+1i*e)*I)));
ew  = ew(1);
ewr = ewr(1);
ewi = ewi(1);

f = real(ew);
g = [real(ewr-ew); imag(ewi-ew)]/e;
return

% sort by real part
function xx = sortreal(xx)
[~,indx] = sort(real(xx));
xx = xx(indx);
return
