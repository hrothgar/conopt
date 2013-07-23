% constrained optimization
% hrothgar, 23 july 2013
function [f, x] = conopt(func, pars, c, gam, x0, tol)

xx = x0(:);
k = 1;
lastdiff = Inf;

% should be a while loop but I am just testing
while norm(lastdiff) > tol,
    k = k + 1;
    x = xx(:,k-1);
    [f,g]  = func(x, pars);

    % FIXME: still not general code, but it's closer
    lambda = norm(c)/sqrt(2*gam*(pars.ep^2-f+1/2/gam*norm(g)^2));

    % FIXME: we actually have to check whether to do +/- lambda
    xx(:,k) = x + (c/lambda - g)/gam;

    lastdiff = xx(:,k) - xx(:,k-1);
end

x = xx(:,k);
f = func(x, pars);

return
