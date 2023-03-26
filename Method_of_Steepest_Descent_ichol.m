function [ x, niters ] = Method_of_Steepest_Descent_ichol( A, b, x0 )
    L = sparse(ichol(sparse(A)));
    invL = L \ eye(size(b,1));
    invM = transpose(invL) * invL;
    x = x0;
    r = b - (A*x0);
    niters = 0;
    while not (norm(r) < eps(1.0) * norm(invM*b))
        p = invM*r;
        q = A*p;
        a = (transpose(p) * r) / (transpose(p) * q);
        x = x + (a*p);
        r = r - (a*q);
        niters = niters + 1;
    end
end