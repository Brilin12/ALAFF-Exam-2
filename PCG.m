function [ x, niters ] = PCG( A, b, x0 )
    L = sparse(ichol(sparse(A)));
    invL = L \ eye(size(b,1));
    invM = transpose(invL) * invL;
    x = x0;
    r = b - (A*x0);
    niters = 0;
    z = invM*r;
    while not (norm(r) < eps(1.0) * norm(invM*b))
        z_prev = z;
        z = invM*r;
        if niters == 0
            p = z;
        else
            g = (transpose(r)*z) / (transpose(r_prev)*z_prev);
            p = z + (g * p);
        end
        q = A*p;
        a = (transpose(r) * z) / (transpose(p) * q);
        x = x + (a*p);
        r_prev = r;
        r = r - (a*q);
        niters = niters + 1;
    end
end