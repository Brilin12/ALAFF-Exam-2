function [ x, niters ] = CG( A, b, x0 )
    x = x0;
    r = b - (A*x0);
    niters = 0;
    while not (norm(r) < eps(1.0) * norm(b))
        if niters == 0
            p = r;
        else
            g = (-(transpose(p)*(A*r))) / (transpose(p)*(A*p));
            p = r + (g * p);
        end
        q = A*p;
        a = (transpose(p) * r) / (transpose(p) * q);
        x = x + (a*p);
        r = r - (a*q);
        niters = niters + 1;
    end
end