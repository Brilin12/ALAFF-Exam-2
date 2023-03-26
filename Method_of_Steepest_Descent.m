function [ x, niters ] = Method_of_Steepest_Descent( A, b, x0)
    x = x0;
    r = b - (A*x0);
    niters = 0;
    while not (norm(r) < eps(1.0) * norm(b))
        p = r;
        q = A*p;
        a = (transpose(p) * r) / (transpose(p) * q);
        x = x + (a*p);
        r = r - (a*q);
        niters = niters + 1;
    end
end