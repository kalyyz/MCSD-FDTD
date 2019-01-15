% Matlab example with bi-complex number

[j1 j2] = initialize;
f  = @(x) 4*x*x*x + 2/x -1;
F = f(1+j1+j2);

%to find real and imaginary parts

real(F)

im1(F)

im2(F)

im12(F)

%second-order derivative at f(1)

h = 1e-10;
dd_f = im12( f(x+j1h+j2h) )/(h*h)
