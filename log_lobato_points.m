function [x, x_lin] = log_lobato_points(x1, x2, n_points)
%[x, x_lin] = log_lobato_points(x1, x2, n_points)
%   given x1 and x2 it returns x, the n_points chebyshev nodes spaced on
%   a log scale. x_lin are the linear chebyshev nodes usefull for
%   interpolation.

x_lin = cos((2*(n_points:-1:1)-1)/(2*n_points)*pi);

a = log10(x1);
b = log10(x2);

log_x = 0.5*(a+b) + 0.5*(b-a)*x_lin;

x = 10.^(log_x);

end