function basis_set = chebyshev_polynomials(x, number_of_polynomials)
%basis_set = chebyshev_polynomials(x, number_of_polynomials)
% it returns the first number_of_polynomials basis_set of chebyshev 
% polynomials calculated on x. x should be -1<=x<=1.

x = x(:);
basis_set = cos(acos(x)*(0:number_of_polynomials-1));

end
