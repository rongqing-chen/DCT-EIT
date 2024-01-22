function basisSet = chebyshev_polynomials(x, number_of_polynomials)


basisSet = cos(acos(x)*(0:number_of_polynomials-1));

end
