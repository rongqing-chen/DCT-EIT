function basis_set = legendre_polynomials(x, number_of_polynomials)
%basis_set = legendre_polynomials(x, number_of_polynomials)
% it returns the first number_of_polynomials basis_set of legendre 
% polynomials calculated on x. x should be -1<=x<=1.

order = number_of_polynomials-1;

% compute the Legendre polynomial coefficients matrix coeff
% coeff(i,j) gives the polynomial coefficient for term x^{j-1} in P_{i-1}(x)
if order > 1
    coeff = zeros(order+1);
    coeff([1 order+3]) = 1; % set coefficients of P_0(x) and P_1(x)
    % now compute for higher order: nP_n(x) = (2n-1)xP_{n-1}(x) - (n-1)P_{n-2}(x)
    for ii = 3:order+1
        coeff(ii,:) = (2-1/(ii-1))*coeff(ii-1,[end 1:end-1]) - (1-1/(ii-1))*coeff(ii-2,:);
    end
else
    % simple case
    coeff = eye(order+1);
end

x = x(:); % make it a column

% Evaluate the polynomials for every element in X
basis_set = cumprod([ones(size(x)) x(:,ones(1,order))], 2) * coeff.';
% or D = cumprod([ones(m,1) repmat(X,[1 N])], 2) * coeff.';
end
