function [basisSet, varargout] = legendre_polynomials(x, number_of_polynomials)
% newer version where the number of coefficients/polynomials is passed
% instead of the order
% arguments
%     L {mustBeInteger}
%     number_of_polynomials {mustBeInteger}
%     sflag  {mustBeMember(sflag,{'periodic','symmetric'})} = 'periodic'
% end
% 
% if number_of_polynomials < 1
%     error('number_of_polynomials must be higher than 1, but it was: %d', number_of_polynomials)
% end

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

% %     x = -1:2/(L-1):1; % Legendre polynomials are supported for |x|<=1
% if strcmp(sflag, 'periodic') 
%     x = 2*((0:L-1)/L)-1;
% elseif strcmp(sflag, 'symmetric')
%     x = linspace(-1,1,L);
% end
x = x(:); % make it a column

varargout{1} = x;

% Evaluate the polynomials for every element in X
basisSet = cumprod([ones(size(x)) x(:,ones(1,order))], 2) * coeff.';
% or D = cumprod([ones(m,1) repmat(X,[1 N])], 2) * coeff.';
end
