function [M] = zigzag_matrix(input_matrix)
%[M] = zigzag_matrix(input_matrix)
% return the values of input_matrix in a zigzag order. It is useful for
% having a more or less constant order traversing of the 2D coefficients of
% a tensor product. 
% Found on 
% https://stackoverflow.com/questions/3024939/matrix-zigzag-reordering

M = input_matrix;

ind = reshape(1:numel(M), size(M));         %# indices of elements
ind = fliplr( spdiags( fliplr(ind) ) );     %# get the anti-diagonals
ind(:,1:2:end) = flipud( ind(:,1:2:end) );  %# reverse order of odd columns
ind(ind==0) = [];                           %# keep non-zero indices

M = M(ind);
end