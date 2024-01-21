function [M] = zigzag_matrix(input_matrix)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
% https://stackoverflow.com/questions/3024939/matrix-zigzag-reordering

M = input_matrix;

ind = reshape(1:numel(M), size(M));         %# indices of elements
ind = fliplr( spdiags( fliplr(ind) ) );     %# get the anti-diagonals
ind(:,1:2:end) = flipud( ind(:,1:2:end) );  %# reverse order of odd columns
ind(ind==0) = [];                           %# keep non-zero indices

M = M(ind);
end