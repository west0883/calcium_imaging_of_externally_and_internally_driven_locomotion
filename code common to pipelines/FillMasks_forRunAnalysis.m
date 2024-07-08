% FillMasks.m
% Sarah West
% 9/1/21

% Takes a vector of mask indices and a matrix of data and fills in the mask
% indices (with NaNs) that have been removed by masking.Returns the images
% in the 2D format (because using it in 2D is the only time you'd want to
% fill in the mask anyway).

function [parameters]=FillMasks_forRunAnalysis(parameters)
 % Inputs:
 % data_matrix-- must be a 2D matrix with the list of pixels in the 1st
 % dimension (each row is a different pixel). 
 % indices_of_mask -- a vector of indices where data SHOULD be (a positive
 % mask)
   
   parameters.data_matrix_filled = FillMasks(parameters.data, parameters.indices_of_mask...
       , parameters.pixels(2), parameters.pixels(1)); 
end 