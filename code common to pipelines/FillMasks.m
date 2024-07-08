% FillMasks.m
% Sarah West
% 9/1/21

% Takes a vector of mask indices and a matrix of data and fills in the mask
% indices (with NaNs) that have been removed by masking.Returns the images
% in the 2D format (because using it in 2D is the only time you'd want to
% fill in the mask anyway).

function [data_matrix_filled]=FillMasks(data_matrix, indices_of_mask, yDim, xDim)
 % Inputs:
 % data_matrix-- must be a 2D matrix with the list of pixels in the 1st
 % dimension (each row is a different pixel). 
 % indices_of_mask -- a vector of indices where data SHOULD be (a positive
 % mask)
 
    % Create a matrix of NaNs 
    data_matrix_filled=NaN(yDim*xDim, size(data_matrix,2)); 
            
    % Fill in missing pixels with NaNs 
   data_matrix_filled(indices_of_mask,:)=data_matrix;
   
   % Reshape into the desired stack of 2D images
   data_matrix_filled=reshape(data_matrix_filled, yDim, xDim, size(data_matrix,2)); 
end 