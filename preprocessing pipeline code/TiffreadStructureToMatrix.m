% TiffreadStructureToMatrix.m
% Sarah West
% 8/26/21

% Takes the output structure from tiffread functions and turns the .data
% fields into one big matrix. 

function [data_matrix]=TiffreadStructureToMatrix(data_structure, frames_list)

% Inputs:
% data_structure -- a structure of images outputted by tiffread
% frames_list -- the indices/frames of the structure you want to include in the
    % data matrix. Is a vector. 
    
% Outputs:
% data_matrix -- a 3D matrix of an image stack (pixels, pixels, frames)

    % Use the first frame of the structure to determine image dimensions 
    [yDim, xDim]=size(data_structure(1).data); 
    
    data_matrix = cell2mat({data_structure(frames_list).data});
    data_matrix = reshape(data_matrix, yDim, xDim, []);

end 