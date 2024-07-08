% HemoCorrection.m
% Sarah West
% 8/26/21

% Corrects the hemodynamics of a stack. Requires input of two matrices--
% the blue and violet data from a single stack. 

% Inputs:
% bData -- the blue-channel only stack, in a 3D matrix (pixels, pixels,
    % frames)
% vData -- the violet-channel only stack, in a 3D matrix (pixels, pixels,
    % frames). Must be same size in all dimensions as the blue channel
    % matrix. 
    
% Outputs:
% hData -- The hemodynamic-corrected stack. 

function [hData]=HemoCorrection(bData, vData)
    
    n_dims = ndims(bData);

    if n_dims > 2
        % Find blue matrix sizes
        yDim=size(bData,1);
        xDim=size(bData,2);
        frames=size(bData, 3); 
        
        % Reshape matrices so the "detrend" function below will work
        % properly. Needs to be 2D (frames, pixels*pixels).
        bData=[reshape(bData, yDim*xDim, frames)]; 
        vData=[reshape(vData, yDim*xDim, frames)]; 
    
    end 

    % Flipe matrices so "detrend" function below works properly. Needs to be 
    % (frames, pixels*pixels).
    bData = bData';
    vData = vData';

    % Hemodynamic correction: 
    % Needs to be (frames, pixels*pixels). For each channel, detrend each pixel along the time dimension and 
    % divide by that pixel's mean. Then subtract the violet data from the 
    % blue data. 
    hData=detrend(bData)./mean(bData)-detrend(vData)./mean(vData);
    
    % Flip hData back.
    hData = hData';

    % If data was put in with 3 dimensions, reshape it back.
    if n_dims > 2
        hData=reshape(hData, yDim, xDim, frames); 
    end
    
end 