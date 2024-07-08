% ApplyMasks.m
% Sarah West
% 5/25/22

% Applies brain masks to data using RunAnalysis.

function [parameters] = ApplyMasks(parameters)

    % Message 
    MessageToUser('Masking ', parameters);

    % If there are masks to apply 
    if isfield(parameters, 'indices_of_mask')
    
        % Apply mask to data (for right now assumes different
        % sources in 3rd dimension)
        holder = reshape(parameters.data, parameters.yDim * parameters.xDim, []);
        parameters.data_masked = holder(parameters.indices_of_mask, :);
    
    else
        % 
        disp('No indices of mask given.');

        % Assign masked sources just as the inputted sources
        parameters.data_masked = parameters.data; 
    end

end