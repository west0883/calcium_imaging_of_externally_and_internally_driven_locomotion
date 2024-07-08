% FisherTransform.m
% Sarah West
% 5/9/22

% Called by RunAnalysis, Applies the FisherTransform (arctanh(x)) to 
% Pearson correlation coefficients.

function parameters = FisherTransform(parameters)
    
    % Tell user what's happening.
    MessageToUser('Transforming', parameters);
   
    % Apply transformation. Use cellfun if it's a cell. 
    if iscell(parameters)
        parameters.data_transformed = cellfun(@atanh, parameters.data, 'UniformOutput', false);
    else
         parameters.data_transformed = atanh(parameters.data);
    end
end 