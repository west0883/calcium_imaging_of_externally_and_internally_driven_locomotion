% PermuteData
% Sarah West
% 5/11/22

% Permutes data for use with RunAnalysis.m. 

function [parameters] = PermuteData(parameters)
    
    % Display progress message to user.
    MessageToUser('Permuting ', parameters); 
    
    % Permute
    parameters.data_permuted = permute(parameters.data, parameters.DimOrder); 
end 