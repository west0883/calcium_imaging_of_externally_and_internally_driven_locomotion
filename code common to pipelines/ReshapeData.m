% ReshapeData.m
% Sarah West
% 4/13/22

% Reshapes data for use with RunAnalysis.m. Uses a string cell array as
% input for dimensions.

function [parameters] = ReshapeData(parameters)
    
    % Display progress message to user.
    MessageToUser('Reshaping ', parameters); 
    
    % Put in reshaping directions to make this as flexible as possible.
    toReshape_string = CreateStrings(parameters.toReshape, parameters.keywords, parameters.values);
    eval(['toReshape = ' toReshape_string ';']);
    
    dimensions_string = CreateStrings(parameters.reshapeDims, parameters.keywords, parameters.values);
    eval(['dimensions = ' dimensions_string ';']);
    
    parameters.data_reshaped = reshape(toReshape, dimensions{:});

end 