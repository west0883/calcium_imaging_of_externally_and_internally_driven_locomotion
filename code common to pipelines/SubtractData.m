% SubtractData.m
% Sarah West
% 4/22/22

function [parameters] = SubtractData(parameters)

    % Display progress message to user.
    MessageToUser('Subtracting ', parameters);

    parameters.data_subtracted = parameters.subtract_from_this - parameters.subtract_this;

end 