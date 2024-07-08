% MessageToUser.m
% Sarah West
% 4/26/22

% Is called within functions that are themselves called by RunAnalysis.m
% Displays a message to user in the command window reporting the progress
% of the function by saying what iterators it's on.

% Inputs: 
% parameters.values -- created by RunAnalysis.m
% message -- a string array with the base message to be reported to user,
% should reflect the action of the calling function.

function [] = MessageToUser(message, parameters)
   
    % If there's a flag saying this function was called through RunAnalysis
    if  isfield(parameters, 'RunAnalysis_flag') && parameters.RunAnalysis_flag

        % If there's a "values" field from RunAnalysis, print updating message
        % for user. 
        if isfield(parameters, 'values')
            holder = strjoin(parameters.values(1:numel(parameters.values)/2), ', ');
            disp([message holder]); 
        end
    end
end