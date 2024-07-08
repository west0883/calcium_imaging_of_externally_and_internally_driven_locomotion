% EvaluateOnData.m
% Sarah West
% 5/10/22

% A function that allows you to apply any expression to data using
% RunAnalysis.

function [parameters] = EvaluateOnData(parameters)
    
    % Display progress message to user.
    MessageToUser('Evaluating on ', parameters); 

    % Make evaluation string
    evaluation_string = CreateStrings(parameters.evaluation_instructions{parameters.functioni}, parameters.keywords, parameters.values);
    
    % Calculate
    eval(evaluation_string);

    % Put into output.
    parameters.data_evaluated = data_evaluated;

end 