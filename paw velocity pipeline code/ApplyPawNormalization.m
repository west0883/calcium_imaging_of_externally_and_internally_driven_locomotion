% ApplyPawNormalization.m
% Sarah West
% 8/30/23

% A place to edit paw velocity normalizations with more control (was
% previously doing everything with EvaluateOnData.m)

function [parameters] =ApplyPawNormalization(parameters)

    MessageToUser('Applying normalization on ', parameters);

    % Inputs
    data = parameters.data; % The velocity timeseries for the stack
    average = parameters.average; % The average velocity for that recording day
    std_dev = parameters.std_dev; % The std of velocity for that recording day
    std_threshold = parameters.std_threshold;
    threshold_body_parts = parameters.threshold_body_parts;

    % find the body part
    body_part = parameters.values{strcmp(parameters.keywords, 'body_part')};

    % if the std_dev is below the threshold
    if std_dev < std_threshold

        % & this is a body part that the threshold applies to 
        if any(strcmp(body_part, threshold_body_parts))

            % change the std_dev of this day to be equal to the minimum
            % threshold
            std_dev = std_threshold;
        end
    end 

    % Divide data by the standard deviation
    data_normalized = data./std_dev;

    % Put normalized data into output structure.
    parameters.data_normalized = data_normalized;


end 