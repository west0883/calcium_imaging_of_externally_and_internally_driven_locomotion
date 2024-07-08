% ShortenFluorescenceWarningPeriods.m
% Sarah West
% 7/31/23

function [parameters] = ShortenFluorescenceWarningPeriods(parameters)

    MessageToUser('Shortening ', parameters);

    % Pull out input values
    indices_to_shorten = parameters.indices_to_shorten;
    data = parameters.data;
    shorten_dimensions = parameters.shorten_dimensions;

    % Get period iterator
    period_iterator = parameters.values{strcmp(parameters.keywords, 'period_iterator')};

    % If the period iterator matches an index to shorten, shorten it
    if any(period_iterator == indices_to_shorten)

        % Shorten with instructions 
        eval(['data_shortened = data(' shorten_dimensions ');']);
    else 
        data_shortened = data;

    end 


    parameters.data_shortened = data_shortened;

end 