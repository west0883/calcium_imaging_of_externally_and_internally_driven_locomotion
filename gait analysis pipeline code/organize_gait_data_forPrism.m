% organize_gait_data_forPrism.m
% Sarah West
% 8/29/23

function [parameters] = organize_gait_data_forPrism(parameters)

    % within each group:
    % # mice x mean, SD, N
    % 5 groups: spon, 1600, 2000, 2400, 2800 concatenated horizontally
    
    MessageToUser('Organizing for Prism ', parameters);

    % Inputs
    data = parameters.data;
    number_of_mice = parameters.number_of_mice;
    number_of_periods = parameters.number_of_periods;
    instancesDim = parameters.instancesDim;

    % if no data_organized_out field, make a new data_organized
    if ~isfield(parameters, 'data_organized')
        
        data_organized = repmat({NaN(number_of_mice,3)}, 1, number_of_periods); 

    % Otherwise, pull out existing field
    else
        data_organized = parameters.data_organized;
    end 

    % Get mouse iterator
    mouse_iterator = parameters.values{strcmp(parameters.keywords, 'mouse_iterator')};

    % Get period iterator 
    period_iterator = parameters.values{strcmp(parameters.keywords, 'period_iterator')};

    % Get N, mean, and std
    N = size(data, instancesDim);
    average = mean(data, instancesDim, 'omitnan');
    std_dev = std(data, [], instancesDim, 'omitnan');

    % Put into appropriate place
    data_organized{period_iterator}(mouse_iterator, 1) = average;
    data_organized{period_iterator}(mouse_iterator, 2) = std_dev;
    if isfield(parameters, 'NisMice') && parameters.NisMice
        data_organized{period_iterator}(mouse_iterator, 3) = 1;
    else
        data_organized{period_iterator}(mouse_iterator, 3) = N;
    end

    % Concatenate horizontally
    data_organized_out = horzcat(data_organized{:});

    % Put into output structure
    parameters.data_organized = data_organized;
    parameters.data_organized_out = data_organized_out;

end