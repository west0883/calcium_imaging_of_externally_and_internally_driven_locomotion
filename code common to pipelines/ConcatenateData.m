% ConcatenateData.m
% Sarah West
% 4/9/22
% Concatenates data together. 

% Input: parameters.data, 
%       parameters.concatenated_data - previously concatenated data that
%           you're adding this data to. 
%       parameters.concatDim -- the dimension the data should be
%       concatenated across.
% Output: 

% parameters.concatenated_data - the concatenated data.
% parameters.concatenated_origin  - a cell array listing where this data came from. 
  

function [parameters] = ConcatenateData(parameters)
    
    % Display progress message to user.
    MessageToUser('Concatenating ', parameters);

    origin = parameters.values; 

    % If the user gave a concatenation level value field
    if isfield(parameters,'concatenation_level')

        % Get the current iterator value for that level
        iterator_level = find(strcmp(parameters.loop_list.iterators(:,1), parameters.concatenation_level));
        current_iterator = parameters.values{numel(parameters.values)/2 + iterator_level};

        % If the current iterator is 1, that means you're starting a new
        % concatenation, clear any previously concatenated data.
        if current_iterator == 1

            % If there are iterators below this iterator 
            if size(parameters.loop_list.iterators, 1) > iterator_level
    
                all_lowest_iterators = [parameters.values{numel(parameters.values)/2 + iterator_level : numel(parameters.values)}];
                new_concatenation_flag = all(all_lowest_iterators == 1);

            else 
                new_concatenation_flag = true; 
            end
    
            % Clear any previously concatenated data
            if new_concatenation_flag && isfield(parameters, 'concatenated_data')
                parameters = rmfield(parameters, 'concatenated_data'); 
            end 
        end 
    end 
   
    % If parameters.data has only one entry (no cells/not a cell array), or if user said to
    % concatenated across (instead of within) cells. 
    if ~iscell(parameters.data)
        % Make an empty variable for celli to determine display
        % message (just puts an empty place in the message)
        celli = []; 

        % If concatenated_data doesn't exist yet, create it as empty
        % array.
        if ~isfield(parameters, 'concatenated_data')
            parameters.concatenated_data = [];
        end

        % If a cell array with info about the origin of the entries doesn't
        % exist yet, create it as empty cell array.
        if ~isfield(parameters, 'concatenated_origin')
            parameters.concatenated_origin = {};
        end

        % If data to concatenate isn't empty, add it to the
        % concatenated data 
        if ~isempty(parameters.data)
            [parameters.concatenated_data, parameters.concatenated_origin] = SubConcatenateData(parameters.concatenated_data, parameters.data, parameters.concatDim, celli, parameters.concatenated_origin, origin); 
        end

    % If parameters.data has more than one entry, in cell form....
  
    % ... and you want to concatenate across cells
    elseif (isfield(parameters, 'concatenate_across_cells') && parameters.concatenate_across_cells)
        
        % If concatenated_data doesn't exist yet, create it as empty
        if ~isfield(parameters, 'concatenated_data')
            parameters.concatenated_data = [];
        end

        % If a cell array with info about the origin of the entries doesn't
        % exist yet, create it as empty cell array.
        if ~isfield(parameters, 'concatenated_origin')
            parameters.concatenated_origin = {};
        end

        [parameters.concatenated_data, parameters.concatenated_origin] = SubConcatenateData(parameters.concatenated_data, parameters.data, parameters.concatDim, [], parameters.concatenated_origin, origin); 
  
    % ... and you want to concatenate within cells across other levels 
    else    
        % If concatenated_data doesn't exist yet, create it as empty
        % cell array with same size as number of cells in data.
        if ~isfield(parameters, 'concatenated_data')
            parameters.concatenated_data = cell(size(parameters.data));
        end

        % If a cell array with info about the origin of the entries doesn't
        % exist yet, create it as empty cell array.
        if ~isfield(parameters, 'concatenated_origin')
            parameters.concatenated_origin = cell(size(parameters.data));
        end

        % Could do cellfun, but I want to disp where the errors occur.
        for celli = 1:numel(parameters.data)
            if ~isempty(parameters.data{celli})

                [parameters.concatenated_data{celli}, parameters.concatenated_origin{celli}] = SubConcatenateData(parameters.concatenated_data{celli}, parameters.data{celli}, parameters.concatDim, celli, parameters.concatenated_origin{celli}, origin); 

            end
        end
    end 
end

function [concatenated_data, concatenated_origin] = SubConcatenateData(concatenated_data, data, concatDim, celli, concatenated_origin, origin) 
    
    % Update origin with number of instances included in data.
    origin = [size(data, concatDim); origin];
    % Initialize empty output
    try
        concatenated_data = cat(concatDim, concatenated_data, data);
        concatenated_origin = cat(concatDim, concatenated_origin, {origin});
    catch 
        if isempty(celli)
            disp(['Dimension error.']);
        else
           disp(['Dimension error in cell ' num2str(celli) '.']);
        end
    end
end 
