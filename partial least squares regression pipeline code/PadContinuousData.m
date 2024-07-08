% PadContinuousData.m
% Sarah West
% 8/4/22

% Pads continuous data to have speed, accel, duration, pupil diameter in
% first - 4th entries of data, respectively. 

function [parameters] = PadContinuousData(parameters)

    MessageToUser('Padding ', parameters);

    % Make default empty data padded output = data.
    parameters.data_padded = parameters.data;

    % Make default comparison type.
    comparison_type = 'default';

    % Find a comparison type field.
    if isfield(parameters, 'comparison_type') 
        comparison_type = parameters.comparison_type;
    end

    % Find a comparison type iterator (supercedes the field above).
    if any(strcmp(parameters.keywords, 'comparison_type'))
        comparison_type = parameters.values{strcmp(parameters.keywords, 'comparison_type')};
    end 

    % If comparison_type is continuous
    if strcmp(comparison_type, 'continuous')

        % Find the comparison name
        comparison = parameters.values{strcmp(parameters.keywords, 'comparison')};

        % Find the comparison in the continuous comparison structure, get
        % variables to use.
        variablesToUse= parameters.comparisons_continuous(strcmp({parameters.comparisons_continuous(:).name}, comparison)).variablesToUse;

        % If not useing all variables, 
        if numel(variablesToUse) < numel(parameters.continuous_variable_names)

            % Make new matrix of data.
            sizes = size(parameters.data);
            sizes(parameters.padDim) =  numel(parameters.continuous_variable_names);
            data_padded = NaN(sizes);

            % For each of the total continuous variables being used,
            for variablei = 1:numel(parameters.continuous_variable_names)

                variable = [parameters.continuous_variable_names{variablei} '_vector'];

                Cin= repmat({':'}, 1, ndims(parameters.data));
                Cout = repmat({':'}, 1, ndims(parameters.data));

                % Find if this variable is in "variables to use."
                index = find(strcmp(variablesToUse, variable));

                if ~isempty(index)
                    % put into dimensions
                    Cin{parameters.padDim} = index;
                    Cout{parameters.padDim} = variablei;
    
                    % Put into padded form
                    data_padded(Cout{:}) = parameters.data(Cin{:});

                end
            end

            % Put into output.
            parameters.data_padded = data_padded;

        end
    end
end 