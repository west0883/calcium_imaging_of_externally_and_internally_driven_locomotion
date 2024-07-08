% ReshapeContinuousData.m
% Sarah West
% 1/23/23

% Reshapes continuous data from PLSR so each contnuous data type is in its
% own entry in a new, last dimension.

% Needed inputs:
% parameters.data 
% parameters.variablesDimIn 
% parameters.this_comparison_set

% Outputs:
% parameters.data_reshaped

function [parameters] = ReshapeContinuousData(parameters)

    MessageToUser('Reshaping ', parameters);

    variablesDimIn = parameters.variablesDimIn;
    data = parameters.data;
    continuous_variable_names = parameters.continuous_variable_names;

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
        comparisons_continuous = parameters.this_comparison_set; 
        variablesToUse = comparisons_continuous(strcmp({comparisons_continuous(:).name}, comparison)).variablesToUse;

        % Get default sizes for data_reshaped, adding a dimension to the
        % end for each variable name
        sizes = [size(parameters.data) numel(parameters.continuous_variable_names)];
        
        % Make new matrix of data.
        new_size = sizes(variablesDimIn) ./ numel(variablesToUse);
        sizes(variablesDimIn) = new_size;
        data_reshaped = NaN(sizes);

        % Create matrices for variying dimension indexing. 
        Cin= repmat({':'}, 1, numel(sizes) - 1);
        Cout = repmat({':'}, 1, numel(sizes));

        % For each continuous variable,
        for variablei = 1:numel(parameters.continuous_variable_names)
            
            variable = [continuous_variable_names{variablei} '_vector'];

            % Find if this variable is in "variables to use."
            index = find(strcmp(variablesToUse, variable));
            
            % If this variable is in variablesToUse
            if ~isempty(index)
                
                % Put into new shaped data
                Cin{variablesDimIn} = (index - 1) * new_size + 1 :  index * new_size; 
                Cout{end} = variablei; 

                data_reshaped(Cout{:}) = data(Cin{:}); 

            end 
        end 

        parameters.data_reshaped = data_reshaped;
        
end