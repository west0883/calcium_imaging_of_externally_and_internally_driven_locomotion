% FluorescenceContinuousAverageHemispheres.m
% Sarah West
% 8/23/23

function [parameters] = FluorescenceContinuousAverageHemispheres_Inverted(parameters)

    MessageToUser('Hemisphere averageing ', parameters);
    data = parameters.data; 
    comparisons_continuous = parameters.comparisons_continuous;
    variables_to_average_across_hemispheres = parameters.variables_to_average_across_hemispheres;
    nodeDim = parameters.nodeDim;
    variableDim = parameters.variableDim; 

    % Find the name of this comparison
    this_comparison = parameters.values{strcmp(parameters.keywords, 'comparison')};

    % Find the index in the set of comparisons & the variables to use
    index = strcmp({comparisons_continuous(:).name}, this_comparison);
    variables_to_use = comparisons_continuous(index).variablesToUse;

    indices = [];

    % Find which variables in this comparison should be averaged
    for i = 1:numel(variables_to_average_across_hemispheres)

        variable = variables_to_average_across_hemispheres{i};
        indices_holder = find(strcmp(variable, variables_to_use));
    
        indices = [indices indices_holder];
    end 
    
    % Get dimensions
    C = repmat({':'}, 1, ndims(data));
    
    % If output type is being used, separate by Covs or betas

    % If Betas,
    if any(strcmp(parameters.keywords, 'output_type')) && strcmp(parameters.values{strcmp(parameters.keywords, 'output_type')}, 'BETA') || strcmp(parameters.output_type, 'BETA')
        
        % For the columns marked as "indices", average across hemispheres.

        % Skip the first variable dimension
        C1 = C; 
        C1{variableDim} = 2:size(data, variableDim);
        data_no_intercepts= data(C1{:});
        dims1 = 1:2:size(data_no_intercepts, nodeDim);
        dims2 = 2:2:size(data_no_intercepts, nodeDim);

        % Put intercepts back in (are removed later)
        C2 = C;
        C2{variableDim} = 1;
        data_intercepts = data(C2{:});
       
        output = data_no_intercepts;

    % Otherwise, assume Covs
    else

        % For the columns marked as "indices", average across hemisphere
        dims1 = 1:2:size(data, nodeDim);
        dims2 = 2:2:size(data, nodeDim);
        
        output = data;
    end 

    C{variableDim} = indices;
        
    C{nodeDim} = dims1;
    a = data_no_intercepts(C{:});
    
    C{nodeDim} = dims2;
    b = data_no_intercepts(C{:});

    c = mean(cat(ndims(data_no_intercepts) + 1, a, b), ndims(data_no_intercepts) + 1, 'omitnan'); % concatenate and average across one extra dimension
    
    % Put into output
    C{nodeDim} = dims1;
    output(C{:}) = c;

    C{nodeDim} = dims2;
    output(C{:}) = c;

    % Put intercepts back in (will be removed later)
    output = cat(variableDim, data_intercepts, output);

    % Put into output structure
    parameters.output = output;

end 