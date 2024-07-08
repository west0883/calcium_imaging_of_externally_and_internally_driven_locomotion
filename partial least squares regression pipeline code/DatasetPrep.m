
function [parameters] = DatasetPrep(parameters)

    MessageToUser('Prepping ' , parameters);

    % Find the location & value of the comparison iterator in
    % parameters.values
    iterator_location = find(cellfun(@strcmp, parameters.keywords, repmat({'comparison_iterator'}, size(parameters.keywords)))); 
    comparison_iterator = parameters.values{iterator_location};

    % Grab comparison name, indices, variables to use.
    %comparison_name = parameters.this_comparison_set(comparison_iterator).name; % Might not need (is really only for saving files)
    comparison_indices = parameters.this_comparison_set(comparison_iterator).indices;
    comparison_variablesToUse = parameters.this_comparison_set(comparison_iterator).variablesToUse; 

    % Get the explanatory data period indices you're interested in. 
    explanatoryVariables = parameters.explanatory(comparison_indices);

    % Get the response variables period indices & variables you're
    % interested in. 
    % responseVariables is now an array
    holder = cell(1, numel(comparison_variablesToUse));
    for i = 1:numel(comparison_variablesToUse)
        
        holder{i} = parameters.response{comparison_indices, comparison_variablesToUse{i}}; 

    end
    responseVariables = horzcat(holder{:}); 

    % Reshape explanatoryVariables to make each roll into an instance, inside each period cell. 
    % The 2 reshape inputs 
    input1 =  repmat({size(explanatoryVariables{1}, 1)}, size(explanatoryVariables));
    input2 = repmat({[]}, size(explanatoryVariables));
    explanatoryVariables = cellfun(@reshape, explanatoryVariables, input1, input2, 'UniformOutput', false);
 
    % For each response variable (each column of responseVariables),
    % reshape to make each roll into an instance, inside each period cell.
    for variablei = 1:numel(comparison_variablesToUse)
        
        % The 2 reshape inputs
        input1 =  repmat({size(responseVariables{1, variablei}, 1)}, size(explanatoryVariables));
        input2 = repmat({[]}, size(explanatoryVariables));
        
        % Reshape
        responseVariables(:, variablei) =  cellfun(@reshape, responseVariables(:, variablei), input1, input2, 'UniformOutput', false);

    end 

    % **Transpose
    explanatoryVariables = cellfun(@transpose, explanatoryVariables, 'UniformOutput', false);
    responseVariables = cellfun(@transpose, responseVariables, 'UniformOutput', false);

    % **Concatenate explanatory data & response variables vertically across periods. 

    % explanatoryVariables has only one cell array column, so can be done all at once.
    explanatoryVariables_precat = explanatoryVariables;
    explanatoryVariables = vertcat(explanatoryVariables{:});

    % If it's still a cell (can happen if there's a cell with '[]' in it),
    % do again.
    if iscell(explanatoryVariables)

       explanatoryVariables = vertcat(explanatoryVariables{:});
    end

    % Each variable in responseVariables is own cell array column, so have
    % to do each of those first before horizontally concatenating the
    % different variables into  same matrix. Keep the separated variables
    % for later.
    responseVariables_separateVariables = cell(1, numel(comparison_variablesToUse));
    for variablei =  1:numel(comparison_variablesToUse)
        responseVariables_separateVariables{variablei} = vertcat(responseVariables{:, variablei});
    end
   
    % **Remove unnecessary columns of variables.
    for variablei =  1:numel(comparison_variablesToUse)
        
        columns_to_remove = ~any(responseVariables_separateVariables{variablei});

         % If there are columns to remove, remove them.
        if ~isempty(columns_to_remove)
            responseVariables_separateVariables{variablei}(:, columns_to_remove) = [];  
        end
    end

    % Get number of remaining response columns per category, for permuting
    % independently later. 
    variable_category_column_numbers = cellfun('size', responseVariables_separateVariables, 2);

    columns_to_use = cell(1, numel(comparison_variablesToUse)); 
    column_counter = 0;
    for variablei = 1:numel(comparison_variablesToUse)
        columns_to_use{variablei} = column_counter + 1:variable_category_column_numbers(variablei); % Defined this separately for clarity.
        column_counter = column_counter + variable_category_column_numbers(variablei);
    end
    % Put number of columns into dataset.
    dataset.columns_to_use = columns_to_use;

    % Horizontally concatenate the different response variable categories.
    responseVariables = horzcat(responseVariables_separateVariables{:}); 

    % Remove any rows/observations that are all NaNs in the explanatory
    % variables. (Comes from any previous outlier removal). 
    removed_observations = all(isnan(explanatoryVariables), 2);
%     if any(removed_observations)
%         explanatoryVariables(removed_observations, :) = [];
%         responseVariables(removed_observations, :) = [];
%     end

    % *** Deal with explanatory variable outliers, replace with NaNs ***
    % Will be imputed int next step.
    % If user says to (default is not to), remove outliers values
    if isfield(parameters, 'removeOutliers') && parameters.removeOutliers
       
        disp('Removing outliers.');

        explanatoryVariables_old = explanatoryVariables;
       
        % Run outlier removal code. Uses a PCA + Trimmed score regression
        % method.
        %[explanatoryVariables, outliers] = OUTLIERS(explanatoryVariables_old);
        % Find outliers per column.
        if strcmp(parameters.comparison_type, 'categorical')

            change = find(diff(responseVariables(:, 1)) ~= 0);
            outliers1 = isoutlier(explanatoryVariables_old(1:change, :));
            outliers2 = isoutlier(explanatoryVariables_old(change + 1:end, :));

            outliers = [outliers1; outliers2];

        elseif strcmp(parameters.comparison_type, 'continuous')

            outliers_cell = cellfun(@isoutlier, explanatoryVariables_precat, 'UniformOutput', false);
            outliers = vertcat(outliers_cell{:});

        else
            error('No comparison_type defined! Line 134');
            
        end

        % Replace outliers with NaNs.
        explanatoryVariables(outliers) = NaN;

        % Put important info into output structures.
        dataset.outliers.tf = outliers;
        dataset.outliers.indices = find(outliers);
        if ~isempty(outliers)
            [row, column] =  ind2sub([size(explanatoryVariables)], find(outliers));
            dataset.outliers.indices_rowcolumn = [row, column];
        else
            dataset.outliers.indices = [];
        end
        dataset.outliers.values_old = explanatoryVariables_old(dataset.outliers.indices);

        % If any rows (observations) have more than 1/10 the explanatory values
        % as outliers, remove that observation.
%         summation = sum(outliers,2);
%         holder = find(summation > size(explanatoryVariables,2) * 1/10);
%         dataset.outliers.removed_indices = holder;
% 
%         % If there are observations with this, remove. 
%         if ~isempty(holder)
%             explanatoryVariables(holder, :) = [];
%             responseVariables(holder, :) = [];
%             
%         end
    end

   
    % *** Deal with missing data values ***

    % Calculate number of missing values for each variable.
    dataset.NaN_ratios.responseVariables = sum(isnan(responseVariables), 1)/size(responseVariables,1);
    dataset.NaN_ratios.explanatoryVariables = sum(isnan(explanatoryVariables), 1)/size(explanatoryVariables,1);

    % If user says to (default is not to), impute missing values. MUST
    % impute missing if outliers were removed.
    if (isfield(parameters, 'imputeMissing') && parameters.imputeMissing) || (isfield(parameters, 'removeOutliers') && parameters.removeOutliers)

        disp('Imputing missing data.')

        % Put indices of missing data into output structure.
        dataset.missing_data_imputation.indices = find(isnan(responseVariables));
        [row, column] = ind2sub(size(responseVariables), dataset.missing_data_imputation.indices);
        dataset.missing_data_imputation.induces_rowcolumn = [row, column];

        % Run the modified code for trimmed square regression (TSR) for PLS
        [explanatoryVariables, responseVariables, iterations_needed, tolerance_reached, components_needed] = plsmbtsr1_TSRonly(explanatoryVariables, responseVariables, parameters.imputation_components_variance_explained, parameters.imputation_max_components);% parameters.imputation_ncomponents); 

        % Put iterations needed and tolerance reached into dataset
        % structure.
        dataset.missing_data_imputation.iterations_needed = iterations_needed;
        dataset.missing_data_imputation.tolerance_reached = tolerance_reached;
        dataset.missing_data_imputation.components_needed = components_needed;
        dataset.missing_data_imputation.values_new = responseVariables(dataset.missing_data_imputation.indices);

        % If you were calculating outliers, too, 
        if isfield(parameters, 'removeOutliers') && parameters.removeOutliers



            % Get new calculated values.
            %dataset.outlier.values_new = explanatoryVariables(dataset.outliers.indices);

        end 
    end

    % Zscore both variable sets. Keep mu & sigmas for better interprebility
    % of betas later.
    [explanatoryVariables, mu_explanatory, sigma_explanatory] = zscore(explanatoryVariables);
    [responseVariables, mu_response, sigma_response] = zscore(responseVariables);
    dataset.zscoring.explanatoryVariables.mu = mu_explanatory;
    dataset.zscoring.explanatoryVariables.sigma = sigma_explanatory;
    dataset.zscoring.responseVariables.mu = mu_response;
    dataset.zscoring.responseVariables.sigma = sigma_response;

    % For convenience, put both variable sets for this comparison into a stucture for saving. 
    dataset.explanatoryVariables = explanatoryVariables;
    dataset.responseVariables = responseVariables;
    %dataset.responseVariables_separateVariables = responseVariables_separateVariables;
    dataset.variable_category_column_numbers = variable_category_column_numbers; % For telling which columns belong to which category later

    % Put dataset into output structure.
    parameters.dataset = dataset;
    
end