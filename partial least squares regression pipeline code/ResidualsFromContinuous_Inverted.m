% SubtractContinuousVariables.m
% Sarah West
% 6/13/22

% Removes Beta * continuous Y response variable from X explanatory variables
% after continuous variable regressions.

function [parameters] = ResidualsFromContinuous_Inverted(parameters)

    MessageToUser('Residuals from ', parameters);

    dataset_info = parameters.dataset;
    PLSR_results = parameters.PLSR_results;

    % Find the location & value of the type iterator in
    % parameters.values
    iterator_location = logical(cellfun(@strcmp, parameters.keywords, repmat({'comparison_iterator'}, size(parameters.keywords)))); 
    comparison_iterator = parameters.values{iterator_location};

    % Get the comparison & variables used for regression from comparison 
    indicesToUse = parameters.comparisons_continuous(comparison_iterator).indices;

    % *** INVERT HERE ***
    % Grab the X residuals from the continuous variables regression.
    %Xnew = parameters.PLSR_results.stats.Xresiduals;
    Xnew = parameters.PLSR_results.stats.Yresiduals;

    % Get out response variables (for outlier replacement imputation).
    responseVariables = parameters.dataset.responseVariables;

    % Make a version of Xnew that won't have any outliers removed (for
    % counting where to find instances that were entirely removed as
    % outliers). Same with responseVariables.
    Xnew_old = Xnew; 


    % If going to behavior intercepts, 
    if parameters.useIntercepts
       
        % Grab the behavior intercepts
        none  = parameters.intercepts; 

        % Find which intercepts to keep/use
        % get the current comparison 
        comparison = parameters.values{strcmp(parameters.keywords, 'comparison')};

        % find this comparison in continous comparisons
        comparison_iterator = find(strcmp({parameters.comparisons_continuous(:).name}, comparison));

        variables_to_use = strrep(parameters.comparisons_continuous(comparison_iterator).variablesToUse, '_vector', '');
        variable_indices = NaN(1, numel(parameters.continuous_variable_names));
        for i = 1:numel(parameters.continuous_variable_names)
            variable_indices(i) = any(strcmp(variables_to_use, parameters.continuous_variable_names{i}));
        end

        none = none(logical(variable_indices));

        % convert behavior intercept to zscore space (is inverted, so use
        % the responseVariables field from dataset_info)
        none_zscored = (none - dataset_info.zscoring.responseVariables.mu)./ dataset_info.zscoring.responseVariables.sigma;
        
        % find expected brain values
        brain_expected_zscored = none_zscored * PLSR_results.BETA(2:end, :);

        % remove expected brain values from zscore space (is inverted, so use
        % the explanatoryVariables field from dataset_info)
        brain_expected = brain_expected_zscored .* dataset_info.zscoring.explanatoryVariables.sigma + dataset_info.zscoring.explanatoryVariables.mu;
        
        % Pull out original behavior variabels, un z-scored.
        brain_original = dataset_info.explanatoryVariables .* dataset_info.zscoring.explanatoryVariables.sigma + dataset_info.zscoring.explanatoryVariables.mu;
        
        % Un z-score the brain residuals
        brain_residuals =  Xnew .* dataset_info.zscoring.explanatoryVariables.sigma + dataset_info.zscoring.explanatoryVariables.mu;
        
        % Add the residuals to the expected brain values at 0.
        Xnew =  brain_residuals + brain_expected;

        % put these variables into dataset_out
        dataset_out.intercepts.brain_original = brain_original;
        dataset_out.intercepts.brain_expected = brain_expected;
        dataset_out.intercepts.brain_residuals= brain_residuals;
        dataset_out.intercepts.brain_at_intercepts = Xnew;

    else 

        % Un-normalize (because will need to be normalized with other behavior
        % types in across-category comparisons later).
        Xnew = Xnew .* parameters.dataset.zscoring.explanatoryVariables.sigma + parameters.dataset.zscoring.explanatoryVariables.mu;
    end 

    % Remove any outliers, if user says so.
    if isfield(parameters, 'removeOutliers') && parameters.removeOutliers
        outliers = isoutlier(Xnew, 1);

        % Replace outliers with NaN
        Xnew(outliers) = NaN;

        % Put important info into output structures.
        dataset_out.outliers.indices = find(outliers);
        if ~isempty(outliers)
            [row, column] =  ind2sub([size(Xnew)], find(outliers));
            dataset_out.outliers.indices_rowcolumn = [row, column];
        else
            dataset_out.outliers.indices = [];
        end
        dataset_out.outliers.values_old = Xnew_old(dataset_out.outliers.indices);

        % If any rows (observations) have more than 1/10 the explanatory values
        % as outliers, remove that observation.
        summation = sum(outliers,2);
        holder = find(summation > size(Xnew,2) * 1/10);
        dataset_out.outliers.removed_indices = holder;

        % If there are observations with this, remove. 
        if ~isempty(holder)
            Xnew(holder, :) = [];
            responseVariables(holder, :) = [];
            Xnew_old(holder, :) = NaN;
        end
    else 
        dataset_out.outliers = [];
    end

    % If user says to (default is not to), impute missing values. MUST
    % impute missing if outliers were removed.
    if (isfield(parameters, 'imputeMissing') && parameters.imputeMissing) || (isfield(parameters, 'removeOutliers') && parameters.removeOutliers)

        disp('Imputing missing data.')

        % Run the modified code for trimmed square regression (TSR) for PLS
        [Xnew, ~, iterations_needed, tolerance_reached, components_needed] = plsmbtsr1_TSRonly(Xnew, responseVariables, parameters.imputation_components_variance_explained, parameters.imputation_max_components);% parameters.imputation_ncomponents); 

        % Put iterations needed and tolerance reached into dataset_out
        % structure.
        dataset_out.missing_data_imputation.iterations_needed = iterations_needed;
        dataset_out.missing_data_imputation.tolerance_reached = tolerance_reached;
        dataset_out.missing_data_imputation.components_needed = components_needed;
    end

    % Put the removed instances back in. 
    if (isfield(parameters, 'removeOutliers') && parameters.removeOutliers) && (isfield(parameters, 'imputeMissing') && parameters.imputeMissing)

        Xholder = Xnew_old;

        % find indices of Xnew_old that aren't NaNs, put in Xnew.
        Xholder(~isnan(Xnew_old)) = Xnew;
        
        % Make Xnew into Xholder. 
        Xnew = Xholder;
    end 

    % Put in instances removed during continuous variable regression prep, too. 
    if isfield(parameters.dataset, 'outliers')
%         if ~isempty(parameters.dataset.outliers.removed_indices)
% 
%             removed_indices = parameters.dataset.outliers.removed_indices; % Pull out for clarity
%             Xholder = NaN(size(Xnew,1) + numel(removed_indices), size(Xnew, 2));
% 
%             % Handle first entry
%             Xholder(1:removed_indices(1) - 1, :) = Xnew(1:removed_indices(1) - 1, :);
%             counter = 1 + removed_indices(1);
%             
%             % If more than one removed instance, need to skip lines
%             if numel(removed_indices) > 1
%                 
%                 for instancei = 2:numel(removed_indices)
%     
%                     % In Xnew, need to remove each instance already accounted
%                     % for when looking at the range.
%                     Xholder(counter:removed_indices(instancei) - 1, :) = Xnew((counter - instancei + 1): ...
%                         (removed_indices(instancei) - instancei), :);
%                     
%                     % Add removed instances to counter.
%                     counter = removed_indices(instancei) + 1;
%                 end
%             end
%             
%             % Now for everything after the last removed index.
%             Xholder(counter:end, :) = Xnew((counter - numel(removed_indices)):end, :);
% 
%             % Rename Xholder to Xnew.
%             Xnew = Xholder;

%        end
    end

    % Start putting Xnew back into the big correlations matrix

    % Transpose 
    Xnew = Xnew'; 

    % If there isn't a field for values_new yet, make it & fill it with the
    % old explanatory values matrix.
    if ~isfield(parameters, 'values_new')
        parameters.values_new = parameters.values_old;
    end

    % Set up a counter to hold your place in new values. Can do this
    % because you concatenated only these indices the first time in
    % DatasetPrep.m
    counter = 1;
   
    % Find each cell entry of old values that corresponds to this behavior
    % type.
    for indexi = 1:numel(indicesToUse)
        index = indicesToUse(indexi);

        % Get the roll number & instances from this 
        rollnumber = size(parameters.values_old{index}, 2);
        instancesnumber = size(parameters.values_old{index}, 3);

        % Get a vector of the instances of Xnew that correspond to this
        % period/index of values_old.
        subset_instances_vector = counter : (counter + rollnumber * instancesnumber - 1);

        % Pull out the matrices that belong to this cell/period
        subset = Xnew(:, subset_instances_vector);

        % Reshape to match this cell/period, put into place.
        parameters.values_new{index} = reshape(subset, size(subset,1), rollnumber, instancesnumber);

        % Update counter by number of matrices included in that cell
        counter = counter + rollnumber * instancesnumber;
    end

    % Counter should be 1 more than the size of Xnew by the end
    if counter ~= size(Xnew, 2) + 1
        error('counter should equal 1 more than size of the new explanatory value');
    end

    % Put dataset_out info into output structure
    parameters.dataset_out = dataset_out;
end 
