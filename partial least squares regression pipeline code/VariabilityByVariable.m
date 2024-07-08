% VariabilityByVariable.m
% Sarah West
% 10/3/23

% Is called by RunAnalysis.m
function [parameters] = VariabilityByVariable(parameters)

    MessageToUser('Variability of ', parameters);

    % Inputs
    responseVariables = parameters.dataset.responseVariables;
    comparisons_continuous = parameters.comparisons_continuous;
    continuous_variable_names = parameters.continuous_variable_names;
    runSingleVariables = parameters.runSingleVariables;
    omitSingleVariables = parameters.omitSingleVariables;
    shuffleOthers = parameters.shuffleOthers;
    shuffleNumber = parameters.shuffleNumber;

    % Get the current continuous comparison 
    comparison = parameters.values{strcmp(parameters.keywords, 'comparison')};

    % Find this comparison in the list of comparisons
    comparison_index = find(strcmp({comparisons_continuous(:).name}, comparison));

    % Get the variables_to_use for this comparison
    variables_to_use = comparisons_continuous(comparison_index).variablesToUse;
    
    % Cycle through each continuous variable 
    for variablei = 1:numel(continuous_variable_names)
        
        % Get the current continuous response variable
        variable = continuous_variable_names{variablei};

        disp(variable);

        % Find if this variable is used in this comparison
        variable_index = find(strcmp(variables_to_use, [variable '_vector'])); 
    
        % If this variable isn't used in this comparison, continue to next
        % variable
        if isempty(variable_index)

            % Tell user
            disp(['no ' variable ' for this comparison']);
            
            % Go to next variable
            continue 
        end 

        % Use only this variable in response variables

        % If user said to run each variable alone, 
        if runSingleVariables

            % Tell user
            disp('running single variable');

            % If not shuffling others (only using 1 variable)
            if ~shuffleOthers

                % Make a new version of responseVariables, keeping only that variable 
                responseVariables_runSingle = responseVariables(:, variable_index);
    
    %             % Run with only 1 PLSR component 
    %             parameters.ncomponents_max = 1;
    %             parameters.findBestNComponents = false;
    
                parameters.ncomponents_max =  1;  %parameters.ncomponents_max_default;
    
                % Put back into parameters
                parameters.dataset.responseVariables = responseVariables_runSingle;
    
                % Run PLSR
                parameters = PLSR_forRunAnalysis_Inverted(parameters);
    
                % Put results into output structure
                all_variability.(variable).singleVariable = parameters.results; 
            
            % If shuffling others
            else


                % make a vector of 1:each column
                indices = 1:size(responseVariables, 2);
    
                % make all other indices true, this index false
                indices_tf = indices ~= variable_index; 

                all_variability.(variable).singleVariable = NaN(2, shuffleNumber);

                % Create shuffles for each
                for shufflei = 1:shuffleNumber

                    permutation_indices = NaN(size(responseVariables));
                    % For each variable 
                    for var_sub_i = 1:size(responseVariables, 2)

                        permutation_indices(:, var_sub_i) = randperm(size(responseVariables, 1));

                    end 
                    
                    permutation_indices(:, variable_index) = 1:size(responseVariables, 1);

                    responseVariables_runSingle = NaN(size(responseVariables));

                    % Make a new version of responseVariables, keeping only that variable 
                    for var_sub_i = 1:size(responseVariables, 2)
                        responseVariables_runSingle(:, var_sub_i) = responseVariables(permutation_indices(:, var_sub_i), var_sub_i);
                    end 
                    % Put back into parameters
                    parameters.dataset.responseVariables = responseVariables_runSingle;
        
                    % Run PLSR
                    parameters = PLSR_forRunAnalysis_Inverted(parameters);
       
                    all_variability.(variable).singleVariable(:, shufflei) = sum(parameters.results.PCTVAR, 2); 
 
                end 

            end 

        end 

        % Run all other variables in response variables
         parameters.findBestNComponents = true;
        
        % If user said to run all except this variable, 
        if omitSingleVariables

            % Tell user
            disp('omitting variable');

            % Get true-false indices, one for each column of response
            % variables.

            % make a vector of 1:each column
            indices = 1:size(responseVariables, 2);

            % make all other indices true, this index false
            indices_tf = indices ~= variable_index; 

            % Make a new version of responseVariables, removing that variable 
            responseVariables_omitSingle = responseVariables(:, indices_tf);

            % Put back into parameters
            parameters.dataset.responseVariables = responseVariables_omitSingle;

            % Run with the inputted "default" number of components
            parameters.ncomponents_max = parameters.ncomponents_max_default;

            % Run PLSR
            parameters = PLSR_forRunAnalysis_Inverted(parameters);

            % Put results into output structure
            all_variability.(variable).omittedVariable = parameters.results; 

        end 
    end 

    % Put all results into output structure 
    parameters.all_variability = all_variability; 
end 