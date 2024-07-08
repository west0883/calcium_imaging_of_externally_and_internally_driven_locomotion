% PLSR_forRunAnalysis.m
% Sarah West
% 6/6/22

% A wrapper that runs the Matlab function plsregress_fullcode via RunAnalysis.m

function [parameters] = PLSR_forRunAnalysis_Inverted(parameters)

    % *** Inputs *** 
    % parameters.dataset.explanatoryVariables --
    % parameters.dataset.responseVariables -- 
    % parameters.ncomponents_max = 10; 

    % parameters.comparison_type -- string; either 'continuous' or
    %   'categorical'; denotes if response variables are continuous or
    %   categorical

        % ** Validation parameters **

        % parameters.findBestComponents -- true/false; whether or not to
            % find optimal number of components to use using cross-validation
        % parameters.kFold -- number of 'folds' to do for k-fold cross
            % validataion
        % parameters.MonteCarloReps
        % parameters.stratify 
        % parameters.contiguous_partitions

        % ** null distribution generation parameters ** 

        % parameters.permutationGeneration
        % parameters.n_permutations
        % parameters.onPermutations
        % parameters.useBootstrapping

    % *** Outputs ***


    MessageToUser('PLS regressing ', parameters);

    % Tell user all the parameters they sent in (so they can cancel before 
    % it runs for a long time if they made a mistake).
    if isfield(parameters, 'findBestNComponents') && parameters.findBestNComponents

        if isfield(parameters, 'contiguous_partitions') && ~parameters.contiguous_partitions 
        else

            message = ['Finding best number of components with maximum ' num2str(parameters.ncomponents_max) ...
                ' components, ' num2str(parameters.kFolds) ' -folds.'];
            disp(message);
        end
       
    end

    if isfield(parameters, 'permutationGeneration') && parameters.permutationGeneration

        disp(['Will run ' num2str(parameters.n_permutations) ' permutations.'])

    end
    
    if isfield(parameters, 'useBootstrapping') && parameters.useBootstrapping

        disp(['Will run ' num2str(parameters.n_bootstraps) ' Bootstrap samples.'])

    end
    comparison_type = parameters.comparison_type;
    
    % Run plsregress_fullcode to find the optimal number of components, using a maximal number of components (somewhat
    % arbitrary)
    % ***** INVERT HERE***
    explanatoryVariables = parameters.dataset.responseVariables;
    responseVariables = parameters.dataset.explanatoryVariables;

    % If there are fewer observations than asked-for components (happens in
    % 2nd level analyses), lower ncomponents to 1 less than number of rows.
    if size(responseVariables, 1) <= parameters.ncomponents_max
        ncomponents_max = size(responseVariables, 1) - 1;

    else
        ncomponents_max = parameters.ncomponents_max;
    end

    % If user says so
    if isfield(parameters, 'findBestNComponents') && parameters.findBestNComponents

        if isfield(parameters, 'contiguous_partitions') && parameters.contiguous_partitions

            % Adjust k-folds size to response variable observation number.
            if size(responseVariables, 1) <= parameters.kFolds
                kFolds = size(responseVariables, 1) - 1;
        
            else
                kFolds = parameters.kFolds;
            end
       
            stratify = parameters.stratify;

            % Run cross-validation with contiguous partitions.
    
            % Divide data into parameters.kFolds number of partitions.
           
            % If the data is categorical, you have to stratify it within each
            % category (partition within each category).--> only when
            % there's just one "class" of category being compared.
            if strcmp(comparison_type, 'categorical') && stratify 
    
                % Make a holder cell array of indices to use for each
                % partition. Rows are different partitions, columns
                % different variables, dimension 3 is different monte carlo
                % divisions.
                partition_indices_holder = cell(parameters.kFolds, parameters.MonteCarloReps, size(responseVariables, 2));
    

                % Have to switch to explanatory for inversion
                % For each explanatory variable 
                for variablei = 1:size(explanatoryVariables, 2)
    
                    % Get out variable indices

                    % These are normalized category labels. Assume the
                    % relevant positive category is the max in this column.
                    % (1 vs 0 or even -1).
                    
                    variable_indices = find(explanatoryVariables(:,variablei) > 0); 
                    
                    % Find number of observations that will go into each fold.
                    % Remainder will go into the last fold.
                    nobservations = floor(numel(variable_indices)/kFolds);
                    remainder = rem(numel(variable_indices), kFolds);

                    % Make a list of offsets to generate different
                    % partitions. Ranging from 1:nobservations. 
                    if nobservations <= parameters.MonteCarloReps
                       disp(['Number of observations allows for only ' num2str(nobservations) ' MonteCarlo repitions.']);
                       offset_vector = 1:nobservations;
                    else
                       offset_vector = randperm(nobservations, parameters.MonteCarloReps);
                    end

                    % For each monteCarlo repetition, 
                    for repititioni = 1:numel(offset_vector)

                        % Get the offset for this repition. 
                        offset = offset_vector(repititioni);

                        % For each fold/partition.
                        parfor foldi = 1:kFolds
        
                            % Make a vector of indices to use fot this partition.
                            vector_indices = offset + [1:nobservations] + (foldi - 1) * nobservations;

                            % If vector_indices go past the number of
                            % total observations for this variable, adjust them. 
                            if any(vector_indices > numel(variable_indices))
                               
                                overs = find(vector_indices > numel(variable_indices));
                                vector_indices(overs) = vector_indices(overs) - numel(variable_indices);
                            end 

                            partition_indices_holder{foldi, repititioni, variablei} = variable_indices(vector_indices);
                        end
        
                        % Adjust remainders.
                        remainders = offset + variable_indices((end - remainder + 1) : end);

                        % If vector_indices go past the number of total observations for this variable, adjust them. 
                        if any(remainders > numel(variable_indices))
                            overs = find(remainders > numel(variable_indices));
                            remainders(overs) = remainders(overs) - numel(variable_indices);
                        end 

                        % Put the remainder observations in the last partition. 
                        partition_indices_holder{end, repititioni, variablei} =  [partition_indices_holder{end, repititioni, variablei}; remainders];
               
                    end
                end
    
                % Concatenate indices across different variables (so you have
                % different variables in different partitions).
                partition_indices = cellfun(@vertcat, partition_indices_holder(:, :, 1), ...
                     partition_indices_holder(:,:,2), 'UniformOutput', false);
    
            % If continuous, don't need to stratify.
            else
                % Make a holder cell array of indices to use for each partition. Rows
                % are different partitions. Dimension 2 is different monte carlo
                % divisions.
                partition_indices = cell(kFolds, parameters.MonteCarloReps);
    
                % Get out variable indices (for continuous, are just 1: total number
                % of observations)
                variable_indices = 1:size(responseVariables, 1);
                
                % Find number of observations that will go into each fold.
                % Remainder will go into the last fold.
                nobservations = floor(numel(variable_indices)/kFolds);
                remainder = rem(numel(variable_indices), kFolds);
                
                % Make a list of offsets to generate different
                % partitions. Ranging from 1:nobservations. 
                if nobservations <= parameters.MonteCarloReps
                   disp(['Number of observations allows for only ' num2str(nobservations) ' MonteCarlo repitions.']);
                   offset_vector = 1:nobservations;
                else
                   offset_vector = randperm(nobservations, parameters.MonteCarloReps);
                end
                

                % For each monteCarlo repetition, 
                for repititioni = 1:numel(offset_vector)

                    % Get the offset for this repition. 
                    offset = offset_vector(repititioni);
                    
                    % For each fold/partition.
                    parfor foldi = 1:kFolds
                    
                        % Put the number of observations into the partition
                        % indices holder. Make a vector of indices for easier
                        % reading.
                        vector_indices = offset + [1:nobservations] + (foldi - 1) * nobservations; % Do separately here for readability.
                        
                        % If vector_indices go past the number of
                        % total observations for this variable, adjust them. 
                        if any(vector_indices > numel(variable_indices))
                           
                            overs = find(vector_indices > numel(variable_indices));
                            vector_indices(overs) = vector_indices(overs) - numel(variable_indices);
                        end 
                        
                        partition_indices{foldi, repititioni} = variable_indices(vector_indices);
                    end
                
                    % Put the remainder observations in the last partition.
                    remainders = offset + variable_indices((end - remainder + 1) : end);
                    if any(remainders > numel(variable_indices))
                        overs = find(remainders > numel(variable_indices));
                        remainders(overs) = remainders(overs) - numel(variable_indices);
                    end 

                    partition_indices{end, repititioni} =  [partition_indices{end, repititioni} remainders];
                end
            end

            % Run calculations of sum squared error of PLSR model with each fold.
    
            % Make holder for sum squared errors. Rows are each fold, dimension 2 are 
            % explanatory vs response, dimension 3 are each number of components 
            % (plus the 0 component null condition).
           
            SSEs_byrepitition = NaN(parameters.MonteCarloReps, 2, ncomponents_max + 1);
            
            % For each repitition.
            for repititioni = 1:parameters.MonteCarloReps

                SSEs_byfold = NaN(kFolds, 2, ncomponents_max + 1);
        
                % For each fold, 
                for foldi = 1:kFolds
        
                    % Make a vector of the fold numbers. Make new on each fold iteration.
                    fold_numbers_vector = 1:kFolds;
        
                    % Remove foldi from the vector, leaving only the indices to be
                    % used for training.
                    fold_numbers_vector(foldi) = [];
        
                    % Set up testing data. (One of the folds).
                    Xtest = explanatoryVariables(partition_indices{foldi, repititioni}, :);
                    Ytest = responseVariables(partition_indices{foldi, repititioni}, :);
        
                    % Set up training data. (The rest of the folds). 
        
                    % Concatenate indices of all other folds.
                    if strcmp(comparison_type, 'categorical') && stratify
                        train_indices = vertcat(partition_indices{fold_numbers_vector, repititioni});
                    else
                        train_indices = horzcat(partition_indices{fold_numbers_vector, repititioni});
                    end
                    Xtrain = explanatoryVariables(train_indices, :);
                    Ytrain = responseVariables(train_indices, :);
        
                    % Calculate model & sum squared error. 
                    SSEs_byfold(foldi, :, :) = SSEFunction(Xtrain,Ytrain,Xtest,Ytest,ncomponents_max); 
        
                end
                SSEs_byrepitition(repititioni, :, :) = squeeze(mean(SSEs_byfold, 1));
            end

            % Find mean squared error by taking mean across folds of SSEs.
            MSEP_original = squeeze(sum(SSEs_byrepitition, 1)./(parameters.MonteCarloReps * size(responseVariables, 2)));
        
            % Calculate AIC and BIC. Don't use first entry (is null model with 0 components).
            [aicy, bicy] = aicbic(-MSEP_original(2,2:end), 3:size(MSEP_original,2) - 1 + 2, size(responseVariables,1));

            % Use knee point if user says so. 
            if isfield(parameters , 'useKneePoint') && parameters.useKneePoint
                 [~, ncomponents] = knee_pt(bicy,1:numel(bicy),true);
                 ncomponents = ncomponents + 1; 
                 
            elseif isfield(parameters, 'useDifference') && parameters.useDifference

                bicy_difference = diff(bicy);
                ncomponents = find(bicy_difference > parameters.minDifference, 1);

            else
                % Get number of components to use from reponse variable BIC
                % minimum.
                [~ , ncomponents] = min(bicy);
            end 

            % Also calculate for explanatory, for completeness.
            [aicx, bicx] = aicbic(-MSEP_original(1,2:end), 3:size(MSEP_original) -1 + 2, size(responseVariables,1));

            % Concatenate aics & bics to match style of MSEP
            aic = [aicx; aicy];
            bic = [bicx; bicy];
           
            % Put MSE_original, ncomponents,aic, bic, & W_original into the results.
            results.maximal_components.MSEP = MSEP_original;
            results.ncomponents_used = ncomponents;
            results.maximal_components.AIC = aic;
            results.maximal_components.BIC = bic;

            % If user said to run with max components (to get out loadings and
            % weights and stuff), run it.
            if isfield(parameters, 'run_with_max_components') && parameters.run_with_max_components
                disp('Running with max components.');
                 [results.maximal_components.XL, results.maximal_components.YL, ~, ~, results.maximal_components.BETA, ~, ...
                       results.maximal_components.MSE, results.maximal_components.stats] = plsregress_fullcode(explanatoryVariables, responseVariables, ncomponents_max); 
            end

        % If not contguous partitions, run with random paritions.
        else 
                 [~, ~, ~, ~, ~, ~, MSEP_original, stats_original, MSEP_byVars_original] ...
               = plsregress_fullcode(explanatoryVariables, responseVariables, ncomponents_max, 'cv', parameters.crossValidationReps, 'mcreps', parameters.MonteCarloReps, 'Options', statset('UseParallel',true) );
            
            % Save the original weights of Y for later (in case you want to look at
            % what those components look like later)
            W_original = stats_original.W; 
        
            % Find component with minimum response-variable MSEP. Don't use
            % first entry (is null model with 0 components).
            [~ , ncomponents] = min(MSEP_original(2,2:end));

            % Put MSE_original, ncomponents, & W_original into the results.
            results.maximal_components.MSEP = MSEP_original;
            results.maximal_components.MSEP_byVars = MSEP_byVars_original;
            results.maximal_components.W = W_original;
            results.ncomponents_used = ncomponents;
        end 

    % Otherwise, just run with ncomponents as ncomponents_max.
    else
        ncomponents = ncomponents_max;
        results.ncomponents_used = ncomponents;
        
    end

    % Now run with optimal (or user given) number of components.

    disp(['Running PLSR with ' num2str(ncomponents) ' components.']);

    % If this isn't a second-level comparison being run on previously
    % calculated random permutations, run just the one asked-for
    % regression.
    if ~isfield(parameters, 'onPermutations') || (isfield(parameters, 'onPermutations') && ~parameters.onPermutations)
   
        [results.XL, results.YL, results.XS, results.YS, results.BETA, results.PCTVAR, results.MSEP, results.stats, results.MSEP_byVars,  results.Tnotnormal, results.pctVar_byVar] ...
          = plsregress_fullcode(explanatoryVariables, responseVariables, ncomponents); 
       
         % Calculate covariance matrix, divide by n - 1 observations to
         % make correlation matrix. (Covariance matrix was calculated on
         % the normalized variables, so you can do this). 
         results.Cov = (results.XL * results.YL')./ (size(explanatoryVariables, 1) - 1); 

         % Put results into output structure.
         parameters.results = results; 

    % If this IS a second-level comparison being run on previously
    % calculated random permutations, 
    else

        disp('Running on permutations.')

        % Set up a holder that will hold the newly generated covariance matrices (number
        % of mice + 1 for intercept x number of correlations x number of permutations)
        Covs_randomPermutations_2ndlevel = NaN(size(responseVariables,1) + 1, size(responseVariables,2), size(responseVariables,3));

        % For each permutation,
        parfor repi = 1:size(responseVariables, 3)

            % Run PLSR regression.
            [XL, YL,] = plsregress_fullcode(explanatoryVariables, responseVariables(:,:, repi), ncomponents);

            % Calculate normalize covariance matrix
            Cov = XL * YL' ./ (size(explanatoryVariables, 1) - 1);

            % Put into holder.
            Covs_randomPermutations_2ndlevel(:, :, repi) = Cov;

        end

        % Put all Covs into output structure.
        parameters.Covs_randomPermutations_2ndlevel = Covs_randomPermutations_2ndlevel;
    end

    % Run iterative permutations for permutation significance testing. Randomly
    % permute the order of the response variables. 

    % If user says so
    if isfield(parameters, 'permutationGeneration') && parameters.permutationGeneration

        disp('Running permutations'); 
        
        if isfield(parameters, 'comparison_type') && strcmp(comparison_type, 'categorical')
            columns_to_use = 1;
        else 
            columns_to_use = 1:size(explanatoryVariables, 2);
        end

        % Make a holding matrix for Cov & BETA permutations.
        % if categorical, only use one row
        if isfield(parameters, 'comparison_type') && strcmp(comparison_type, 'categorical')
            BETAs_permutations = NaN(numel(columns_to_use), size(results.BETA, 2), parameters.n_permutations);  %size(results.BETA, 2)
            Covs_permutations = NaN(numel(columns_to_use),size(results.Cov, 2), parameters.n_permutations);  %size(results.BETA, 2)
        else
            BETAs_permutations = NaN(numel(columns_to_use) + 1, size(results.Cov,2), parameters.n_permutations);  %size(results.BETA, 2)
            Covs_permutations = NaN(numel(columns_to_use), size(results.Cov, 2), parameters.n_permutations);  %size(results.BETA, 2)
        end 

        parfor repi = 1:parameters.n_permutations % parfor 

            % For each response variable being looked at (want the diffent 
            % categories to vary independently), make a mixing vector that's
            % made up of a random permutation of the number of periods
            % included.

%             % Make a holder for mixed response variables
%             responseVariables_mixed = NaN(size(responseVariables, 1), numel(columns_to_use));
% 
%             for variablei = 1:numel(columns_to_use)     %(responseVariables, 2)
%                 column = columns_to_use(variablei);
% 
%                 % Randomize order
%                 vect_mix = randperm(size(responseVariables, 1));
%     
%                 % Mix/permute response varables for this category's columns
%                 % only. 
%                 responseVariables_mixed(:, column) = responseVariables(vect_mix, column);
% 
%             end
             % Run the plsregress_fullcode on the mixed/permuted data.
           % [XL,  YL, ~, ~, BETA] = plsregress_fullcode(explanatoryVariables, responseVariables_mixed, ncomponents);

             % Make a holder for mixed explanatory variables
            explanatoryVariables_mixed = NaN(size(explanatoryVariables, 1), numel(columns_to_use));

            for variablei = 1:numel(columns_to_use)     %(explanatoryVariables, 2)
                column = columns_to_use(variablei);

                % Randomize order
                vect_mix = randperm(size(explanatoryVariables, 1));
    
                % Mix/permute explanatory varables for this category's columns
                % only. 
                explanatoryVariables_mixed(:, column) = explanatoryVariables(vect_mix, column);

            end
            
            % If categorical comparison, duplicate explanatoryVariables_mixed
            % with a negative value
            if isfield(parameters, 'comparison_type') && strcmp(comparison_type, 'categorical')
                explanatoryVariables_mixed = [explanatoryVariables_mixed -explanatoryVariables_mixed];
            end 

            % Run the plsregress_fullcode on the mixed/permuted data.
           [XL,  YL, ~, ~, BETA] = plsregress_fullcode(explanatoryVariables_mixed, responseVariables, ncomponents);
            
            % Calculate normalized covariance matrix 
            Cov = XL * YL' ./ (size(explanatoryVariables,1) - 1);
       
           permutated_output_BETA = BETA; 
           permutated_output_Cov = Cov;

            % If comparison type is categorical, keep just the second row
            % of Cov & BETA. (In inverted version, first row is intercepts
            if isfield(parameters, 'comparison_type') && strcmp(comparison_type, 'categorical')
                permutated_output_Cov = permutated_output_Cov(2, :);
                permutated_output_BETA = permutated_output_BETA(2, :);
            end 

            % Put into holding matrix.
            Covs_permutations(:, :, repi) = permutated_output_Cov; 
            BETAs_permutations(:, :, repi) = permutated_output_BETA; 

        end 

        % Put Covs_permutations into output structure, as single precision.
        parameters.Covs_randomPermutations = single(Covs_permutations);
        parameters.BETAs_randomPermutations = single(BETAs_permutations);

    end 

    % Run Bootstrappipng resampling & calculations

     % If user says so
    if isfield(parameters, 'useBootstrapping') && parameters.useBootstrapping

        disp('Running bootstrapping.'); 
%         
%         if isfield(parameters, 'comparison_type') && strcmp(comparison_type, 'categorical')
%             columns_to_use = 1;
%         else 
            columns_to_use = 1:size(responseVariables, 2);
     %   end

       % If stratifying,
       if parameters.stratify
           % For each response variable,
           for variablei = 1:size(responseVariables, 2)

                % Get out variable indices
    
                % These are normalized category labels. Assume the
                % relevant positive category is the max in this column.
                % (1 vs 0 or even -1).
                cat_value = max(responseVariables(:,variablei));
                variable_indices = find(responseVariables(:,variablei) == cat_value);
                
                [~,bootsam] = bootstrp(parameters.n_bootstraps, [], responseVariables(variable_indices));
                bootstrap_indices_holder{variablei} = variable_indices(bootsam); 
               
           end

           % Concatenate vertically. 
           bootstrap_indices = [bootstrap_indices_holder{1}; bootstrap_indices_holder{2}];
           %cellfun(@horzcat, bootstrap_indices_holder(1) bootstrap_indices_holder(2), 'UniformOutput', true);

       % If not stratifying,
       else
            [~, bootstrap_indices] = bootstrp(parameters.n_bootstraps, [], responseVariables);
           
       end 

       % Make a holding matrix for Cov permutations.
       Covs_bootstrap = NaN(size(results.Cov,1), numel(columns_to_use), size(bootstrap_indices,2));  %size(results.BETA, 2)
        
       % Now run bootstraps
       parfor repi = 1:size(bootstrap_indices,2) 
           
            indices = bootstrap_indices(:,repi);

            % Run the plsregress_fullcode on the mixed/permuted data.
            [XL, YL] = plsregress_fullcode(explanatoryVariables(indices, :), responseVariables(indices,:), ncomponents);
            
            % Calculate normalized covariance matrix.
            Cov = XL * YL' ./ (size(explanatoryVariables, 1) - 1);

            % Put into holding matrix.
            Covs_bootstrap(:, :, repi) = Cov; 
            
        end 

        % If categorical, 
        if isfield(parameters, 'comparison_type') && strcmp(comparison_type, 'categorical')
            
            % Keep only the first variable
            Covs_bootstrap = Covs_bootstrap(:, 1, :);
        end 
        
        % Convert to single precision to take up less space. 
        % Put Covs_bootstraps into output structure
        parameters.Covs_bootstrap = single(Covs_bootstrap);

    end 
end 