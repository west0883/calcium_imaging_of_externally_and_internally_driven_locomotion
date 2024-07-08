% PLSR_Simplified.m
% Sarah West
% 6/6/22

% A wrapper that runs the Matlab function plsregress_fullcode.m 

function [parameters] = PLSR_Simplified(parameters)

    % for notes purposes:
    % n = number of observations
    % m = number of explanatory variables
    % j = number of response variables

    % *** Inputs *** 

    % parameters.dataset.explanatoryVariables -- n x m matrix of explanatory variables;
            % Columns should be normalized (z scored);
    % parameters.dataset.responseVariables -- n x j matrix of response 
            % variables; Columns should be normalized (z scored);
    % parameters.ncomponents_max -- scalar; (if no cross validation) number
            % of latent components you want PLSR to find; (if using cross validation)
            % the maximum number of latent components you want PLSR to try
    % parameters.comparison_type -- string; either 'continuous' or
            % 'categorical'; denotes if response variables are continuous or
            % categorical

        % ** Validation parameters **

        % parameters.findBestComponents -- true/false; whether or not to
                % find optimal number of components (whether or not to
                % cross-validate)
        % parameters.kFold -- number of 'folds' to do for k-fold cross
                % validataion
        % parameters.MonteCarloReps -- number of times you repeat the
                % k-fold cross validation
        % parameters.stratify -- true/false; if a categorical comparison,
                % need to make the number of each category is consistent across
                % partitions
        % parameters.contiguous_partitions -- true/false; if you want to do
                % cross-validation on partitions with adjacent observations instead
                % of randomly selected (should be "true" when adjacent data belongs
                % to same/similar time periods, mice, etc)

        % ** null distribution generation parameters ** 

        % parameters.permutationGeneration -- true/false; whether or not to
                % run PLSR on random permutations; Results from random
                % permutations are used for statistical validation.
        % parameters.n_permutations -- scalar; number of random
                % permutations you want to run for a null distribution (at least
                % 500 or so)

    % *** Outputs ***
    
    % parametesr.results -- structure with the following fields:
            % fields with same meaning as Matlab plsregress function output:
                % XL -- > these are the "subnetwork" latent components
                % YL
                % XS --> "scores"
                % YS
                % BETA --> normalized Betas (isn't as useful as the
                          % covariance, below)
                % PCTVAR 
                % MSEP
                % stats 

            % additional fields:
                % Cov -- m x j matrix; covariances-- important!; should
                        % use this as the change-in-value instead of the normalized BETA variable. 
                % MSEP_byVars  -- calculated MSE contribution of each
                        % variable (can only get if you use
                        % plsregress_fullcode.m instead of plsregress.m)
                % pctVar_byVar -- calculated pctVar contribution of each
                        % variable (can only get if you use
                        % plsregress_fullcode.m instead of plsregress.m)
                % n_components_used -- the number of components determined to
                        % be significant from cross-validation and used in final output
                % maximal_components == a structure with the following fields
                        % (are the same as above, but for the maximum number of
                        % components asked for if you ran cross
                        % -validataion):
                            % MSEP  
                            % MSEP_byVars 
                            % W -- (same as stats.W above)

    % parameters.Covs_randomPermutations -- m x j x number of permutations; Covariances from randomn
                % permutatoins (a null distribution of covariances)


    % **** Begin function ****

    % Tell user all the parameters they sent in (so they can cancel before 
    % it runs for a long time if they made a mistake).
    if isfield(parameters, 'findBestNComponents') && parameters.findBestNComponents

            message = ['Finding best number of components with maximum ' num2str(parameters.ncomponents_max) ...
                ' components, ' num2str(parameters.kFolds) ' -folds.'];
            disp(message);
    end

    if isfield(parameters, 'permutationGeneration') && parameters.permutationGeneration

        disp(['Will run ' num2str(parameters.n_permutations) ' permutations.'])

    end
    
    % Pull some variables out of 'parameters' structure
    comparison_type = parameters.comparison_type;
    explanatoryVariables = parameters.dataset.explanatoryVariables;
    responseVariables = parameters.dataset.responseVariables;

    % If there are fewer observations than asked-for components, lower
    % ncomponents to 1 less than number of observation
    if size(responseVariables, 1) <= parameters.ncomponents_max
        ncomponents_max = size(responseVariables, 1) - 1;
    else
        ncomponents_max = parameters.ncomponents_max;
    end


    % ********************************************************************
    % If user says so, run cross-validation
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
    
                % For each response variable,
                for variablei = 1:size(responseVariables, 2)
    
                    % Get out variable indices

                    % These are normalized category labels. Assume the
                    % relevant positive category is the max in this column.
                    % (1 vs 0 or even -1).
                    
                    variable_indices = find(responseVariables(:,variablei) > 0); 
                    
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
            parfor repititioni = 1:parameters.MonteCarloReps

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

            % Get number of components to use from reponse variable BIC
            % minimum.
            [~ , ncomponents] = min(bicy);

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


    % *********************************************************************
    % Now run with optimal (or user given) number of components.

    disp(['Running PLSR with ' num2str(ncomponents) ' components.']);

    % Run PLSR
    [results.XL, results.YL, results.XS, results.YS, results.BETA, results.PCTVAR, results.MSEP, results.stats, results.MSEP_byVars,  results.Tnotnormal, results.pctVar_byVar] ...
      = plsregress_fullcode(explanatoryVariables, responseVariables, ncomponents); 

     % Calculate covariance matrix, divide by n - 1 observations to
     % make correlation matrix. (Covariance matrix was calculated on
     % the normalized variables, so you can do this). 
     results.Cov = (results.XL * results.YL')./ (size(explanatoryVariables, 1) - 1); 

     % Put results into output structure.
     parameters.results = results; 


    % *********************************************************************
    % If user says so, run iterative permutations for permutation significance testing. Randomly
    % permute the order of the response variables. 
    if isfield(parameters, 'permutationGeneration') && parameters.permutationGeneration

        disp('Running permutations'); 
        
        % Make a holding matrix for Cov permutations.
        Covs_permutations = NaN(size(results.Cov,1), size(reponseVariables, 2), parameters.n_permutations);  %size(results.BETA, 2)

        parfor repi = 1:parameters.n_permutations 

            % For each response variable being looked at (want the diffent 
            % categories to vary independently), make a mixing vector that's
            % made up of a random permutation of the number of periods
            % included.

            % Make a holder for mixed response variables
            responseVariables_mixed = NaN(size(responseVariables, 1), size(responseVariables,2));

            for variablei = 1:size(responseVariables, 2)

                % Randomize order
                vect_mix = randperm(size(responseVariables, 1));
    
                % Mix/permute response varables for this category's columns
                % only. 
                responseVariables_mixed(:, variablei) = responseVariables(vect_mix, variablei);

            end

            % Run the plsregress_fullcode on the mixed/permuted data.
            [XL,  YL] = plsregress_fullcode(explanatoryVariables, responseVariables_mixed, ncomponents);
            
            % Calculate normalized covariance matrix 
            Cov = XL * YL' ./ (size(explanatoryVariables,1) - 1);

            % Put into holding matrix.
            Covs_permutations(:, :, repi) = Cov; 
            
        end 

        % Put Covs_permutations into output structure, as single precision.
        parameters.Covs_randomPermutations = single(Covs_permutations);

    end 
end 