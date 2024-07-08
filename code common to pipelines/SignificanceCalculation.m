% SignificanceCalculation.m
% Sarah West
% 6/20/22

% Takes the "real" sample values and compares it to a null distribution of
% values to calculate p-values. Run with RunAnalysis.
% Inputs:
% parameters.test_values
% parameters.null_distribution
% parameters.shufflesDim
% parameters.alphaValue -- the statistical alpha value.
% parameters.twoTailed -- true/false-- if the statistical significance should be checked
% as a two-tailed test. Default is "true"

% Outputs:
% parameters.significance -- .pvals_raw, .increase,
% .decrease, .all. 


function [parameters] = SignificanceCalculation(parameters)

    % Give user a progress message
    MessageToUser('Testing significance of ', parameters);

    % Find a p value by fitting a normal distibution to the null data & 
    % calculating the cumulative distribution function at the "real" sample values. 

    % if not using bootstrapping
    if ~isfield(parameters, 'useBootstrapping') || ~parameters.useBootstrapping 
        
        % If user doesn't want to use a normal distribution,
        if isfield(parameters, 'useNormalDistribution')  && ~parameters.useNormalDistribution
    
            % Divide alpha value in half (to get increase or decrease)
            alpha = parameters.alphaValue ./ 2;
    
            % (can't get out p-values from this, exactly)
            % Sort null distribution data
            null_distribution_sorted = sort(parameters.null_distribution, 2); 
    
            increase_threshold_element_index = round(size(null_distribution_sorted,2) * (1 - alpha));
            decrease_threshold_element_index = round(size(null_distribution_sorted,2) * alpha); 
    
            % threshold values. 
            increase_threshold_values = null_distribution_sorted(:, increase_threshold_element_index);
            decrease_threshold_values = null_distribution_sorted(:, decrease_threshold_element_index);
    
            % If the sizes of parameters.test_values & the threshold matrices
            % don't match, transpose the threshold matrices. 
            if ~isequal([size(parameters.test_values)], [size(increase_threshold_values)])
    
                increase_threshold_values = increase_threshold_values';
                decrease_threshold_values = decrease_threshold_values';
            end
            
            % Calculate significance for increase & decrease.
            significance.increase = parameters.test_values > increase_threshold_values;
            significance.decrease = parameters.test_values < decrease_threshold_values;
    
            % All significant results.
            significance.all = significance.increase | significance.decrease;

    
        % Otherwise, use a normal distribution
        else
    
           % Set up holders for means, standard deviations, and results (ztest won't take
           % vectors).
           mus = mean(parameters.null_distribution, parameters.shufflesDim, 'omitnan');
           sigmas =  std(parameters.null_distribution,[], parameters.shufflesDim, 'omitnan');
           pvals_raw = NaN(size(mus));
           test_values = parameters.test_values;
           for variablei = 1:numel(mus)
                mu = mus(variablei);
                sigma = sigmas(variablei);
                value = test_values(variablei);
    
                [~, pval] = ztest(value, mu, sigma);
                pvals_raw(variablei) = pval;
                
           end 
            
           significance.pvals_raw = pvals_raw;
    
           % Default is 2-tailed
           % If  2-tailed
           % if ~isfield(parameters, 'twoTailed') || parameters.twoTailed
     
                % Divide alpha value in half (to get increase or decrease)
                alpha = parameters.alphaValue ./ 2;
           % end
                
           % If using false discovery rate,
           if isfield(parameters, 'useFDR') && parameters.useFDR 

               decreases = 1.000000 - significance.pvals_raw; 
               increases = significance.pvals_raw;

               pvals = min(decreases, increases); 

               [significance.all, significance.crit_p] = fdr_bh(pvals, alpha, 'dep','no');

               significance.increase = significance.pvals_raw > 0.5 & significance.all;
               significance.decrease = significance.pvals_raw < 0.5 & significance.all;

           % If no false discovery rate,
           else

                % Calculate significance for increase & decrease.
                significance.increase = significance.pvals_raw > (1 - alpha); 
                significance.decrease = significance.pvals_raw < alpha; 

                % All significant results.
                significance.all = significance.increase | significance.decrease;
               

           end
    

        end

    % Use bootstrapping
    else
        % Divide alpha value in half (to get increase or decrease)
        alpha = parameters.alphaValue ./ 2;

        % (can't get out p-values from this, exactly)
        % Sort null distribution data
        null_distribution_sorted = sort(parameters.null_distribution, 2); 

        % If not using normal distribution,
        if ~isfield(parameters, 'useNormalDistribution') ||  ~parameters.useNormalDistribution
            upper_threshold_element_index = round(size(null_distribution_sorted,2) * (1 - alpha));
            lower_threshold_element_index = round(size(null_distribution_sorted,2) * alpha); 
    
            % threshold values. 
            upper_threshold_values = null_distribution_sorted(:, upper_threshold_element_index);
            lower_threshold_values = null_distribution_sorted(:, lower_threshold_element_index);
    
            % Calculate significance for upper & lower.
            significance.increase = lower_threshold_values > 0;
            significance.decrease = upper_threshold_values < 0;
    
            % All significant results.
            significance.all = significance.increase | significance.decrease;

        % If using normal distribution,
        else

           % Set up holders for means, standard deviations, and results (ztest won't take
           % vectors).
           mus = mean(parameters.null_distribution, parameters.shufflesDim, 'omitnan');
           sigmas =  std(parameters.null_distribution,[], parameters.shufflesDim, 'omitnan');
           pvals_raw = NaN(numel(mus),1);
           test_value = 0;
           for variablei = 1:numel(mus)
                mu = mus(variablei);
                sigma = sigmas(variablei);
     
                [~, pval] = ztest(test_value, mu, sigma);
                pvals_raw(variablei) = pval;
                
           end 
            
           significance.pvals_raw = pvals_raw;

            % Divide alpha value in half (to get increase or decrease)
            alpha = parameters.alphaValue ./ 2;

            % Calculate significance for increase & decrease.
            significance.increase = pvals_raw > (1 - alpha); 
            significance.decrease = pvals_raw < alpha; 

            % All significant results.
            significance.all = significance.increase | significance.decrease;

        end
      
    end
    
    % Put significance into parameters structure.
    parameters.significance = significance;

end 