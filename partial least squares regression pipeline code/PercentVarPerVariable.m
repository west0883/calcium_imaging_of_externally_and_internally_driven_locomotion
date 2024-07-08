% PercentvarPerVariable
% Sarah West
% 8/4/22

% Calculates the percent variance of brain data (explanatory or "X" data) 
% explained by each behavior variable (response or "Y" data) in a PLSR
% regression. Is run with RunAnalysis.

% Inputs: 
% parameters.XL -- the loadings of each component/latent variable for the explanatory
% variables. # of X variables x # components.
% parameters.YL -- is the loadings of each component/latent variable for the response
% variables. # of Y variables x # components.
% parameters.X0 -- the normalized explanatory variable dataset. #
% observations x # of X variables.

% Outputs: 
% parameters.variance_perVar -- variance of X0 explained by each Y
% variable. 1 x # of Y variables. Is a ratio of 1, not actually percents.

% parameters.comparison_type -- if this is a categorical or continuous
% variable response. 
function [parameters] = PercentVarPerVariable(parameters)

    MessageToUser('Percent var for ', parameters);

    % Calculate variance of X0 explained per response variable. 
    % Is the variance of the covariance matrix (XL * YL')
    % divided by the total variance of the X0 data. 

    var_of_cov = abs(parameters.XL * parameters.YL');
   

    var_of_X0 = sum(parameters.X0.^2,1);
    variance_allVars= var_of_cov./var_of_X0'; % ./ var_of_X0'; % Total proportion of var of each X variable, per Y variable.
    variance_perVar = mean(variance_allVars, 1); 
   
    % Put into output structure.
    parameters.variance_perVar = variance_perVar;
    parameters.variance_allVars = variance_allVars;

end