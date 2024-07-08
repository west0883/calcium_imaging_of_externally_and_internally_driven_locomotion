% ZScoreData.m
% Sarah West
% 5/27/22
% Zscore data (across all dimensions) using RunAnalysis. 

function [parameters] = ZScoreData(parameters)

    MessageToUser('Zscoring ', parameters);

    [parameters.data_zscored, mu, sigma] = zscore(parameters.data, 0, parameters.zscoreDim);

    parameters.normal_values = [mu, sigma]; 

end 