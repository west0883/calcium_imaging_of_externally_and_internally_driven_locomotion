% ExtractFlorescenceTimeseries.m
% Sarah West
% 4/5/22

% Extracts mean flourescence timeseries from source masks. 

function [parameters] = ExtractFluorescenceTimeseries(parameters)

     % Announce what extraction you're on.
     MessageToUser('Extracting ', parameters);

     % Make sure data & source dimensions agree.
     if size(parameters.data, 1) ~= size(parameters.sources, 1)
        error('Make sure your data and sources are inputted with pixels in first dimensions.');
     end 

     % Perform extraction

     % Done by matrix multiplication of stack & mask weights (weighted) or 
     % indices (not weighted), then divided by the sum of the mask (weights
     % or indices, weighted vs not weighted, respectively).
      
     % Weighted
     if isfield(parameters, 'weightedMean') && parameters.weightedMean
        parameters.timeseries = parameters.data'*parameters.sources ./sum(parameters.sources,1);
     
     % Not weighted.
     else
         parameters.timeseries = (parameters.data' * (parameters.sources > 0))./sum(parameters.sources, 1); 
     end

end