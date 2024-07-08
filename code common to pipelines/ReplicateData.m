% ReplicateData.m
% Sarah West
% 4/13/22

% Replicates data for use with RunAnalysis.m. Uses a string cell array as
% input for dimensions.

function [parameters] = ReplicateData(parameters)
    
    % Display progress message to user.
    MessageToUser('Replicating ', parameters);

    toReplicate_string = CreateStrings(parameters.toReplicate, parameters.keywords, parameters.values);
    eval(['toReplicate = ' toReplicate_string ';']);
    
    replicateDims_string = CreateStrings(parameters.replicateDims, parameters.keywords, parameters.values);
    eval(['replicateDims = ' replicateDims_string ';']);

    parameters.data_replicated = repmat(toReplicate, replicateDims);

end 