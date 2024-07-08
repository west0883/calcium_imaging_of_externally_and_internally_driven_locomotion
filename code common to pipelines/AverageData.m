% AverageData.m
% Sarah West
% 4/11/22

function [parameters] = AverageData(parameters)

    % Display progress message to user.
    MessageToUser('Averaging ', parameters);

    % Pull out data.
    data = parameters.data; 

    % If user says to remove outliers
    if isfield(parameters, 'removeOutliers') && parameters.removeOutliers

        % Remove outliers along averageDim, replace with NaNs.
        outliers = isoutlier(data, parameters.averageDim);
        data(outliers) = NaN;
    end 

    % if user says not to squeeze the data, don't
    if isfield(parameters, 'useSqueeze') && ~parameters.useSqueeze
        %  Take the mean
        average = mean(data, parameters.averageDim, 'omitnan'); 
    
        % Take the standard deviation
        std_dev = std(data, 0, parameters.averageDim, 'omitnan'); 
    else
         %  Take the mean
        average = squeeze(mean(data, parameters.averageDim, 'omitnan')); 
    
        % Take the standard deviation
        std_dev = squeeze(std(data, 0, parameters.averageDim, 'omitnan')); 
    end

    % If user says to put the average & std_dev in the same file,
    if isfield(parameters, 'average_and_std_together') && parameters.average_and_std_together
       parameters.average = [average, std_dev];
    else
        % Put them in different places/files
       parameters.average = average;
       parameters.std_dev = std_dev;
    end

end 