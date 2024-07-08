% FitCircles.m
% Sarah West
% 6/14/22

% Adjusted from puplitrack.m,  original author Justin Aronson 10/19
% Updated by Angela Nietz 6/13/2022. 
% Uses functon CircleFitByPratt.m to fit the circles.

function [parameters] = FitCircles(parameters)

    MessageToUser('Fitting ', parameters);

    % Pull out of parameters to make code easier to read/write.
    DLC_data = parameters.DLC_data;

    % If DLC_data is in wrong dimension order, permute so frames increase
    % along rows, different data types change across columns.
    if parameters.framesDim ~= 1
        DLC_data = permute(DLC_data, [parameters.framesDim parameters.dataDim]);
    end

    % Skip the first frames, indicated by parameters.skip

    % If user wants to skip frames,
    if isfield(parameters, 'skip')
        
        % Remove those frames
        DLC_data = DLC_data(parameters.skip/parameters.channelNumber + 1 :end, :);

    end 

    % **Find & remove outliers & low condfidence points from data. Replace with NaN. 

    % Find outliers in each column, using the generalized extreme
    % Studentized deviate test for outliers. 
    outliersTF = isoutlier(DLC_data,'gesd');

    % Make any data point that is an outlier & is also in the position
    % column an NaN.
    positionColumns_logical = zeros(size(DLC_data));
    positionColumns_logical(:, [parameters.xPositionColumns parameters.yPositionColumns]) = 1;
    positionColumns_logical = logical(positionColumns_logical); 

    DLC_data(outliersTF & positionColumns_logical) = NaN; 

    % Find low-confidence positions.
    low_confidence_positions = DLC_data(:, parameters.confidenceColumns) < parameters.confidenceThreshold; 

    % Turn low-confidence points into NaNs.
    for pointi = 1:parameters.numberOfPoints
        rows = find(low_confidence_positions(:, pointi));
        DLC_data(rows, [parameters.xPositionColumns(pointi) parameters.yPositionColumns(pointi)]) = NaN;
    end

    
    % ** Begin circle fitting process.

    % Create holder to hold the output of CircleFitByPratt (centroid x,
    % centroid y, radius). 
    circData = NaN(size(DLC_data, 1), 3);
 
    % Pull out certain fields of parameters so you don't have to send it
    % all to the parfor loops.
    numberOfPoints = parameters.numberOfPoints;
    xPositionColumns =  parameters.xPositionColumns;
    yPositionColumns =  parameters.yPositionColumns;

    % For each frame.
    parfor framei = 1:size(DLC_data,1) 

        % Make a holder for points. (formatting for entry to
        % CircleFitByPratt)
        points = NaN(numberOfPoints, 2);

        % For each point
        for pointi = 1:numberOfPoints

            % Put the corresponding x and y position data into their place.
            points(pointi,1) = DLC_data(framei, xPositionColumns(pointi));
            points(pointi,2) = DLC_data(framei, yPositionColumns(pointi)); 
        end

        % Remove points (rows) that have NaN values. 
        points(any(isnan(points),2), :) = [];

        % If fewer than 2 points are still good, leave all entries for this frame NaN
        if isempty(points) || size(points, 1) < 2

        else 
            % Fit circle to tracked points at each video frame.
            try
            circData(framei,:) = CircleFitByPratt(points);
            catch
            end
        end
    end

    % Calculate diameter, min diameter, max diameter
    diameter = circData(:, 3) * 2;

    % Remove outliers from the calculated diameters. 
    diameter(isoutlier(diameter)) = NaN;

    % Calculate min & max diameters.
    parameters.circle_info.max_diameter = max(diameter);
    parameters.circle_info.min_diameter = min(diameter);

    % Put diameter & circData into output structure.
    parameters.circle_info.diameter = diameter;
    parameters.circle_info.circData = circData;

end
    