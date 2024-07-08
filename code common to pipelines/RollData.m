% RollData.m
% Sarah West
% 4/12/22

function [parameters] = RollData(parameters)

    % Display progress message to user.
    MessageToUser('Rolling ', parameters);

    % If the roll dimension is greater than the roll window plus 1 step
    % (can get at least 1 roll out of it) 
    % Figure out how many rolls you can get out of it (should be a
    % whole number). 
    duration = size(parameters.data, parameters.rollDim);
    number_of_rolls = CountRolls(duration, parameters.windowSize, parameters.stepSize);

    if number_of_rolls > 1
        % ***Set up new holder matrix with abstracted dimensions. Make 1
        % extra*** dimension for adding new rolls.
        dimensions_holder = repmat({':'},1, ndims(parameters.data) + 1);
       
        % Put in dimension sizes to match with data.
        for i = 1:ndims(parameters.data)
            dimensions_holder{i} = size(parameters.data, i);
        end
        
        % Replace holder dimensions with needed rolled dimensions
        dimensions_holder{parameters.rollDim} = parameters.windowSize;
        dimensions_holder{end} = number_of_rolls; 
        data_rolled = NaN(dimensions_holder{:});

        % Put new rolls in.

        % Make new holder dimensions for extraction & insertion
        dimensions_insert =  repmat({':'},1, ndims(data_rolled));
        dimensions_extract = repmat({':'},1, ndims(parameters.data));

        % For each new roll
        for i = 1:number_of_rolls
            
            dimensions_insert{end} = i;
            startframe = (i - 1) * parameters.stepSize + 1;
            endframe = startframe + parameters.windowSize - 1;
            dimensions_extract{parameters.rollDim} = [startframe:endframe];

            data_rolled(dimensions_insert{:}) = parameters.data(dimensions_extract{:}); 

        end

    % Otherwise, data_rolled is just the original data
    else 
        data_rolled = parameters.data;
    end

    % Put data_rolled into parameters.
    parameters.data_rolled = data_rolled;
    parameters.roll_number = number_of_rolls;

end