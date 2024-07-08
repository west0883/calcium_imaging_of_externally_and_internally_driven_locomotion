% DurationPlaceForFluorescencePLSR.m
% Sarah West
% 7/27/23

% Takes duration_place for cont. rest & walk averaged across 1 second
% You don't have to worry about rolling overlap, becuase you don't do that 


function [parameters] = DurationPlaceForFluorescencePLSR(parameters)

    MessageToUser('Duration for ', parameters);

    data = parameters.data; 
    window_length = parameters.window_length;
    fps = parameters.fps;

    % If data is empty, skip it (all non-continued periods)
    if isempty(data)
        
        parameters.duration_place_new = [];
        
        return
    end 


    % Each 1-s duration_place is the END of the 20 time-point window 
    duration_place_new = NaN(size(data, 1), window_length);

    holder = data - 1;

    for i = 1:size(data, 1)

        holder2 = holder(i) : 1/fps : (data(i) - 1/fps) ; 
    
        duration_place_new(i, :) = holder2;
    end


    % Reshape so it's all in one row, transpose to column
    duration_place_new_reshaped = reshape(duration_place_new', 1,  size(data, 1) * window_length)'; 

    parameters.duration_place_new = duration_place_new_reshaped;


end 