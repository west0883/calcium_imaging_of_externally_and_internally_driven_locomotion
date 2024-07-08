% extractMotorData.m
% Sarah West
% 11/16/21

% Function that takes the output files from the Arduino serial monitor
% during motorized treadmill experiements and saves data as .mat files.

function [parameters] = extractMotorData(parameters)

    putty_flag = parameters.putty_flag;

    % Display progress message to user.
    MessageToUser('Extracting ', parameters);

    % If not using PUTTY
    if ~putty_flag 
      
        % Just change the name and return, you've already read in the text at the RunAnalysis
        % step.
        parameters.trial = parameters.log; 
        
        return
                
    % If using Putty,          
    else

        % Make default NaN values for start_point and end_point 
        % so you get an error if the real start point wasn't found,
        % instead of using the wrong data.
        start_point = NaN;
        end_point = NaN;

        % Find trial starting point by looking for 'Trial ...' string at
        % beginning of each trial.
        for i = 1:size(parameters.log, 1)

            % Look for the trial number, have to convert to a
            % number and back to remove the leading 0s. Will always be the
            % 'stack' keyword-value pair 
            stack_number = CreateStrings({'stack'}, parameters.keywords, parameters.values);

            if strcmp(parameters.log{i,1}, ['Trial ' num2str(str2num(stack_number))])

                % mark the start point, break the loop
                start_point = i; 
               
                break
            end 
        end 
        
        % If start point didn't get a value, skip to next
        % trial
        if isnan(start_point)
           disp(['Start of trial ' stack_number ' not found']);
            
           parameters.trial = [];
           parameters.dont_save = true; 
           return 
        end
        
        % Now look for the end point
        
        % Get the next possible stacknumber. (Needs to
        % include next stack that wasn't recorded, because
        % a new trial is always loaded up instantly by
        % Arduino).
        next_stacknum = str2num(stack_number) + 1; 
        
        for i = start_point:size(parameters.log,1)
            if strcmp(parameters.log{i,1}, ['Trial ' num2str(next_stacknum)])

                % mark the end point as the last row before the new Trial label, break the loop
                end_point = i-1; 
                break
            end 
        end

        % If end point not found, assign the last value of
        % matrix and skip
        if isnan(end_point)
           disp(['End of trial ' stack_number ' not found, assigning to end of log.']);
           end_point = size(parameters.log,1);
        end
        
        % Assign range of this trial.
        parameters.trial = parameters.log(start_point:end_point, :);
    end 
      
end 