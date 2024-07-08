% MotorBehaviorPeriodsTable.m
% Sarah West
% 4/8/22

function [parameters] = MotorBehaviorPeriodsTable(parameters)

    % Call this for each stack/trial
    
    % Get list of period names in this all_periods
    period_names = fieldnames(parameters.all_periods);

    % Establish all_periods_table as Conditions plus a field for time
    % ranges, a field for duration_place
    parameters.all_periods_table = parameters.Conditions;
    parameters.all_periods_table.time_ranges(:) = cell(size(parameters.Conditions,1), 1); 
    parameters.all_periods_table.duration_place(:) = cell(size(parameters.Conditions,1), 1);
   
    % List field names needed for information about each time range. 
    information_fields = {'speed', 'accel', 'previous_speed', 'two_speeds_ago'};

    % Display progress message to user.
    MessageToUser('Finding ', parameters);
    
    % For each period/condition name in all_periods,
    for periodi = 1:numel(period_names)

        % Get the information about individual occurances
        occurances_structure = getfield(parameters.all_periods, period_names{periodi});
        
        % For each instance inside occurances_structure (skip first entry,
        % is always empty
        for instancei =2:size(occurances_structure,1)
            
            % Establish logical array holder;
           % holder_logical = NaN(size(parameters.all_periods_table, 1), numel(information_fields)+1); 
            
            % Put in first logical, for the condition name;
            holder_logical = (parameters.Conditions.condition) == convertCharsToStrings(period_names{periodi});
            
            % Grab first instance of this condition (for finding NaNs)
            ind = find(holder_logical, 1);

            % Assign each information field
            for infoi = 1:numel(information_fields)
                
                 % If this condition doesn't care about this value, make the search criteria be NaN 
                 if strcmp(parameters.Conditions.(information_fields{infoi})(ind), "NaN")
                
                     holder_logical(:, infoi + 1) =  (parameters.Conditions.(information_fields{infoi})) == "NaN"; 
                    

                 else
                     % Get out relevant information.
                     value = getfield(occurances_structure, {instancei}, information_fields{infoi}); 

                     % Search and put into logical array.
                     holder_logical(:, infoi + 1) =  (parameters.Conditions.(information_fields{infoi})) == convertCharsToStrings(num2str(value)); 
                 
                    % Concatenate to holder_logical 
                 end
                 
            end

            % Get location where this instance belongs. 
            index = find(all(holder_logical,2));
            
            % If there is no matching position (can happen if it's an old data collection setup and we don't want to analyze that) 
            if isempty(index)
                
                % Skip to next instance
                continue
            end

            % Check if this instance is empty. If so, skip. 
            if isempty(occurances_structure(instancei).time_range)
                continue
            else
                % Check that the duration of this occurance matches the
                % duration of the matched index. Otherwise, skip.
                
                % Calculate duration, adding 1 for comparison to the periods
                % table duration.s
                duration_instance = occurances_structure(instancei).time_range(2) - occurances_structure(instancei).time_range(1) +1;
                duration_desired = parameters.all_periods_table.duration{index};
                if duration_instance ~= duration_desired
    
    %                 % If it's in the "finished" or "maintining" phases, only
    %                 % keep if between the continued chunk length and desired duration
    %                 if (strcmp(parameters.all_periods_table.condition{index}(1:3), 'm_f') ...
    %                     || strcmp(periods.condition{periodi}, 'm_maint') ...
    %                     || strcmp(periods.condition{periodi}, 'm_p_nochange') ...
    %                     || strcmp(periods.condition{periodi}, 'm_p_nowarn_maint'))
    %                     & (duration_instance < duration_desired) & duration_instance  
    % 
    
    
                    % Skip to next instance
                    disp(['Dimension error in ' num2str(index) ', ' parameters.all_periods_table.condition{index}]);
                    continue
                end
    
                % Put these in the correct time range entry, concatenated with
                % any previous instances.
                parameters.all_periods_table.time_ranges{index} = [parameters.all_periods_table.time_ranges{index}; occurances_structure(instancei).time_range];
            
                % do the same with duration_place
                 parameters.all_periods_table.duration_place{index} = [parameters.all_periods_table.duration_place{index}; occurances_structure(instancei).duration_place];
            end
        end 
    end 

end 