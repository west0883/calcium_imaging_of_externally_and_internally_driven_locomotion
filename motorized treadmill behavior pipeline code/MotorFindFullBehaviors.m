% MotorFindFullBehaviors.m
% Sarah West
% 5/22/22

% Runs with RunAnalysis.m. Uses behavior ranges from
% MotorBehaviorPeriodsTable.m to put together full trainsition behaviors
% (from "warning" to "finished" portions).

function [parameters] = MotorFindFullBehaviors(parameters)


  % Load & define full_periods_table.

   no_warns_m = {
        'm_p_nowarn_start';
        'm_p_nowarn_stop';
        'm_p_nowarn_accel';
        'm_p_nowarn_decel';
        'm_p_nowarn_maint'}; 

   no_warns_f = {'m_fstart'; 'm_fstop';'m_faccel'; 'm_fdecel'; NaN};

   all_warnings = {'w_start', 'm_start', 'm_fstart'; 
        'w_stop', 'm_stop', 'm_fstop'; 
        'w_accel', 'm_accel', 'm_faccel';
        'w_decel', 'm_decel', 'm_fdecel';
        'w_maint', 'm_maint', NaN;
        'w_p_nochange_start', 'm_p_nochange', NaN;
        'w_p_nochange_stop', 'm_p_nochange', NaN;
        'w_p_nochange_accel', 'm_p_nochange', NaN;
        'w_p_nochange_decel', 'm_p_nochange', NaN;
        'w_p_nochange_maint', 'm_p_nochange', NaN;
        'w_p_nowarn', no_warns_m, no_warns_f };

    full_period_conditions = { 'start'; 'stop'; 'accel'; 'decel'; 'maint'; 
        'nochange_start'; 'nochange_stop'; 'nochange_accel'; 'nochange_decel'; 'nochange_maint';
        'nowarn_start'; 'nowarn_stop'; 'nowarn_accel'; 'nowarn_decel'; 'nowarn_maint'};

    % For each warning type, 
    for warningi = 1:numel(all_warnings)
        warning = all_warnings{warningi}; 

        % Get the time ranges for this warning type from this all_periods table.
        warning_table = all_periods(string(periods_table.condition)== warning);
        warning_timeranges = warning_table.time_ranges; 

        % Look at only non-empty entries.
        nonempty_indices = find(~cellfun(@isempty, warning_timeranges));

        % For each non-empty entry,
        for instancei = 1:numel(nonempty_indices)
            instance = nonempty_indices(instancei);

            % Make full time range anew, put in this warning range.
            full_timerange = warning_timeranges{instance};

            % Get the end number of the time range.
            ending = warning_timeranges{instance}(2); 

            % If this warning is w_p_nowarn, have to do something
            % different

            % If not no warning,
            if ~strcmp(warning, 'w_p_nowarn')

                % Get all the motor possibilites.
                motor_table = all_periods(string(periods_table.condition) == all_warnings{warningi, 2});
                motor_timeranges =  motor_table.time_ranges; 
                
                % Look at only non-empty motor entries.
                nonempty_indices = find(~cellfun(@isempty, warning_timeranges));

                % Look for number at beginning of time range that aligns. 
                motor_index = find(cellfun(@(x) x(1) == ending + 1, motor_timeranges(nonempty_indices))); 
                
                % Check that motor_index didn't contain more than 1 element
                % or if it's empty. 
                if numel(motor_index) > 1 
                    error('Found more than one entry.');
                elseif isempty(motor_index)
                    error('Did not find an entry.');
                end

                % Add this to full time range. 
                full_timerange = [full_timerange motor_timeranges{motor_index}];

                % Re-define ending based on this new timerange.
                ending = motor_timeranges{motor_index}(2);

                % *******Repeat for finished periods. **** 

                % Don't do if this entry is an NaN.
                if ~isnan(all_warnings{instance, 3})

                    % Get all the finished possibilites.
                    finished_table = all_periods(string(periods_table.condition) == all_warnings{warningi, 3});
                    finished_timeranges =  finished_table.time_ranges; 
                    
                    % Look at only non-empty finished entries.
                    nonempty_indices = find(~cellfun(@isempty, warning_timeranges));
    
                    % Look for number at beginning of time range that aligns. 
                    finished_index = find(cellfun(@(x) x(1) == ending + 1, finished_timeranges(nonempty_indices))); 
                    
                    % Check that finished_index didn't contain more than 1 element
                    % or if it's empty. 
                    if numel(finished_index) > 1 
                        error('Found more than one entry.');
                    elseif isempty(finished_index)
                        error('Did not find an entry.');
                    end
    
                    % Add this to full time range. 
                    full_timerange = [full_timerange finished_timeranges{finished_index}];

                end 

                % Figure out beginning speed, accel, and ending speed from
                % motor segment. 


                % Find the indices of conditions names
                full_condition_index = full_periods_table(string(full_periods_table.condition) == full_period_conditions{warningi} & string(full_periods_table.previous_speed) == previous_speed & string(full_periods_table.accel) == accel & string(full_periods_table.speed) == speed).index;
               
                % Put full_timerange into proper place in periods table. 
                full_periods.time_ranges(full_condition_index) = [full_periods.time_ranges(full_condition_index); full_timerange];
           
            
            % If it is no warning,
            else 

                % Go through each possible motor type, 
                for motori = 1:numel(no_warns_m)

                    % Get all the motor possibilites.
                    motor_table = all_periods(string(periods_table.condition) == all_warnings{warningi, 2}{motori});
                    motor_timeranges =  motor_table.time_ranges; 
                    
                    % Look at only non-empty motor entries.
                    nonempty_indices = find(~cellfun(@isempty, warning_timeranges));
    
                    % Look for number at beginning of time range that aligns. 
                    motor_index = find(cellfun(@(x) x(1) == ending + 1, motor_timeranges(nonempty_indices))); 
                    
                    % Check that motor_index didn't contain more than 1 element
                    % or if it's empty. 
                    if numel(motor_index) > 1 
                        error('Found more than one entry.');
                    elseif isempty(motor_index)
                        error('Did not find an entry.');
                    end
    
                    % Add this to full time range. 
                    full_timerange = [full_timerange motor_timeranges{motor_index}];
                    
                   
                    % *****Repeat for corresponding finished entry.****
                    
                    % Get all the finished possibilites.
                    finished_table = all_periods(string(periods_table.condition) == all_warnings{warningi, 3}{motori});
                    finished_timeranges =  finished_table.time_ranges; 
                    
                    % Look at only non-empty finished entries.
                    nonempty_indices = find(~cellfun(@isempty, warning_timeranges));
    
                    % Look for number at beginning of time range that aligns. 
                    finished_index = find(cellfun(@(x) x(1) == ending + 1, finished_timeranges(nonempty_indices))); 
                    
                    % Check that finished_index didn't contain more than 1 element
                    % or if it's empty. 
                    if numel(finished_index) > 1 
                        error('Found more than one entry.');
                    elseif isempty(finished_index)
                        error('Did not find an entry.');
                    end
    
                    % Add this to full time range. 
                    full_timerange = [full_timerange finished_timeranges{finished_index}];

                    % Find the indices of conditions names
                    full_condition_index = full_periods_table(string(full_periods_table.condition) == full_period_conditions{warningi + motori} & string(full_periods_table.previous_speed) == previous_speed & string(full_periods_table.accel) == accel & string(full_periods_table.speed) == speed).index;
               
                    % Put full_timerange into proper place in periods table. 
                    full_periods.time_ranges(full_condition_index) = [full_periods.time_ranges(full_condition_index); full_timerange];
                    
                end 
            end 

        end 
    end
end