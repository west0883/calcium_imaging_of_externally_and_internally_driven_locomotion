% motorFindBehaviorPeriods.m
% Sarah West
% 11/16/21

% Makes time ranges of different behavior periods output-ed by the Arduino
% program controlling the random motor stages. 

% Activity tag list of possible periods 
%  *  1 = motor accelerating
%  *  2 = motor decelerating
%  *  3 = motor maintaining
%  *  4 = motor stopping
%  *  5 = motor finished stopping // Keeping this because it quickly finds the start of a rest period.
%  *  6 = motor reached faster speed
%  *  7 = motor reached slower speed
%  *  8 = tone: starting
%  *  9 = tone: stopping
%  *  10 = tone: accelerating
%  *  11 = tone: decellerating
%  *  12 = tone: maintaining
%  *  13 = Warning cue: probe, no warning.
%  *  14 = Motor: probe, no warning tone, accelerating.
%  *  15 = Motor: probe, no warning tone, decelerating.
%  *  16 = Motor: probe, no warning tone, stopping.
%  *  17 = Motor: probe, no warning tone, maintaining.
%  *  18 = Warning cue: probe, starting cue, no change in motor.
%  *  19 = Warning cue: probe, stopping cue, no change in motor.
%  *  20 = Warning cue: probe, accerlerating cue, no change in motor.
%  *  21 = Warning cue: probe, decelerating cue, no change in motor.
%  *  22 = Warning cue: probe, maintaining cue, no change in motor.
%  *  23 = Motor probe: no change.
%  *  24 = motor starting
%  *  25 = motor finished starting    
%  *  26 = Continued walking.
%  *  27 = Continued rest.
%  *  28 = Motor: probe, no warning tone, starting

% Input: 
% useAccel -- number of accelerations used. If only one accel used, the columns will be
% different.
% trial -- an nx 4 or n x 5 cell created by extractMotorData.m. Time, speed, [accel,] activity tag, message 
% Output: 
% all_periods -- a structure that holds all the behavior periods.
function [parameters] = motorFindBehaviorPeriods(parameters)
    
    % Display progress message to user.
    MessageToUser('Finding ', parameters);
   
    % For convenience
    trial = parameters.trial; 
    duration_place_maximum_default = parameters.duration_place_maximum_default;

    % Put all info per stage into a structure called behavior_period, then 
    % concatenate each field per activity tag. 
    
    % Find putty and accel flags, based on inputed location.
    putty_flag_location_string = CreateStrings(parameters.location_putty_flag, parameters.keywords, parameters.values);
    eval(['parameters.putty_flag = ' putty_flag_location_string ';']); 
  
    accel_flag_location_string = CreateStrings(parameters.location_accel_flag, parameters.keywords, parameters.values);
    eval(['parameters.useAccel = ' accel_flag_location_string ';']); 
    

    % Set columns based on useAccel, putty_flag
    
    % If accels used
    if parameters.useAccel > 1 
        
        % And no putty
        if ~strcmp(parameters.putty_flag, 'yes')
            time_column = 1;
            speed_column = 2;
            accel_column = 3; 
            activity_column = 4; 
        
        % And with putty.
        else
            time_column = 2;
            speed_column = 3;
            accel_column = 4; 
            activity_column = 5; 
        end
    else
        % & No putty, 
        if ~strcmp(parameters.putty_flag, 'yes')
            time_column = 1;
            speed_column = 2;
            activity_column = 3;
        
        % & with putty,
        else
            time_column = 2;
            speed_column = 3;
            activity_column = 4;
        end
        % Set a defaul accel that was used. 
        default_accel = 800; 
    end 
  
    % Initialize empty behavior condition variables with fields
    for condi = 1:size(parameters.Conditions,2)
        all_periods.(parameters.Conditions(condi).short).time_range = [];
        all_periods.(parameters.Conditions(condi).short).speed = []; % The current speed. Is the ending speed for speed transitions
        all_periods.(parameters.Conditions(condi).short).accel = []; % Needed for speed transitions
        all_periods.(parameters.Conditions(condi).short).previous_speed = []; % Needed for immediately post-transition. Is the starting speed for transition periods 
        all_periods.(parameters.Conditions(condi).short).previous_accel = []; % Needed for immediately post-transition
        all_periods.(parameters.Conditions(condi).short).two_speeds_ago = []; % May need for immediately post-transition analysis (to know what kind of transition just happened)
        all_periods.(parameters.Conditions(condi).short).duration_place = []; % Needed for continued rest and walk
    end
    continued_rest_long = []; 
    continued_walk_long_1600 = []; 
    continued_walk_long_2000 = []; 
    continued_walk_long_2400 = []; 
    continued_walk_long_2800 = []; 

    % Find the start of the trial by finding the string 'Starting Mouse
    % Runner' in the first colomn. 
    % Set up a holder for start_point
    start_point = NaN; 
    % Check each entry of trial.
    for i = 1:size(trial, 1)
        if strcmp(trial{i,1}, 'Starting Mouse Runner')
            
            % When you find it, get the index and break out of the loop
            start_point = i; % find(strcmp(trial{:,1}, 'Starting Mouse Runner'));
            
            break
        end 
    end 

    % If no start point found, display to user, tell RunAnalysis not to save, and skip to next iteration.
    if isnan(start_point)

        disp('No start point found.')
        parameters.dont_save = true;
        return

    end
    
    % Go through every entry after the start point, not including the last
    % 'finished stopping' or 'Done' entries
    for i =start_point + 1  : size(trial,1) - 2
        
        % Pull out the time range of this stage. Need to subtract out the
        % first reported time.
        behavior_period.time_range = [trial{i,  time_column} - trial{start_point +1 , time_column},  trial{i + 1,  time_column} - trial{start_point +1 , time_column} - 1]; 
        
        % Convert the time range to hemo-corrected frame range. (Convert to 
        % seconds, then multiply by brain sampling rate & round to nearest integer.) 
        behavior_period.time_range = round(behavior_period.time_range ./ parameters.wheel_Hz .* parameters.fps);
     
        % If the END of the time range falls within the skipped time range, don't
        % record it.
        if behavior_period.time_range(2) <= parameters.skip/parameters.channelNumber
            continue 
        % If only the beginning of the time range falls before the skipped 
        % time range, truncate that period to not include the skip, make the activity tag be the continued rest    
        elseif behavior_period.time_range(1) <= parameters.skip/parameters.channelNumber
             behavior_period.time_range(1) = parameters.skip/parameters.channelNumber + 1; 
             trial{i, activity_column} = 27;
        end 
         
        % If two consecutive times are the same (happens when motor is
        % already stopped), skip this stage and continue next iteration
        if trial{i, time_column} == trial{i + 1, time_column}
            continue
        end 
        
        % Round to nearest multiple of parameters.smallest_time_increment.
        behavior_period.time_range = parameters.smallest_time_increment .* round(behavior_period.time_range ./ parameters.smallest_time_increment); 
        
        % Add + 1 to start of time range, so each time point belongs to exactly one period 
        behavior_period.time_range(1) = behavior_period.time_range(1) + 1;
        
        % Subtract out the skip from the time range.
        behavior_period.time_range = behavior_period.time_range - parameters.skip/parameters.channelNumber; 

        % Make sure the time range doesn't excede the max frame number
        % (already converted to fps, skip subtracted). If it doers, skip to
        % next stage. (This can happen sometimes with mis-timings from
        % Arduino/Arduino told to run an extra few seconds for safety).
        if any(behavior_period.time_range > parameters.frames)
            continue
        end
        
        % Pull out the speed, previous speed, and the  "activity tag" for
        % categorizing.
        behavior_period.speed = trial{i, speed_column}; 
        if parameters.useAccel > 1
            behavior_period.accel = trial{i, accel_column}; 
        else
            behavior_period.accel = default_accel;
        end

        activity_tag = trial{i, activity_column};
    
        % If it's the first stage, set previous speed to 0, just in case
        if i == start_point +1 
            behavior_period.previous_speed = 0;
            behavior_period.previous_accel = 0;
        else 
            % Otherwise, get speed and accel from previous stage.
            behavior_period.previous_speed = trial{i - 1, speed_column};
            if parameters.useAccel > 1 
                behavior_period.previous_accel = trial{i - 1, accel_column};
            else
                behavior_period.previous_accel = default_accel;
            end
        end
        
        % If the stage isn't one of the first two stages, get the 2 speeds
        % ago value.
        if i > start_point + 2
            behavior_period.two_speeds_ago = trial{i - 2, speed_column};
        else
            % Make it empty so you don't get yelled at by Matlab field
            % concatenation sizes.
            behavior_period.two_speeds_ago = []; 
        end 
        
        % Divide acceleration tag into actual accel and starting
        
        % Motor accelerating
        if activity_tag == 1 

            if behavior_period.previous_speed == 0

                    % Starting
                    activity_tag = 24;

            end
        
        % Motor finished accelerating     
        elseif activity_tag == 6     
             
            % Have to use the speed from 2 stages ago. If 0, is a finished
            % starting. 
            if behavior_period.two_speeds_ago == 0
                
                activity_tag = 25; 
                
            end 
        
        % No warning, motor accelerating    
        elseif activity_tag == 14
            if behavior_period.previous_speed == 0

                    % Starting
                    activity_tag = 28;

            end
            
        % Motor maintaining at rest is reported as motor finished stopping
        % (skips the motor stopping because that's always entered as the
        % same time occurance)
        elseif activity_tag == 5
            
            % If two speeds ago(skip the motor stopping report), then this
            % is actually a motor maintaining.
            if behavior_period.two_speeds_ago == 0
                activity_tag = 3;
            end
         
        % For motor probe no change, change to "maintaining" if the
        % previous cue was a warning tone for maintaining.
        elseif activity_tag == 23 && trial{i-1, activity_column} == 12
            activity_tag = 3;
        end

        % Now concatenate fields based on activity_tag. All concatenated at
        % once. Must divide the "motor: finished ... " periods into recently
        % finished (< 3 seconds) and the "continued" periods. Also divide 
        % motor:maintaining into 3 s of maintaining and appropriate type of 
        % continued. Have to do this here because otherwise you're increasing 
        % the size of your trial list.
        switch activity_tag
            
            % Rest, walking, maintaining, and motor probe no change periods
            case num2cell([3, 5, 6, 7, 23, 25, 27]) 
                
                % If the available time range is above the specified
                % continued window, divide them and concatenat separately
                if (behavior_period.time_range(2) - behavior_period.time_range(1)) > (parameters.continued_window * parameters.fps -1)
                   
                   % Get all the associated parameters. 
                   finished_behavior_period = behavior_period;
                   continued_behavior_period = behavior_period;
                   
                   % Divide
                   finished_behavior_period.time_range = [behavior_period.time_range(1), [behavior_period.time_range(1)+parameters.continued_window*parameters.fps-1]];
                   % If activity tag is already 27 (continued rest), also
                   % divide the finished period and put it into the
                   % continued periods, because this only occurs in the
                   % very first period. 
                   if activity_tag == 27
                      continued_behavior_period.time_range = behavior_period.time_range;
                   else
                      continued_behavior_period.time_range = [[behavior_period.time_range(1) + parameters.continued_window * parameters.fps], behavior_period.time_range(2)]; 
                   end
                   
                   % Further divide REMAINING continued periods into the desired
                   % time-chunks. Size of chunkcs given by parameters.continued_chunk_length . 
                   period_length = continued_behavior_period.time_range(2) - (continued_behavior_period.time_range(1)-1); 


                   % If this instances is greater than the desired chunk
                   % length
                   if period_length > parameters.continued_chunk_length * parameters.fps 
                      
                      % find how many chunks the instance can make 
                      quotient=floor(period_length/(parameters.continued_chunk_length * parameters.fps));

                      % make the first chunk start where the instance starts
                      new_chunk_start=continued_behavior_period.time_range(1); 

                      % make the first chunk end 1 time window length after the
                      % start of the instance
                      new_chunk_end=continued_behavior_period.time_range(1)+parameters.continued_chunk_length * parameters.fps-1; 

                      % Make list of broken down chunks for the stack.
                      brokendown=[new_chunk_start, new_chunk_end]; 

                      if continued_behavior_period.time_range(1) == 1
                          duration_place = duration_place_maximum_default;
                      else 
                          duration_place = 1;
                      end 

                      % if the instance can create more than 1 3-second chunk
                      if quotient>1
                          for quotienti=1:(quotient -1) % don't use the last one because you're cycling through the *start* of each chunk

                              % find the start of the given chunk
                              new_chunk_start= continued_behavior_period.time_range(1) + parameters.continued_chunk_length * parameters.fps*quotienti;

                              % find the end of the given chunk
                              new_chunk_end = continued_behavior_period.time_range(1) + parameters.continued_chunk_length * parameters.fps*(quotienti +1)-1;

                              % concatenate the chunk into your list of 
                              brokendown=[brokendown; new_chunk_start, new_chunk_end ] ;

                              % if the continued period is the very start
                              % of the stack, make duration place the
                              % maximum default
                              if continued_behavior_period.time_range(1) == 1
                                  duration_place = [duration_place; duration_place_maximum_default];
                              else 
                                  duration_place = [duration_place; quotienti + 1];
                              end 
                          end 
                      end
                   else
                       %if it isn't too long (can only happen if exactly the length of the time window)
                       %make only 1 chunk using the start and stop of the
                       %instancee
                       brokendown=[continued_behavior_period.time_range];   
                       duration_place = [1];
                   end

                   continued_behavior_long = continued_behavior_period.time_range;

                   % Replace continued instances with broken-down versions,
                   % have to do each individually. 
                   for chunki = 1:size(brokendown,1)
                        continued_behavior_period(chunki,1).time_range = brokendown(chunki, :);
                        continued_behavior_period(chunki,1).speed = continued_behavior_period(1).speed;
                        continued_behavior_period(chunki,1).accel = continued_behavior_period(1).accel;
                        continued_behavior_period(chunki,1).previous_speed = continued_behavior_period(1).previous_speed;
                        continued_behavior_period(chunki,1).previous_accel = continued_behavior_period(1).previous_accel;
                        continued_behavior_period(chunki,1).two_speeds_ago = continued_behavior_period(1).two_speeds_ago;
                        continued_behavior_period(chunki,1).duration_place = duration_place(chunki);
                   end 

                   % Concatenate finished period. Don't include if the
                   % activity tag is already 27 (continued rest), because
                   % that only occurs in the very first rest period, and
                   % you don't want the first 3 seconds double-counted.
                   if activity_tag ~=27
                       eval(['all_periods.' parameters.Conditions(activity_tag).short '= [all_periods.' parameters.Conditions(activity_tag).short '; finished_behavior_period];']);
                   end 

                   % Concatenate continued period, depending on if rest,
                   % maintaining, motor probe no change, or other (walk). 
                   if activity_tag == 5
                       % Put in continued rest
                       continued_rest_long = [continued_rest_long; continued_behavior_long];
                       eval(['all_periods.' parameters.Conditions(27).short '= [all_periods.' parameters.Conditions(27).short '; continued_behavior_period];']);
                  
                   elseif activity_tag == 3 || activity_tag == 23 || activity_tag == 27
                       % If maintaining or motor probe no change, need to see if this is at rest or
                       % at walking.

                       % If at rest, put in continued rest.
                       if behavior_period.speed == 0
                           continued_rest_long = [continued_rest_long; continued_behavior_long];
                           eval(['all_periods.' parameters.Conditions(27).short '= [all_periods.' parameters.Conditions(27).short '; continued_behavior_period];']);

                       % Otherwise, put in continued walking
                       else
                           switch behavior_period.speed 
                               case 1600
                                    continued_walk_long_1600 = [continued_walk_long_1600; continued_behavior_long];
                               case 2000
                                    continued_walk_long_2000 = [continued_walk_long_2000; continued_behavior_long];
                               case 2400
                                    continued_walk_long_2400 = [continued_walk_long_2400; continued_behavior_long];
                               case 2800
                                    continued_walk_long_2800 = [continued_walk_long_2800; continued_behavior_long];
                           end 
                           
                           eval(['all_periods.' parameters.Conditions(26).short '= [all_periods.' parameters.Conditions(26).short '; continued_behavior_period];']);
                       end

                   else
                       % Put in continued walking
                       switch behavior_period.speed 
                           case 1600
                                continued_walk_long_1600 = [continued_walk_long_1600; continued_behavior_long];
                           case 2000
                                continued_walk_long_2000 = [continued_walk_long_2000; continued_behavior_long];
                           case 2400
                                continued_walk_long_2400 = [continued_walk_long_2400; continued_behavior_long];
                           case 2800
                                continued_walk_long_2800 = [continued_walk_long_2800; continued_behavior_long];
                       end 

                       eval(['all_periods.' parameters.Conditions(26).short '= [all_periods.' parameters.Conditions(26).short '; continued_behavior_period];']);
                   end
                   
                % If time range isn't long enough, concatenate just as the "finished"
                else    
                     eval(['all_periods.' parameters.Conditions(activity_tag).short '= [all_periods.' parameters.Conditions(activity_tag).short '; behavior_period];']);
                end 
            
            % If periods don't need to be divided, concatenate normally     
            otherwise     
                duration = behavior_period.time_range(2) - behavior_period.time_range(1);
                behavior_period.duration_place = [];
                eval(['all_periods.' parameters.Conditions(activity_tag).short '= [all_periods.' parameters.Conditions(activity_tag).short '; behavior_period];']);
        end 
   
    end

    % Put all_periods into parameters.
    parameters.all_periods = all_periods;
    parameters.long_periods.rest = continued_rest_long;
    parameters.long_periods.walk_1600 = continued_walk_long_1600;
    parameters.long_periods.walk_2000 = continued_walk_long_2000;
    parameters.long_periods.walk_2400 = continued_walk_long_2400;
    parameters.long_periods.walk_2800 = continued_walk_long_2800;
end

