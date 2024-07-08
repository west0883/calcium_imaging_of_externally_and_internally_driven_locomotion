% motor_behavior_period_structures.m
% Sarah West
% 12/9/21

% Splits the all_periods structure into a structure where each possible
% type of period is split into its own field. Is harder to read for humans,
% but will be easier for using with established code. 

function [] = motor_behavior_period_structures(parameters)
  
   % Establish base input directory
    dir_in_base=[parameters.dir_exper 'behavior\period instances\'];

    % Establish base output directory
    dir_out_base=[parameters.dir_exper 'behavior\period instances\all structure format\'];
   
    clear holding_structure; 
    
    % Re-arrange all possible periods as structures
    for condi =1:size(parameters.Conditions,2)
        short_name = parameters.Conditions(condi).short;
        
        eval(['current_type = all_periods.' short_name ';']); 
      
        % Create a field for the current type
        eval(['holding_structure.' short_name ' = [];']);
        
        % Only make fields for things we actually care about
        % per condition type
        
        switch condi
            
            % For all but start & maintaining warning cues, the two 
            % motor maintaining types, the continued walking
            % types, and the no change in motor probe,we only 
            % care about the current speed subdivisions
            case num2cell([9:11 19:21 26])
                
                % Initialize empty fields 
                for speedi = 1: size (parameters.speeds,2)
                    eval(['holding_structure.' short_name '.x' num2str(parameters.speeds(speedi)) '= [];'])
                end 
                
                for instancei = 1:size(current_type,1)
                    time_range = current_type(instancei).time_range; 
                    speed = current_type(instancei).speed;
                    
                    if ~isempty(time_range)
                        struc_name = ['holding_structure.' short_name '.x' num2str(speed)];
                        
                        % Sometimes Arduino did something weird
                        % with ending of trials, so ignore
                        % things that don't fit.
                        try 
                            eval([struc_name '= [' struc_name '; time_range];']);  
                    
                        end 
                    end 
                end
                
            % For accel, decel, faccel, & fdecel we care about the
            % previous speed, current speed, and acceleration
            case num2cell([1,2, 14, 15])
                
                % Initialize empty fields 
                for speedi = 1: size (parameters.speeds,2)
                    speed = parameters.speeds(speedi);
                    for acceli = 1:size(parameters.accels_acceldecel,2)
                        accel = parameters.accels_acceldecel(acceli);
                        
                        for speed2i = 1:size(parameters.speeds,2)
                            speed2 = parameters.speeds(speed2i);
                            eval(['holding_structure.' short_name '.x' num2str(speed) '.x' num2str(accel) '.x' num2str(speed2) '= [];'])
                        end 
                    end 
                end
                
                for instancei = 1:size(current_type,1)
                    time_range = current_type(instancei).time_range; 
                    speed = current_type(instancei).speed;
                    accel = current_type(instancei).accel;
                    speed2 = current_type(instancei).previous_speed;
                    if ~isempty(time_range) 
                        struc_name = ['holding_structure.' short_name '.x' num2str(speed) '.x' num2str(accel) '.x' num2str(speed2)];
                        eval([struc_name '= [' struc_name '; time_range];']);  
                    end 
                end 
            % For finished accel & finished decel, we care about 
            % current speed, 2 speeds ago, and acceleration
            case num2cell([6, 7])   
                % Initialize empty fields 
                for speedi = 1: size (parameters.speeds,2)
                    speed = parameters.speeds(speedi);
                    for acceli = 1:size(parameters.accels_acceldecel,2)
                        accel = parameters.accels_acceldecel(acceli);
                        
                        for speed2i = 1:size(parameters.speeds,2)
                            speed2 = parameters.speeds(speed2i);
                            eval(['holding_structure.' short_name '.x' num2str(speed) '.x' num2str(accel) '.x' num2str(speed2) '= [];'])
                        end 
                    end 
                end
                
                for instancei = 1:size(current_type,1)
                    time_range = current_type(instancei).time_range; 
                    speed = current_type(instancei).speed;
                    accel = current_type(instancei).accel;
                    speed2 = current_type(instancei).two_speeds_ago;
                    if ~isempty(time_range) 
                        struc_name = ['holding_structure.' short_name '.x' num2str(speed) '.x' num2str(accel) '.x' num2str(speed2)];
                        eval([struc_name '= [' struc_name '; time_range];']);  
                    end 
                end 
                
                
            % For start, no warning start, & finished start, care only about current speed and accel
            case num2cell([24, 25, 28])
                % Initialize empty fields 
                for speedi = 1: size (parameters.speeds,2)
                    speed = parameters.speeds(speedi);
                    for acceli = 1:size(parameters.accels_startstop, 2)
                        accel = parameters.accels_startstop(acceli);
                        eval(['holding_structure.' short_name '.x' num2str(speed) '.x' num2str(accel) '= [];'])
                
                    end 
                end
                
                for instancei = 1:size(current_type,1)
                    time_range = current_type(instancei).time_range; 
                    speed = current_type(instancei).speed;
                    accel = current_type(instancei).accel;
                    if ~isempty(time_range) 
                        struc_name = ['holding_structure.' short_name '.x' num2str(speed) '.x' num2str(accel)];
                        eval([struc_name '= [' struc_name '; time_range];']);  
                    end 
                end 
                
            % For stop & no warning stop, care only about
            % previous speed and accel. Include 0 in speeds
            % becuase sometimes it's weird.
            case num2cell([4, 16])
                % Initialize empty fields 
                all_speeds = [0 parameters.speeds];
                for speedi = 1: size (all_speeds,2)
                    speed = all_speeds(speedi);
                    for acceli = 1:size(parameters.accels_startstop,2)
                        accel = parameters.accels_startstop(acceli);
                        eval(['holding_structure.' short_name '.x' num2str(speed) '.x' num2str(accel) '= [];'])
                
                    end 
                end
                
                for instancei = 1:size(current_type,1)
                    time_range = current_type(instancei).time_range; 
                    speed = current_type(instancei).previous_speed;
                    accel = current_type(instancei).accel;
                    if ~isempty(time_range) 
                        struc_name = ['holding_structure.' short_name '.x' num2str(speed) '.x' num2str(accel)];
                        eval([struc_name '= [' struc_name '; time_range];']);  
                    end 
                end 
                
            % For finished stop, care only about 2 speeds ago
            % and accel. Need to include 0.
            case 5 
                % Initialize empty fields
                all_speeds = [0 parameters.speeds];
                for speedi = 1: size (all_speeds,2)
                    speed = all_speeds(speedi);
                    for acceli = 1:size(parameters.accels_startstop, 2)
                        accel = parameters.accels_startstop(acceli);
                        eval(['holding_structure.' short_name '.x' num2str(speed) '.x' num2str(accel) '= [];'])
                
                    end 
                end
                
                for instancei = 1:size(current_type,1)
                    time_range = current_type(instancei).time_range; 
                    speed = current_type(instancei).two_speeds_ago;
                    accel = current_type(instancei).accel;
                    if ~isempty(time_range) 
                        struc_name = ['holding_structure.' short_name '.x' num2str(speed) '.x' num2str(accel)];
                        eval([struc_name '= [' struc_name '; time_range];']);  
                    end 
                end 
                
            % For warning start, warning start probe, and continued rest don't
            % care about anything
            case num2cell([8,18, 27])
                 for instancei = 1:size(current_type,1)
                    time_range = current_type(instancei).time_range; 
                    
                    if ~isempty(time_range)
                        struc_name = ['holding_structure.' short_name];
                        eval([struc_name '= [' struc_name '; time_range];']);  
                    end 
                 end 
                
            % Motor maintaining, warning maintaining, probe warning maintaining,and motor no change
            % care about all speeds, including 0; and no
            % warning
            case num2cell([3,12,13,17,22, 23])   
                 
                all_speeds = [0 parameters.speeds];
                % Initialize empty fields 
                for speedi = 1: size (all_speeds,2)
                    eval(['holding_structure.' short_name '.x' num2str(all_speeds(speedi)) '= [];'])
                end 
                
                for instancei = 1:size(current_type,1)
                    time_range = current_type(instancei).time_range; 
                    speed = current_type(instancei).speed;
                    
                    if ~isempty(time_range)
                        struc_name = ['holding_structure.' short_name '.x' num2str(speed)];
                        eval([struc_name '= [' struc_name '; time_range];']);  
                    end 
                end 

        end 
    end 
% Rename holding structure
parameters.all_periods_structure = holding_structure;
          

end