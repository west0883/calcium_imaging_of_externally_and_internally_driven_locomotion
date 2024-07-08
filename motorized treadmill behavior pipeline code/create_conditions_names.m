% create_conditions_names.m
% Sarah West

%% Parameters for directories
clear all;

% Time for continued time chunks, in seconds.
continued_chunk_length = 1;

experiment_name=['Random Motorized Treadmill\'];

dir_base='Y:\Sarah\Analysis\Experiments\';
dir_exper=[dir_base experiment_name '\']; 

dir_out=dir_exper; 
mkdir(dir_out);

% Frame rate of brain recordings-- the frame rate all calculations will be
% done in.
fps = 20; 


% Time for finished xyz chunkcs, in seconds
finished_chunk_length = 3; 

% Warning period length (seconds)
warning_length = 5; 

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
%  *  24 = Motor: starting
%  *  25 = Motor: finished starting
%  *  26 = Continued walking.
%  *  27 = Continued rest.
%  *  28 = Motor: provbe, no warning tone, starting

%%
Conditions(1).name = 'motor accelerating';
Conditions(1).short = 'm_accel';
Conditions(2).name = 'motor decelerating';
Conditions(2).short = 'm_decel';
Conditions(3).name = 'motor maintaining';
Conditions(3).short = 'm_maint';
Conditions(4).name = 'motor stopping';
Conditions(4).short = 'm_stop';
Conditions(5).name = 'motor finished stopping';
Conditions(5).short = 'm_fstop';
Conditions(6).name = 'motor finished accelerating';
Conditions(6).short = 'm_faccel';
Conditions(7).name = 'motor finished decelerating';
Conditions(7).short = 'm_fdecel';
Conditions(8).name = 'warning cue: starting';
Conditions(8).short = 'w_start';
Conditions(9).name = 'warning cue: stopping';
Conditions(9).short = 'w_stop';
Conditions(10).name = 'warning cue: accelerating';
Conditions(10).short = 'w_accel';
Conditions(11).name = 'warning cue: decelerating';
Conditions(11).short = 'w_decel';
Conditions(12).name = 'warning cue: maintaining';
Conditions(12).short = 'w_maint';
Conditions(13).name = 'warning cue: probe, no warning';
Conditions(13).short = 'w_p_nowarn';
Conditions(14).name = 'motor: probe, no warning, accelerating';
Conditions(14).short = 'm_p_nowarn_accel';
Conditions(15).name = 'motor: probe, no warning, decelerating';
Conditions(15).short = 'm_p_nowarn_decel';
Conditions(16).name = 'motor: probe, no warning, stopping';
Conditions(16).short = 'm_p_nowarn_stop';
Conditions(17).name = 'motor: probe, no warning, maintaining';
Conditions(17).short = 'm_p_nowarn_maint';
Conditions(18).name = 'warning cue: probe, starting cue, no change in motor.';
Conditions(18).short = 'w_p_nochange_start';
Conditions(19).name = 'warning cue: probe, stopping cue, no change in motor.';
Conditions(19).short = 'w_p_nochange_stop';
Conditions(20).name = 'warning cue: probe, accelerating cue, no change in motor.';
Conditions(20).short = 'w_p_nochange_accel';
Conditions(21).name = 'warning cue: probe, decelerating cue, no change in motor.';
Conditions(21).short = 'w_p_nochange_decel';
Conditions(22).name = 'warning cue: probe, maintaining cue, no change in motor.';
Conditions(22).short = 'w_p_nochange_maint';
Conditions(23).name = 'motor: probe, no change';
Conditions(23).short = 'm_p_nochange';

% Extras to sort.
Conditions(24).name = 'motor: starting';
Conditions(24).short = 'm_start'; 
Conditions(25).name = 'motor: finished starting';
Conditions(25).short = 'm_fstart'; 
Conditions(26).name = 'continued walking';
Conditions(26).short = 'c_walk'; 
Conditions(27).name = 'continued rest';
Conditions(27).short = 'c_rest'; 
Conditions(28).name = 'motor: probe, no warning, starting';
Conditions(28).short = 'm_p_nowarn_start'; 

save([dir_out 'Behavior_Conditions.mat'], 'Conditions');

%% FROM BEFORE RUNANALYSIS IMPLEMENTATION
% Create a periods structure
% speeds = [1600, 2000, 2400, 2800]; 
% accels_startstop = [400, 800];
% accels_acceldecel = [200, 800]; 
% 
% periods = cell(1,1); 
% for condi =1:size(Conditions,2)
%     short_name = Conditions(condi).short;
% 
%     % Only make fields for things we actually care about
%     % per condition type
% 
%     switch condi
% 
%         % For all but start & maintaining warning cues, the two 
%         % motor maintaining types, the continued walking
%         % types, and the no change in motor probe,we only 
%         % care about the current speed subdivisions
%         case num2cell([9:11 19:21 26])
% 
%             % Initialize empty fields 
%             for speedi = 1: size (speeds,2)
%                 new_name = {[short_name '.x' num2str(speeds(speedi))]};
%                 periods =[periods; new_name];
%             end 
% 
%         % For accel, decel, we care about the
%         % previous speed, current speed, and acceleration
%         
%         % For accel and probe accel, previous speed has to be lower than
%         % current speed 
%         case num2cell([1,14])
% 
%             % Initialize empty fields 
%             for speedi = 1: size (speeds,2)
%                 speed =speeds(speedi);
%                 for acceli = 1:size(accels_acceldecel,2)
%                     accel =accels_acceldecel(acceli);
% 
%                     for speed2i = 1: (speedi-1)
%                         speed2 =speeds(speed2i);
%                         
%                         if speed ~=speed2
%                             new_name = {[ short_name '.x' num2str(speed) '.x' num2str(accel) '.x' num2str(speed2)]};
%                             periods =[periods; new_name];
%                         end 
%                     end 
%                 end 
%             end
%             
%             % For decel and probe decel, previous speed has to be higher than
%             % current speed 
%             
%             case num2cell([2,15])
% 
%                 % Initialize empty fields 
%                 for speedi = 1: size (speeds,2)
%                     speed =speeds(speedi);
%                     for acceli = 1:size(accels_acceldecel,2)
%                         accel =accels_acceldecel(acceli);
% 
%                         for speed2i = speedi + 1 : size(speeds,2)
%                             speed2 =speeds(speed2i);
% 
%                             if speed ~=speed2
%                                 new_name = {[ short_name '.x' num2str(speed) '.x' num2str(accel) '.x' num2str(speed2)]};
%                                 periods =[periods; new_name];
%                             end 
%                         end 
%                     end 
%                 end
% 
%         % For finished accel & finished decel, we care about 
%         % current speed, 2 speeds ago, and acceleration.
%         
%         % faccel 's second speed can't be higher than first
%         case num2cell([6])   
%             % Initialize empty fields 
%             for speedi = 1: size (speeds,2)
%                 speed =speeds(speedi);
%                 for acceli = 1:size(accels_acceldecel,2)
%                     accel =accels_acceldecel(acceli);
% 
%                     for speed2i = 1:speedi-1
%                         speed2 =speeds(speed2i);
%                         
%                         if speed ~=speed2
%                             new_name = {[ short_name '.x' num2str(speed) '.x' num2str(accel) '.x' num2str(speed2)]};
%                             periods =[periods; new_name];
%                         end
%                     end 
%                 end 
%             end
%             
%         % fdecel 's second speed can't be lower than first
%         case num2cell(7)   
%             % Initialize empty fields 
%             for speedi = 1: size (speeds,2)
%                 speed =speeds(speedi);
%                 for acceli = 1:size(accels_acceldecel,2)
%                     accel =accels_acceldecel(acceli);
% 
%                     for speed2i = speedi+1:size(speeds,2)
%                         speed2 =speeds(speed2i);
%                         
%                         if speed ~=speed2
%                             new_name = {[ short_name '.x' num2str(speed) '.x' num2str(accel) '.x' num2str(speed2)]};
%                             periods =[periods; new_name];
%                         end
%                     end 
%                 end 
%             end    
% 
%         % For start, no warning start, & finished start, care only about current speed and accel
%         case num2cell([24, 25, 28])
%             % Initialize empty fields 
%             for speedi = 1: size (speeds,2)
%                 speed =speeds(speedi);
%                 for acceli = 1:size(accels_startstop, 2)
%                     accel =accels_startstop(acceli);
%                     new_name = {[ short_name '.x' num2str(speed) '.x' num2str(accel)]};
%                     periods =[periods; new_name];
% 
%                 end 
%             end
% 
%         % For stop & no warning stop, care only about
%         % previous speed and accel.
%         case num2cell([4, 16])
%             % Initialize empty fields 
%             all_speeds = [speeds];
%             for speedi = 1: size (all_speeds,2)
%                 speed = all_speeds(speedi);
%                 for acceli = 1:size(accels_startstop,2)
%                     accel =accels_startstop(acceli);
%                     new_name = {[ short_name '.x' num2str(speed) '.x' num2str(accel)]};
%                     periods =[periods; new_name];
%                 end 
%             end
% 
%         % For finished stop, care only about 2 speeds ago
%         % and accel. 
%         case 5 
%             % Initialize empty fields
%             all_speeds = [ speeds];
%             for speedi = 1: size (all_speeds,2)
%                 speed = all_speeds(speedi);
%                 for acceli = 1:size(accels_startstop, 2)
%                     accel =accels_startstop(acceli);
%                     new_name = {[ short_name '.x' num2str(speed) '.x' num2str(accel)]};
%                     periods =[periods; new_name];
%                 end 
%             end
% 
%         % For warning start, warning start probe, and continued rest don't
%         % care about anything
%         case num2cell([8,18, 27])
%             periods =[ periods; {short_name}];
% 
%         % Motor maintaining, warning maintaining, probe warning maintaining,and motor no change
%         % care about all speeds, including 0; and no
%         % warning
%         case num2cell([3,12,13,17,22, 23])   
% 
%             all_speeds = [0 speeds];
%             % add names
%             for speedi = 1: size (all_speeds,2)
%                 new_name = {[short_name '.x' num2str(all_speeds(speedi))]};
%                 periods =[periods; new_name];
%             end 
% 
%     end
% end
% % Remove the empty first entry
% periods = periods(2:end);
% 
% % Save
% save([dir_out 'periods_namestructure.mat'], 'periods');


%% Make an iterations structure for each period type. 
% Create a periods structure (make them all strings because they're names
% of speeds/accels).
speeds = {'0', '1600', '2000', '2400', '2800'}; 
accels_startstop = {'400', '800'};
accels_acceldecel = {'200', '800'}; 

% Load conditions structure.
load([dir_exper 'Behavior_Conditions.mat']);
 
for condi =1:size(Conditions,2)

    short_name = Conditions(condi).short;
    periods(condi).name = short_name;

    % Only make fields for things we actually care about
    % per condition type

    switch condi

        % For all but start & maintaining warning cues, the two 
        % motor maintaining types, the continued walking
        % types, and the no change in motor probe,we only care about the 
        % current speed subdivisions.  Motor maintaining, warning maintaining, probe warning maintaining,and motor no change
        % care about speeds, including 0; and no warning
        % Warning for stopping, accel, decel
        % Probes for those,
        % Continued walking
        % Motor maintaining
        % Maintaining probes
        % Motor no change probe.
        case num2cell([9:11 19:21 26 3,12,13,17,22, 23])
             periods(condi).speed = speeds;
             periods(condi).accel = NaN; 
             periods(condi).previous_speed = NaN;
             periods(condi).two_speeds_ago = NaN; 

        % For accel, decel, & corresponding probes, we care about the
        % previous speed (no 0), current speed (no 0), and acceleration
        case num2cell([1, 2, 14, 15])
            
            periods(condi).speed = speeds(2:end);
            periods(condi).accel = accels_acceldecel; 
            periods(condi).previous_speed = speeds(2:end);
            periods(condi).two_speeds_ago = NaN; 
            

        % For finished accel & finished decel, we care about 
        % current speed (no 0), 2 speeds ago (no 00, and acceleration.
        case num2cell([6, 7]) 
            periods(condi).speed = speeds(2:end);
            periods(condi).accel = accels_acceldecel; 
            periods(condi).previous_speed = NaN;
            periods(condi).two_speeds_ago = speeds(2:end); 

        % For start, no warning start, & finished start, care only about current speed (no 0)and accel
        case num2cell([24, 25, 28])
            periods(condi).speed = speeds(2:end);
            periods(condi).accel = accels_startstop; 
            periods(condi).previous_speed = NaN; 
            periods(condi).two_speeds_ago = NaN; 

        % For stop & no warning stop, care only about
        % previous speed (no 0) and accel.
        case num2cell([4, 16])
            periods(condi).speed = NaN;
            periods(condi).accel = accels_startstop; 
            periods(condi).previous_speed = speeds(2:end);
            periods(condi).two_speeds_ago = NaN; 
           
        % For finished stop, care only about 2 speeds ago (no 0).
        % and accel. 
        case 5 
            periods(condi).speed = NaN;
            periods(condi).accel = accels_startstop; 
            periods(condi).previous_speed = NaN;
            periods(condi).two_speeds_ago = speeds(2:end); 

        % For warning start, warning start probe, and continued rest don't
        % care about anything (there's only one possible value of each).
        case num2cell([8,18, 27])
            periods(condi).speed = NaN;
            periods(condi).accel = NaN; 
            periods(condi).previous_speed = NaN;
            periods(condi).two_speeds_ago = NaN; 

    end
end

% Save
save([dir_out 'periods.mat'], 'periods');

%% Now make something that works better with code

load([dir_out 'periods.mat']);

loop_list.iterators = {
               'condition', {'loop_variables.periods(:).name'}, 'condition_iterator';
               'speed', {'loop_variables.periods(', 'condition_iterator', ').speed'}, 'speed_iterator';
               'accel', {'loop_variables.periods(', 'condition_iterator', ').accel'}, 'accel_iterator';
               'previous_speed', {'loop_variables.periods(', 'condition_iterator', ').previous_speed'}, 'previous_speed_iterator';
               'two_speeds_ago', {'loop_variables.periods(', 'condition_iterator', ').two_speeds_ago'}, 'two_speeds_ago_iterator'
               };

 loop_variables.periods = periods; 

 looping_output_list = LoopGenerator(loop_list, loop_variables);

holder = NaN(size(looping_output_list,1), size(loop_list.iterators,1));
for iteratori = 1:size(loop_list.iterators,1)

    eval(['holder_sub = [looping_output_list.' loop_list.iterators{iteratori,3} '];']);
    holder(:, iteratori) = holder_sub;
    
end     

% Keep only unique entries (may not be necessary if LoopGenerator
% worked correctly
[~, ia, ~] = unique(holder, 'rows','stable');
looping_output_list = looping_output_list(ia);
periods = looping_output_list;

% Remove unecessary "iterator" fields.
fields = {'condition_iterator', 'speed_iterator', 'accel_iterator', 'previous_speed_iterator', 'two_speeds_ago_iterator'};
periods = rmfield(periods, fields);

% Make into table.
periods = struct2table(periods);

indices_to_remove = []; 
% Add a column for duration, in frames;
for periodi = 1:size(periods,1)

    % For all finished periods & maintaining
    if strcmp(periods.condition{periodi}(1:3), 'm_f') ...
            || strcmp(periods.condition{periodi}, 'm_maint') ...
            || strcmp(periods.condition{periodi}, 'm_p_nochange') ...
            || strcmp(periods.condition{periodi}, 'm_p_nowarn_maint')
        periods.duration{periodi} = finished_chunk_length * fps;  
        

        % Still calculate the duration for removing not-real ones, though
        final = str2double(periods.speed{periodi});
        first = str2double(periods.two_speeds_ago{periodi}); 
        accel = str2double(periods.accel{periodi}); 

        duration = (final - first)/accel * fps;

        
        % If this is an accel or its probe & the duration is negative, you
        % can remove it because it doesn't actually exist.
        if duration < 0 & strcmp(periods.condition{periodi}, 'm_faccel') 
            indices_to_remove = [indices_to_remove; periodi];
        
        elseif duration > 0 & strcmp(periods.condition{periodi}, 'm_fdecel') 
            indices_to_remove = [indices_to_remove; periodi];
        
        % If duration is 0, remove.
        elseif duration == 0
            
            indices_to_remove = [indices_to_remove; periodi];
       
        end 

    % For continued periods. 
    elseif strcmp(periods.condition{periodi}, 'c_walk') || strcmp(periods.condition{periodi}, 'c_rest')
           

       periods.duration{periodi} = continued_chunk_length * fps; 
    
    % For warning periods
    elseif strcmp(periods.condition{periodi}(1:2), 'w_')
        periods.duration{periodi} = warning_length * fps ;

    % For starts, stops, accels, decels, & their probes
    else

        % Grab current & previous speeds. If either are "NaN", set to 0.
        if strcmp(periods.speed{periodi}, "NaN")
            final = 0; 
        else 
            final = str2double(periods.speed{periodi});  
        end

        if strcmp(periods.previous_speed{periodi}, "NaN")
            first = 0; 
        else 
            first = str2double(periods.previous_speed{periodi});  
        end
        
        accel = str2double(periods.accel{periodi});

        duration = (final - first)/accel * fps; 

        % If this is an accel or its probe & the duration is negative, you
        % can remove it because it doesn't actually exist.
        if duration < 0 & (strcmp(periods.condition{periodi}, 'm_accel') || strcmp(periods.condition{periodi}, 'm_p_nowarn_accel')) 
            indices_to_remove = [indices_to_remove; periodi];
        
        % If this is an decel or its probe & the duration is positive, you
        % can remove it because it doesn't actually exist.
        elseif duration > 0 & (strcmp(periods.condition{periodi}, 'm_decel') || strcmp(periods.condition{periodi}, 'm_p_nowarn_decel')) 
            indices_to_remove = [indices_to_remove; periodi];
  
        % If duration is 0, remove.
        elseif duration == 0
            indices_to_remove = [indices_to_remove; periodi];
        end 

        % Now assign the absolute value
        periods.duration{periodi} = abs(duration);

    end

end 

% Remove periods (they don't actually exist)
periods(indices_to_remove, :) = [];

% Add indices to the end, for searching later.
index = [1:size(periods,1)]';
a = table(index);
periods = [periods a];

% Save
save([dir_out 'periods_nametable.mat'], 'periods');


% Can search tables like this: 
% g = periods_table(string(periods_table.condition)=="m_accel" & string(periods_table.speed)=="2000", :);
%% Create structure of groupings of periods that are similar.--For timeseries of different durations

% clear all;

% % Directories
% experiment_name='Random Motorized Treadmill';
% dir_base='Y:\Sarah\Analysis\Experiments\';
% dir_exper=[dir_base experiment_name '\']; 
% dir_out=dir_exper; 
% mkdir(dir_out);
% 
% % Possible speeds and accels. 
% speeds = [1600, 2000, 2400, 2800]; 
% accels_startstop = [400, 800];
% accels_acceldecel = [200, 800]; 
% 
% periods_grouped_varyingduration = cell(1,1);
% 
% % Accelerations -- same magnitude, accel 
% 
% % For each possible acceleration.
% for acceli = 1:numel(accels_acceldecel) 
%     accel = accels_acceldecel(acceli);
%     
%     % For each possible magnitude
%     for magi = 1:numel(speeds)-1
%     
%         % Find the magnitude of the change.
%         magnitude = magi * 400;
%    
%         new_name = {[ 'm_accel.x' num2str(accel) '.' num2str(magnitude)]};
%            
%         periods_grouped_varyingduration = [periods_grouped_varyingduration; new_name];
%         
%     end
% end 
% 
% % Decelerations -- same magnitude, accel
% 
% % For each possible acceleration.
% for acceli = 1:numel(accels_acceldecel) 
%     accel = accels_acceldecel(acceli);
%     
%     % For each possible magnitude
%     for magi = 1:numel(speeds)-1
%     
%         % Find the magnitude of the change.
%         magnitude = magi * 400;
%    
%         new_name = {[ 'm_decel.x' num2str(accel) '.' num2str(magnitude)]};
%            
%         periods_grouped_varyingduration = [periods_grouped_varyingduration; new_name];
%         
%     end
% end 
% 
% % Starts
% 
% % Stops
% 
% % Finish of accels, same accel
% for acceli = 1:numel(accels_acceldecel) 
%     accel = accels_acceldecel(acceli);
% 
%     new_name = {[ 'm_faccel.x' num2str(accel)]};
%     periods_grouped_varyingduration = [periods_grouped_varyingduration; new_name];
% end 
% 
% % Finish of decels, same accel
% for acceli = 1:numel(accels_acceldecel) 
%     accel = accels_acceldecel(acceli);
%    
%     new_name = {[ 'm_fdecel.x' num2str(accel)]};
%     periods_grouped_varyingduration = [periods_grouped_varyingduration; new_name];
% end 
% 
% % Finish of starts, same accels, all speeds
% for acceli = 1:numel(accels_startstop) 
%     accel = accels_acceldecel(acceli);
%    
%     new_name = {[ 'm_fstart.x' num2str(accel)]};
%     periods_grouped_varyingduration = [periods_grouped_varyingduration; new_name];   
% end 
% 
% % Finish of stops, same accels, all speeds.
% for acceli = 1:numel(accels_startstop) 
%     accel = accels_acceldecel(acceli);
%     new_name = {[ 'm_fstop.x' num2str(accel)]};
%     periods_grouped_varyingduration = [periods_grouped_varyingduration; new_name];
% end 
% 
% % Maintaining -- rest and all speeds.
% periods_grouped_varyingduration = [periods_grouped_varyingduration; 'm_maint.rest'];
% periods_grouped_varyingduration = [periods_grouped_varyingduration; 'm_maint.walk'];
% 
% % Warning of accels -- all speeds
% periods_grouped_varyingduration = [periods_grouped_varyingduration; 'w_accel'];
% 
% % Warning of decels -- all speeds
% periods_grouped_varyingduration = [periods_grouped_varyingduration; 'w_decel'];
% 
% % Warning of stops -- all speeds
% periods_grouped_varyingduration = [periods_grouped_varyingduration; 'w_stop'];
% 
% % Warning of maintaining -- rest and all speeds
% periods_grouped_varyingduration = [periods_grouped_varyingduration; 'w_maint.rest'];
% periods_grouped_varyingduration = [periods_grouped_varyingduration; 'w_maint.walk'];
% 
% % All continued walk (all speeds)
% periods_grouped_varyingduration = [periods_grouped_varyingduration; 'c_walk'];
% 
% % Continued rest. 
% periods_grouped_varyingduration = [periods_grouped_varyingduration; 'c_rest'];
% 
% % Probes.
% 
% % Save

%% Create structure of groupings of periods that are similar.--For summary correlations or other across-duration summarys

%clear all;
% 
% % Directories
% dir_exper=[dir_base experiment_name '\']; 
% dir_out=dir_exper; 
% mkdir(dir_out);
% 
% % Possible speeds and accels. 
% speeds = [1600, 2000, 2400, 2800]; 
% accels_startstop = [400, 800];
% accels_acceldecel = [200, 800]; 
% 
% % Initialize
% periods_grouped_allduration = cell(1,1);
% 
% % Accelerations
% % For each possible acceleration .
% periods_grouped_allduration = [periods_grouped_allduration; 'm_accel'];
% 
% % Decelerations 
% periods_grouped_allduration = [periods_grouped_allduration; 'm_decel'];
% 
% % Starts
% periods_grouped_allduration = [periods_grouped_allduration; 'm_start'];
% 
% % Stops
% periods_grouped_allduration = [periods_grouped_allduration; 'm_stop'];
% 
% % Finish of accel
% periods_grouped_allduration = [periods_grouped_allduration; 'm_faccel'];
% 
% % Finish of decels
% periods_grouped_allduration = [periods_grouped_allduration; 'm_faccel'];
% 
% % Finish of starts.
% periods_grouped_allduration = [periods_grouped_allduration; 'm_fstart'];
% 
% % Finish of stops 
% periods_grouped_allduration = [periods_grouped_allduration; 'm_fstop'];
% 
% % Maintaining -- rest and all speeds.
% periods_grouped_allduration = [periods_grouped_allduration; 'm_maint.rest'];
% periods_grouped_allduration = [periods_grouped_allduration; 'm_maint.walk'];
% 
% % Warning of accels
% periods_grouped_allduration = [periods_grouped_allduration; 'w_accel'];
% 
% % Warning of decels
% periods_grouped_allduration = [periods_grouped_allduration; 'w_decel'];
% 
% % Warning of stops
% periods_grouped_allduration = [periods_grouped_allduration; 'w_stop'];
% 
% % Warning of maintaining 
% periods_grouped_allduration = [periods_grouped_allduration; 'w_maint.rest'];
% periods_grouped_allduration = [periods_grouped_allduration; 'w_maint.walk'];
% 
% % All continued walk (all speeds)
% periods_grouped_allduration = [periods_grouped_allduration; 'c_walk'];
% 
% % Continued rest. 
% periods_grouped_allduration = [periods_grouped_allduration; 'c_rest'];
% 
% % Probes.
% 
% % Save
% save([dir_out 'periods_grouped_allduration.mat'], 'periods_grouped_allduration');
