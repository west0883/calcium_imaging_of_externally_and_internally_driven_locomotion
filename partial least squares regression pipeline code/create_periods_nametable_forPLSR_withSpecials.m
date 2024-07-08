% create_periods_nametable_forPLSR.m
% Sarah West
% 6/4/22

% Script that organizes periods_nametable for easy creation of response
% variables for partial least squares regression for Random Motorized
% Treadmill experiment.

clear all; 

% Create the experiment name.
parameters.experiment_name='Random Motorized Treadmill';

% Output directory name bases
parameters.dir_base='Y:\Sarah\Analysis\Experiments\';
parameters.dir_exper=[parameters.dir_base parameters.experiment_name '\']; 

% Load mice_all, pass into parameters structure
load([parameters.dir_exper '\mice_all.mat']);
parameters.mice_all = mice_all;

% ****Change here if there are specific mice, days, and/or stacks you want to work with**** 
parameters.mice_all = parameters.mice_all;

% Load names of motorized periods
load([parameters.dir_exper 'periods_nametable.mat']);
periods_motorized = periods;

% Load names of spontaneous periods
load([parameters.dir_exper 'periods_nametable_spontaneous.mat']);
periods_spontaneous = periods(1:6, :);
clear periods; 

% Create a shared motorized & spontaneous list to edit as your new
% nametable.
periods = [periods_motorized; periods_spontaneous]; 

fps = 20;
window_size = 1;
window_step_size = 1; 

% Make a partial least squares regression folder. 
if ~isdir([parameters.dir_exper 'PLSR\'])
    mkdir([parameters.dir_exper 'PLSR\']);
end

%% Make a list of dummy variable categories for using later. 
% What to put the creation of this here in this script to keep it with the
% other organizing period_nametable code.
% No prep for now.
categories.motorized_vs_spon = {'motorized', 'spontaneous'};
categories.type = {'rest', 'walk', 'start', 'stop', 'accel', 'decel', 'finished_start','finished_stop', 'finished_accel', 'finished_decel', ...
    'w_start', 'w_maint_rest', 'prewalk'}; % 'prep'

save([parameters.dir_exper 'PLSR\response_categories.mat'], 'categories');

%% Make values for speed, accel into numerics. 
speeds = periods{:,'speed'};
periods.speed = cellfun(@str2num, speeds, 'UniformOutput', false); 

accels = periods{:,'accel'};
periods.accel = cellfun(@str2num, accels, 'UniformOutput', false); 

clear speeds accels; 

%% Make all speeds & accels into cm/s 
% Only values in here so far are motorized, so can apply to all speed
% values.
original_speeds = [1600, 2000, 2400, 2800];
cms_speeds = [2.25, 2.77, 3.33, 3.88];
conversion_factor = round(mean(original_speeds./cms_speeds),-1);

speeds = periods{:,'speed'};
periods.speed = cellfun(@(x) x/conversion_factor , speeds, 'UniformOutput', false); 

accels = periods{:,'accel'};
periods.accel = cellfun(@(x) x/conversion_factor, accels, 'UniformOutput', false); 

%% Remove behavior fields you don't need.
% previous speed, two speeds ago, previous accel.
columns_to_remove = {'previous_speed', 'two_speeds_ago'}; % 'previous_accel'
for coli = 1:numel(columns_to_remove)
    periods.(columns_to_remove{coli}) = [];
end

clear columns_to_remove; 

%% Create labels for "spontaneous" vs "motorized"
% Do this before removing periods you don't need so you can use size of
% periods_motorized and periods_spontaneous.
motorized_vs_spon = cell(size(periods,1),1);
motorized_vs_spon([1:size(periods_motorized,1)]) = repmat({'motorized'}, size(periods_motorized,1) , 1);
motorized_vs_spon([size(periods_motorized,1) + 1:end]) = repmat({'spontaneous'}, size(periods_spontaneous,1) , 1);

periods.motorized_vs_spon = motorized_vs_spon;
clear motorized_vs_spon;

%% Keep only the set of periods you want to use for PLS regression 
% For first-pass of PLSR, don't use probe trials, maintaining, or full onset/offset. 
% Save the indices you removed, too, so you can easily remove them from
% correlation data variables later. 
% DON'T re-number the index field--> that would make it very hard to find
% the corresponding velocities and accels later. 
% No prep periods for now.
conditions_to_remove = {'m_maint', 'w_decel', 'w_accel', 'm_p', 'w_p', 'w_stop', 'full_onset', 'full_offset'};

% Also include some of the meaningless ones.
indices_to_remove = [71; 76; 80; 81; 82; 175]; 

for i = 1:size(periods,1)

    if contains(cell2mat(periods{i, 'condition'}), conditions_to_remove)
        indices_to_remove = [indices_to_remove; i];
    end
end 

periods(indices_to_remove, :) = [];
save([parameters.dir_exper 'PLSR\indices_to_remove_Specials.mat'], 'indices_to_remove');

%% Label by "super-type"
% start, stop, accel, decel, rest, walk, finished, or finished stop. For 
% spontaneous: prewalk = prep, startwalk = start, stopwalk = stop, postwalk 
% = finished.
% No prep for now.
look_for = {'m_start','startwalk', 'm_stop', 'stopwalk', 'm_accel', 'm_decel',...
          'c_rest', 'rest', 'c_walk', 'walk', 'postwalk', 'm_fstop', 'm_fstart', 'm_faccel', 'm_fdecel',...
          'prewalk', 'w_start', 'w_maint'}; 
           
types = cell(size(periods,1), 1);

for i = 1:size(periods, 1)
   
    this_condition = repmat(periods{i, 'condition'}, 1, numel(look_for));
    type = find(cellfun(@strcmp, this_condition, look_for));

    if isempty(type)
        disp('empty')
    end 

    switch type

        case {1, 2}   % start
            types{i} = 'start';
    
        case {3, 4} % stop
            types{i} = 'stop';

        case 5 % accel
             types{i} = 'accel';

        case 6 % decel
            types{i} = 'decel';     

        case {7, 8} 
            types{i} = 'rest';

        case {9, 10} % walk
            types{i} = 'walk';

        case {11,12}  % postwalk, m_fstop --> finished_stop
            types{i} = 'finished_stop'; 

        case 13  % finished start
            types{i} = 'finished_start'; 

        case 14 % m_faccel
            types{i} = 'finished_accel';

        case 15 % m_fdecel
            types{i} = 'finished_decel';

        case 16 % prewalk
            types{i} = 'prewalk';

        case 17 % w_start
            types{i} = 'w_start';

        case 18 % the w_maint s 

            % Separate into rest w_maint and walking w_maintain
            if periods{i, 'speed'}{1} == 0
                types{i} = 'w_maint_rest';
            else 
                types{i} = 'w_maint_walk';
            end 

    end
end
periods.type = types; 
clear types type type1 type2 look_for_pattern look_for this_condition;

%% Make column for pupil diameter, tail, nose, FL, HL, x
% All periods. Are all NaN for now, will be put in at "populate response
% variables" step in main pipline.
extra_variables = {'pupil_diameter',  'tail', 'nose', 'FL', 'HL', 'x'};
for i = 1:numel(extra_variables)
    holder = repmat({NaN}, size(periods,1),1);
    periods.(extra_variables{i}) = holder;
end 

clear holder;

%% Put in accel = 0 for all finisheds, walk, & rest periods
% Leave spontaneous blank (Nan) for now
% No prep periods for now.
for i = 1:size(periods,1)

    if strcmp(periods{i, 'motorized_vs_spon'}, 'motorized') 
        
       if strcmp(periods{i, 'type'}, 'finished_start') || strcmp(periods{i, 'type'}, 'finished_stop') ...
               || strcmp(periods{i, 'type'}, 'finished_accel') || strcmp(periods{i, 'type'}, 'finished_decel') ...
               || strcmp(periods{i, 'type'}, 'walk') || strcmp(periods{i, 'type'}, 'rest')  % strcmp(periods{i, 'type'}, 'prep') ||
       
           periods{i, 'accel'} = {0} ;      
       end
    end
end

clear i;
%% Make speed = 0 for motorized stop, finished stop, warning start, c_rest
% Leave spontaneous blank (Nan) for now. No prep periods for now.
look_for = {'m_stop', 'm_fstop', 'c_rest', 'w_start', 'w_maint_rest', 'w_maint_walk'};
for i = 1:size(periods,1)
    this_condition = repmat(periods{i, 'condition'}, 1, numel(look_for));

    change = any(cellfun(@strcmp, this_condition, look_for));

    if change
           periods{i, 'speed'} = {0} ;    
    end
end

clear look_for i this_condition change;

%% Calculate new "roll numbers" for periods & duration you have left. 
% Figure out how many rolls you can get out of it (should be a whole number). 
number_of_rolls = cell(size(periods,1), 1);
for i = 1:size(periods,1)
   number_of_rolls{i} = CountRolls(periods{i, 'duration'}{1}, window_size, window_step_size);
end
periods.number_of_rolls = number_of_rolls;

%% Create "instantaneous" duration vectors 
% Based on roll numbers for all periods.
% The time of a given roll is the center of the roll window.
duration_vector = cellfun(@(x) ([1:x] + 1) * (window_step_size/fps) , number_of_rolls, 'UniformOutput', false);
periods.duration_vector = duration_vector;

%% Convert the rest & walk periods' duration vectors to NaN
% Will be replaced by duration_place vectors later

for i = 1:size(periods,1)
    if strcmp(periods{i, 'type'}{1}, 'walk') || strcmp(periods{i, 'type'}{1}, 'rest')
        periods{i, 'duration_vector'} = {NaN}; 
    end
end 

%% For motorized transitions, calculate "instantaneous" speed 
% For certain periods, based on speed, previous speed, accel. (or duration,
% roll number, speed, accel)/ Leave spontaneous blank for now. 
 speed_vector = cell(size(periods,1),1);

% The current speeds are where the motor is GOING to. 
for i = 1:size(periods,1)
    if strcmp(periods{i, 'type'}{1}, 'start') || strcmp(periods{i, 'type'}{1}, 'accel') 

        % The time is the center of the roll window
        duration_vector = periods{i, "duration_vector"}{1};
        speed_vector{i} = fliplr(periods{i, 'speed'}{1} - duration_vector * periods{i, 'accel'}{1});


    elseif strcmp(periods{i, 'type'}{1}, 'stop') || strcmp(periods{i, 'type'}{1}, 'decel')
        % The time is the center of the roll window
        duration_vector = periods{i, "duration_vector"}{1};
        speed_vector{i} = fliplr(periods{i, 'speed'}{1} + duration_vector * periods{i, 'accel'}{1});

    else
        % If not one of those special transitions, replicate by
        % roll_number.
        speed_vector(i) = {repmat(periods{i, 'speed'}{1}, 1, periods{i,'number_of_rolls'}{1})};
        
    end
end 
periods.speed_vector = speed_vector;

%% Replicate accels & pupil diameter & other variables by roll number.

variables = [{'accel'} extra_variables];

% set up holders
for ii = 1:numel(variables)
    holder.([variables{ii} '_vector']) =  cell(size(periods,1),1);
end 

for i = 1:size(periods,1)
    for ii = 1:numel(variables)

        holder.([variables{ii} '_vector'])(i) =  {repmat(periods{i, variables{ii}}{1}, 1, periods{i,'number_of_rolls'}{1})};
    end 
end 
for ii = 1:numel(variables)
    periods.([variables{ii} '_vector'])= holder.([variables{ii} '_vector']);
end

%% Make a column that says if this period is one of the "active" transitions 
% Is for later when you want to compare all transitions to other things. 
transition_or_not = cell(size(periods, 1), 1); 
for i = 1:size(periods,1)

    if strcmp(periods{i, 'type'}{1}, 'start') || strcmp(periods{i, 'type'}{1}, 'accel') || strcmp(periods{i, 'type'}{1}, 'decel') || strcmp(periods{i, 'type'}{1}, 'stop')
        transition_or_not{i} = 1;
    else 
        transition_or_not{i} = 0;
    end
    
end
periods.transition_or_not = transition_or_not; 


%% Create dummy variables for categoricals -- motorized_vs_spon, type, & transition or not
%categories.motorized_vs_spon = {'motorized', 'spontaneous'};
%categories.type = {'rest', 'walk', 'prep', 'start', 'stop', 'accel', 'decel', 'finished'};
motorized_vs_spon_dummyvars = cell(size(periods,1),1);
type_dummyvars = cell(size(periods,1),1);
transition_or_not_dummyvars = cell(size(periods,1),1);

for i = 1:size(periods,1)

    % Motorized vs spon
    motorized_vs_spon =  strcmp(repmat(periods{i, 'motorized_vs_spon'}, size(categories.motorized_vs_spon)), categories.motorized_vs_spon);
    motorized_vs_spon_dummyvars{i} = double(motorized_vs_spon)';

    % Type 
    type = strcmp(repmat(periods{i, 'type'}, size(categories.type)), categories.type);
    type_dummyvars{i} = double(type)';

    % Transition or not. 
    if periods{i, 'transition_or_not'}{1}
        transition_or_not_dummyvars{i} = [1; 0];
    else
       transition_or_not_dummyvars{i} = [0; 1];
    end

end 
periods.motorized_vs_spon_dummyvars = motorized_vs_spon_dummyvars;
periods.type_dummyvars = type_dummyvars;
periods.transition_or_not_dummyvars = transition_or_not_dummyvars;

%% Replicate dummy vars by roll number. 
motorized_vs_spon_dummyvars_vector = cell(size(periods,1),1);
type_dummyvars_vector = cell(size(periods,1),1);
transition_or_not_dummyvars_vector = cell(size(periods,1),1);
for i = 1:size(periods,1)

    % Motorized vs spon
    motorized_vs_spon_dummyvars_vector{i} =  repmat(periods{i, 'motorized_vs_spon_dummyvars'}{1}, [1 number_of_rolls{i}]);

    % Type 
    type_dummyvars_vector{i} = repmat(periods{i, 'type_dummyvars'}{1}, [1 number_of_rolls{i}]);

    % Transition or not 
    transition_or_not_dummyvars_vector{i} = repmat(periods{i, 'transition_or_not_dummyvars'}{1}, [1 number_of_rolls{i}]);

end 
periods.motorized_vs_spon_dummyvars_vector = motorized_vs_spon_dummyvars_vector;
periods.type_dummyvars_vector = type_dummyvars_vector;
periods.transition_or_not_dummyvars_vector = transition_or_not_dummyvars_vector;


%% Save 
save([parameters.dir_exper 'PLSR\periods_nametable_forPLSR_specials.mat'], 'periods', '-v7.3');

%clear all;