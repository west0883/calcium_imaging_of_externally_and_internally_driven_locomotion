% create_PLSR_comparisons_level1_categorical.m
% Sarah West
% 6/9/22

% A script that creates a structure of information to run comparisons
% between different subsets of data for the Random Motorized Treadmill
% experiments. Each comparison will be run per mouse, eventually will have
% a second-level comparison across mice.

% (This is for first-level PLSR comparisons);

% 3 main pieces for each comparison:
% name --> the working name of the comparison, for file names.
% variablesToUse --> the response variable categories you're interested in
% indices --> the period indices from periods_nametable that are relevant
% to this comparison & will find the related brain data & response
% variable pairs. Is important for when comparison dataset is first
% removecd from larger dataset. 
% Will also have a "figure_type" identifier, which is for determining the
% color range of final figures. Can be "continued", "startstop", or
% "acceldecel".

%% Initial setup
clear all; 

% Create the experiment name.
parameters.experiment_name='Random Motorized Treadmill';

% Output directory name basis
parameters.dir_base='Y:\Sarah\Analysis\Experiments\';
parameters.dir_exper=[parameters.dir_base parameters.experiment_name '\']; 

% Load the periods_nametable for this experiment/analysis.
% (For PLSR, was created with create_PLSR_comparisons.m).
load([parameters.dir_exper 'PLSR\periods_nametable_forPLSR_specials.mat'], 'periods');

% Continuous variable names
continuous_variable_names = {'speed_vector', 'accel_vector', 'duration_vector'};

% Pull out the relevant columns/info from periods
period_motorized_vs_spon = periods.motorized_vs_spon;
period_types = periods.type;
indices_motorized = strcmp(period_motorized_vs_spon, 'motorized');
indices_spontaneous = strcmp(period_motorized_vs_spon, 'spontaneous');

% Output filename
filename_out = [parameters.dir_exper 'PLSR\comparisons_special_categorical.mat'];
% can search tables like this:
% g = periods_table(string(periods_table.condition)=="m_accel" & string(periods_table.speed)=="2000", :);

counter = 1;


%% Rests, motorized vs spon.
comparisons(counter).name = 'rest_motorizedvsspon_categorical';
comparisons(counter).variablesToUse = {'motorized_vs_spon_dummyvars_vector'};
comparisons(counter).mice_not_to_use = {'1100'};
comparisons(counter).plotMultiplier = 1;
comparisons(counter).figure_type = 'continued';
comparisons(counter).fromContinuous = 'all_rest_continuousVars';

% Get relevent indices for this type.
indices_type = strcmp(period_types, 'rest');

comparisons(counter).indices = find(indices_type);

counter = counter + 1;

%% Walks, motorized vs spon.
comparisons(counter).name = 'walk_motorizedvsspon_categorical';
comparisons(counter).variablesToUse = {'motorized_vs_spon_dummyvars_vector'};
comparisons(counter).mice_not_to_use = {'1100'};
comparisons(counter).plotMultiplier = 1;
comparisons(counter).figure_type = 'continued';
comparisons(counter).fromContinuous = 'all_walk_continuousVars';

% Get relevent indices for this type.
indices_type = strcmp(period_types, 'walk');

comparisons(counter).indices = find(indices_type);

counter = counter + 1;

%% Starts, motorized vs spon
comparisons(counter).name = 'start_motorizedvsspon_categorical';
comparisons(counter).variablesToUse = {'motorized_vs_spon_dummyvars_vector'};
comparisons(counter).mice_not_to_use = {'1100'};
comparisons(counter).plotMultiplier = 1;
comparisons(counter).figure_type = 'startstop';
comparisons(counter).fromContinuous = 'all_start_continuousVars';

% Get relevent indices for this type.
indices_type = strcmp(period_types, 'start');

comparisons(counter).indices = find(indices_type);

counter = counter + 1;

%% Stops, motorized vs spon
comparisons(counter).name = 'stop_motorizedvsspon_categorical';
comparisons(counter).variablesToUse = {'motorized_vs_spon_dummyvars_vector'};
comparisons(counter).mice_not_to_use = {'1100'};
comparisons(counter).plotMultiplier = 1;
comparisons(counter).figure_type = 'startstop';
comparisons(counter).fromContinuous = 'all_stop_continuousVars';

% Get relevent indices for this type.
indices_type = strcmp(period_types, 'stop');

comparisons(counter).indices = find(indices_type);

counter = counter + 1;

%% Finished stops, motorized vs spon
comparisons(counter).name = 'finished_stop_motorizedvsspon_categorical';
comparisons(counter).variablesToUse = {'motorized_vs_spon_dummyvars_vector'};
comparisons(counter).mice_not_to_use = {'1100'};
comparisons(counter).plotMultiplier = 1;
comparisons(counter).figure_type = 'startstop';
comparisons(counter).fromContinuous = 'all_finished_stop_continuousVars';

% Get relevent indices for this type.
indices_type = strcmp(period_types, 'finished_stop');

comparisons(counter).indices = find(indices_type);

counter = counter + 1;

%% Prewalk vs wstart
comparisons(counter).name = 'prewalkvswstart';
comparisons(counter).variablesToUse = {'motorized_vs_spon_dummyvars_vector'};
comparisons(counter).mice_not_to_use = {'1100'};
comparisons(counter).plotMultiplier = 1;
comparisons(counter).figure_type = 'warningPeriods';
comparisons(counter).fromContinuous = 'pre_and_post';

% Get relevent indices for this type.
indices_type = strcmp(period_types, 'prewalk') | strcmp(period_types, 'w_start');

comparisons(counter).indices = find(indices_type);

counter = counter + 1;

%% Rest vs walk 
% Motorized
comparisons(counter).name = 'motorized_restvswalk';
comparisons(counter).variablesToUse = {'type_dummyvars_vector'};
comparisons(counter).mice_not_to_use = {'1100'};
comparisons(counter).plotMultiplier = -1;
comparisons(counter).figure_type = 'continued';

% Get relevent indices for this type.
indices_type = strcmp(period_types, 'walk') | strcmp(period_types, 'rest');

% Get intersection of motorized & type. 
comparisons(counter).indices = find(indices_motorized & indices_type);

counter = counter + 1; 

% Spontaneous
comparisons(counter).name = 'spontaneous_restvswalk';
comparisons(counter).variablesToUse = {'type_dummyvars_vector'};
comparisons(counter).mice_not_to_use = {'1100'};
comparisons(counter).plotMultiplier = -1;
comparisons(counter).figure_type = 'continued';

% Get relevent indices for this type.
indices_type = strcmp(period_types, 'walk') | strcmp(period_types, 'rest');

% Get intersection of spontaneous & type. 
comparisons(counter).indices = find(indices_spontaneous & indices_type);

counter = counter + 1; 

%% rest vs start 
% motorized
comparisons(counter).name = 'motorized_restvsstart_categorical';
comparisons(counter).variablesToUse = {'type_dummyvars_vector'};
comparisons(counter).mice_not_to_use = {'1100'};
comparisons(counter).plotMultiplier = -1;
comparisons(counter).figure_type = 'startstop';

% Get relevent indices for this type.
indices_type = strcmp(period_types, 'rest') | strcmp(period_types, 'start');

% Get intersection of motorized & type. 
comparisons(counter).indices = find(indices_motorized & indices_type);

counter = counter + 1;

% spontaneous
comparisons(counter).name = 'spontaneous_restvsstart_categorical';
comparisons(counter).variablesToUse = {'type_dummyvars_vector'};
comparisons(counter).mice_not_to_use = {'1100'};
comparisons(counter).plotMultiplier = -1;
comparisons(counter).figure_type = 'startstop';

% Get relevent indices for this type.
indices_type = strcmp(period_types, 'rest') | strcmp(period_types, 'start');

% Get intersection of spontaneous & type. 
comparisons(counter).indices = find(indices_spontaneous & indices_type);

counter = counter + 1; 

%% rest vs stop 
% motorized
comparisons(counter).name = 'motorized_restvsstop_categorical';
comparisons(counter).variablesToUse = {'type_dummyvars_vector'};
comparisons(counter).mice_not_to_use = {'1100'};
comparisons(counter).plotMultiplier = -1;
comparisons(counter).figure_type = 'startstop';

% Get relevent indices for this type.
indices_type = strcmp(period_types, 'rest') | strcmp(period_types, 'stop');

% Get intersection of motorized & type. 
comparisons(counter).indices = find(indices_motorized & indices_type);

counter = counter + 1;

% spontaneous
comparisons(counter).name = 'spontaneous_restvstop_categorical';
comparisons(counter).variablesToUse = {'type_dummyvars_vector'};
comparisons(counter).mice_not_to_use = {'1100'};
comparisons(counter).plotMultiplier = -1;
comparisons(counter).figure_type = 'startstop';

% Get relevent indices for this type.
indices_type = strcmp(period_types, 'rest') | strcmp(period_types, 'stop');

% Get intersection of spontaneous & type. 
comparisons(counter).indices = find(indices_spontaneous & indices_type);

counter = counter + 1; 

%% Rest vs finished stop categorical 

% Motorized
comparisons(counter).name = 'motorized_restvsfinishedstop_categorical';
comparisons(counter).variablesToUse = {'type_dummyvars_vector'};
comparisons(counter).mice_not_to_use = {'1100'};
comparisons(counter).plotMultiplier = -1;
comparisons(counter).figure_type = 'startstop';

% Get relevent indices for this type.
indices_type = strcmp(period_types, 'rest') | strcmp(period_types, 'finished_stop');

% Get intersection of motorized & type. 
comparisons(counter).indices = find(indices_motorized & indices_type);

counter = counter + 1; 

% Spontaneous
comparisons(counter).name = 'spontaneous_restvsfinishedstop_categorical';
comparisons(counter).variablesToUse = {'type_dummyvars_vector'};
comparisons(counter).mice_not_to_use = {'1100'};
comparisons(counter).plotMultiplier = -1;
comparisons(counter).figure_type = 'startstop';

% Get relevent indices for this type.
indices_type = strcmp(period_types, 'rest') | strcmp(period_types, 'finished_stop');

% Get intersection of motorized & type. 
comparisons(counter).indices = find(indices_spontaneous & indices_type);

counter = counter + 1;

%%  walk vs start
% motorized
comparisons(counter).name = 'motorized_walkvsstart_categorical';
comparisons(counter).variablesToUse = {'type_dummyvars_vector'};
comparisons(counter).mice_not_to_use = {'1100'};
comparisons(counter).plotMultiplier = 1;
comparisons(counter).figure_type = 'startstop';

% Get relevent indices for this type.
indices_type = strcmp(period_types, 'walk') | strcmp(period_types, 'start');

% Get intersection of motorized & type. 
comparisons(counter).indices = find(indices_motorized & indices_type);

counter = counter + 1;

% spontaneous
comparisons(counter).name = 'spontaneous_walkvsstart_categorical';
comparisons(counter).variablesToUse = {'type_dummyvars_vector'};
comparisons(counter).mice_not_to_use = {'1100'};
comparisons(counter).plotMultiplier = 1;
comparisons(counter).figure_type = 'startstop';

% Get relevent indices for this type.
indices_type = strcmp(period_types, 'walk') | strcmp(period_types, 'start');

% Get intersection of spontaneous & type. 
comparisons(counter).indices = find(indices_spontaneous & indices_type);

counter = counter + 1; 

%% walk vs stop 
% motorized
comparisons(counter).name = 'motorized_walkvsstop_categorical';
comparisons(counter).variablesToUse = {'type_dummyvars_vector'};
comparisons(counter).mice_not_to_use = {'1100'};
comparisons(counter).plotMultiplier = -1;
comparisons(counter).figure_type = 'startstop';

% Get relevent indices for this type.
indices_type = strcmp(period_types, 'walk') | strcmp(period_types, 'stop');

% Get intersection of motorized & type. 
comparisons(counter).indices = find(indices_motorized & indices_type);

counter = counter + 1;

% spontaneous
comparisons(counter).name = 'spontaneous_walkvsstop_categorical';
comparisons(counter).variablesToUse = {'type_dummyvars_vector'};
comparisons(counter).mice_not_to_use = {'1100'};
comparisons(counter).plotMultiplier = -1;
comparisons(counter).figure_type = 'startstop';

% Get relevent indices for this type.
indices_type = strcmp(period_types, 'walk') | strcmp(period_types, 'stop');

% Get intersection of spontaneous & type. 
comparisons(counter).indices = find(indices_spontaneous & indices_type);

counter = counter + 1; 


%% Stop vs finished stop
% motorized
comparisons(counter).name = 'motorized_stopvsfstop_categorical';
comparisons(counter).variablesToUse = {'type_dummyvars_vector'};
comparisons(counter).mice_not_to_use = {'1100'};
comparisons(counter).plotMultiplier = -1;
comparisons(counter).figure_type = 'startstop';

% Get relevent indices for this type.
indices_type = strcmp(period_types, 'stop') | strcmp(period_types, 'finished_stop');

% Get intersection of motorized & type. 
comparisons(counter).indices = find(indices_motorized & indices_type);

counter = counter + 1;

% Spontaneous
comparisons(counter).name = 'spontaneous_stopvsfstop_categorical';
comparisons(counter).variablesToUse = {'type_dummyvars_vector'};
comparisons(counter).mice_not_to_use = {'1100'};
comparisons(counter).plotMultiplier = -1;
comparisons(counter).figure_type = 'startstop';

% Get relevent indices for this type.
indices_type = strcmp(period_types, 'stop') | strcmp(period_types, 'finished_stop');

% Get intersection of motorized & type. 
comparisons(counter).indices = find(indices_spontaneous & indices_type);

counter = counter + 1;

%% prewalk vs spontaneous rest
comparisons(counter).name = 'prewalkvsrest';
comparisons(counter).variablesToUse = {'type_dummyvars_vector'};
comparisons(counter).mice_not_to_use = {'1100'};
comparisons(counter).plotMultiplier = -1;
comparisons(counter).figure_type = 'warningPeriods';

% Get relevent indices for this type.
indices_type = strcmp(period_types, 'prewalk') | (strcmp(period_types, 'rest')  & strcmp(periods.motorized_vs_spon(:), 'spontaneous'));
comparisons(counter).indices = find(indices_type);

counter = counter + 1; 

%% wstart vs wmaint
comparisons(counter).name = 'wstartvswmaint';
comparisons(counter).variablesToUse = {'type_dummyvars_vector'};
comparisons(counter).mice_not_to_use = {'1100'};
comparisons(counter).plotMultiplier = 1;
comparisons(counter).figure_type = 'warningPeriods';

% Get relevent indices for this type.
indices_type = strcmp(period_types, 'w_start') | (strcmp(period_types, 'w_maint_rest')  & strcmp(periods.motorized_vs_spon(:), 'motorized'));
comparisons(counter).indices = find(indices_type);

counter = counter + 1; 

%% wmaint vs motorized rest
comparisons(counter).name = 'wmaintvsrest';
comparisons(counter).variablesToUse = {'type_dummyvars_vector'};
comparisons(counter).mice_not_to_use = {'1100'};
comparisons(counter).plotMultiplier = -1;
comparisons(counter).figure_type = 'warningPeriods';

% Get relevent indices for this type.
indices_type = (strcmp(period_types, 'rest') & strcmp(periods.motorized_vs_spon(:), 'motorized')) | (strcmp(period_types, 'w_maint_rest')  & strcmp(periods.motorized_vs_spon(:), 'motorized'));
comparisons(counter).indices = find(indices_type);

counter = counter + 1; 

%% wstart vs motorized rest (for completeness)
comparisons(counter).name = 'wstartvsrest';
comparisons(counter).variablesToUse = {'type_dummyvars_vector'};
comparisons(counter).mice_not_to_use = {'1100'};
comparisons(counter).plotMultiplier = -1;
comparisons(counter).figure_type = 'warningPeriods';

% Get relevent indices for this type.
indices_type = strcmp(period_types, 'w_start') | (strcmp(period_types, 'rest')  & strcmp(periods.motorized_vs_spon(:), 'motorized'));
comparisons(counter).indices = find(indices_type);

counter = counter + 1; 

%% Motorized walk vs accel
comparisons(counter).name = 'motorized_walkvsaccel_categorical';
comparisons(counter).variablesToUse = {'type_dummyvars_vector'};
comparisons(counter).mice_not_to_use = {}; %{'1100'};
comparisons(counter).plotMultiplier = -1;
comparisons(counter).figure_type = 'acceldecel';
comparisons(counter).fromContinuous = 'accel_and_decel';

% Get relevent indices for this type.
indices_type = strcmp(period_types, 'walk') | strcmp(period_types, 'accel');

% Get intersection of motorized & type. 
comparisons(counter).indices = find(indices_motorized & indices_type);

counter = counter + 1;

%% Motorized walk vs decel 
comparisons(counter).name = 'motorized_walkvsdecel_categorical';
comparisons(counter).variablesToUse = {'type_dummyvars_vector'};
comparisons(counter).mice_not_to_use =  {}; % {'1100'};
comparisons(counter).plotMultiplier = -1;
comparisons(counter).figure_type = 'acceldecel';
comparisons(counter).fromContinuous = 'accel_and_decel';

% Get relevent indices for this type.
indices_type = strcmp(period_types, 'walk') | strcmp(period_types, 'decel');

% Get intersection of motorized & type. 
comparisons(counter).indices = find(indices_motorized & indices_type);

counter = counter + 1;

%% Accel vs decel categorical.
% Motorized
comparisons(counter).name = 'motorized_accelvsdecel_categorical';
comparisons(counter).variablesToUse = {'type_dummyvars_vector'};
comparisons(counter).mice_not_to_use =  {}; %{'1100'};
comparisons(counter).plotMultiplier = 1;
comparisons(counter).figure_type = 'acceldecel';
comparisons(counter).fromContinuous = 'accel_and_decel';

% Get relevent indices for this type.
indices_type = strcmp(period_types, 'accel') | strcmp(period_types, 'decel');

% Get intersection of motorized & type. 
comparisons(counter).indices = find(indices_motorized & indices_type);

counter = counter + 1; 

%% Accel vs finished accel 
% motorized
comparisons(counter).name = 'motorized_accelvsfaccel_categorical';
comparisons(counter).variablesToUse = {'type_dummyvars_vector'};
comparisons(counter).mice_not_to_use =  {}; % {'1100'};
comparisons(counter).plotMultiplier = -1;
comparisons(counter).figure_type = 'acceldecel';
comparisons(counter).fromContinuous = 'accel_and_decel';

% Get relevent indices for this type.
indices_type = strcmp(period_types, 'accel') | strcmp(period_types, 'finished_accel');

% Get intersection of motorized & type. 
comparisons(counter).indices = find(indices_motorized & indices_type);

counter = counter + 1;

%% Decel vs finished decel    

comparisons(counter).name = 'motorized_decelvsfdecel_categorical';
comparisons(counter).variablesToUse = {'type_dummyvars_vector'};
comparisons(counter).mice_not_to_use =  {}; % {'1100'};
comparisons(counter).plotMultiplier = -1;
comparisons(counter).figure_type = 'acceldecel';
comparisons(counter).fromContinuous = 'accel_and_decel';

% Get relevent indices for this type.
indices_type = strcmp(period_types, 'decel') | strcmp(period_types, 'finished_decel');

% Get intersection of motorized & type. 
comparisons(counter).indices = find(indices_motorized & indices_type);

counter = counter + 1;

%% Walk vs finished accel
% motorized
comparisons(counter).name = 'motorized_walkvsfaccel_categorical';
comparisons(counter).variablesToUse = {'type_dummyvars_vector'};
comparisons(counter).mice_not_to_use =  {}; %{'1100'};
comparisons(counter).plotMultiplier = -1;
comparisons(counter).figure_type = 'acceldecel';
comparisons(counter).fromContinuous = 'accel_and_decel';

% Get relevent indices for this type.
indices_type = strcmp(period_types, 'walk') | strcmp(period_types, 'finished_accel');

% Get intersection of motorized & type. 
comparisons(counter).indices = find(indices_motorized & indices_type);

counter = counter + 1;

%% Walk vs finished decel 
comparisons(counter).name = 'motorized_walkvsfdecel_categorical';
comparisons(counter).variablesToUse = {'type_dummyvars_vector'};
comparisons(counter).mice_not_to_use =  {}; %{'1100'};
comparisons(counter).plotMultiplier = -1;
comparisons(counter).figure_type = 'acceldecel';
comparisons(counter).fromContinuous = 'accel_and_decel';

% Get relevent indices for this type.
indices_type = strcmp(period_types, 'walk') | strcmp(period_types, 'finished_decel');

% Get intersection of motorized & type. 
comparisons(counter).indices = find(indices_motorized & indices_type);

counter = counter + 1;


%% finsihed accel vs finished decel
comparisons(counter).name = 'motorized_faccelvsfdecel_categorical';
comparisons(counter).variablesToUse = {'type_dummyvars_vector'};
comparisons(counter).mice_not_to_use =  {}; %{'1100'};
comparisons(counter).plotMultiplier = -1;
comparisons(counter).figure_type = 'acceldecel';
comparisons(counter).fromContinuous = 'accel_and_decel';

% Get relevent indices for this type.
indices_type = strcmp(period_types, 'finished_accel') | strcmp(period_types, 'finished_decel');

% Get intersection of motorized & type. 
comparisons(counter).indices = find(indices_motorized & indices_type);

counter = counter + 1;
%% prewalk vs rest
% comparisons(counter).name = 'prewalkvsrest';
% comparisons(counter).variablesToUse = {'type_dummyvars_vector'};
% comparisons(counter).mice_not_to_use = {'1100'};
% comparisons(counter).plotMultiplier = 1;
% comparisons(counter).figure_type = 'warningPeriods';
% comparisons(counter).fromContinuous = 'pre_and_post';
% 
% % Get relevent indices for this type.
% indices_type = strcmp(period_types, 'prewalk') | strcmp(period_types, 'rest');
% 
% % Get intersection of spontaneous & type. 
% comparisons(counter).indices = find(indices_spontaneous & indices_type);
% 
% counter = counter + 1;
% 
% %% w_start vs w_maint_rest
% comparisons(counter).name = 'wstartvswmaint';
% comparisons(counter).variablesToUse = {'type_dummyvars_vector'};
% comparisons(counter).mice_not_to_use = {'1100'};
% comparisons(counter).plotMultiplier = 1;
% comparisons(counter).figure_type = 'warningPeriods';
% comparisons(counter).fromContinuous = 'pre_and_post';
% 
% % Get relevent indices for this type.
% indices_type = strcmp(period_types, 'w_start') | strcmp(period_types, 'w_maint_rest');
% 
% % Get intersection of spontaneous & type. 
% comparisons(counter).indices = find(indices_motorized & indices_type);
% 
% counter = counter + 1;
% 


%% Save 
save(filename_out, 'comparisons');