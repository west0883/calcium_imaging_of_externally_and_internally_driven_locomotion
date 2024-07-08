% create_PLSR_comparisons_special_continuous.m
% Sarah West
% 9/12/23

% A script that creates a structure of information to run comparisons
% between different subsets of data for the Random Motorized Treadmill
% experiments. Each comparison will be run per mouse, eventually run a
% second-level analysis across mice. 

% (This is for first-level PLSR comparisons);

% 3 main pieces for each comparison:
% name --> the working name of the comparison, for file names.
% variablesToUse --> the response variable categories you're interested in
% indices --> the period indices from periods_nametable that are relevant
% to this comparison & will find the related brain data & response
% variable pairs. Need these indices when you're first getting the comparison's
% dataset from the larger 
% type --> just the name of the behavior type in its own entry, to make it
% easier to manipulate later.
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
continuous_variable_names = {'speed_vector', 'accel_vector', 'duration_vector', 'pupil_diameter_vector', 'tail_vector', 'nose_vector', 'FL_vector', 'HL_vector', 'x_vector'};

% Pull out the relevant columns/info from periods
period_motorized_vs_spon = periods.motorized_vs_spon;
period_types = periods.type;
indices_motorized = strcmp(period_motorized_vs_spon, 'motorized');
indices_spontaneous = strcmp(period_motorized_vs_spon, 'spontaneous');

% can search tables like this:
% g = periods_table(string(periods_table.condition)=="m_accel" & string(periods_table.speed)=="2000", :);

counter = 1; 


%% Walk: motorized & Spontaneous 
% Include accel as well.
comparisons(counter).name = 'all_walk_continuousVars';
comparisons(counter).variablesToUse = {'speed_vector', 'accel_vector', 'duration_vector', 'pupil_diameter_vector', 'tail_vector', 'nose_vector', 'FL_vector', 'HL_vector', 'x_vector'};
comparisons(counter).type = 'walk';
comparisons(counter).mice_not_to_use = {'1100'};
comparisons(counter).figure_type = 'continued';

% Get relevent indices for this type.
indices_type = strcmp(period_types, 'walk');

% Get intersection of spontaneous & type. 
comparisons(counter).indices = find(indices_type);

counter = counter + 1;

%% rest: motorized & spon
comparisons(counter).name = 'all_rest_continuousVars';
comparisons(counter).variablesToUse = {'duration_vector', 'pupil_diameter_vector', 'tail_vector', 'nose_vector', 'FL_vector', 'HL_vector', 'x_vector' };
comparisons(counter).type = 'rest';
comparisons(counter).mice_not_to_use = {'1100'};
comparisons(counter).figure_type = 'continued';

% Get relevent indices for this type.
indices_type = strcmp(period_types, 'rest');

% Get intersection of motorized & type. 
comparisons(counter).indices = find(indices_type);

counter = counter + 1; 

%% start: motorized & spon
comparisons(counter).name = 'all_start_continuousVars';
comparisons(counter).variablesToUse = {'speed_vector', 'accel_vector', 'duration_vector', 'pupil_diameter_vector', 'tail_vector', 'nose_vector', 'FL_vector', 'HL_vector', 'x_vector' };
comparisons(counter).type = 'start';
comparisons(counter).mice_not_to_use = {'1100'};
comparisons(counter).figure_type = 'startstop';

% Get relevent indices for this type.
indices_type = strcmp(period_types, 'start');

% Get intersection of motorized & type. 
comparisons(counter).indices = find(indices_type);

counter = counter + 1; 

%% stop: motorized & spon
comparisons(counter).name = 'all_stop_continuousVars';
comparisons(counter).variablesToUse = {'speed_vector', 'accel_vector', 'duration_vector', 'pupil_diameter_vector', 'tail_vector', 'nose_vector', 'FL_vector', 'HL_vector', 'x_vector' };
comparisons(counter).type = 'stop';
comparisons(counter).mice_not_to_use = {'1100'};
comparisons(counter).figure_type = 'startstop';

% Get relevent indices for this type.
indices_type = strcmp(period_types, 'stop');

% Get intersection of motorized & type. 
comparisons(counter).indices = find(indices_type);

counter = counter + 1; 

%% finished stop periods: Motorized & spon
comparisons(counter).name = 'all_finished_stop_continuousVars';
comparisons(counter).variablesToUse = {'speed_vector', 'accel_vector', 'duration_vector', 'pupil_diameter_vector', 'tail_vector', 'nose_vector', 'FL_vector', 'HL_vector', 'x_vector'};
comparisons(counter).type = 'finished_stop';
comparisons(counter).mice_not_to_use = {'1100'};
comparisons(counter).figure_type = 'startstop';

% Get relevent indices for this type.
indices_type = strcmp(period_types, 'finished_stop');

% Get type. 
comparisons(counter).indices = find(indices_type);

counter = counter + 1; 

%% Prewalks: motorized vs spon
types = {'w_start', 'prewalk'};
comparisons(counter).name = 'all_prewalks';
comparisons(counter).variablesToUse = {'speed_vector', 'accel_vector', 'duration_vector', 'pupil_diameter_vector', 'tail_vector', 'nose_vector', 'FL_vector', 'HL_vector', 'x_vector'};
comparisons(counter).type = '';
comparisons(counter).mice_not_to_use = {'1100'};
comparisons(counter).figure_type = 'warningPeriods';

indices_type = [];

% Get relevent indices for this type.
for typei = 1:numel(types)
    indices = find(strcmp(period_types, types{typei}));
    indices_type = [indices_type; indices];
end  

% Get type. 
comparisons(counter).indices = indices_type;

counter = counter + 1; 

%% Prewalks 2: include wmaint
types = {'w_start', 'prewalk', 'w_maint_rest'};
comparisons(counter).name = 'all_prewalks_2';
comparisons(counter).variablesToUse = {'speed_vector', 'accel_vector', 'duration_vector', 'pupil_diameter_vector', 'tail_vector', 'nose_vector', 'FL_vector', 'HL_vector', 'x_vector'};
comparisons(counter).type = '';
comparisons(counter).mice_not_to_use = {'1100'};
comparisons(counter).figure_type = 'warningPeriods';

indices_type = [];

% Get relevent indices for this type.
for typei = 1:numel(types)
    indices = find(strcmp(period_types, types{typei}));
    indices_type = [indices_type; indices];
end  

% Get type. 
comparisons(counter).indices = indices_type;

counter = counter + 1; 


%% Accelerating

% Only relevant for motorized.
comparisons(counter).name = 'motorized_accel_continuousVars';
comparisons(counter).variablesToUse = {'speed_vector', 'accel_vector', 'duration_vector', 'pupil_diameter_vector', 'tail_vector', 'nose_vector', 'FL_vector', 'HL_vector', 'x_vector'};
comparisons(counter).type = 'accel';
comparisons(counter).mice_not_to_use = {'1100'};
comparisons(counter).figure_type = 'acceldecel';

% Get relevent indices for this type.
indices_type = strcmp(period_types, 'accel');

% Get intersection of motorized & type. 
comparisons(counter).indices = find(indices_motorized & indices_type);

counter = counter + 1;

%% Decelerating
% Only relevant for motorized.
comparisons(counter).name = 'motorized_decel_continuousVars';
comparisons(counter).variablesToUse = {'speed_vector', 'accel_vector', 'duration_vector', 'pupil_diameter_vector', 'tail_vector', 'nose_vector', 'FL_vector', 'HL_vector', 'x_vector'};
comparisons(counter).type = 'decel';
comparisons(counter).mice_not_to_use = {'1100'};
comparisons(counter).figure_type = 'acceldecel';

% Get relevent indices for this type.
indices_type = strcmp(period_types, 'decel');

% Get intersection of motorized & type. 
comparisons(counter).indices = find(indices_motorized & indices_type);

counter = counter + 1;
%% Continuous variables, finished_accel
% Use speed & duration. 

% Only relevant for motorized.
comparisons(counter).name = 'motorized_finished_accel_continuousVars';
comparisons(counter).variablesToUse = {'speed_vector', 'duration_vector', 'pupil_diameter_vector', 'tail_vector', 'nose_vector', 'FL_vector', 'HL_vector', 'x_vector'};
comparisons(counter).type = 'finished_accel';
comparisons(counter).mice_not_to_use = {'1100'};
comparisons(counter).figure_type = 'acceldecel';

% Get relevent indices for this type.
indices_type = strcmp(period_types, 'finished_accel');

% Get intersection of motorized & type. 
comparisons(counter).indices = find(indices_motorized & indices_type);

counter = counter + 1;

%% Continuous variables, finished_decel
% Use speed & duration. 

% Only relevant for motorized.
comparisons(counter).name = 'motorized_finished_decel_continuousVars';
comparisons(counter).variablesToUse = {'speed_vector', 'duration_vector', 'pupil_diameter_vector', 'tail_vector', 'nose_vector', 'FL_vector', 'HL_vector', 'x_vector'};
comparisons(counter).type = 'finished_decel';
comparisons(counter).mice_not_to_use = {'1100'};
comparisons(counter).figure_type = 'acceldecel';

% Get relevent indices for this ty pe.
indices_type = strcmp(period_types, 'finished_decel');

% Get intersection of motorized & type. 
comparisons(counter).indices = find(indices_motorized & indices_type);

counter = counter + 1;

%% pre,post,and rest
% types = {'w_start', 'prewalk', 'w_maint_rest', 'c_rest', 'rest', 'finished_stop', 'postwalk'};
% comparisons(counter).name = 'pre_and_post';
% comparisons(counter).variablesToUse = {'speed_vector', 'accel_vector', 'duration_vector', 'pupil_diameter_vector', 'tail_vector', 'nose_vector', 'FL_vector', 'HL_vector', 'x_vector'};
% comparisons(counter).type = '';
% comparisons(counter).mice_not_to_use = {'1100'};
% comparisons(counter).figure_type = 'warningPeriods';
% 
% indices_type = [];
% 
% % Get relevent indices for this type.
% for typei = 1:numel(types)
%     indices = find(strcmp(period_types, types{typei}));
%     indices_type = [indices_type; indices];
% end  
% 
% % Get type. 
% comparisons(counter).indices = indices_type;
% 
% counter = counter + 1; 

%% accels and decels 
% types = {'walk', 'accel', 'decel', 'finished_accel', 'finished_decel'};
% comparisons(counter).name = 'accel_and_decel';
% comparisons(counter).variablesToUse = {'speed_vector', 'accel_vector', 'duration_vector', 'pupil_diameter_vector', 'tail_vector', 'nose_vector', 'FL_vector', 'HL_vector', 'x_vector'};
% comparisons(counter).type = '';
% comparisons(counter).mice_not_to_use = {};
% comparisons(counter).figure_type = 'acceldecel';
% 
% indices_type = [];
% 
% % Get relevent indices for this type.
% for typei = 1:numel(types)
%     indices = strcmp(period_types, types{typei});
% 
%     % Make sure you only grab motorized
%     indices = find(indices_motorized & indices);
% 
%     indices_type = [indices_type; indices];
% end  
% 
% % Get type. 
% comparisons(counter).indices = indices_type;
% 
% counter = counter + 1; 

%% prewalk periods 
% types = {'c_rest', 'rest', 'w_maint_rest', 'w_start', 'prewalk'};
% comparisons(counter).name = 'all_prewalks';
% comparisons(counter).variablesToUse = {'speed_vector', 'accel_vector', 'duration_vector', 'pupil_diameter_vector', 'tail_vector', 'nose_vector', 'FL_vector', 'HL_vector', 'x_vector'};
% comparisons(counter).type = '';
% comparisons(counter).mice_not_to_use = {'1100'};
% comparisons(counter).figure_type = 'startstop';
% 
% indices_type = [];
% 
% % Get relevent indices for this type.
% for typei = 1:numel(types)
%     indices = find(strcmp(period_types, types{typei}));
%     indices_type = [indices_type; indices];
% end  
% 
% % Get type. 
% comparisons(counter).indices = indices_type;
% 
% counter = counter + 1; 
%% Save 
save([parameters.dir_exper 'PLSR\comparisons_special_continuous.mat'], 'comparisons');