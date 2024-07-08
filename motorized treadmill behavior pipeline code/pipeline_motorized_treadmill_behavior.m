% pipeline_motorized_treadmill_behavior.m
% Sarah West
% 11/16/21 

% Use "create_mice_all.m" and "create_conditions_names.m" before using this.

%% Initial setup
% Put all needed parameters in a structure called "parameters", which you
% can then easily feed into your functions. 
clear all; 

% Output Directories

% Create the experiment name. This is used to name the output folder. 
parameters.experiment_name='Random Motorized Treadmill';

% Output directory name bases
parameters.dir_base='Y:\Sarah\Analysis\Experiments\';
parameters.dir_exper=[parameters.dir_base parameters.experiment_name '\']; 

% *********************************************************
% Data to process

% (DON'T EDIT). Load the "mice_all" variable you've created with "create_mice_all.m"
load([parameters.dir_exper 'mice_all.mat']);

% Add mice_all to parameters structure.
parameters.mice_all = mice_all; 

% Load the names of conditions from create_conditions_names.m
load([parameters.dir_exper 'Behavior_Conditions.mat']); 

% Add Conditions to parameters structure.
parameters.Conditions = Conditions; 

% ****Change here if there are specific mice, days, and/or stacks you want to work with**** 
parameters.mice_all = parameters.mice_all;
%parameters.mice_all(1).days = parameters.mice_all(1).days(13:end);


% **********************************************************************8
% Input Directories

% Establish the format of the daily/per mouse directory and file names of 
% the collected data. Will be assembled with CreateStrings.m Each piece 
% needs to be a separate entry in a cell  array. Put the string 'mouse 
% number', 'day', or 'stack number' where the mouse, day, or stack number 
% will be. If you concatenated this as a sigle string, it should create a 
% file name, with the correct mouse/day/stack name inserted accordingly. 
parameters.dir_dataset_name = {'Y:\Sarah\Data\Random Motorized Treadmill\', 'day', '\', 'mouse number', '\Arduino output\'};
%parameters.input_data_name = {'0', 'stack number', '.txt'};
parameters.input_data_name={'Motor_ArduinoOutput*.log'}; 

% Give the number of digits that should be included in each stack number.
parameters.digitNumber=2; 

% *************************************************************************
% Parameters

% Sampling frequency of serial monitor, in Hz. Arduino displays time in ms,
% so is almost always 1000 Hz.
parameters.wheel_Hz = 1000;

% Sampling frequency of collected brain data (per channel), in Hz or frames per
% second.
parameters.fps= 20; 

% Amount of time after a transition that you want to start counting as a
% "continued" behavior, in seconds. (In spontanuous locomotion paper,
% fluorescence took ~2.5 s to return to baseline after the mice stopped.)
parameters.continued_window = 3; 

% Length of time (in seconds) you want continued behaviors to be divided
% into.
parameters.continued_chunk_length = 1; 

% Number of channels from brain data (need this to calculate correct
% "skip" time length).
parameters.channelNumber = 2;

% Number of frames you recorded from brain and want to keep (don't make chunks longer than this)  
parameters.frames=6000; 

% Number of initial brain frames to skip, allows for brightness/image
% stabilization of camera. Need this to know how much to skip in the
% behavior.
parameters.skip = 1200; 

% Using variable accelerations?
parameters.useAccel = true;

% Radius of wheel, in cm.
parameters.wheel_radius = 8.5;
                                  
% Smallest time increment. All accelerations and stage times should be an
% exact whole second or half second in duration. 
parameters.smallest_time_increment = 0.5 * parameters.fps;

% FOR 3S CONTINUED BEHAVIOR
% parameters.experiment_name= ['Random Motorized Treadmill\' num2str(parameters.continued_chunk_length) 's continued'];
% parameters.dir_exper=[parameters.dir_base parameters.experiment_name '\']; 

%% Extract data and save as .mat file -- with Putty

% Go back and reset mice_all each time! 
parameters.mice_all = mice_all(4:6);
% parameters.mice_all(1).days = parameters.mice_all(1).days(12);
% parameters.mice_all(1).days(1).stacks = {'16'};

% Always clear loop list first. 
if isfield(parameters, 'loop_list')
parameters = rmfield(parameters,'loop_list');
end

% Was PUTTY used for the recording? 
parameters.putty_flag = true;

% Get only mice_all days with putty_for_motor = 'yes'.
for i = 1:numel(parameters.mice_all)
    putty_flags = NaN(numel(parameters.mice_all(i).days),1);
    for dayi = 1: numel(parameters.mice_all(i).days)
        % Get the putty value, remove any spaces
        value = strrep(parameters.mice_all(i).days(dayi).putty_for_motor, ' ', '');
        % Compare to 'yes'
        putty_flags(dayi) = strcmp(value, 'yes');
    end
    parameters.mice_all(i).days(~putty_flags) = [];
end 


% Iterations.
parameters.loop_list.iterators = {'mouse', {'loop_variables.mice_all(:).name'}, 'mouse_iterator'; 
               'day', {'loop_variables.mice_all(', 'mouse_iterator', ').days(:).name'}, 'day_iterator';
               'log', { 'dir("Y:\Sarah\Data\Random Motorized Treadmill\', 'day', '\', 'mouse', '\Arduino Output\Motor_ArduinoOutput*.log").name'}, 'log_iterator'; 
               'stack', {'loop_variables.mice_all(',  'mouse_iterator', ').days(', 'day_iterator', ').stacks'}, 'stack_iterator'};

parameters.loop_variables.mice_all = parameters.mice_all;

% Input values
parameters.loop_list.things_to_load.log.dir = {'Y:\Sarah\Data\Random Motorized Treadmill\', 'day', '\', 'mouse', '\Arduino Output\'};
parameters.loop_list.things_to_load.log.filename= {'log'}; 
parameters.loop_list.things_to_load.log.variable= {}; 
parameters.loop_list.things_to_load.log.level = 'log';
parameters.loop_list.things_to_load.log.load_function = @readtext;

% Output values. 
parameters.loop_list.things_to_save.trial.dir = {[parameters.dir_exper 'behavior\motorized\extracted motor data\'], 'mouse', '\', 'day', '\'};
parameters.loop_list.things_to_save.trial.filename= {'trial', 'stack', '.mat'};
parameters.loop_list.things_to_save.trial.variable= {'trial'}; 
parameters.loop_list.things_to_save.trial.level = 'stack';

RunAnalysis({@extractMotorData}, parameters);

%% Extract data and save as .mat file -- with NO Putty

% Go back and reset mice_all each time! 
parameters.mice_all = mice_all;

% Always clear loop list first. 
if isfield(parameters, 'loop_list')
parameters = rmfield(parameters,'loop_list');
end

% Get only mice_all days with putty_for_motor = 'no'.
for i = 1:numel(parameters.mice_all)
    putty_flags = NaN(numel(parameters.mice_all(i).days),1);
    for dayi = 1: numel(parameters.mice_all(i).days)
        % Get the putty value, remove any spaces
        value = strrep(parameters.mice_all(i).days(dayi).putty_for_motor, ' ', '');
        % Compare to 'no'
        putty_flags(dayi) = strcmp(value, 'no');
    end
    parameters.mice_all(i).days(~putty_flags) = [];
end 


% Was PUTTY used for the recording? 
parameters.putty_flag = false;

% Iterations.
parameters.loop_list.iterators = {'mouse', {'loop_variables.mice_all(:).name'}, 'mouse_iterator'; 
               'day', {'loop_variables.mice_all(', 'mouse_iterator', ').days(:).name'}, 'day_iterator';
               'stack', {'loop_variables.mice_all(',  'mouse_iterator', ').days(', 'day_iterator', ').stacks'}, 'stack_iterator'};

parameters.loop_variables.mice_all = parameters.mice_all;

% Input values
parameters.loop_list.things_to_load.log.dir = {'Y:\Sarah\Data\Random Motorized Treadmill\', 'day', '\', 'mouse', '\Arduino Output\'};
parameters.loop_list.things_to_load.log.filename= {'0', 'stack', '.txt'}; 
parameters.loop_list.things_to_load.log.variable= {}; 
parameters.loop_list.things_to_load.log.level = 'stack';
parameters.loop_list.things_to_load.log.load_function = @readtext;

% Output values. 
parameters.loop_list.things_to_save.trial.dir = {[parameters.dir_exper 'behavior\motorized\extracted motor data\'], 'mouse', '\', 'day', '\'};
parameters.loop_list.things_to_save.trial.filename= {'trial', 'stack', '.mat'};
parameters.loop_list.things_to_save.trial.variable= {'trial'}; 
parameters.loop_list.things_to_save.trial.level = 'stack';

RunAnalysis({@extractMotorData}, parameters);

%% Find behavior periods. 
% (Don't need to separate Putty vs non-Putty)

% Always clear loop list first. 
if isfield(parameters, 'loop_list')
parameters = rmfield(parameters,'loop_list');
end

% Skip any files that don't exist. 
parameters.load_abort_flag = true; 

% Don't remove empty iterations. 
parameters.removeEmptyIterations = false;  

% Reset mice_all (because you're changing it for the Putty only stacks)
% parameters.mice_all = mice_all; 
% parameters.mice_all = parameters.mice_all(1);

% Iterators
parameters.loop_list.iterators = {'mouse', {'loop_variables.mice_all(:).name'}, 'mouse_iterator'; 
               'day', {'loop_variables.mice_all(', 'mouse_iterator', ').days(:).name'}, 'day_iterator';
               'stack', {'loop_variables.mice_all(',  'mouse_iterator', ').days(', 'day_iterator', ').stacks'}, 'stack_iterator'};

parameters.loop_variables.mice_all = parameters.mice_all;
parameters.duration_place_maximum_default = 25; 

% Tell code where to find Putty and Accel flags as strings in CreateStrings
% format
parameters.location_putty_flag = {'parameters.mice_all(', 'mouse_iterator', ').days(', 'day_iterator', ').putty_for_motor'};
parameters.location_accel_flag = {'parameters.mice_all(', 'mouse_iterator', ').days(', 'day_iterator', ').number_accels_used'};

% Input
parameters.loop_list.things_to_load.trial.dir = {[parameters.dir_exper 'behavior\motorized\extracted motor data\'], 'mouse', '\', 'day', '\'};
parameters.loop_list.things_to_load.trial.filename= {'trial', 'stack', '.mat'};
parameters.loop_list.things_to_load.trial.variable= {'trial'}; 
parameters.loop_list.things_to_load.trial.level = 'stack';

% Output
% all periods
parameters.loop_list.things_to_save.all_periods.dir = {[parameters.dir_exper 'behavior\motorized\period instances\'], 'mouse', '\', 'day', '\'};
parameters.loop_list.things_to_save.all_periods.filename= {'all_periods_', 'stack', '.mat'};
parameters.loop_list.things_to_save.all_periods.variable= {'all_periods'}; 
parameters.loop_list.things_to_save.all_periods.level = 'stack';
% time ranges of not-brokendown rest and walk
parameters.loop_list.things_to_save.long_periods.dir = {[parameters.dir_exper 'behavior\motorized\period instances\'], 'mouse', '\', 'day', '\'};
parameters.loop_list.things_to_save.long_periods.filename= {'long_periods_', 'stack', '.mat'};
parameters.loop_list.things_to_save.long_periods.variable= {'long_periods'}; 
parameters.loop_list.things_to_save.long_periods.level = 'stack';

RunAnalysis({@motorFindBehaviorPeriods}, parameters);

%% Make into more later code-readable table format

% Use the periods.mat to help you tell which fields to care about, use the
% periods_nametable.mat to collect the the time ranges. 

% Load periods.mat for conditions names 
load([parameters.dir_exper 'periods_nametable.mat']); 
parameters.Conditions = periods; 

% Always clear loop list first. 
if isfield(parameters, 'loop_list')
parameters = rmfield(parameters,'loop_list');
end

% Skip any files that don't exist. 
parameters.load_abort_flag = true; 

% Iterators
parameters.loop_list.iterators = {'mouse', {'loop_variables.mice_all(:).name'}, 'mouse_iterator'; 
               'day', {'loop_variables.mice_all(', 'mouse_iterator', ').days(:).name'}, 'day_iterator';
               'stack', {'loop_variables.mice_all(',  'mouse_iterator', ').days(', 'day_iterator', ').stacks'}, 'stack_iterator'};

parameters.loop_variables.mice_all = parameters.mice_all;

% Input
parameters.loop_list.things_to_load.all_periods.dir = {[parameters.dir_exper 'behavior\motorized\period instances\'], 'mouse', '\', 'day', '\'};
parameters.loop_list.things_to_load.all_periods.filename= {'all_periods_', 'stack', '.mat'};
parameters.loop_list.things_to_load.all_periods.variable= {'all_periods'}; 
parameters.loop_list.things_to_load.all_periods.level = 'stack';

% Output 
parameters.loop_list.things_to_save.all_periods_table.dir = {[parameters.dir_exper 'behavior\motorized\period instances table format\'], 'mouse', '\', 'day', '\'};
parameters.loop_list.things_to_save.all_periods_table.filename= {'all_periods_', 'stack', '.mat'};
parameters.loop_list.things_to_save.all_periods_table.variable= {'all_periods'}; 
parameters.loop_list.things_to_save.all_periods_table.level = 'stack';
 
RunAnalysis({@MotorBehaviorPeriodsTable}, parameters); 

