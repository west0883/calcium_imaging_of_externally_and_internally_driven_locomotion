% pipeline_duration_place_concatenation.m
% Sarah West 
% 6/12/23

%% Initial setup
% Put all needed paramters in a structure called "parameters", which you
% can then easily feed into your functions. 
clear all; 

% Output Directories

% Create the experiment name. This is used to name the output folder. 
parameters.experiment_name='Random Motorized Treadmill';

% Output directory name bases
parameters.dir_base='Y:\Sarah\Analysis\Experiments\';
parameters.dir_exper=[parameters.dir_base parameters.experiment_name '\']; 

% *********************************************************

% (DON'T EDIT). Load the "mice_all" variable you've created with "create_mice_all.m"
load([parameters.dir_exper 'mice_all.mat']);

% Add mice_all to parameters structure.
parameters.mice_all = mice_all; 

% ****Change here if there are specific mice, days, and/or stacks you want to work with****
parameters.mice_all = parameters.mice_all;
% 

% Load names of motorized periods
load([parameters.dir_exper 'periods_nametable.mat']);
periods_motorized = periods;

% Load names of spontaneous periods
load([parameters.dir_exper 'periods_nametable_spontaneous.mat']);
periods_spontaneous = periods(1:6, :);
clear periods; 

% Create a shared motorized & spontaneous list.
periods_bothConditions = [periods_motorized; periods_spontaneous]; 
parameters.periods_bothConditions = periods_bothConditions; 

% Give the number of digits that should be included in each stack number.
parameters.digitNumber=2; 

parameters.loop_variables.periods = {'rest', 'walk'}; % spontaneous continued periods 
parameters.loop_variables.mice_all = parameters.mice_all;
parameters.loop_variables.conditions =   {'motorized'; 'spontaneous'};
parameters.loop_variables.conditions_stack_locations =  {'stacks'; 'spontaneous'};
parameters.loop_variables.periods_bothConditions = periods_bothConditions.condition;

%% Motorized: Look for stacks when number of instances don't match those of fluorescence.
% (Just do motorized rest, that's the one you're having problems with).
% Don't need to load, so don't use RunAnalysis.
% Only need to run 1 body part & velocity direction (FL, total magnitude).

parameters.loop_list.iterators = {
               'mouse', {'loop_variables.mice_all(:).name'}, 'mouse_iterator'; 
               'day', {'loop_variables.mice_all(', 'mouse_iterator', ').days(:).name'}, 'day_iterator';
               'stack', {'loop_variables.mice_all(',  'mouse_iterator', ').days(', 'day_iterator', ').stacks'}, 'stack_iterator';
               };

% If using a sub-structure, need to use regular loading
parameters.use_substructure = true; 

parameters.checkingDim = 1;
parameters.check_againstDim = 3;

% Input values
parameters.loop_list.things_to_check.dir = {[parameters.dir_exper 'behavior\motorized\period instances table format\'], 'mouse', '\', 'day', '\'};
parameters.loop_list.things_to_check.filename= {'all_periods_', 'stack', '.mat'};
parameters.loop_list.things_to_check.variable= {'all_periods.duration_place(176:180)'}; 
parameters.loop_list.things_to_check.level = 'stack';

parameters.loop_list.check_against.dir = {[parameters.dir_exper 'fluorescence analysis\segmented timeseries\motorized\'],  'mouse', '\', 'day', '\'};
parameters.loop_list.check_against.filename= {'segmented_timeseries_', 'stack', '.mat'};  
parameters.loop_list.check_against.variable = {'segmented_timeseries(176:180)'};

% Output
parameters.loop_list.mismatched_data.dir = {[parameters.dir_exper 'behavior\duration place concatenated\']};
parameters.loop_list.mismatched_data.filename= {'mismatched_data_motorized.mat'};

CheckSizes2(parameters);

%% Spontaneous: Look for stacks when number of instances don't match those of fluorescence.
% (Just do motorized rest, that's the one you're having problems with).
% Don't need to load, so don't use RunAnalysis.
% Only need to run 1 body part & velocity direction (FL, total magnitude).

parameters.loop_list.iterators = {
               'mouse', {'loop_variables.mice_all(:).name'}, 'mouse_iterator'; 
               'day', {'loop_variables.mice_all(', 'mouse_iterator', ').days(:).name'}, 'day_iterator';
               'stack', {'loop_variables.mice_all(',  'mouse_iterator', ').days(', 'day_iterator', ').spontaneous'}, 'stack_iterator';
               'period', {'loop_variables.periods'}, 'period_iterator';
               };

% If using a sub-structure, need to use regular loading
parameters.use_substructure = true; 

parameters.checkingDim = 1;
parameters.check_againstDim = 3;

% Input values
parameters.loop_list.things_to_check.dir = {[parameters.dir_exper 'behavior\spontaneous\segmented behavior periods\'], 'mouse', '\', 'day', '\'};
parameters.loop_list.things_to_check.filename= {'duration_places_', 'stack', '.mat'};
parameters.loop_list.things_to_check.variable= {'duration_places.', 'period'}; 
parameters.loop_list.things_to_check.level = 'stack';

parameters.loop_list.check_against.dir = {[parameters.dir_exper 'fluorescence analysis\segmented timeseries\spontaneous\'],  'mouse', '\', 'day', '\'};
parameters.loop_list.check_against.filename= {'segmented_timeseries_', 'stack', '.mat'};  
parameters.loop_list.check_against.variable = {'segmented_timeseries(', 'period_iterator', ')'};

% Output
parameters.loop_list.mismatched_data.dir = {[parameters.dir_exper 'behavior\duration place concatenated\']};
parameters.loop_list.mismatched_data.filename= {'mismatched_data_spontaneous.mat'};

CheckSizes2(parameters);
%% Notes for removal.

% For these, are okay to just take all but the end instances 

% motorized 

% load mismatched data 
load('Y:\Sarah\Analysis\Experiments\Random Motorized Treadmill\behavior\duration place concatenated\mismatched_data_motorized.mat')

% For each item in mismatched data, 
for itemi = 2:size(mismatched_data, 1)
   
    % load the duration data
    load(['Y:\Sarah\Analysis\Experiments\Random Motorized Treadmill\behavior\motorized\period instances table format\' mismatched_data{itemi, 1} '\' mismatched_data{itemi, 2} '\all_periods_' mismatched_data{itemi, 3} '.mat'])
    duration = all_periods.duration_place(176:180);
    
    % load the flourescence data 
    load(['Y:\Sarah\Analysis\Experiments\Random Motorized Treadmill\fluorescence analysis\segmented timeseries\motorized\' mismatched_data{itemi, 1} '\' mismatched_data{itemi, 2} '\segmented_timeseries_' mismatched_data{itemi, 3} '.mat'])
    fluor = segmented_timeseries;

    % find number of instances you should have
    % if empty, make it set = 0 (because 3rd dim will be 1)
    if  ~isempty(fluor{mismatched_data{itemi, 4} + 175})
        correct_num = size(fluor{mismatched_data{itemi, 4} + 175}, 3);
    else
        correct_num = 0;
    end
    duration_new = duration;
    duration_new{mismatched_data{itemi, 4}} = duration{mismatched_data{itemi, 4}}(1:correct_num);

    % Save over previous 
    all_periods.duration_place(176:180) = duration_new;
    save(['Y:\Sarah\Analysis\Experiments\Random Motorized Treadmill\behavior\motorized\period instances table format\' mismatched_data{itemi, 1} '\' mismatched_data{itemi, 2} '\all_periods_' mismatched_data{itemi, 3} '.mat'], 'all_periods');

end 


%% correct mismatches in spontaneous 

% load mismatched data 
load('Y:\Sarah\Analysis\Experiments\Random Motorized Treadmill\behavior\duration place concatenated\mismatched_data_spontaneous.mat')

% For each item in mismatched data, 
for itemi = 2:size(mismatched_data, 1)
   
    % load the duration data
    load(['Y:\Sarah\Analysis\Experiments\Random Motorized Treadmill\behavior\spontaneous\segmented behavior periods\' mismatched_data{itemi, 1} '\' mismatched_data{itemi, 2} '\duration_places_' mismatched_data{itemi, 3} '.mat'])
    duration = duration_places;
    
    % load the flourescence data 
    load(['Y:\Sarah\Analysis\Experiments\Random Motorized Treadmill\fluorescence analysis\segmented timeseries\spontaneous\' mismatched_data{itemi, 1} '\' mismatched_data{itemi, 2} '\segmented_timeseries_' mismatched_data{itemi, 3} '.mat'])
    fluor = segmented_timeseries;

    % find number of instances you should have
    % if empty, make it set = 0 (because 3rd dim will be 1)
    if  ~isempty(fluor{mismatched_data{itemi, 5}})
        correct_num = size(fluor{mismatched_data{itemi, 5}}, 3);
    else
        correct_num = 0;
    end
    duration_new = duration;
    duration_new.(mismatched_data{itemi, 4}) = duration.(mismatched_data{itemi, 4})(1:correct_num);

    % Save over previous 
    duration_places = duration_new;
    save(['Y:\Sarah\Analysis\Experiments\Random Motorized Treadmill\behavior\spontaneous\segmented behavior periods\' mismatched_data{itemi, 1} '\' mismatched_data{itemi, 2} '\duration_places_' mismatched_data{itemi, 3} '.mat'], 'duration_places')

end 

%From first round:
% '1087'	'011222'	'motorized'	'11'	180
% --> get rid of the last instance, the fluoresence was short by ~20 frames
% load('Y:\Sarah\Analysis\Experiments\Random Motorized Treadmill\behavior\motorized\period instances table format\1087\011222\all_periods_11.mat');
% all_periods.duration_place{180} = all_periods.duration_place{180}(1:end-1);
% save('Y:\Sarah\Analysis\Experiments\Random Motorized Treadmill\behavior\motorized\period instances table format\1087\011222\all_periods_11.mat', 'all_periods');
% 
% %'1088'	'010522'	'motorized'	'11'	180 
% % --> same as above
% load('Y:\Sarah\Analysis\Experiments\Random Motorized Treadmill\behavior\motorized\period instances table format\1088\010522\all_periods_11.mat');
% all_periods.duration_place{180} = all_periods.duration_place{180}(1:end-1);
% save('Y:\Sarah\Analysis\Experiments\Random Motorized Treadmill\behavior\motorized\period instances table format\1088\010522\all_periods_11.mat', 'all_periods');
% 
% % '1088'	'011322'	'motorized'	'06'	180
% % short by 2 frames; remove last 2 instances. 
% load('Y:\Sarah\Analysis\Experiments\Random Motorized Treadmill\behavior\motorized\period instances table format\1088\011322\all_periods_06.mat');
% all_periods.duration_place{180} = all_periods.duration_place{180}(1:end-2);
% load('Y:\Sarah\Analysis\Experiments\Random Motorized Treadmill\behavior\motorized\period instances table format\1088\011322\all_periods_06.mat', 'all_periods');


%% Motorized 

% Always clear loop list first. 
if isfield(parameters, 'loop_list')
parameters = rmfield(parameters,'loop_list');
end

% Iterators
parameters.loop_list.iterators = {'mouse', {'loop_variables.mice_all(:).name'}, 'mouse_iterator'; 
               'day', {'loop_variables.mice_all(', 'mouse_iterator', ').days(:).name'}, 'day_iterator';
               'stack', {'loop_variables.mice_all(',  'mouse_iterator', ').days(', 'day_iterator', ').stacks'}, 'stack_iterator'};

parameters.concatenation_level = 'day';
parameters.concatDim = 1;
parameters.concatenate_across_cells = false;

% Inputs
parameters.loop_list.things_to_load.data.dir = {[parameters.dir_exper 'behavior\motorized\period instances table format\'], 'mouse', '\', 'day', '\'};
parameters.loop_list.things_to_load.data.filename= {'all_periods_', 'stack', '.mat'};
parameters.loop_list.things_to_load.data.variable= {'all_periods.duration_place'}; 
parameters.loop_list.things_to_load.data.level = 'stack';

% Outputs
parameters.loop_list.things_to_save.concatenated_data.dir = {[parameters.dir_exper 'behavior\duration place concatenated\motorized\'], 'mouse', '\'};
parameters.loop_list.things_to_save.concatenated_data.filename= {'all_duration_place.mat'};
parameters.loop_list.things_to_save.concatenated_data.variable= {'all_duration_place'}; 
parameters.loop_list.things_to_save.concatenated_data.level = 'mouse';

RunAnalysis({@ConcatenateData}, parameters);

%% Spontaneous 
% Always clear loop list first. 
if isfield(parameters, 'loop_list')
parameters = rmfield(parameters,'loop_list');
end

% Iterators
parameters.loop_list.iterators = {'mouse', {'loop_variables.mice_all(:).name'}, 'mouse_iterator'; 
               'period', {'loop_variables.periods'}, 'period_iterator';
               'day', {'loop_variables.mice_all(', 'mouse_iterator', ').days(:).name'}, 'day_iterator';
               'stack', {'loop_variables.mice_all(',  'mouse_iterator', ').days(', 'day_iterator', ').spontaneous'}, 'stack_iterator';
               };

parameters.concatenation_level = 'day';
parameters.concatDim = 1;
parameters.concatenate_across_cells = false;

% Inputs
parameters.loop_list.things_to_load.data.dir = {[parameters.dir_exper 'behavior\spontaneous\segmented behavior periods\'], 'mouse', '\', 'day', '\'};
parameters.loop_list.things_to_load.data.filename= {'duration_places_', 'stack', '.mat'};
parameters.loop_list.things_to_load.data.variable= {'duration_places.', 'period'}; 
parameters.loop_list.things_to_load.data.level = 'stack';

% Outputs
parameters.loop_list.things_to_save.concatenated_data.dir = {[parameters.dir_exper 'behavior\duration place concatenated\spontaneous\'], 'mouse', '\'};
parameters.loop_list.things_to_save.concatenated_data.filename= {'all_duration_place_', 'period', '.mat'};
parameters.loop_list.things_to_save.concatenated_data.variable= {'all_duration_place'}; 
parameters.loop_list.things_to_save.concatenated_data.level = 'period';

RunAnalysis({@ConcatenateData}, parameters);

%% Convert spontaneous into cell form 
% with empty placeholder cells for the other 4 periods

% Always clear loop list first. 
if isfield(parameters, 'loop_list')
parameters = rmfield(parameters,'loop_list');
end

% Iterators
parameters.loop_list.iterators = {'mouse', {'loop_variables.mice_all(:).name'}, 'mouse_iterator'; 
               };

parameters.evaluation_instructions = {{'holder = [{parameters.rest} ; {parameters.walk}];' ...
                                        'data_evaluated = [holder; cell(4, 1)];'}};
        
% Inputs
% rest
parameters.loop_list.things_to_load.rest.dir = {[parameters.dir_exper 'behavior\duration place concatenated\spontaneous\'], 'mouse', '\'};
parameters.loop_list.things_to_load.rest.filename= {'all_duration_place_rest.mat'};
parameters.loop_list.things_to_load.rest.variable= {'all_duration_place'}; 
parameters.loop_list.things_to_load.rest.level = 'mouse';
% walk
parameters.loop_list.things_to_load.walk.dir = {[parameters.dir_exper 'behavior\duration place concatenated\spontaneous\'], 'mouse', '\'};
parameters.loop_list.things_to_load.walk.filename= {'all_duration_place_walk.mat'};
parameters.loop_list.things_to_load.walk.variable= {'all_duration_place'}; 
parameters.loop_list.things_to_load.walk.level = 'mouse';

% Outputs
parameters.loop_list.things_to_save.data_evaluated.dir = {[parameters.dir_exper 'behavior\duration place concatenated\spontaneous\'], 'mouse', '\'};
parameters.loop_list.things_to_save.data_evaluated.filename= {'all_duration_place.mat'};
parameters.loop_list.things_to_save.data_evaluated.variable= {'all_duration_place'}; 
parameters.loop_list.things_to_save.data_evaluated.level = 'mouse';

RunAnalysis({@EvaluateOnData}, parameters);


%% Put both conditions together in same cell array 

% Always clear loop list first. 
if isfield(parameters, 'loop_list')
parameters = rmfield(parameters,'loop_list');
end

parameters.loop_list.iterators = {'mouse', {'loop_variables.mice_all(:).name'}, 'mouse_iterator'; 
               'condition', 'loop_variables.conditions', 'condition_iterator';
                };

parameters.concatenation_level = 'condition';
parameters.concatDim = 1;
parameters.concatenate_across_cells = true;

% Inputs
parameters.loop_list.things_to_load.data.dir = {[parameters.dir_exper 'behavior\duration place concatenated\'], '\', 'condition', '\', 'mouse', '\'};
parameters.loop_list.things_to_load.data.filename= {'all_duration_place.mat'};
parameters.loop_list.things_to_load.data.variable= {'all_duration_place'}; 
parameters.loop_list.things_to_load.data.level = 'condition';

% Outputs
parameters.loop_list.things_to_save.concatenated_data.dir = {[parameters.dir_exper 'behavior\duration place concatenated\both conditions\'], 'mouse', '\'};
parameters.loop_list.things_to_save.concatenated_data.filename= {'all_duration_place.mat'};
parameters.loop_list.things_to_save.concatenated_data.variable= {'all_duration_place'}; 
parameters.loop_list.things_to_save.concatenated_data.level = 'mouse';

RunAnalysis({@ConcatenateData}, parameters);

parameters.concatenate_across_cells = false;

%%  Make duration_place for fluorescence PLSR 
% (Single time point)
% Always clear loop list first. 
if isfield(parameters, 'loop_list')
parameters = rmfield(parameters,'loop_list');
end

% Iterators
parameters.loop_list.iterators = {'mouse', {'loop_variables.mice_all(:).name'}, 'mouse_iterator'; 
               'period', {'loop_variables.periods_bothConditions'}, 'period_iterator';
               };

parameters.window_length = 20; % window length in time points
parameters.fps = 20; % Frame per second of the fluorescence

% Input
parameters.loop_list.things_to_load.data.dir = {[parameters.dir_exper 'behavior\duration place concatenated\both conditions\'], 'mouse', '\'};
parameters.loop_list.things_to_load.data.filename= {'all_duration_place.mat'};
parameters.loop_list.things_to_load.data.variable= {'all_duration_place{', 'period_iterator', '}'}; 
parameters.loop_list.things_to_load.data.level = 'mouse';

% Output
parameters.loop_list.things_to_save.duration_place_new.dir = {[parameters.dir_exper 'behavior\duration place concatenated\duration place for fluorescence PLSR\'], 'mouse', '\'};
parameters.loop_list.things_to_save.duration_place_new.filename= {'duration_place_forFluorescence.mat'};
parameters.loop_list.things_to_save.duration_place_new.variable= {'duration_place_forFluorescence{', 'period_iterator', ', 1}'}; 
parameters.loop_list.things_to_save.duration_place_new.level = 'mouse';

RunAnalysis({@DurationPlaceForFluorescencePLSR}, parameters);