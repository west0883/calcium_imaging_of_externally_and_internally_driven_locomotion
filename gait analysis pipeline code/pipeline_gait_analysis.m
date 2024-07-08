% pipeline_gait_analysis.m
% Sarah West
% 7/3/23

% Uses extracted and behavior-segmented traces of paw/body part velocities from
% DeepLabCut to find locomotion strides. Then uses those strides to run
% gait analysis.

% NOTE: From DLC, positive x is RIGHT, positive y is DOWN

%% Initial Setup  
% Put all needed paramters in a structure called "parameters", which you
% can then easily feed into your functions. 
% Use correlations, Fisher transformed, mean removed within mice (mean
% removed for at least the cases when you aren't using mice as response
% variables).

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
parameters.mice_all(6).days = parameters.mice_all(6).days([1:7 9:10]); %:7); %(6); %[1 3:end]); % Remove day 2/13/24   (2/14/23) for mouse 1107

% Other parameters
parameters.digitNumber = 2;
parameters.yDim = 256;
parameters.xDim = 256;
parameters.number_of_sources = 32; 
parameters.indices = find(tril(ones(parameters.number_of_sources), -1));
parameters.fps = 20; % converted fps

% Load periods_nametable.m for motorized & spontaneous. Concatenate
% together

load([parameters.dir_exper 'periods_nametable.mat']);
periods_motorized = periods;
  
load([parameters.dir_exper 'periods_nametable_spontaneous.mat']);
periods_spontaneous = periods;

parameters.periods = [periods_motorized; periods_spontaneous];

clear periods periods_motorized periods_spontaneous;

% Names of all continuous variables.
parameters.continuous_variable_names = {'speed', 'accel', 'duration', 'pupil_diameter', 'tail', 'nose', 'FL', 'HL'};


parameters.load_abort_flag = false;

% list of body parts to use
body_parts_all(1).name = 'FL';
body_parts_all(1).directions(1).name = 'total_magnitude';
body_parts_all(1).directions(2).name  = 'x';
body_parts_all(2).name = 'HL';
body_parts_all(2).directions(2).name  = 'x';
body_parts_all(2).directions(1).name  = 'total_magnitude';
body_parts_all(3).name = 'tail';
body_parts_all(3).directions(2).name  = 'y';
body_parts_all(3).directions(1).name  = 'total_magnitude';
body_parts_all(4).name = 'nose';
body_parts_all(4).directions(1).name  = 'total_magnitude';

% list of body parts to use
body_parts_all_position(1).name = 'FL';
body_parts_all_position(1).directions(1).name  = 'x';
body_parts_all_position(1).directions(2).name  = 'y';
body_parts_all_position(2).name = 'HL';
body_parts_all_position(2).directions(1).name  = 'x';
body_parts_all_position(2).directions(2).name  = 'y';
body_parts_all_position(3).name = 'tail';
body_parts_all_position(3).directions(1).name  = 'x';
body_parts_all_position(3).directions(2).name  = 'y';

parameters.body_parts_all = body_parts_all(1:3);
parameters.body_parts_all_position = body_parts_all_position;
parameters.paws = {'FL', 'HL', 'tail'};
% Put relevant variables into loop_variables.
parameters.loop_variables.mice_all = parameters.mice_all;
parameters.loop_variables.periods = parameters.periods.condition(1:194); % Don't include full_onset & full_offset 
parameters.loop_variables.conditions = {'motorized'; 'spontaneous'};
parameters.loop_variables.conditions_stack_locations = {'stacks'; 'spontaneous'};
parameters.loop_variables.conditions_stack_locations_long = {'spontaneous', 'stacks', 'stacks', 'stacks', 'stacks'}; 
parameters.loop_variables.variable_type = {'response variables', 'correlations'};
parameters.loop_variables.paws = parameters.paws;
parameters.loop_variables.body_parts = {'FL', 'HL', 'tail', 'nose'}; % {'FR', 'FL', 'HL', 'tail', 'nose', 'eye'};
parameters.loop_variables.body_parts_all = parameters.body_parts_all;
parameters.loop_variables.body_parts_all_position = parameters.body_parts_all_position;
parameters.loop_variables.velocity_directions = {'x', 'y', 'total_magnitude'}; %'total_angle'};
parameters.loop_variables.type_tags = {'longWalk_spontaneous', 'longWalk_motorized1600', 'longWalk_motorized2000', 'longWalk_motorized2400', 'longWalk_motorized2800'}; % No "allPeriods" anymore; For concatenated majority of peiods and the long versions of motorized & spontaneous walk
parameters.loop_variables.type_tags2 = {'longWalk_spontaneous', 'longWalk_motorized_1600', 'longWalk_motorized_2000', 'longWalk_motorized_2400', 'longWalk_motorized2800'};
parameters.loop_variables.motorSpeeds = {'1600', '2000', '2400', '2800'};
parameters.loop_variables.periods_withLongs = [parameters.periods.condition(1:194); {'walkLong_spon'}; {'walkLong_1600'}; {'walkLong_2000'}; {'walkLong_2400'}; {'walkLong_2800'}];
parameters.loop_variables.periods_longsOnly = [{'walkLong_spon'}; {'walkLong_1600'}; {'walkLong_2000'}; {'walkLong_2400'}; {'walkLong_2800'}];
parameters.loop_variables.peak_depression = {'depression'}; % {'peak', 'depression'};
parameters.loop_variables.velocity_directions_sublist = {'total_magnitude'};%{ 'y', 'total_magnitude'}; % 'y', % 'total_angle'
parameters.loop_variables.paws_sublist =  {'HL', 'tail'};
parameters.loop_variables.segmentation_types =  { 'stride segmentations from own x depressions', 'stride segmentations from FL x depressions'};
parameters.loop_variables.position_directions = {'x', 'y'};

parameters.average_and_std_together = false;

%% Segment body velocities based on long periods: Motorized
% (similar to paw_velocity_pipeline_code.m)
if isfield(parameters, 'loop_list')
parameters = rmfield(parameters,'loop_list');
end

% Iterators
parameters.loop_list.iterators = {'mouse', {'loop_variables.mice_all(:).name'}, 'mouse_iterator'; 
               'day', {'loop_variables.mice_all(', 'mouse_iterator', ').days(:).name'}, 'day_iterator';
                   'stack', {'loop_variables.mice_all(',  'mouse_iterator', ').days(', 'day_iterator', ').stacks'}, 'stack_iterator';
                   'body_part', {'loop_variables.body_parts_all(:).name'}, 'body_part_iterator';
                   'velocity_direction', {'loop_variables.body_parts_all(', 'body_part_iterator', ').directions(:).name'}, 'direction_iterator';
                   'motorSpeed', {'loop_variables.motorSpeeds'}, 'motorSpeed_iterator';
                   };

% Skip any files that don't exist (spontaneous or problem files)
parameters.load_abort_flag = true; 

% Dimension of different time range pairs.
parameters.rangePairs = 1; 

% 
parameters.segmentDim = 1;
parameters.concatDim = 2;

% Are the segmentations the same length? (If false, will put outputs into
% cell array)
parameters.uniformSegments = false;

% Input values. 
% Extracted timeseries.
parameters.loop_list.things_to_load.timeseries.dir = {[parameters.dir_exper 'behavior\body\normalized with std only\paw velocity normalized\'], 'mouse', '\', 'day', '\'};
parameters.loop_list.things_to_load.timeseries.filename= {'velocity', 'stack', '.mat'};
parameters.loop_list.things_to_load.timeseries.variable= {'velocity.', 'body_part', '.', 'velocity_direction'}; 
parameters.loop_list.things_to_load.timeseries.level = 'stack';
% Time ranges
parameters.loop_list.things_to_load.time_ranges.dir = {[parameters.dir_exper 'behavior\motorized\period instances\'], 'mouse', '\', 'day', '\'};
parameters.loop_list.things_to_load.time_ranges.filename= {'long_periods_', 'stack', '.mat'};
parameters.loop_list.things_to_load.time_ranges.variable= {'long_periods.walk_', 'motorSpeed'}; 
parameters.loop_list.things_to_load.time_ranges.level = 'stack';

% Output Values
parameters.loop_list.things_to_save.segmented_timeseries.dir = {[parameters.dir_exper 'behavior\body\normalized with std only\segmented velocities\'], 'body_part', '\', 'velocity_direction', '\motorized\', 'mouse', '\', 'day', '\'};
parameters.loop_list.things_to_save.segmented_timeseries.filename= {'segmented_timeseries_longPeriods_walk_', 'motorSpeed', '_', 'stack', '.mat'};
parameters.loop_list.things_to_save.segmented_timeseries.variable= {'segmented_timeseries'}; 
parameters.loop_list.things_to_save.segmented_timeseries.level = 'motorSpeed';

RunAnalysis({@SegmentTimeseriesData}, parameters);

%% Segment body velocities based on long periods: Spontaneous
% (similar to paw_velocity_pipeline_code.m)
if isfield(parameters, 'loop_list')
parameters = rmfield(parameters,'loop_list');
end

% Iterators
parameters.loop_list.iterators = {'mouse', {'loop_variables.mice_all(:).name'}, 'mouse_iterator'; 
                    'day', {'loop_variables.mice_all(', 'mouse_iterator', ').days(:).name'}, 'day_iterator';
                   'stack', {'loop_variables.mice_all(',  'mouse_iterator', ').days(', 'day_iterator', ').spontaneous'}, 'stack_iterator';
                    'body_part', {'loop_variables.body_parts_all(:).name'}, 'body_part_iterator';
                    'velocity_direction', {'loop_variables.body_parts_all(', 'body_part_iterator', ').directions(:).name'}, 'velocity_direction_iterator';
                    };

% Skip any files that don't exist (spontaneous or problem files)
parameters.load_abort_flag = true; 

% Dimension of different time range pairs.
parameters.rangePairs = 1; 

% 
parameters.segmentDim = 1;
parameters.concatDim = 2;

% Are the segmentations the same length? (If false, will put outputs into
% cell array)
parameters.uniformSegments = false;

% Input values. 
% Extracted timeseries.
parameters.loop_list.things_to_load.timeseries.dir = {[parameters.dir_exper 'behavior\body\normalized with std only\paw velocity normalized\'], 'mouse', '\', 'day', '\'};
parameters.loop_list.things_to_load.timeseries.filename= {'velocity', 'stack', '.mat'};
parameters.loop_list.things_to_load.timeseries.variable= {'velocity.', 'body_part', '.', 'velocity_direction'}; 
parameters.loop_list.things_to_load.timeseries.level = 'stack';
% Time ranges
parameters.loop_list.things_to_load.time_ranges.dir = {[parameters.dir_exper 'behavior\spontaneous\segmented behavior periods\'], 'mouse', '\', 'day', '\'};
parameters.loop_list.things_to_load.time_ranges.filename= {'long_periods_', 'stack', '.mat'};
parameters.loop_list.things_to_load.time_ranges.variable= {'long_periods.walk'}; 
parameters.loop_list.things_to_load.time_ranges.level = 'stack';

% Output Values
parameters.loop_list.things_to_save.segmented_timeseries.dir = {[parameters.dir_exper 'behavior\body\normalized with std only\segmented velocities\'], 'body_part', '\', 'velocity_direction', '\spontaneous\', 'mouse', '\', 'day', '\'};
parameters.loop_list.things_to_save.segmented_timeseries.filename= {'segmented_timeseries_longPeriods_walk', '_', 'stack', '.mat'};
parameters.loop_list.things_to_save.segmented_timeseries.variable= {'segmented_timeseries'}; 
parameters.loop_list.things_to_save.segmented_timeseries.level = 'velocity_direction';

RunAnalysis({@SegmentTimeseriesData}, parameters);

%% Concatenate long periods -- motorized

if isfield(parameters, 'loop_list')
parameters = rmfield(parameters,'loop_list');
end

% Iterators
parameters.loop_list.iterators = {
               'body_part', {'loop_variables.body_parts_all(:).name'}, 'body_part_iterator';
                'velocity_direction', {'loop_variables.body_parts_all(', 'body_part_iterator', ').directions(:).name'}, 'velocity_direction_iterator';
               'mouse', {'loop_variables.mice_all(:).name'}, 'mouse_iterator'; 
               'motorSpeed', {'loop_variables.motorSpeeds'}, 'motorSpeed_iterator';
               'day', {'loop_variables.mice_all(', 'mouse_iterator', ').days(:).name'}, 'day_iterator';
               'stack', {'loop_variables.mice_all(',  'mouse_iterator', ').days(', 'day_iterator', ').stacks'}, 'stack_iterator';
                    };

parameters.concatDim = 1;
%parameters.concatenation_level = 'day';
parameters.concatenate_across_cells = true;

% Input
parameters.loop_list.things_to_load.data.dir = {[parameters.dir_exper 'behavior\body\normalized with std only\segmented velocities\'], 'body_part', '\', 'velocity_direction', '\motorized',  '\' 'mouse', '\', 'day', '\'};
parameters.loop_list.things_to_load.data.filename = {'segmented_timeseries_longPeriods_walk_', 'motorSpeed', '_', 'stack', '.mat'};
parameters.loop_list.things_to_load.data.variable = {'segmented_timeseries'}; 
parameters.loop_list.things_to_load.data.level = 'stack';

% Output
parameters.loop_list.things_to_save.concatenated_data.dir = {[parameters.dir_exper 'behavior\body\normalized with std only\concatenated velocity\'], 'body_part', '\', 'velocity_direction', '\motorized\', 'mouse', '\'};
parameters.loop_list.things_to_save.concatenated_data.filename = {'concatenated_velocity_longPeriods_walk_', 'motorSpeed', '.mat'};
parameters.loop_list.things_to_save.concatenated_data.variable = {'velocity_all'}; 
parameters.loop_list.things_to_save.concatenated_data.level = 'motorSpeed';

RunAnalysis({@ConcatenateData}, parameters);

%% Concatenate long periods -- spontaneous 

if isfield(parameters, 'loop_list')
parameters = rmfield(parameters,'loop_list');
end

% Iterators
parameters.loop_list.iterators = {
               'body_part', {'loop_variables.body_parts_all(:).name'}, 'body_part_iterator';
               'velocity_direction', {'loop_variables.body_parts_all(', 'body_part_iterator', ').directions(:).name'}, 'velocity_direction_iterator';
               'mouse', {'loop_variables.mice_all(:).name'}, 'mouse_iterator'; 
               'day', {'loop_variables.mice_all(', 'mouse_iterator', ').days(:).name'}, 'day_iterator';
               'stack', {'loop_variables.mice_all(',  'mouse_iterator', ').days(', 'day_iterator', ').spontaneous'}, 'stack_iterator';
                    };

parameters.concatDim = 1;
%parameters.concatenation_level = 'day';
parameters.concatenate_across_cells = true;

% Input
parameters.loop_list.things_to_load.data.dir = {[parameters.dir_exper 'behavior\body\normalized with std only\segmented velocities\'], 'body_part', '\', 'velocity_direction', '\spontaneous',  '\' 'mouse', '\', 'day', '\'};
parameters.loop_list.things_to_load.data.filename = {'segmented_timeseries_longPeriods_walk_', 'stack', '.mat'};
parameters.loop_list.things_to_load.data.variable = {'segmented_timeseries'}; 
parameters.loop_list.things_to_load.data.level = 'stack';

% Output
parameters.loop_list.things_to_save.concatenated_data.dir = {[parameters.dir_exper 'behavior\body\normalized with std only\concatenated velocity\'], 'body_part', '\', 'velocity_direction', '\spontaneous\', 'mouse', '\'};
parameters.loop_list.things_to_save.concatenated_data.filename = {'concatenated_velocity_longPeriods_walk.mat'};
parameters.loop_list.things_to_save.concatenated_data.variable = {'velocity_all'}; 
parameters.loop_list.things_to_save.concatenated_data.level = 'mouse';

RunAnalysis({@ConcatenateData}, parameters);

%% Find peaks & depressions in x velocity traces for both left paws: 
% % All periods besides long rest & walk, motorized and spontaneous
% 
% % Using peakdet.m 
% % this worked well; a is a paw velocity trace:
% % b = a - mean(a);
% % [pks,dep,pid,did] = peakdet(b, 0.05 , 'zero', 2);
% 
% % Should use x-velocity only to find strides
% 
% if isfield(parameters, 'loop_list')
% parameters = rmfield(parameters,'loop_list');
% end
% 
% % Is so you can use a single loop for calculations. 
% parameters.loop_list.iterators = {
%                'paw', {'loop_variables.paws'}, 'paw_iterator';
%                'mouse', {'loop_variables.mice_all(:).name'}, 'mouse_iterator';
%                'period', {'loop_variables.periods'}, 'period_iterator';
%                 };
% 
% parameters.timeDim = 1;
% parameters.instanceDim = 2;
% % Minnimum height from mean to count as a peak
% parameters.peakMinHeight = 0.1;
% % Minnimum time point separation between peaks to count as different peaks
% % (4 = 5 Hz stride)
% parameters.peakMinSeparation = 5; 
% parameters.instancesAsCells = false;
% 
% % Inputs
% parameters.loop_list.things_to_load.data.dir = {[parameters.dir_exper 'behavior\body\normalized with std only\concatenated velocity\'], 'paw', '\', 'x', '\both conditions\', 'mouse', '\'};
% parameters.loop_list.things_to_load.data.filename = {'concatenated_velocity_all_periods.mat'};
% parameters.loop_list.things_to_load.data.variable = {'velocity_all{', 'period_iterator', '}'}; 
% parameters.loop_list.things_to_load.data.level = 'mouse';
% 
% % Outputs
% % peaks
% parameters.loop_list.things_to_save.peaks.dir = {[parameters.dir_exper 'behavior\gait analysis\normalized with std only\x peaks\all periods\'],'paw', '\', 'mouse', '\'};
% parameters.loop_list.things_to_save.peaks.filename = {'x_peaks.mat'};
% parameters.loop_list.things_to_save.peaks.variable = {'x_peaks{', 'period_iterator', ', 1}'}; 
% parameters.loop_list.things_to_save.peaks.level = 'mouse';
% % velocities segmented from depressions
% parameters.loop_list.things_to_save.segmentations_peak.dir = {[parameters.dir_exper 'behavior\gait analysis\normalized with std only\stride segmentations\from peaks\all periods\'],'paw', '\', 'mouse', '\'};
% parameters.loop_list.things_to_save.segmentations_peak.filename = {'stride_segmentations_allPeriods.mat'};
% parameters.loop_list.things_to_save.segmentations_peak.variable = {'stride_segmentations_peak{', 'period_iterator', ', 1}'}; 
% parameters.loop_list.things_to_save.segmentations_peak.level = 'mouse';
% % segmented from depressions
% parameters.loop_list.things_to_save.segmentations_depression.dir = {[parameters.dir_exper 'behavior\gait analysis\normalized with std only\stride segmentations\from depressions\all periods\'],'paw', '\', 'mouse', '\'};
% parameters.loop_list.things_to_save.segmentations_depression.filename = {'stride_segmentations_allPeriods.mat'};
% parameters.loop_list.things_to_save.segmentations_depression.variable = {'stride_segmentations_depression{', 'period_iterator', ', 1}'}; 
% parameters.loop_list.things_to_save.segmentations_depression.level = 'mouse';
% 
% RunAnalysis({@FindStrides}, parameters);

%% Spontaneous long walk: find peaks & depressions in x or y velocity traces for FLx, HLx, tail y
% FL only
if isfield(parameters, 'loop_list')
parameters = rmfield(parameters,'loop_list');
end

% make body_parts all specific to x or y directoins
parameters.loop_variables.body_parts_all(1).directions = parameters.loop_variables.body_parts_all(1).directions(2);
parameters.loop_variables.body_parts_all(2).directions = parameters.loop_variables.body_parts_all(2).directions(2);
parameters.loop_variables.body_parts_all(3).directions = parameters.loop_variables.body_parts_all(3).directions(2);

% Is so you can use a single loop for calculations. 
parameters.loop_list.iterators = {
               'body_part', {'loop_variables.body_parts_all(:).name'}, 'body_part_iterator';
               'velocity_direction', {'loop_variables.body_parts_all(', 'body_part_iterator', ').directions(:).name'}, 'velocity_direction_iterator';
               'mouse', {'loop_variables.mice_all(:).name'}, 'mouse_iterator';
                };

parameters.timeDim = 1;
parameters.instanceDim = 2;
parameters.peakMinHeight = 0.1;
% Minnimum time point separation between peaks to count as different peaks
% (5 = 4 Hz stride)
parameters.peakMinSeparation = 5; 
parameters.instancesAsCells = true; 

% Inputs
parameters.loop_list.things_to_load.data.dir = {[parameters.dir_exper 'behavior\body\normalized with std only\concatenated velocity\'], 'body_part', '\', 'velocity_direction', '\spontaneous\', 'mouse', '\'};
parameters.loop_list.things_to_load.data.filename = {'concatenated_velocity_longPeriods_walk.mat'};
parameters.loop_list.things_to_load.data.variable = {'velocity_all'}; 
parameters.loop_list.things_to_load.data.level = 'mouse';

% Outputs
% peaks
parameters.loop_list.things_to_save.peaks.dir = {[parameters.dir_exper 'behavior\gait analysis\normalized with std only\x peaks\all periods\'],'body_part', '\', 'mouse', '\'};
parameters.loop_list.things_to_save.peaks.filename = {'x_peaks_longWalk_spontaneous.mat'};
parameters.loop_list.things_to_save.peaks.variable = {'x_peaks'}; 
parameters.loop_list.things_to_save.peaks.level = 'mouse';
% velocities segmented into strides
parameters.loop_list.things_to_save.segmentations_peak.dir = {[parameters.dir_exper 'behavior\gait analysis\normalized with std only\stride segmentations\from peaks\all periods\'],'body_part', '\', 'mouse', '\'};
parameters.loop_list.things_to_save.segmentations_peak.filename = {'stride_segmentations_longWalk_spontaneous.mat'};
parameters.loop_list.things_to_save.segmentations_peak.variable = {'stride_segmentations_peak'}; 
parameters.loop_list.things_to_save.segmentations_peak.level = 'mouse';

parameters.loop_list.things_to_save.segmentations_depression.dir = {[parameters.dir_exper 'behavior\gait analysis\normalized with std only\stride segmentations\from depressions\all periods\'],'body_part', '\', 'mouse', '\'};
parameters.loop_list.things_to_save.segmentations_depression.filename = {'stride_segmentations_longWalk_spontaneous.mat'};
parameters.loop_list.things_to_save.segmentations_depression.variable = {'stride_segmentations_depression'}; 
parameters.loop_list.things_to_save.segmentations_depression.level = 'mouse';

RunAnalysis({@FindStrides}, parameters);

% reset body_parts_all
parameters.loop_variables.body_parts_all = parameters.body_parts_all;

%% Motorized long walk: find peaks & depressions in x or y velocity traces for FLx, HLx, tail y

if isfield(parameters, 'loop_list')
parameters = rmfield(parameters,'loop_list');
end

% make body_parts all specific to x or y directoins
parameters.loop_variables.body_parts_all(1).directions = parameters.loop_variables.body_parts_all(1).directions(2);
parameters.loop_variables.body_parts_all(2).directions = parameters.loop_variables.body_parts_all(2).directions(2);
parameters.loop_variables.body_parts_all(3).directions = parameters.loop_variables.body_parts_all(3).directions(2);

% Is so you can use a single loop for calculations. 
parameters.loop_list.iterators = {
               'body_part', {'loop_variables.body_parts_all(:).name'}, 'body_part_iterator';
               'velocity_direction', {'loop_variables.body_parts_all(', 'body_part_iterator', ').directions(:).name'}, 'velocity_direction_iterator';
               'mouse', {'loop_variables.mice_all(:).name'}, 'mouse_iterator';
               'motorSpeed', {'loop_variables.motorSpeeds'}, 'motorSpeed_iterator';
                };

parameters.timeDim = 1;
parameters.instanceDim = 2;
parameters.peakMinHeight = 0.1;
% Minnimum time point separation between peaks to count as different peaks
% (4 = 5 Hz stride)
parameters.peakMinSeparation = 5; 
parameters.instancesAsCells = true; 

% Inputs
parameters.loop_list.things_to_load.data.dir = {[parameters.dir_exper 'behavior\body\normalized with std only\concatenated velocity\'], 'body_part', '\', 'velocity_direction', '\motorized\', 'mouse', '\'};
parameters.loop_list.things_to_load.data.filename = {'concatenated_velocity_longPeriods_walk_', 'motorSpeed', '.mat'};
parameters.loop_list.things_to_load.data.variable = {'velocity_all'}; 
parameters.loop_list.things_to_load.data.level = 'motorSpeed';

% Outputs
% peaks
parameters.loop_list.things_to_save.peaks.dir = {[parameters.dir_exper 'behavior\gait analysis\normalized with std only\x peaks\all periods\'],'body_part', '\', 'mouse', '\'};
parameters.loop_list.things_to_save.peaks.filename = {'x_peaks_longWalk_motorized', 'motorSpeed', '.mat'};
parameters.loop_list.things_to_save.peaks.variable = {'x_peaks'}; 
parameters.loop_list.things_to_save.peaks.level = 'motorSpeed';
% velocities segmented into strides from peaks
parameters.loop_list.things_to_save.segmentations_peak.dir = {[parameters.dir_exper 'behavior\gait analysis\normalized with std only\stride segmentations\from peaks\all periods\'],'body_part', '\', 'mouse', '\'};
parameters.loop_list.things_to_save.segmentations_peak.filename = {'stride_segmentations_longWalk_motorized', 'motorSpeed', '.mat'};
parameters.loop_list.things_to_save.segmentations_peak.variable = {'stride_segmentations_peak'}; 
parameters.loop_list.things_to_save.segmentations_peak.level = 'motorSpeed';
% velocities segmented into strides from depressions
parameters.loop_list.things_to_save.segmentations_depression.dir = {[parameters.dir_exper 'behavior\gait analysis\normalized with std only\stride segmentations\from depressions\all periods\'],'body_part', '\', 'mouse', '\'};
parameters.loop_list.things_to_save.segmentations_depression.filename = {'stride_segmentations_longWalk_motorized', 'motorSpeed', '.mat'};
parameters.loop_list.things_to_save.segmentations_depression.variable = {'stride_segmentations_depression'}; 
parameters.loop_list.things_to_save.segmentations_depression.level = 'motorSpeed';

RunAnalysis({@FindStrides}, parameters);

% reset body_parts_all
parameters.loop_variables.body_parts_all = parameters.body_parts_all;


%% Concatenate segmentations of long walk motorized, spontaneous to all periods
% So you can run next steps all together

if isfield(parameters, 'loop_list')
parameters = rmfield(parameters,'loop_list');
end

% Is so you can use a single loop for calculations. 
parameters.loop_list.iterators = {
               'paw', {'loop_variables.paws'}, 'paw_iterator';
               'mouse', {'loop_variables.mice_all(:).name'}, 'mouse_iterator';
               'peak_depression', {'loop_variables.peak_depression'}, 'peak_depression_iterator';
               'type_tag',  {'loop_variables.type_tags'}, 'type_tag_iterator';
               };

parameters.concatDim = 1;
parameters.concatenation_level = 'type_tag';
parameters.concatenate_across_cells = true;

% Inputs
parameters.loop_list.things_to_load.data.dir = {[parameters.dir_exper 'behavior\gait analysis\normalized with std only\stride segmentations\from '], 'peak_depression', 's\all periods\','paw', '\', 'mouse', '\'};
parameters.loop_list.things_to_load.data.filename = {'stride_segmentations_', 'type_tag', '.mat'};
parameters.loop_list.things_to_load.data.variable = {'stride_segmentations_' 'peak_depression'}; 
parameters.loop_list.things_to_load.data.level = 'type_tag';
% Outputs
parameters.loop_list.things_to_save.concatenated_data.dir = {[parameters.dir_exper 'behavior\gait analysis\normalized with std only\stride segmentations\from '], 'peak_depression', 's\concatenated periods\','paw', '\', 'mouse', '\'};
parameters.loop_list.things_to_save.concatenated_data.filename = {'stride_segmentations', '.mat'};
parameters.loop_list.things_to_save.concatenated_data.variable = {'stride_segmentations'}; 
parameters.loop_list.things_to_save.concatenated_data.level = 'peak_depression';

RunAnalysis({@ConcatenateData}, parameters);

%% For each period, plot, resample, take means & standard deviations 

if isfield(parameters, 'loop_list')
parameters = rmfield(parameters,'loop_list');
end

parameters.loop_variables.paws = parameters.paws;
% Is so you can use a single loop for calculations. 
parameters.loop_list.iterators = {
               'body_part', {'loop_variables.body_parts_all(:).name'}, 'body_part_iterator';
               'peak_depression', {'loop_variables.peak_depression'}, 'peak_depression_iterator';
               'mouse', {'loop_variables.mice_all(:).name'}, 'mouse_iterator';
               'period', {'loop_variables.periods_longsOnly'}, 'period_iterator';
               }; 

% close each figure after saving it
parameters.closeFigures = true;
% the number of timepoints to resample each stride velocity segment to.
parameters.resampleLength = 10; % to 0.5 s = 10 time points

% Inputs
parameters.loop_list.things_to_load.data.dir = {[parameters.dir_exper 'behavior\gait analysis\normalized with std only\stride segmentations\from '], 'peak_depression', 's\concatenated periods\','body_part', '\', 'mouse', '\'};
parameters.loop_list.things_to_load.data.filename = {'stride_segmentations', '.mat'};
parameters.loop_list.things_to_load.data.variable = {'stride_segmentations{', 'period_iterator', '}'}; 
parameters.loop_list.things_to_load.data.level = 'period';

% Outputs
% stride segmentations reformatted so strides are all on same cell level
parameters.loop_list.things_to_save.segmentations_together.dir = {[parameters.dir_exper 'behavior\gait analysis\normalized with std only\stride segmentations\from '], 'peak_depression', 's\concatenated periods\','body_part', '\', 'mouse', '\'};
parameters.loop_list.things_to_save.segmentations_together.filename = {'stride_segmentations_together', '.mat'};
parameters.loop_list.things_to_save.segmentations_together.variable = {'stride_segmentations_together{', 'period_iterator', ', 1}'}; 
parameters.loop_list.things_to_save.segmentations_together.level = 'mouse';
% resampled segmentations
parameters.loop_list.things_to_save.resampled.dir = {[parameters.dir_exper 'behavior\gait analysis\normalized with std only\stride segmentations\from '], 'peak_depression', 's\concatenated periods\','body_part', '\', 'mouse', '\'};
parameters.loop_list.things_to_save.resampled.filename = {'stride_segmentations_resampled', '.mat'};
parameters.loop_list.things_to_save.resampled.variable = {'stride_segmentations_resampled{', 'period_iterator', ', 1}'}; 
parameters.loop_list.things_to_save.resampled.level = 'mouse';
% mean 
parameters.loop_list.things_to_save.average.dir = {[parameters.dir_exper 'behavior\gait analysis\normalized with std only\stride segmentations\from '], 'peak_depression', 's\concatenated periods\','body_part', '\', 'mouse', '\'};
parameters.loop_list.things_to_save.average.filename = {'stride_segmentations_average', '.mat'};
parameters.loop_list.things_to_save.average.variable = {'stride_segmentations_average{', 'period_iterator', ', 1}'}; 
parameters.loop_list.things_to_save.average.level = 'mouse';
% std
parameters.loop_list.things_to_save.std_dev.dir = {[parameters.dir_exper 'behavior\gait analysis\normalized with std only\stride segmentations\from '], 'peak_depression', 's\concatenated periods\','body_part', '\', 'mouse', '\'};
parameters.loop_list.things_to_save.std_dev.filename = {'stride_segmentations_std_dev', '.mat'};
parameters.loop_list.things_to_save.std_dev.variable = {'stride_segmentations_std_dev{', 'period_iterator', ', 1}'}; 
parameters.loop_list.things_to_save.std_dev.level = 'mouse';
% padded segmentations
parameters.loop_list.things_to_save.padded.dir = {[parameters.dir_exper 'behavior\gait analysis\normalized with std only\stride segmentations\from '], 'peak_depression', 's\concatenated periods\','body_part', '\', 'mouse', '\'};
parameters.loop_list.things_to_save.padded.filename = {'stride_segmentations_padded', '.mat'};
parameters.loop_list.things_to_save.padded.variable = {'stride_segmentations_padded{', 'period_iterator', ', 1}'}; 
parameters.loop_list.things_to_save.padded.level = 'mouse';
% padded mean 
parameters.loop_list.things_to_save.average_padded.dir = {[parameters.dir_exper 'behavior\gait analysis\normalized with std only\stride segmentations\from '], 'peak_depression', 's\concatenated periods\','body_part', '\', 'mouse', '\'};
parameters.loop_list.things_to_save.average_padded.filename = {'stride_segmentations_average_padded', '.mat'};
parameters.loop_list.things_to_save.average_padded.variable = {'stride_segmentations_average_padded{', 'period_iterator', ', 1}'}; 
parameters.loop_list.things_to_save.average_padded.level = 'mouse';
% padded std
parameters.loop_list.things_to_save.std_dev_padded.dir = {[parameters.dir_exper 'behavior\gait analysis\normalized with std only\stride segmentations\from '], 'peak_depression', 's\concatenated periods\','body_part', '\', 'mouse', '\'};
parameters.loop_list.things_to_save.std_dev_padded.filename = {'stride_segmentations_std_dev_padded', '.mat'};
parameters.loop_list.things_to_save.std_dev_padded.variable = {'stride_segmentations_std_dev_padded{', 'period_iterator', ', 1}'}; 
parameters.loop_list.things_to_save.std_dev_padded.level = 'mouse';

% figure: not-resampled segmentations 
parameters.loop_list.things_to_save.fig_segmentations_together.dir = {[parameters.dir_exper 'behavior\gait analysis\normalized with std only\stride segmentations\from '], 'peak_depression', 's\concatenated periods\','body_part', '\', 'mouse', '\'};
parameters.loop_list.things_to_save.fig_segmentations_together.filename = {'stride_segmentations_together_','period_iterator', '.fig'};
parameters.loop_list.things_to_save.fig_segmentations_together.variable = {'fig'}; 
parameters.loop_list.things_to_save.fig_segmentations_together.level = 'period';
% figure: resampled segmentations
parameters.loop_list.things_to_save.fig_resampled.dir = {[parameters.dir_exper 'behavior\gait analysis\normalized with std only\stride segmentations\from '], 'peak_depression', 's\concatenated periods\','body_part', '\', 'mouse', '\'};
parameters.loop_list.things_to_save.fig_resampled.filename = {'stride_segmentations_resampled_','period_iterator', '.fig'};
parameters.loop_list.things_to_save.fig_resampled.variable = {'fig'}; 
parameters.loop_list.things_to_save.fig_resampled.level = 'period';
% figure: mean and std
parameters.loop_list.things_to_save.fig_average.dir = {[parameters.dir_exper 'behavior\gait analysis\normalized with std only\stride segmentations\from '], 'peak_depression', 's\concatenated periods\','body_part', '\', 'mouse', '\'};
parameters.loop_list.things_to_save.fig_average.filename = {'stride_segmentations_average_','period_iterator', '.fig'};
parameters.loop_list.things_to_save.fig_average.variable = {'fig'}; 
parameters.loop_list.things_to_save.fig_average.level = 'period';
% figure: padded segmentations
parameters.loop_list.things_to_save.fig_padded.dir = {[parameters.dir_exper 'behavior\gait analysis\normalized with std only\stride segmentations\from '], 'peak_depression', 's\concatenated periods\','body_part', '\', 'mouse', '\'};
parameters.loop_list.things_to_save.fig_padded.filename = {'stride_segmentations_padded_','period_iterator', '.fig'};
parameters.loop_list.things_to_save.fig_padded.variable = {'fig'}; 
parameters.loop_list.things_to_save.fig_padded.level = 'period';
% figure: padded mean and std
parameters.loop_list.things_to_save.fig_average_padded.dir = {[parameters.dir_exper 'behavior\gait analysis\normalized with std only\stride segmentations\from '], 'peak_depression', 's\concatenated periods\','body_part', '\', 'mouse', '\'};
parameters.loop_list.things_to_save.fig_average_padded.filename = {'stride_segmentations_average_padded_','period_iterator', '.fig'};
parameters.loop_list.things_to_save.fig_average_padded.variable = {'fig'}; 
parameters.loop_list.things_to_save.fig_average_padded.level = 'period';

RunAnalysis({@GaitResampling}, parameters);
parameters.closeFigures = false;

% reset paws variable
parameters.loop_variables.paws = parameters.paws;

%% Plot all mices' averages of each period together

if isfield(parameters, 'loop_list')
parameters = rmfield(parameters,'loop_list');
end

% Is so you can use a single loop for calculations. 
parameters.loop_list.iterators = {
               'paw', {'loop_variables.paws'}, 'paw_iterator';
               'period', {'loop_variables.periods_longsOnly'}, 'period_iterator';
               'mouse', {'loop_variables.mice_all(:).name'}, 'mouse_iterator';
               }; 

parameters.instancesDim = 1;
parameters.ylimits = [-10 10];
parameters.mymap = flipud(hsv(7));

% Inputs
% mean 
parameters.loop_list.things_to_load.average.dir = {[parameters.dir_exper 'behavior\gait analysis\normalized with std only\stride segmentations\from depressions\concatenated periods\'],'paw', '\', 'mouse', '\'};
parameters.loop_list.things_to_load.average.filename = {'stride_segmentations_average', '.mat'};
parameters.loop_list.things_to_load.average.variable = {'stride_segmentations_average{', 'period_iterator', ', 1}'}; 
parameters.loop_list.things_to_load.average.level = 'mouse';
% std
parameters.loop_list.things_to_load.std_dev.dir = {[parameters.dir_exper 'behavior\gait analysis\normalized with std only\stride segmentations\from depressions\concatenated periods\'],'paw', '\', 'mouse', '\'};
parameters.loop_list.things_to_load.std_dev.filename = {'stride_segmentations_std_dev', '.mat'};
parameters.loop_list.things_to_load.std_dev.variable = {'stride_segmentations_std_dev{', 'period_iterator', ', 1}'}; 
parameters.loop_list.things_to_load.std_dev.level = 'mouse';
% resampled segmentations (to get the number of instances for standard error of the mean calculation)
parameters.loop_list.things_to_load.resampled.dir = {[parameters.dir_exper 'behavior\gait analysis\normalized with std only\stride segmentations\from depressions\concatenated periods\'],'paw', '\', 'mouse', '\'};
parameters.loop_list.things_to_load.resampled.filename = {'stride_segmentations_resampled', '.mat'};
parameters.loop_list.things_to_load.resampled.variable = {'stride_segmentations_resampled{', 'period_iterator', ', 1}'}; 
parameters.loop_list.things_to_load.resampled.level = 'mouse';

% Outputs
% figure
parameters.loop_list.things_to_save.fig.dir = {[parameters.dir_exper 'behavior\gait analysis\normalized with std only\stride segmentations\from depressions\all mice\overlays\'],'paw', '\'};
parameters.loop_list.things_to_save.fig.filename = {'overlay_','period_iterator','_', 'period', '.fig'};
parameters.loop_list.things_to_save.fig.variable = {'fig'}; 
parameters.loop_list.things_to_save.fig.level = 'period';

RunAnalysis({@PlotMiceStrideOverlays}, parameters);

close all;
%% Average each period across mice
% include m1100, aren't comparing across spon and motorized yet
% concatenate & average

if isfield(parameters, 'loop_list')
parameters = rmfield(parameters,'loop_list');
end

% Is so you can use a single loop for calculations. 
parameters.loop_list.iterators = {
               'paw', {'loop_variables.paws'}, 'paw_iterator';
               'period', {'loop_variables.periods_longsOnly'}, 'period_iterator';
               'mouse', {'loop_variables.mice_all(:).name'}, 'mouse_iterator';
               }; 

parameters.concatDim = 1;
parameters.concatenation_level = 'mouse';
parameters.averageDim = 1;

% Don't include mouse 1100 in spontaneous averages
parameters.evaluation_instructions = {{'period_iterator = parameters.values{strcmp(parameters.keywords, "period_iterator")};'...;
                                     'mouse = parameters.values{strcmp(parameters.keywords, "mouse")};' ...
                                     'if  (period_iterator == 1) && strcmp(mouse, "1100");'...
                                     'data_evaluated = [];'...
                                     'else;'...
                                     'data_evaluated = parameters.data;'...
                                     'end'}};
% Inputs
% each mouse
parameters.loop_list.things_to_load.data.dir = {[parameters.dir_exper 'behavior\gait analysis\normalized with std only\stride segmentations\from depressions\concatenated periods\'],'paw', '\', 'mouse', '\'};
parameters.loop_list.things_to_load.data.filename = {'stride_segmentations_average', '.mat'};
parameters.loop_list.things_to_load.data.variable = {'stride_segmentations_average{', 'period_iterator', ', 1}'}; 
parameters.loop_list.things_to_load.data.level = 'mouse';

% Outputs
% concatenated data
parameters.loop_list.things_to_save.concatenated_data.dir = {[parameters.dir_exper 'behavior\gait analysis\normalized with std only\stride segmentations\from depressions\all mice\data\'],'paw', '\'};
parameters.loop_list.things_to_save.concatenated_data.filename = {'data_all_','period_iterator','_', 'period', '.mat'};
parameters.loop_list.things_to_save.concatenated_data.variable = {'data_all'}; 
parameters.loop_list.things_to_save.concatenated_data.level = 'period';
% average
parameters.loop_list.things_to_save.average.dir = {[parameters.dir_exper 'behavior\gait analysis\normalized with std only\stride segmentations\from depressions\all mice\data\'],'paw', '\'};
parameters.loop_list.things_to_save.average.filename = {'average_','period_iterator','_', 'period', '.mat'};
parameters.loop_list.things_to_save.average.variable = {'average'}; 
parameters.loop_list.things_to_save.average.level = 'period';
% std_dev
parameters.loop_list.things_to_save.std_dev.dir = {[parameters.dir_exper 'behavior\gait analysis\normalized with std only\stride segmentations\from depressions\all mice\data\'],'paw', '\'};
parameters.loop_list.things_to_save.std_dev.filename = {'std_dev_','period_iterator','_', 'period', '.mat'};
parameters.loop_list.things_to_save.std_dev.variable = {'std_dev'}; 
parameters.loop_list.things_to_save.std_dev.level = 'period';

parameters.loop_list.things_to_rename = {   {'data_evaluated', 'data'}
                                            {'concatenated_data', 'data'}};

RunAnalysis({@EvaluateOnData, @ConcatenateData, @AverageData}, parameters);

%% Plot average of each period across mice
% use standard error of the mean (SEM) as the errors
% use consitent axes limits
if isfield(parameters, 'loop_list')
parameters = rmfield(parameters,'loop_list');
end

% Is so you can use a single loop for calculations. 
parameters.loop_list.iterators = {
               'paw', {'loop_variables.paws'}, 'paw_iterator';
               'period', {'loop_variables.periods_longsOnly'}, 'period_iterator';
               }; 
parameters.errorType = 'SEM';
parameters.instancesDim = 1; % for calculating SEM
parameters.ylimits = [-2 2];
parameters.colors = [ 1, 0, 0;  % spontaneous
                     0, 0, 0;  % 1600
                     0, 0, 0.586;  % 2000
                     0, 0.234, 1;  % 2400
                     0, 0.625, 1  % 2800
                     ]; 
parameters.lineTypes = {'-';  % spontaneous
                       '-';  % 1600
                    '-'; %'-.'; % 2000
                  '-';  %':'; % 2400
                   '-'; %'--'; % 2800
                     }; 

parameters.do_allPeriods = true; 


% Inputs
% average
parameters.loop_list.things_to_load.average.dir = {[parameters.dir_exper 'behavior\gait analysis\normalized with std only\stride segmentations\from depressions\all mice\data\'],'paw', '\'};
parameters.loop_list.things_to_load.average.filename = {'average_','period_iterator','_', 'period', '.mat'};
parameters.loop_list.things_to_load.average.variable = {'average'}; 
parameters.loop_list.things_to_load.average.level = 'period';
% std dev
parameters.loop_list.things_to_load.std_dev.dir = {[parameters.dir_exper 'behavior\gait analysis\normalized with std only\stride segmentations\from depressions\all mice\data\'],'paw', '\'};
parameters.loop_list.things_to_load.std_dev.filename = {'std_dev_','period_iterator','_', 'period', '.mat'};
parameters.loop_list.things_to_load.std_dev.variable = {'std_dev'}; 
parameters.loop_list.things_to_load.std_dev.level = 'period';
% data_all (to get the number of mice used for SEM)
parameters.loop_list.things_to_load.concatenated_data.dir = {[parameters.dir_exper 'behavior\gait analysis\normalized with std only\stride segmentations\from depressions\all mice\data\'],'paw', '\'};
parameters.loop_list.things_to_load.concatenated_data.filename = {'data_all_','period_iterator','_', 'period', '.mat'};
parameters.loop_list.things_to_load.concatenated_data.variable = {'data_all'}; 
parameters.loop_list.things_to_load.concatenated_data.level = 'period';

% Outputs
% figure
parameters.loop_list.things_to_save.fig.dir = {[parameters.dir_exper 'behavior\gait analysis\normalized with std only\stride segmentations\from depressions\all mice\average figures\'],'paw', '\'};
parameters.loop_list.things_to_save.fig.filename = {'overlay_','period_iterator','_', 'period', ['_' parameters.errorType]};
parameters.loop_list.things_to_save.fig.variable = {'fig'}; 
parameters.loop_list.things_to_save.fig.level = 'period';
parameters.loop_list.things_to_save.fig.saveas_type = 'svg';
% SEM
parameters.loop_list.things_to_save.SEM.dir = {[parameters.dir_exper 'behavior\gait analysis\normalized with std only\stride segmentations\from depressions\all mice\data\'],'paw', '\'};
parameters.loop_list.things_to_save.SEM.filename = {'SEM_','period_iterator','_', 'period', '.mat'};
parameters.loop_list.things_to_save.SEM.variable = {'SEM'}; 
parameters.loop_list.things_to_save.SEM.level = 'period';
% all together figure;
parameters.loop_list.things_to_save.fig_allPeriods.dir = {[parameters.dir_exper 'behavior\gait analysis\normalized with std only\stride segmentations\from depressions\all mice\average figures\'],'paw', '\'};
parameters.loop_list.things_to_save.fig_allPeriods.filename = {['overlay_allPeriods_' parameters.errorType '_'], 'paw'};
parameters.loop_list.things_to_save.fig_allPeriods.variable = {'fig'}; 
parameters.loop_list.things_to_save.fig_allPeriods.level = 'paw';
parameters.loop_list.things_to_save.fig_allPeriods.saveas_type = 'svg';

RunAnalysis({@PlotMiceStrideAverages}, parameters);
close all;


%% average motorized within animals 
if isfield(parameters, 'loop_list')
parameters = rmfield(parameters,'loop_list');
end

% Is so you can use a single loop for calculations. 
parameters.loop_list.iterators = {
               'paw', {'loop_variables.paws'}, 'paw_iterator';
                'mouse', {'loop_variables.mice_all(:).name'}, 'mouse_iterator';
               'period', {'loop_variables.periods_longsOnly'}, 'period_iterator';      
               }; 


% Don't include spontaneous in this average
parameters.evaluation_instructions = {{'period_iterator = parameters.values{strcmp(parameters.keywords, "period_iterator")};'...
                                     'if  period_iterator == 1;' ...
                                     'data_evaluated = [];'...
                                     'else;'...
                                     'data_evaluated = parameters.data;'...
                                     'end'}};
parameters.concatDim = 1;
parameters.concatenation_level = 'period';
parameters.averageDim = 1;

% Inputs
% each mouse
parameters.loop_list.things_to_load.data.dir = {[parameters.dir_exper 'behavior\gait analysis\normalized with std only\stride segmentations\from depressions\concatenated periods\'],'paw', '\', 'mouse', '\'};
parameters.loop_list.things_to_load.data.filename = {'stride_segmentations_average', '.mat'};
parameters.loop_list.things_to_load.data.variable = {'stride_segmentations_average{', 'period_iterator', ', 1}'}; 
parameters.loop_list.things_to_load.data.level = 'period';

% Outputs
% concatenated data
parameters.loop_list.things_to_save.concatenated_data.dir = {[parameters.dir_exper 'behavior\gait analysis\normalized with std only\stride segmentations\from depressions\average of all motor speeds\'],'mouse', '\'};
parameters.loop_list.things_to_save.concatenated_data.filename = {'data_all_', 'paw', '.mat'};
parameters.loop_list.things_to_save.concatenated_data.variable = {'data_all'}; 
parameters.loop_list.things_to_save.concatenated_data.level = 'mouse';
% average
parameters.loop_list.things_to_save.average.dir = {[parameters.dir_exper 'behavior\gait analysis\normalized with std only\stride segmentations\from depressions\average of all motor speeds\'],'mouse', '\'};
parameters.loop_list.things_to_save.average.filename = {'average_','paw', '.mat'};
parameters.loop_list.things_to_save.average.variable = {'average'}; 
parameters.loop_list.things_to_save.average.level = 'mouse';


parameters.loop_list.things_to_rename = {{'data_evaluated', 'data'}
                                            {'concatenated_data', 'data'}};

RunAnalysis({@EvaluateOnData, @ConcatenateData, @AverageData}, parameters);

%% Average single motor period across mice
if isfield(parameters, 'loop_list')
parameters = rmfield(parameters,'loop_list');
end

% Is so you can use a single loop for calculations. 
parameters.loop_list.iterators = {
               'paw', {'loop_variables.paws'}, 'paw_iterator';
               'mouse', {'loop_variables.mice_all(:).name'}, 'mouse_iterator';
               }; 

parameters.concatDim = 1;
parameters.concatenation_level = 'mouse';
parameters.averageDim = 1;

% Inputs
% each mouse
parameters.loop_list.things_to_load.data.dir = {[parameters.dir_exper 'behavior\gait analysis\normalized with std only\stride segmentations\from depressions\average of all motor speeds\'],'mouse', '\'};
parameters.loop_list.things_to_load.data.filename = {'average_','paw', '.mat'};
parameters.loop_list.things_to_load.data.variable = {'average'}; 
parameters.loop_list.things_to_load.data.level = 'mouse';

% Outputs
% concatenated data
parameters.loop_list.things_to_save.concatenated_data.dir = {[parameters.dir_exper 'behavior\gait analysis\normalized with std only\stride segmentations\from depressions\average of all motor speeds\all mice\']};
parameters.loop_list.things_to_save.concatenated_data.filename = {'data_all_', 'paw', '.mat'};
parameters.loop_list.things_to_save.concatenated_data.variable = {'data_all'}; 
parameters.loop_list.things_to_save.concatenated_data.level = 'paw';
% average
parameters.loop_list.things_to_save.average.dir = {[parameters.dir_exper 'behavior\gait analysis\normalized with std only\stride segmentations\from depressions\average of all motor speeds\all mice\']};
parameters.loop_list.things_to_save.average.filename = {'average_','paw', '.mat'};
parameters.loop_list.things_to_save.average.variable = {'average'}; 
parameters.loop_list.things_to_save.average.level = 'paw';
% std_dev
parameters.loop_list.things_to_save.std_dev.dir = {[parameters.dir_exper 'behavior\gait analysis\normalized with std only\stride segmentations\from depressions\average of all motor speeds\all mice\']};
parameters.loop_list.things_to_save.std_dev.filename = {'std_dev_', 'paw', '.mat'};
parameters.loop_list.things_to_save.std_dev.variable = {'std_dev'}; 
parameters.loop_list.things_to_save.std_dev.level = 'paw';

parameters.loop_list.things_to_rename = {{'concatenated_data', 'data'}};

RunAnalysis({ @ConcatenateData, @AverageData}, parameters);

%% Plot single motorized average across mice
if isfield(parameters, 'loop_list')
parameters = rmfield(parameters,'loop_list');
end

% Is so you can use a single loop for calculations. 
parameters.loop_list.iterators = {
               'paw', {'loop_variables.paws'}, 'paw_iterator';
               }; 
parameters.errorType = 'SEM';
parameters.instancesDim = 1; % for calculating SEM
parameters.ylimits = [-5 5];
parameters.colors = [ 1, 0, 0;  % spontaneous
                     0, 0, 0;  % 1600
                     0, 0, 0.586;  % 2000
                     0, 0.234, 1;  % 2400
                     0, 0.625, 1  % 2800
                     ]; 
parameters.lineTypes = {'-';  % spontaneous
                       '-';  % 1600
                    '-.'; % 2000
                    ':'; % 2400
                    '--'; % 2800
                     }; 

parameters.do_allPeriods = false; 

% Inputs
% concatenated data
parameters.loop_list.things_to_load.concatenated_data.dir = {[parameters.dir_exper 'behavior\gait analysis\normalized with std only\stride segmentations\from depressions\average of all motor speeds\all mice\']};
parameters.loop_list.things_to_load.concatenated_data.filename = {'data_all_', 'paw', '.mat'};
parameters.loop_list.things_to_load.concatenated_data.variable = {'data_all'}; 
parameters.loop_list.things_to_load.concatenated_data.level = 'paw';
% average
parameters.loop_list.things_to_load.average.dir = {[parameters.dir_exper 'behavior\gait analysis\normalized with std only\stride segmentations\from depressions\average of all motor speeds\all mice\']};
parameters.loop_list.things_to_load.average.filename = {'average_','paw', '.mat'};
parameters.loop_list.things_to_load.average.variable = {'average'}; 
parameters.loop_list.things_to_load.average.level = 'paw';
% std_dev
parameters.loop_list.things_to_load.std_dev.dir = {[parameters.dir_exper 'behavior\gait analysis\normalized with std only\stride segmentations\from depressions\average of all motor speeds\all mice\']};
parameters.loop_list.things_to_load.std_dev.filename = {'std_dev_', 'paw', '.mat'};
parameters.loop_list.things_to_load.std_dev.variable = {'std_dev'}; 
parameters.loop_list.things_to_load.std_dev.level = 'paw';

% Outputs
% figure
parameters.loop_list.things_to_save.fig.dir = {[parameters.dir_exper 'behavior\gait analysis\normalized with std only\stride segmentations\from depressions\average of all motor speeds\all mice\']};
parameters.loop_list.things_to_save.fig.filename = {'overlay_','paw', ['_' parameters.errorType]};
parameters.loop_list.things_to_save.fig.variable = {'fig'}; 
parameters.loop_list.things_to_save.fig.level = 'paw';
parameters.loop_list.things_to_save.fig.saveas_type = 'svg';
% SEM
parameters.loop_list.things_to_save.SEM.dir = {[parameters.dir_exper 'behavior\gait analysis\normalized with std only\stride segmentations\from depressions\average of all motor speeds\all mice\']};
parameters.loop_list.things_to_save.SEM.filename = {'SEM_','paw', '.mat'};
parameters.loop_list.things_to_save.SEM.variable = {'SEM'}; 
parameters.loop_list.things_to_save.SEM.level = 'paw';

RunAnalysis({@PlotMiceStrideAverages}, parameters);
close all;

%% Plot single motorized and spontaneous together
if isfield(parameters, 'loop_list')
parameters = rmfield(parameters,'loop_list');
end

% Is so you can use a single loop for calculations. 
parameters.loop_list.iterators = {
               'paw', {'loop_variables.paws'}, 'paw_iterator';
               }; 

parameters.errorType = 'SEM';
parameters.instancesDim = 1; % for calculating SEM
parameters.ylimits = [-5 5];
parameters.colors = [1 0 0; % spontaneous
                    0 0 1]; % motorized
% Inputs
% Motorized
% concatenated data
parameters.loop_list.things_to_load.concatenated_data_motor.dir = {[parameters.dir_exper 'behavior\gait analysis\normalized with std only\stride segmentations\from depressions\average of all motor speeds\all mice\']};
parameters.loop_list.things_to_load.concatenated_data_motor.filename = {'data_all_', 'paw', '.mat'};
parameters.loop_list.things_to_load.concatenated_data_motor.variable = {'data_all'}; 
parameters.loop_list.things_to_load.concatenated_data_motor.level = 'paw';
% average
parameters.loop_list.things_to_load.average_motor.dir = {[parameters.dir_exper 'behavior\gait analysis\normalized with std only\stride segmentations\from depressions\average of all motor speeds\all mice\']};
parameters.loop_list.things_to_load.average_motor.filename = {'average_','paw', '.mat'};
parameters.loop_list.things_to_load.average_motor.variable = {'average'}; 
parameters.loop_list.things_to_load.average_motor.level = 'paw';
% std_dev
parameters.loop_list.things_to_load.std_dev_motor.dir = {[parameters.dir_exper 'behavior\gait analysis\normalized with std only\stride segmentations\from depressions\average of all motor speeds\all mice\']};
parameters.loop_list.things_to_load.std_dev_motor.filename = {'std_dev_', 'paw', '.mat'};
parameters.loop_list.things_to_load.std_dev_motor.variable = {'std_dev'}; 
parameters.loop_list.things_to_load.std_dev_motor.level = 'paw';

% Spontaneous
% average
parameters.loop_list.things_to_load.average_spon.dir = {[parameters.dir_exper 'behavior\gait analysis\normalized with std only\stride segmentations\from depressions\all mice\data\'],'paw', '\'};
parameters.loop_list.things_to_load.average_spon.filename = {'average_1_walkLong_spon.mat'};
parameters.loop_list.things_to_load.average_spon.variable = {'average'}; 
parameters.loop_list.things_to_load.average_spon.level = 'paw';
% std dev
parameters.loop_list.things_to_load.std_dev_spon.dir = {[parameters.dir_exper 'behavior\gait analysis\normalized with std only\stride segmentations\from depressions\all mice\data\'],'paw', '\'};
parameters.loop_list.things_to_load.std_dev_spon.filename = {'std_dev_1_walkLong_spon.mat'};
parameters.loop_list.things_to_load.std_dev_spon.variable = {'std_dev'}; 
parameters.loop_list.things_to_load.std_dev_spon.level = 'paw';
% data_all (to get the number of mice used for SEM)
parameters.loop_list.things_to_load.concatenated_data_spon.dir = {[parameters.dir_exper 'behavior\gait analysis\normalized with std only\stride segmentations\from depressions\all mice\data\'],'paw', '\'};
parameters.loop_list.things_to_load.concatenated_data_spon.filename = {'data_all_1_walkLong_spon.mat'};
parameters.loop_list.things_to_load.concatenated_data_spon.variable = {'data_all'}; 
parameters.loop_list.things_to_load.concatenated_data_spon.level = 'paw';

% Outputs
parameters.loop_list.things_to_save.fig.dir = {[parameters.dir_exper 'behavior\gait analysis\normalized with std only\stride segmentations\from depressions\average of all motor speeds\all mice\']};
parameters.loop_list.things_to_save.fig.filename = {'oneMotor_oneSpon_','paw'};
parameters.loop_list.things_to_save.fig.variable = {'fig'}; 
parameters.loop_list.things_to_save.fig.level = 'paw';
parameters.loop_list.things_to_save.fig.saveas_type = 'svg';

RunAnalysis({@PlotMiceStrideOneMotorAndSpon}, parameters)

close all;

%% Segment FL, HL, and tail with FL x ranges -- Spontaneous long walk
% To find phase differences 
if isfield(parameters, 'loop_list')
parameters = rmfield(parameters,'loop_list');
end

% Is so you can use a single loop for calculations. 
parameters.loop_list.iterators = {
               'body_part', {'loop_variables.body_parts_all(:).name'}, 'body_part_iterator';
               'velocity_direction', {'loop_variables.body_parts_all(', 'body_part_iterator', ').directions(:).name'}, 'direction_iterator';
               'mouse', {'loop_variables.mice_all(:).name'}, 'mouse_iterator';
                };

parameters.segmentDim = 1;
parameters.concatDim = 2;
parameters.instancesAsCells = true; 
parameters.uniformSegments = false;

% Inputs
% timeseries
parameters.loop_list.things_to_load.timeseries.dir = {[parameters.dir_exper 'behavior\body\normalized with std only\concatenated velocity\'], 'body_part', '\', 'velocity_direction', '\spontaneous\', 'mouse', '\'};
parameters.loop_list.things_to_load.timeseries.filename = {'concatenated_velocity_longPeriods_walk.mat'};
parameters.loop_list.things_to_load.timeseries.variable = {'velocity_all'}; 
parameters.loop_list.things_to_load.timeseries.level = 'mouse';
% time ranges
parameters.loop_list.things_to_load.time_ranges.dir = {[parameters.dir_exper 'behavior\gait analysis\normalized with std only\x peaks\all periods\FL\'], 'mouse', '\'};
parameters.loop_list.things_to_load.time_ranges.filename = {'x_peaks_longWalk_spontaneous.mat'};
parameters.loop_list.things_to_load.time_ranges.variable = {'x_peaks.depression_ranges{1}'}; 
parameters.loop_list.things_to_load.time_ranges.level = 'mouse';

% Outputs
% segmented timseries
parameters.loop_list.things_to_save.segmented_timeseries.dir = {[parameters.dir_exper 'behavior\gait analysis\normalized with std only\stride segmentations from FL x depressions\all periods\'],'body_part', '\', 'velocity_direction', '\' 'mouse', '\'};
parameters.loop_list.things_to_save.segmented_timeseries.filename = {'stride_segmentations_longWalk_spontaneous.mat'};
parameters.loop_list.things_to_save.segmented_timeseries.variable = {'stride_segmentations_depression'}; 
parameters.loop_list.things_to_save.segmented_timeseries.level = 'mouse';

RunAnalysis({@StrideSegmentationLooper}, parameters);

%% Segment FL, HL, and tail with FL x ranges -- Motorized long walk
% To find phase differences 
if isfield(parameters, 'loop_list')
parameters = rmfield(parameters,'loop_list');
end

% Is so you can use a single loop for calculations. 
parameters.loop_list.iterators = {
               'body_part', {'loop_variables.body_parts_all(:).name'}, 'body_part_iterator';
               'velocity_direction', {'loop_variables.body_parts_all(', 'body_part_iterator', ').directions(:).name'}, 'direction_iterator';
               'mouse', {'loop_variables.mice_all(:).name'}, 'mouse_iterator';
               'motorSpeed', {'loop_variables.motorSpeeds'}, 'motorSpeed_iterator'
                };

parameters.segmentDim = 1;
parameters.concatDim = 2;
parameters.instancesAsCells = true;

% Inputs
% timeseries
parameters.loop_list.things_to_load.timeseries.dir = {[parameters.dir_exper 'behavior\body\normalized with std only\concatenated velocity\'], 'body_part', '\', 'velocity_direction', '\motorized\', 'mouse', '\'};
parameters.loop_list.things_to_load.timeseries.filename = {'concatenated_velocity_longPeriods_walk_', 'motorSpeed', '.mat'};
parameters.loop_list.things_to_load.timeseries.variable = {'velocity_all'}; 
parameters.loop_list.things_to_load.timeseries.level = 'motorSpeed';
% time ranges
parameters.loop_list.things_to_load.time_ranges.dir = {[parameters.dir_exper 'behavior\gait analysis\normalized with std only\x peaks\all periods\FL\'], 'mouse', '\'};
parameters.loop_list.things_to_load.time_ranges.filename = {'x_peaks_longWalk_motorized', 'motorSpeed', '.mat'};
parameters.loop_list.things_to_load.time_ranges.variable = {'x_peaks.depression_ranges{1}'}; 
parameters.loop_list.things_to_load.time_ranges.level = 'motorSpeed';

% Outputs
% segmented timseries
parameters.loop_list.things_to_save.segmented_timeseries.dir = {[parameters.dir_exper 'behavior\gait analysis\normalized with std only\stride segmentations from FL x depressions\all periods\'],'body_part', '\', 'velocity_direction', '\' 'mouse', '\'};
parameters.loop_list.things_to_save.segmented_timeseries.filename = {'stride_segmentations_longWalk_motorized', 'motorSpeed', '.mat'};
parameters.loop_list.things_to_save.segmented_timeseries.variable = {'stride_segmentations_depression'}; 
parameters.loop_list.things_to_save.segmented_timeseries.level = 'motorSpeed';

RunAnalysis({@StrideSegmentationLooper}, parameters);

%% Body parts segmented with FL x: Concatenate 

if isfield(parameters, 'loop_list')
parameters = rmfield(parameters,'loop_list');
end

% Is so you can use a single loop for calculations. 
parameters.loop_list.iterators = {
               'body_part', {'loop_variables.body_parts_all(:).name'}, 'body_part_iterator';
               'velocity_direction', {'loop_variables.body_parts_all(', 'body_part_iterator', ').directions(:).name'}, 'direction_iterator';
               'mouse', {'loop_variables.mice_all(:).name'}, 'mouse_iterator';
               'type_tag', {'loop_variables.type_tags'}, 'type_tag_iterator'
                };

parameters.concatDim = 1;
parameters.concatenation_level = 'type_tag';
parameters.concatenate_across_cells = true;

% Inputs
% data
parameters.loop_list.things_to_load.data.dir = {[parameters.dir_exper 'behavior\gait analysis\normalized with std only\stride segmentations from FL x depressions\all periods\'],'body_part', '\', 'velocity_direction', '\' 'mouse', '\'};
parameters.loop_list.things_to_load.data.filename = {'stride_segmentations_', 'type_tag', '.mat'};
parameters.loop_list.things_to_load.data.variable = {'stride_segmentations_depression'}; 
parameters.loop_list.things_to_load.data.level = 'type_tag';

% Outputs
parameters.loop_list.things_to_save.concatenated_data.dir = {[parameters.dir_exper 'behavior\gait analysis\normalized with std only\stride segmentations from FL x depressions\concatenated periods\'],'body_part', '\', 'velocity_direction', '\' 'mouse', '\'};
parameters.loop_list.things_to_save.concatenated_data.filename = {'stride_segmentations.mat'};
parameters.loop_list.things_to_save.concatenated_data.variable = {'stride_segmentations'}; 
parameters.loop_list.things_to_save.concatenated_data.level = 'mouse';

RunAnalysis({@ConcatenateData}, parameters);

%% Other paws with FL x depressions: resample and plot

if isfield(parameters, 'loop_list')
parameters = rmfield(parameters,'loop_list');
end

% Is so you can use a single loop for calculations. 
parameters.loop_list.iterators = {
               'body_part', {'loop_variables.body_parts_all(:).name'}, 'body_part_iterator';
               'velocity_direction', {'loop_variables.body_parts_all(', 'body_part_iterator', ').directions(:).name'}, 'direction_iterator';
               'mouse', {'loop_variables.mice_all(:).name'}, 'mouse_iterator';
               'period', {'loop_variables.periods_longsOnly'}, 'period_iterator';
               }; 

% close each figure after saving it
parameters.closeFigures = true;
% the number of timepoints to resample each stride velocity segment to.
parameters.resampleLength = 10; % to 0.5 s = 10 time points

% Inputs
parameters.loop_list.things_to_load.data.dir = {[parameters.dir_exper 'behavior\gait analysis\normalized with std only\stride segmentations from FL x depressions\concatenated periods\'],'body_part', '\', 'velocity_direction', '\' 'mouse', '\'};
parameters.loop_list.things_to_load.data.filename = {'stride_segmentations.mat'};
parameters.loop_list.things_to_load.data.variable = {'stride_segmentations{', 'period_iterator', '}'}; 
parameters.loop_list.things_to_load.data.level = 'mouse';

% Outputs
% stride segmentations reformatted so strides are all on same cell level
parameters.loop_list.things_to_save.segmentations_together.dir =  {[parameters.dir_exper 'behavior\gait analysis\normalized with std only\stride segmentations from FL x depressions\concatenated periods\'], 'body_part', '\', 'velocity_direction', '\', 'mouse', '\'};
parameters.loop_list.things_to_save.segmentations_together.filename = {'stride_segmentations_together', '.mat'};
parameters.loop_list.things_to_save.segmentations_together.variable = {'stride_segmentations_together{', 'period_iterator', ', 1}'}; 
parameters.loop_list.things_to_save.segmentations_together.level = 'mouse';
% resampled segmentations
parameters.loop_list.things_to_save.resampled.dir =  {[parameters.dir_exper 'behavior\gait analysis\normalized with std only\stride segmentations from FL x depressions\concatenated periods\'], 'body_part', '\', 'velocity_direction', '\', 'mouse', '\'};
parameters.loop_list.things_to_save.resampled.filename = {'stride_segmentations_resampled', '.mat'};
parameters.loop_list.things_to_save.resampled.variable = {'stride_segmentations_resampled{', 'period_iterator', ', 1}'}; 
parameters.loop_list.things_to_save.resampled.level = 'mouse';
% mean 
parameters.loop_list.things_to_save.average.dir =  {[parameters.dir_exper 'behavior\gait analysis\normalized with std only\stride segmentations from FL x depressions\concatenated periods\'], 'body_part', '\', 'velocity_direction', '\', 'mouse', '\'};
parameters.loop_list.things_to_save.average.filename = {'stride_segmentations_average', '.mat'};
parameters.loop_list.things_to_save.average.variable = {'stride_segmentations_average{', 'period_iterator', ', 1}'}; 
parameters.loop_list.things_to_save.average.level = 'mouse';
% std
parameters.loop_list.things_to_save.std_dev.dir =  {[parameters.dir_exper 'behavior\gait analysis\normalized with std only\stride segmentations from FL x depressions\concatenated periods\'], 'body_part', '\', 'velocity_direction', '\', 'mouse', '\'};
parameters.loop_list.things_to_save.std_dev.filename = {'stride_segmentations_std_dev', '.mat'};
parameters.loop_list.things_to_save.std_dev.variable = {'stride_segmentations_std_dev{', 'period_iterator', ', 1}'}; 
parameters.loop_list.things_to_save.std_dev.level = 'mouse';
% figure: not-resampled segmentations 
parameters.loop_list.things_to_save.fig_segmentations_together.dir =  {[parameters.dir_exper 'behavior\gait analysis\normalized with std only\stride segmentations from FL x depressions\concatenated periods\'], 'body_part', '\', 'velocity_direction', '\', 'mouse', '\'};
parameters.loop_list.things_to_save.fig_segmentations_together.filename = {'stride_segmentations_together_','period', '_', 'period_iterator', '.fig'};
parameters.loop_list.things_to_save.fig_segmentations_together.variable = {'fig'}; 
parameters.loop_list.things_to_save.fig_segmentations_together.level = 'period';
% figure: resampled segmentations
parameters.loop_list.things_to_save.fig_resampled.dir =  {[parameters.dir_exper 'behavior\gait analysis\normalized with std only\stride segmentations from FL x depressions\concatenated periods\'], 'body_part', '\', 'velocity_direction', '\', 'mouse', '\'};
parameters.loop_list.things_to_save.fig_resampled.filename = {'stride_segmentations_resampled_','period', '_', 'period_iterator', '.fig'};
parameters.loop_list.things_to_save.fig_resampled.variable = {'fig'}; 
parameters.loop_list.things_to_save.fig_resampled.level = 'period';
% figure: mean and std
parameters.loop_list.things_to_save.fig_average.dir = {[parameters.dir_exper 'behavior\gait analysis\normalized with std only\stride segmentations from FL x depressions\concatenated periods\'], 'body_part', '\', 'velocity_direction', '\', 'mouse', '\'};
parameters.loop_list.things_to_save.fig_average.filename = {'stride_segmentations_average_', 'period', '_', 'period_iterator', '.fig'};
parameters.loop_list.things_to_save.fig_average.variable = {'fig'}; 
parameters.loop_list.things_to_save.fig_average.level = 'period';

RunAnalysis({@GaitResampling}, parameters);

parameters.closeFigures = false;

%% segmented with FL x: Plot all mices' averages of each period together

if isfield(parameters, 'loop_list')
parameters = rmfield(parameters,'loop_list');
end

% Is so you can use a single loop for calculations. 
parameters.loop_list.iterators = {
              'body_part', {'loop_variables.body_parts_all(:).name'}, 'body_part_iterator';
               'velocity_direction', {'loop_variables.body_parts_all(', 'body_part_iterator', ').directions(:).name'}, 'direction_iterator';
               'period', {'loop_variables.periods_longsOnly'}, 'period_iterator';
               'mouse', {'loop_variables.mice_all(:).name'}, 'mouse_iterator';
               }; 


parameters.instancesDim = 1;
parameters.ylimits = [-6 6];
parameters.mymap = flipud(hsv(7));

% Inputs
% resampled segmentations
parameters.loop_list.things_to_load.resampled.dir =  {[parameters.dir_exper 'behavior\gait analysis\normalized with std only\stride segmentations from FL x depressions\concatenated periods\'], 'body_part', '\', 'velocity_direction', '\', 'mouse', '\'};
parameters.loop_list.things_to_load.resampled.filename = {'stride_segmentations_resampled', '.mat'};
parameters.loop_list.things_to_load.resampled.variable = {'stride_segmentations_resampled{', 'period_iterator', ', 1}'}; 
parameters.loop_list.things_to_load.resampled.level = 'mouse';
% mean 
parameters.loop_list.things_to_load.average.dir =  {[parameters.dir_exper 'behavior\gait analysis\normalized with std only\stride segmentations from FL x depressions\concatenated periods\'], 'body_part', '\', 'velocity_direction', '\', 'mouse', '\'};
parameters.loop_list.things_to_load.average.filename = {'stride_segmentations_average', '.mat'};
parameters.loop_list.things_to_load.average.variable = {'stride_segmentations_average{', 'period_iterator', ', 1}'}; 
parameters.loop_list.things_to_load.average.level = 'mouse';
% std
parameters.loop_list.things_to_load.std_dev.dir =  {[parameters.dir_exper 'behavior\gait analysis\normalized with std only\stride segmentations from FL x depressions\concatenated periods\'], 'body_part', '\', 'velocity_direction', '\', 'mouse', '\'};
parameters.loop_list.things_to_load.std_dev.filename = {'stride_segmentations_std_dev', '.mat'};
parameters.loop_list.things_to_load.std_dev.variable = {'stride_segmentations_std_dev{', 'period_iterator', ', 1}'}; 
parameters.loop_list.things_to_load.std_dev.level = 'mouse';

% Outputs
% figure
parameters.loop_list.things_to_save.fig.dir = {[parameters.dir_exper 'behavior\gait analysis\normalized with std only\stride segmentations from FL x depressions\all mice\overlays\'], 'body_part', '\' 'velocity_direction', '\'};
parameters.loop_list.things_to_save.fig.filename = {'overlay_','period', '_', 'period_iterator', '.fig'};
parameters.loop_list.things_to_save.fig.variable = {'fig'}; 
parameters.loop_list.things_to_save.fig.level = 'period';

RunAnalysis({@PlotMiceStrideOverlays}, parameters);

close all;

%% segmented with FL x: Average each period across mice
% include m1107, aren't comparing across spon and motorized yet
% concatenate & average

if isfield(parameters, 'loop_list')
parameters = rmfield(parameters,'loop_list');
end

% Is so you can use a single loop for calculations. 
parameters.loop_list.iterators = {
               'body_part', {'loop_variables.body_parts_all(:).name'}, 'body_part_iterator';
               'velocity_direction', {'loop_variables.body_parts_all(', 'body_part_iterator', ').directions(:).name'}, 'direction_iterator';
               'period', {'loop_variables.periods_longsOnly'}, 'period_iterator';
               'mouse', {'loop_variables.mice_all(:).name'}, 'mouse_iterator';
               }; 

parameters.concatDim = 1;
parameters.concatenation_level = 'mouse';
parameters.averageDim = 1;

% Don't include mouse 1100 in spontaneous averages
parameters.evaluation_instructions = {{'period_iterator = parameters.values{strcmp(parameters.keywords, "period_iterator")};'...;
                                     'mouse = parameters.values{strcmp(parameters.keywords, "mouse")};' ...
                                     'if  period_iterator == 1 && strcmp(mouse, "1100");'...
                                     'data_evaluated = [];'...
                                     'else;'...
                                     'data_evaluated = parameters.data;'...
                                     'end'}};
% Inputs
% each mouse
parameters.loop_list.things_to_load.data.dir =  {[parameters.dir_exper 'behavior\gait analysis\normalized with std only\stride segmentations from FL x depressions\concatenated periods\'], 'body_part', '\', 'velocity_direction', '\', 'mouse', '\'};
parameters.loop_list.things_to_load.data.filename = {'stride_segmentations_average', '.mat'};
parameters.loop_list.things_to_load.data.variable = {'stride_segmentations_average{', 'period_iterator', ', 1}'}; 
parameters.loop_list.things_to_load.data.level = 'mouse';

% Outputs
% concatenated data
parameters.loop_list.things_to_save.concatenated_data.dir = {[parameters.dir_exper 'behavior\gait analysis\normalized with std only\stride segmentations from FL x depressions\all mice\data\'],'body_part', '\', 'velocity_direction', '\'};
parameters.loop_list.things_to_save.concatenated_data.filename = {'data_all_','period_iterator','_', 'period', '.mat'};
parameters.loop_list.things_to_save.concatenated_data.variable = {'data_all'}; 
parameters.loop_list.things_to_save.concatenated_data.level = 'period';
% average
parameters.loop_list.things_to_save.average.dir = {[parameters.dir_exper 'behavior\gait analysis\normalized with std only\stride segmentations from FL x depressions\all mice\data\'],'body_part', '\', 'velocity_direction', '\'};
parameters.loop_list.things_to_save.average.filename = {'average_','period_iterator','_', 'period', '.mat'};
parameters.loop_list.things_to_save.average.variable = {'average'}; 
parameters.loop_list.things_to_save.average.level = 'period';
% std_dev
parameters.loop_list.things_to_save.std_dev.dir = {[parameters.dir_exper 'behavior\gait analysis\normalized with std only\stride segmentations from FL x depressions\all mice\data\'],'body_part', '\', 'velocity_direction', '\'};
parameters.loop_list.things_to_save.std_dev.filename = {'std_dev_','period_iterator','_', 'period', '.mat'};
parameters.loop_list.things_to_save.std_dev.variable = {'std_dev'}; 
parameters.loop_list.things_to_save.std_dev.level = 'period';

parameters.loop_list.things_to_rename = {   {'data_evaluated', 'data'}
                                            {'concatenated_data', 'data'}};

RunAnalysis({@EvaluateOnData, @ConcatenateData, @AverageData}, parameters);

%% segmented with FL x: Plot average of each period across mice
% use standard error of the mean (SEM) as the errors
% use consitent axes limits
if isfield(parameters, 'loop_list')
parameters = rmfield(parameters,'loop_list');
end

% Is so you can use a single loop for calculations. 
parameters.loop_list.iterators = {
               'body_part', {'loop_variables.body_parts_all(:).name'}, 'body_part_iterator';
               'velocity_direction', {'loop_variables.body_parts_all(', 'body_part_iterator', ').directions(:).name'}, 'direction_iterator';
               'period', {'loop_variables.periods_longsOnly'}, 'period_iterator';
               }; 

parameters.instancesDim = 1; % for calculating SEM
parameters.ylimits = [-3 3];
parameters.errorType = 'SEM';

% Inputs
% average
parameters.loop_list.things_to_load.average.dir = {[parameters.dir_exper 'behavior\gait analysis\normalized with std only\stride segmentations from FL x depressions\all mice\data\'],'body_part', '\', 'velocity_direction', '\'};
parameters.loop_list.things_to_load.average.filename = {'average_','period_iterator','_', 'period', '.mat'};
parameters.loop_list.things_to_load.average.variable = {'average'}; 
parameters.loop_list.things_to_load.average.level = 'period';
% std dev
parameters.loop_list.things_to_load.std_dev.dir = {[parameters.dir_exper 'behavior\gait analysis\normalized with std only\stride segmentations from FL x depressions\all mice\data\'],'body_part', '\', 'velocity_direction', '\'};
parameters.loop_list.things_to_load.std_dev.filename = {'std_dev_','period_iterator','_', 'period', '.mat'};
parameters.loop_list.things_to_load.std_dev.variable = {'std_dev'}; 
parameters.loop_list.things_to_load.std_dev.level = 'period';
% data_all (to get the number of mice used for SEM)
parameters.loop_list.things_to_load.concatenated_data.dir = {[parameters.dir_exper 'behavior\gait analysis\normalized with std only\stride segmentations from FL x depressions\all mice\data\'],'body_part', '\', 'velocity_direction', '\'};
parameters.loop_list.things_to_load.concatenated_data.filename = {'data_all_','period_iterator','_', 'period', '.mat'};
parameters.loop_list.things_to_load.concatenated_data.variable = {'data_all'}; 
parameters.loop_list.things_to_load.concatenated_data.level = 'period';

% Outputs
% figure
parameters.loop_list.things_to_save.fig.dir = {[parameters.dir_exper 'behavior\gait analysis\normalized with std only\stride segmentations from FL x depressions\all mice\average figures\'],'body_part', '\', 'velocity_direction', '\'};
parameters.loop_list.things_to_save.fig.filename = {'overlay_','period_iterator','_', 'period'};
parameters.loop_list.things_to_save.fig.variable = {'fig'}; 
parameters.loop_list.things_to_save.fig.level = 'period';
parameters.loop_list.things_to_save.fig.saveas_type = 'svg';
% SEM
parameters.loop_list.things_to_save.SEM.dir = {[parameters.dir_exper 'behavior\gait analysis\normalized with std only\stride segmentations from FL x depressions\all mice\data\'],'body_part', '\', 'velocity_direction', '\'};
parameters.loop_list.things_to_save.SEM.filename = {'SEM_','period_iterator','_', 'period', '.mat'};
parameters.loop_list.things_to_save.SEM.variable = {'SEM'}; 
parameters.loop_list.things_to_save.SEM.level = 'period';

RunAnalysis({@PlotMiceStrideAverages}, parameters);
close all;


%% segmented with FLx: average motorized within animals 
if isfield(parameters, 'loop_list')
parameters = rmfield(parameters,'loop_list');
end

% Is so you can use a single loop for calculations. 
parameters.loop_list.iterators = {
               'paw', {'loop_variables.body_parts_all(:).name'}, 'paw_iterator';
               'velocity_direction', {'loop_variables.body_parts_all(', 'paw_iterator', ').directions(:).name'}, 'velocity_direction_iterator';
                'mouse', {'loop_variables.mice_all(:).name'}, 'mouse_iterator';
               'period', {'loop_variables.periods_longsOnly'}, 'period_iterator';      
               }; 


% Don't include spontaneous in this average
parameters.evaluation_instructions = {{'period_iterator = parameters.values{strcmp(parameters.keywords, "period_iterator")};'...
                                     'if  period_iterator == 1;' ...
                                     'data_evaluated = [];'...
                                     'else;'...
                                     'data_evaluated = parameters.data;'...
                                     'end'}};
parameters.concatDim = 1;
parameters.concatenation_level = 'period';
parameters.averageDim = 1;

% Inputs
% each mouse
parameters.loop_list.things_to_load.data.dir = {[parameters.dir_exper 'behavior\gait analysis\normalized with std only\stride segmentations from FL x depressions\concatenated periods\'],'paw', '\' 'velocity_direction', '\', 'mouse', '\'};
parameters.loop_list.things_to_load.data.filename = {'stride_segmentations_average', '.mat'};
parameters.loop_list.things_to_load.data.variable = {'stride_segmentations_average{', 'period_iterator', ', 1}'}; 
parameters.loop_list.things_to_load.data.level = 'period';

% Outputs
% concatenated data
parameters.loop_list.things_to_save.concatenated_data.dir = {[parameters.dir_exper 'behavior\gait analysis\normalized with std only\stride segmentations from FL x depressions\average of all motor speeds\'], 'velocity_direction', '\','mouse', '\'};
parameters.loop_list.things_to_save.concatenated_data.filename = {'data_all_', 'paw', '.mat'};
parameters.loop_list.things_to_save.concatenated_data.variable = {'data_all'}; 
parameters.loop_list.things_to_save.concatenated_data.level = 'mouse';
% average
parameters.loop_list.things_to_save.average.dir = {[parameters.dir_exper 'behavior\gait analysis\normalized with std only\stride segmentations from FL x depressions\average of all motor speeds\'], 'velocity_direction', '\', 'mouse', '\'};
parameters.loop_list.things_to_save.average.filename = {'average_','paw', '.mat'};
parameters.loop_list.things_to_save.average.variable = {'average'}; 
parameters.loop_list.things_to_save.average.level = 'mouse';


parameters.loop_list.things_to_rename = {{'data_evaluated', 'data'}
                                            {'concatenated_data', 'data'}};

RunAnalysis({@EvaluateOnData, @ConcatenateData, @AverageData}, parameters);

%% segmented with FLx: Average single motor period across mice
if isfield(parameters, 'loop_list')
parameters = rmfield(parameters,'loop_list');
end

% Is so you can use a single loop for calculations. 
parameters.loop_list.iterators = {
               'paw', {'loop_variables.body_parts_all(:).name'}, 'paw_iterator';
               'velocity_direction', {'loop_variables.body_parts_all(', 'paw_iterator', ').directions(:).name'}, 'velocity_direction_iterator';
               'mouse', {'loop_variables.mice_all(:).name'}, 'mouse_iterator';
               }; 

parameters.concatDim = 1;
parameters.concatenation_level = 'mouse';
parameters.averageDim = 1;

% Inputs
% each mouse
parameters.loop_list.things_to_load.data.dir = {[parameters.dir_exper 'behavior\gait analysis\normalized with std only\stride segmentations from FL x depressions\average of all motor speeds\'], 'velocity_direction', '\','mouse', '\'};
parameters.loop_list.things_to_load.data.filename = {'average_','paw', '.mat'};
parameters.loop_list.things_to_load.data.variable = {'average'}; 
parameters.loop_list.things_to_load.data.level = 'mouse';

% Outputs
% concatenated data
parameters.loop_list.things_to_save.concatenated_data.dir = {[parameters.dir_exper 'behavior\gait analysis\normalized with std only\stride segmentations from FL x depressions\average of all motor speeds\'], 'velocity_direction', '\all mice\'};
parameters.loop_list.things_to_save.concatenated_data.filename = {'data_all_', 'paw', '.mat'};
parameters.loop_list.things_to_save.concatenated_data.variable = {'data_all'}; 
parameters.loop_list.things_to_save.concatenated_data.level = 'velocity_direction';
% average
parameters.loop_list.things_to_save.average.dir = {[parameters.dir_exper 'behavior\gait analysis\normalized with std only\stride segmentations from FL x depressions\average of all motor speeds\'], 'velocity_direction', '\all mice\'};
parameters.loop_list.things_to_save.average.filename = {'average_','paw', '.mat'};
parameters.loop_list.things_to_save.average.variable = {'average'}; 
parameters.loop_list.things_to_save.average.level = 'velocity_direction';
% std_dev
parameters.loop_list.things_to_save.std_dev.dir = {[parameters.dir_exper 'behavior\gait analysis\normalized with std only\stride segmentations from FL x depressions\average of all motor speeds\'], 'velocity_direction', '\all mice\'};
parameters.loop_list.things_to_save.std_dev.filename = {'std_dev_', 'paw', '.mat'};
parameters.loop_list.things_to_save.std_dev.variable = {'std_dev'}; 
parameters.loop_list.things_to_save.std_dev.level = 'velocity_direction';

parameters.loop_list.things_to_rename = {{'concatenated_data', 'data'}};

RunAnalysis({ @ConcatenateData, @AverageData}, parameters);

%% segmented with FLx: Plot single motorized average across mice
if isfield(parameters, 'loop_list')
parameters = rmfield(parameters,'loop_list');
end

% Is so you can use a single loop for calculations. 
parameters.loop_list.iterators = {
               'paw', {'loop_variables.body_parts_all(:).name'}, 'paw_iterator';
               'velocity_direction', {'loop_variables.body_parts_all(', 'paw_iterator', ').directions(:).name'}, 'velocity_direction_iterator';
               }; 

parameters.errorType = 'SEM';
parameters.instancesDim = 1; % for calculating SEM
parameters.ylimits = [-3 3];
parameters.colors = [ 1, 0, 0;  % spontaneous
                     0, 0, 0;  % 1600
                     0, 0, 0.586;  % 2000
                     0, 0.234, 1;  % 2400
                     0, 0.625, 1  % 2800
                     ]; 
parameters.lineTypes = {'-';  % spontaneous
                       '-';  % 1600
                    '-.'; % 2000
                    ':'; % 2400
                    '--'; % 2800
                     }; 

parameters.do_allPeriods = false; 

% Inputs
% concatenated data
parameters.loop_list.things_to_load.concatenated_data.dir = {[parameters.dir_exper 'behavior\gait analysis\normalized with std only\stride segmentations from FL x depressions\average of all motor speeds\'],'velocity_direction', '\all mice\'};
parameters.loop_list.things_to_load.concatenated_data.filename = {'data_all_', 'paw', '.mat'};
parameters.loop_list.things_to_load.concatenated_data.variable = {'data_all'}; 
parameters.loop_list.things_to_load.concatenated_data.level = 'paw';
% average
parameters.loop_list.things_to_load.average.dir = {[parameters.dir_exper 'behavior\gait analysis\normalized with std only\stride segmentations from FL x depressions\average of all motor speeds\'],'velocity_direction', '\all mice\'};
parameters.loop_list.things_to_load.average.filename = {'average_','paw', '.mat'};
parameters.loop_list.things_to_load.average.variable = {'average'}; 
parameters.loop_list.things_to_load.average.level = 'paw';
% std_dev
parameters.loop_list.things_to_load.std_dev.dir = {[parameters.dir_exper 'behavior\gait analysis\normalized with std only\stride segmentations from FL x depressions\average of all motor speeds\'],'velocity_direction', '\all mice\'};
parameters.loop_list.things_to_load.std_dev.filename = {'std_dev_', 'paw', '.mat'};
parameters.loop_list.things_to_load.std_dev.variable = {'std_dev'}; 
parameters.loop_list.things_to_load.std_dev.level = 'paw';

% Outputs
% figure
parameters.loop_list.things_to_save.fig.dir = {[parameters.dir_exper 'behavior\gait analysis\normalized with std only\stride segmentations from FL x depressions\average of all motor speeds\'],'velocity_direction', '\all mice\'};
parameters.loop_list.things_to_save.fig.filename = {'overlay_','paw', ['_' parameters.errorType]};
parameters.loop_list.things_to_save.fig.variable = {'fig'}; 
parameters.loop_list.things_to_save.fig.level = 'paw';
parameters.loop_list.things_to_save.fig.saveas_type = 'svg';
% SEM
parameters.loop_list.things_to_save.SEM.dir = {[parameters.dir_exper 'behavior\gait analysis\normalized with std only\stride segmentations from FL x depressions\average of all motor speeds\'],'velocity_direction', '\all mice\'};
parameters.loop_list.things_to_save.SEM.filename = {'SEM_','paw', '.mat'};
parameters.loop_list.things_to_save.SEM.variable = {'SEM'}; 
parameters.loop_list.things_to_save.SEM.level = 'paw';

RunAnalysis({@PlotMiceStrideAverages}, parameters);
close all;

%% segmented from FL x:  Plot single motorized and spontaneous together
if isfield(parameters, 'loop_list')
parameters = rmfield(parameters,'loop_list');
end

% Is so you can use a single loop for calculations. 
parameters.loop_list.iterators = {
               'paw', {'loop_variables.body_parts_all(:).name'}, 'paw_iterator';
               'velocity_direction', {'loop_variables.body_parts_all(', 'paw_iterator', ').directions(:).name'}, 'velocity_direction_iterator';
               }; 

parameters.errorType = 'SEM';
parameters.instancesDim = 1; % for calculating SEM
parameters.ylimits = [-5 5];
parameters.colors = [1 0 0; % spontaneous
                    0 0 1]; % motorized
% Inputs
% Motorized
% concatenated data
parameters.loop_list.things_to_load.concatenated_data_motor.dir = {[parameters.dir_exper 'behavior\gait analysis\normalized with std only\stride segmentations from FL x depressions\average of all motor speeds\'], 'velocity_direction', '\all mice\'};
parameters.loop_list.things_to_load.concatenated_data_motor.filename = {'data_all_', 'paw', '.mat'};
parameters.loop_list.things_to_load.concatenated_data_motor.variable = {'data_all'}; 
parameters.loop_list.things_to_load.concatenated_data_motor.level = 'velocity_direction';
% average
parameters.loop_list.things_to_load.average_motor.dir = {[parameters.dir_exper 'behavior\gait analysis\normalized with std only\stride segmentations from FL x depressions\average of all motor speeds\'], 'velocity_direction', '\all mice\'};
parameters.loop_list.things_to_load.average_motor.filename = {'average_','paw', '.mat'};
parameters.loop_list.things_to_load.average_motor.variable = {'average'}; 
parameters.loop_list.things_to_load.average_motor.level = 'velocity_direction';
% std_dev
parameters.loop_list.things_to_load.std_dev_motor.dir = {[parameters.dir_exper 'behavior\gait analysis\normalized with std only\stride segmentations from FL x depressions\average of all motor speeds\'], 'velocity_direction', '\all mice\'};
parameters.loop_list.things_to_load.std_dev_motor.filename = {'std_dev_', 'paw', '.mat'};
parameters.loop_list.things_to_load.std_dev_motor.variable = {'std_dev'}; 
parameters.loop_list.things_to_load.std_dev_motor.level = 'velocity_direction';

% Spontaneous
% average
parameters.loop_list.things_to_load.average_spon.dir = {[parameters.dir_exper 'behavior\gait analysis\normalized with std only\stride segmentations from FL x depressions\all mice\data\'],'paw', '\', 'velocity_direction', '\'};
parameters.loop_list.things_to_load.average_spon.filename = {'average_1_walkLong_spon.mat'};
parameters.loop_list.things_to_load.average_spon.variable = {'average'}; 
parameters.loop_list.things_to_load.average_spon.level = 'velocity_direction';
% std dev
parameters.loop_list.things_to_load.std_dev_spon.dir = {[parameters.dir_exper 'behavior\gait analysis\normalized with std only\stride segmentations from FL x depressions\all mice\data\'],'paw', '\', 'velocity_direction', '\'};
parameters.loop_list.things_to_load.std_dev_spon.filename = {'std_dev_1_walkLong_spon.mat'};
parameters.loop_list.things_to_load.std_dev_spon.variable = {'std_dev'}; 
parameters.loop_list.things_to_load.std_dev_spon.level = 'velocity_direction';
% data_all (to get the number of mice used for SEM)
parameters.loop_list.things_to_load.concatenated_data_spon.dir = {[parameters.dir_exper 'behavior\gait analysis\normalized with std only\stride segmentations from FL x depressions\all mice\data\'],'paw','\', 'velocity_direction', '\'};
parameters.loop_list.things_to_load.concatenated_data_spon.filename = {'data_all_1_walkLong_spon.mat'};
parameters.loop_list.things_to_load.concatenated_data_spon.variable = {'data_all'}; 
parameters.loop_list.things_to_load.concatenated_data_spon.level = 'velocity_direction';

% Outputs
parameters.loop_list.things_to_save.fig.dir = {[parameters.dir_exper 'behavior\gait analysis\normalized with std only\stride segmentations from FL x depressions\average of all motor speeds\'], 'velocity_direction', '\all mice\'};
parameters.loop_list.things_to_save.fig.filename = {'oneMotor_oneSpon_','paw'};
parameters.loop_list.things_to_save.fig.variable = {'fig'}; 
parameters.loop_list.things_to_save.fig.level = 'velocity_direction';
parameters.loop_list.things_to_save.fig.saveas_type = 'svg';

RunAnalysis({@PlotMiceStrideOneMotorAndSpon}, parameters)

close all;



%% segment other velocity directions based on that body part's x-direction: spontaneous

% To find phase differences 
if isfield(parameters, 'loop_list')
parameters = rmfield(parameters,'loop_list');
end

parameters.loop_variables.body_parts_all(1).directions = parameters.body_parts_all(1).directions(1);
parameters.loop_variables.body_parts_all(2).directions = parameters.body_parts_all(2).directions(1);
parameters.loop_variables.body_parts_all(3).directions = parameters.body_parts_all(3).directions(1);

% Is so you can use a single loop for calculations. 
parameters.loop_list.iterators = {
               'body_part', {'loop_variables.body_parts_all(:).name'}, 'body_part_iterator';
               'velocity_direction', {'loop_variables.body_parts_all(', 'body_part_iterator', ').directions(:).name'}, 'direction_iterator';
               'mouse', {'loop_variables.mice_all(:).name'}, 'mouse_iterator';
                };

parameters.segmentDim = 1;
parameters.concatDim = 2;
parameters.instancesAsCells = true; 
parameters.uniformSegments = false;

% Inputs
% timeseries
parameters.loop_list.things_to_load.timeseries.dir = {[parameters.dir_exper 'behavior\body\normalized with std only\concatenated velocity\'], 'body_part', '\', 'velocity_direction', '\spontaneous\', 'mouse', '\'};
parameters.loop_list.things_to_load.timeseries.filename = {'concatenated_velocity_longPeriods_walk.mat'};
parameters.loop_list.things_to_load.timeseries.variable = {'velocity_all'}; 
parameters.loop_list.things_to_load.timeseries.level = 'mouse';
% time ranges
parameters.loop_list.things_to_load.time_ranges.dir = {[parameters.dir_exper 'behavior\gait analysis\normalized with std only\x peaks\all periods\'], 'body_part', '\', 'mouse', '\'};
parameters.loop_list.things_to_load.time_ranges.filename = {'x_peaks_longWalk_spontaneous.mat'};
parameters.loop_list.things_to_load.time_ranges.variable = {'x_peaks.depression_ranges{1}'}; 
parameters.loop_list.things_to_load.time_ranges.level = 'mouse';

% Outputs
% segmented timseries
parameters.loop_list.things_to_save.segmented_timeseries.dir = {[parameters.dir_exper 'behavior\gait analysis\normalized with std only\stride segmentations from own x depressions\all periods\'],'body_part', '\', 'velocity_direction', '\' 'mouse', '\'};
parameters.loop_list.things_to_save.segmented_timeseries.filename = {'stride_segmentations_longWalk_spontaneous.mat'};
parameters.loop_list.things_to_save.segmented_timeseries.variable = {'stride_segmentations_depression'}; 
parameters.loop_list.things_to_save.segmented_timeseries.level = 'mouse';

RunAnalysis({@StrideSegmentationLooper}, parameters);

%% segment other velocity directions based on that body part's x-direction: Motorized long walk
% To find phase differences 
if isfield(parameters, 'loop_list')
parameters = rmfield(parameters,'loop_list');
end

% Is so you can use a single loop for calculations. 
parameters.loop_list.iterators = {
               'body_part', {'loop_variables.body_parts_all(:).name'}, 'body_part_iterator';
               'velocity_direction', {'loop_variables.body_parts_all(', 'body_part_iterator', ').directions(:).name'}, 'direction_iterator';
               'mouse', {'loop_variables.mice_all(:).name'}, 'mouse_iterator';
               'motorSpeed', {'loop_variables.motorSpeeds'}, 'motorSpeed_iterator'
                };

parameters.segmentDim = 1;
parameters.concatDim = 2;
parameters.instancesAsCells = true;

% Inputs
% timeseries
parameters.loop_list.things_to_load.timeseries.dir = {[parameters.dir_exper 'behavior\body\normalized with std only\concatenated velocity\'], 'body_part', '\', 'velocity_direction', '\motorized\', 'mouse', '\'};
parameters.loop_list.things_to_load.timeseries.filename = {'concatenated_velocity_longPeriods_walk_', 'motorSpeed', '.mat'};
parameters.loop_list.things_to_load.timeseries.variable = {'velocity_all'}; 
parameters.loop_list.things_to_load.timeseries.level = 'motorSpeed';
% time ranges
parameters.loop_list.things_to_load.time_ranges.dir = {[parameters.dir_exper 'behavior\gait analysis\normalized with std only\x peaks\all periods\'], 'body_part', '\', 'mouse', '\'};
parameters.loop_list.things_to_load.time_ranges.filename = {'x_peaks_longWalk_motorized', 'motorSpeed', '.mat'};
parameters.loop_list.things_to_load.time_ranges.variable = {'x_peaks.depression_ranges{1}'}; 
parameters.loop_list.things_to_load.time_ranges.level = 'motorSpeed';

% Outputs
% segmented timseries
parameters.loop_list.things_to_save.segmented_timeseries.dir = {[parameters.dir_exper 'behavior\gait analysis\normalized with std only\stride segmentations from own x depressions\all periods\'],'body_part', '\', 'velocity_direction', '\' 'mouse', '\'};
parameters.loop_list.things_to_save.segmented_timeseries.filename = {'stride_segmentations_longWalk_motorized', 'motorSpeed', '.mat'};
parameters.loop_list.things_to_save.segmented_timeseries.variable = {'stride_segmentations_depression'}; 
parameters.loop_list.things_to_save.segmented_timeseries.level = 'motorSpeed';

RunAnalysis({@StrideSegmentationLooper}, parameters);

%% other velocity directions based on that body part's x-direction: Concatenate
if isfield(parameters, 'loop_list')
parameters = rmfield(parameters,'loop_list');
end

% Is so you can use a single loop for calculations. 
parameters.loop_list.iterators = {
              'body_part', {'loop_variables.body_parts_all(:).name'}, 'body_part_iterator';
               'velocity_direction', {'loop_variables.body_parts_all(', 'body_part_iterator', ').directions(:).name'}, 'direction_iterator';
               'mouse', {'loop_variables.mice_all(:).name'}, 'mouse_iterator';
               'type_tag', {'loop_variables.type_tags'}, 'type_tag_iterator'
                };

parameters.concatDim = 1;
parameters.concatenation_level = 'type_tag';
parameters.concatenate_across_cells = true;

% Inputs
% data
parameters.loop_list.things_to_load.data.dir = {[parameters.dir_exper 'behavior\gait analysis\normalized with std only\stride segmentations from own x depressions\all periods\'],'body_part', '\', 'velocity_direction', '\' 'mouse', '\'};
parameters.loop_list.things_to_load.data.filename = {'stride_segmentations_', 'type_tag', '.mat'};
parameters.loop_list.things_to_load.data.variable = {'stride_segmentations_depression'}; 
parameters.loop_list.things_to_load.data.level = 'type_tag';

% Outputs
parameters.loop_list.things_to_save.concatenated_data.dir = {[parameters.dir_exper 'behavior\gait analysis\normalized with std only\stride segmentations from own x depressions\concatenated periods\'],'body_part', '\', 'velocity_direction', '\' 'mouse', '\'};
parameters.loop_list.things_to_save.concatenated_data.filename = {'stride_segmentations.mat'};
parameters.loop_list.things_to_save.concatenated_data.variable = {'stride_segmentations'}; 
parameters.loop_list.things_to_save.concatenated_data.level = 'mouse';

RunAnalysis({@ConcatenateData}, parameters);

%% other velocity directions based on that body part's x-direction: Gait Resampling
if isfield(parameters, 'loop_list')
parameters = rmfield(parameters,'loop_list');
end

% Is so you can use a single loop for calculations. 
parameters.loop_list.iterators = {
              'body_part', {'loop_variables.body_parts_all(:).name'}, 'body_part_iterator';
               'velocity_direction', {'loop_variables.body_parts_all(', 'body_part_iterator', ').directions(:).name'}, 'direction_iterator';
               'mouse', {'loop_variables.mice_all(:).name'}, 'mouse_iterator';
               'period', {'loop_variables.periods_longsOnly'}, 'period_iterator';
               }; 

% close each figure after saving it
parameters.closeFigures = true;
% the number of timepoints to resample each stride velocity segment to.
parameters.resampleLength = 10; % to 0.5 s = 10 time points

% Inputs
parameters.loop_list.things_to_load.data.dir = {[parameters.dir_exper 'behavior\gait analysis\normalized with std only\stride segmentations from own x depressions\concatenated periods\'],'body_part', '\', 'velocity_direction', '\' 'mouse', '\'};
parameters.loop_list.things_to_load.data.filename = {'stride_segmentations.mat'};
parameters.loop_list.things_to_load.data.variable = {'stride_segmentations{', 'period_iterator', '}'}; 
parameters.loop_list.things_to_load.data.level = 'mouse';

% Outputs
% stride segmentations reformatted so strides are all on same cell level
parameters.loop_list.things_to_save.segmentations_together.dir =  {[parameters.dir_exper 'behavior\gait analysis\normalized with std only\stride segmentations from own x depressions\concatenated periods\'], 'body_part', '\', 'velocity_direction', '\', 'mouse', '\'};
parameters.loop_list.things_to_save.segmentations_together.filename = {'stride_segmentations_together', '.mat'};
parameters.loop_list.things_to_save.segmentations_together.variable = {'stride_segmentations_together{', 'period_iterator', ', 1}'}; 
parameters.loop_list.things_to_save.segmentations_together.level = 'mouse';
% resampled segmentations
parameters.loop_list.things_to_save.resampled.dir =  {[parameters.dir_exper 'behavior\gait analysis\normalized with std only\stride segmentations from own x depressions\concatenated periods\'], 'body_part', '\', 'velocity_direction', '\', 'mouse', '\'};
parameters.loop_list.things_to_save.resampled.filename = {'stride_segmentations_resampled', '.mat'};
parameters.loop_list.things_to_save.resampled.variable = {'stride_segmentations_resampled{', 'period_iterator', ', 1}'}; 
parameters.loop_list.things_to_save.resampled.level = 'mouse';
% mean 
parameters.loop_list.things_to_save.average.dir =  {[parameters.dir_exper 'behavior\gait analysis\normalized with std only\stride segmentations from own x depressions\concatenated periods\'], 'body_part', '\', 'velocity_direction', '\', 'mouse', '\'};
parameters.loop_list.things_to_save.average.filename = {'stride_segmentations_average', '.mat'};
parameters.loop_list.things_to_save.average.variable = {'stride_segmentations_average{', 'period_iterator', ', 1}'}; 
parameters.loop_list.things_to_save.average.level = 'mouse';
% std
parameters.loop_list.things_to_save.std_dev.dir =  {[parameters.dir_exper 'behavior\gait analysis\normalized with std only\stride segmentations from own x depressions\concatenated periods\'], 'body_part', '\', 'velocity_direction', '\', 'mouse', '\'};
parameters.loop_list.things_to_save.std_dev.filename = {'stride_segmentations_std_dev', '.mat'};
parameters.loop_list.things_to_save.std_dev.variable = {'stride_segmentations_std_dev{', 'period_iterator', ', 1}'}; 
parameters.loop_list.things_to_save.std_dev.level = 'mouse';
% figure: not-resampled segmentations 
parameters.loop_list.things_to_save.fig_segmentations_together.dir =  {[parameters.dir_exper 'behavior\gait analysis\normalized with std only\stride segmentations from own x depressions\concatenated periods\'], 'body_part', '\', 'velocity_direction', '\', 'mouse', '\'};
parameters.loop_list.things_to_save.fig_segmentations_together.filename = {'stride_segmentations_together_','period', '_', 'period_iterator', '.fig'};
parameters.loop_list.things_to_save.fig_segmentations_together.variable = {'fig'}; 
parameters.loop_list.things_to_save.fig_segmentations_together.level = 'period';
% figure: resampled segmentations
parameters.loop_list.things_to_save.fig_resampled.dir =  {[parameters.dir_exper 'behavior\gait analysis\normalized with std only\stride segmentations from own x depressions\concatenated periods\'], 'body_part', '\', 'velocity_direction', '\', 'mouse', '\'};
parameters.loop_list.things_to_save.fig_resampled.filename = {'stride_segmentations_resampled_','period', '_', 'period_iterator', '.fig'};
parameters.loop_list.things_to_save.fig_resampled.variable = {'fig'}; 
parameters.loop_list.things_to_save.fig_resampled.level = 'period';
% figure: mean and std
parameters.loop_list.things_to_save.fig_average.dir = {[parameters.dir_exper 'behavior\gait analysis\normalized with std only\stride segmentations from own x depressions\concatenated periods\'], 'body_part', '\', 'velocity_direction', '\', 'mouse', '\'};
parameters.loop_list.things_to_save.fig_average.filename = {'stride_segmentations_average_', 'period', '_', 'period_iterator', '.fig'};
parameters.loop_list.things_to_save.fig_average.variable = {'fig'}; 
parameters.loop_list.things_to_save.fig_average.level = 'period';

RunAnalysis({@GaitResampling}, parameters);

parameters.closeFigures = false;

%% Own x segmentations: Plot all mices' averages of each period together

if isfield(parameters, 'loop_list')
parameters = rmfield(parameters,'loop_list');
end

% Is so you can use a single loop for calculations. 
parameters.loop_list.iterators = {
               'body_part', {'loop_variables.body_parts_all(:).name'}, 'body_part_iterator';
               'velocity_direction', {'loop_variables.body_parts_all(', 'body_part_iterator', ').directions(:).name'}, 'direction_iterator';
               'period', {'loop_variables.periods_longsOnly'}, 'period_iterator';
               'mouse', {'loop_variables.mice_all(:).name'}, 'mouse_iterator';
               }; 


parameters.instancesDim = 1;
parameters.ylimits = [-6 6];
parameters.mymap = flipud(hsv(7));

% Inputs
% resampled segmentations
parameters.loop_list.things_to_load.resampled.dir =  {[parameters.dir_exper 'behavior\gait analysis\normalized with std only\stride segmentations from own x depressions\concatenated periods\'], 'body_part', '\', 'velocity_direction', '\', 'mouse', '\'};
parameters.loop_list.things_to_load.resampled.filename = {'stride_segmentations_resampled', '.mat'};
parameters.loop_list.things_to_load.resampled.variable = {'stride_segmentations_resampled{', 'period_iterator', ', 1}'}; 
parameters.loop_list.things_to_load.resampled.level = 'mouse';
% mean 
parameters.loop_list.things_to_load.average.dir =  {[parameters.dir_exper 'behavior\gait analysis\normalized with std only\stride segmentations from own x depressions\concatenated periods\'], 'body_part', '\', 'velocity_direction', '\', 'mouse', '\'};
parameters.loop_list.things_to_load.average.filename = {'stride_segmentations_average', '.mat'};
parameters.loop_list.things_to_load.average.variable = {'stride_segmentations_average{', 'period_iterator', ', 1}'}; 
parameters.loop_list.things_to_load.average.level = 'mouse';
% std
parameters.loop_list.things_to_load.std_dev.dir =  {[parameters.dir_exper 'behavior\gait analysis\normalized with std only\stride segmentations from own x depressions\concatenated periods\'], 'body_part', '\', 'velocity_direction', '\', 'mouse', '\'};
parameters.loop_list.things_to_load.std_dev.filename = {'stride_segmentations_std_dev', '.mat'};
parameters.loop_list.things_to_load.std_dev.variable = {'stride_segmentations_std_dev{', 'period_iterator', ', 1}'}; 
parameters.loop_list.things_to_load.std_dev.level = 'mouse';

% Outputs
% figure
parameters.loop_list.things_to_save.fig.dir = {[parameters.dir_exper 'behavior\gait analysis\normalized with std only\stride segmentations from own x depressions\all mice\overlays\'], 'body_part', '\' 'velocity_direction', '\'};
parameters.loop_list.things_to_save.fig.filename = {'overlay_','period', '_', 'period_iterator', '.fig'};
parameters.loop_list.things_to_save.fig.variable = {'fig'}; 
parameters.loop_list.things_to_save.fig.level = 'period';

RunAnalysis({@PlotMiceStrideOverlays}, parameters);

close all;

%% Own x: Average each period across mice
% include m1107, aren't comparing across spon and motorized yet
% concatenate & average

if isfield(parameters, 'loop_list')
parameters = rmfield(parameters,'loop_list');
end

% Is so you can use a single loop for calculations. 
parameters.loop_list.iterators = {
              'body_part', {'loop_variables.body_parts_all(:).name'}, 'body_part_iterator';
               'velocity_direction', {'loop_variables.body_parts_all(', 'body_part_iterator', ').directions(:).name'}, 'direction_iterator';
               'period', {'loop_variables.periods_longsOnly'}, 'period_iterator';
               'mouse', {'loop_variables.mice_all(:).name'}, 'mouse_iterator';
               }; 

parameters.concatDim = 1;
parameters.concatenation_level = 'mouse';
parameters.averageDim = 1;

% Don't include mouse 1100 in spontaneous averages
parameters.evaluation_instructions = {{'period_iterator = parameters.values{strcmp(parameters.keywords, "period_iterator")};'...;
                                     'mouse = parameters.values{strcmp(parameters.keywords, "mouse")};' ...
                                     'if  period_iterator == 1 && strcmp(mouse, "1100");'...
                                     'data_evaluated = [];'...
                                     'else;'...
                                     'data_evaluated = parameters.data;'...
                                     'end'}};
% Inputs
% each mouse
parameters.loop_list.things_to_load.data.dir =  {[parameters.dir_exper 'behavior\gait analysis\normalized with std only\stride segmentations from own x depressions\concatenated periods\'], 'body_part', '\', 'velocity_direction', '\', 'mouse', '\'};
parameters.loop_list.things_to_load.data.filename = {'stride_segmentations_average', '.mat'};
parameters.loop_list.things_to_load.data.variable = {'stride_segmentations_average{', 'period_iterator', ', 1}'}; 
parameters.loop_list.things_to_load.data.level = 'mouse';

% Outputs
% concatenated data
parameters.loop_list.things_to_save.concatenated_data.dir = {[parameters.dir_exper 'behavior\gait analysis\normalized with std only\stride segmentations from own x depressions\all mice\data\'],'body_part', '\', 'velocity_direction', '\'};
parameters.loop_list.things_to_save.concatenated_data.filename = {'data_all_','period_iterator','_', 'period', '.mat'};
parameters.loop_list.things_to_save.concatenated_data.variable = {'data_all'}; 
parameters.loop_list.things_to_save.concatenated_data.level = 'period';
% average
parameters.loop_list.things_to_save.average.dir = {[parameters.dir_exper 'behavior\gait analysis\normalized with std only\stride segmentations from own x depressions\all mice\data\'],'body_part', '\', 'velocity_direction', '\'};
parameters.loop_list.things_to_save.average.filename = {'average_','period_iterator','_', 'period', '.mat'};
parameters.loop_list.things_to_save.average.variable = {'average'}; 
parameters.loop_list.things_to_save.average.level = 'period';
% std_dev
parameters.loop_list.things_to_save.std_dev.dir = {[parameters.dir_exper 'behavior\gait analysis\normalized with std only\stride segmentations from own x depressions\all mice\data\'],'body_part', '\', 'velocity_direction', '\'};
parameters.loop_list.things_to_save.std_dev.filename = {'std_dev_','period_iterator','_', 'period', '.mat'};
parameters.loop_list.things_to_save.std_dev.variable = {'std_dev'}; 
parameters.loop_list.things_to_save.std_dev.level = 'period';

parameters.loop_list.things_to_rename = {   {'data_evaluated', 'data'}
                                            {'concatenated_data', 'data'}};

RunAnalysis({@EvaluateOnData, @ConcatenateData, @AverageData}, parameters);

%% Own x: Plot average of each period across mice
% use standard error of the mean (SEM) as the errors
% use consitent axes limits
if isfield(parameters, 'loop_list')
parameters = rmfield(parameters,'loop_list');
end

% Is so you can use a single loop for calculations. 
parameters.loop_list.iterators = {
               'body_part', {'loop_variables.body_parts_all(:).name'}, 'body_part_iterator';
               'velocity_direction', {'loop_variables.body_parts_all(', 'body_part_iterator', ').directions(:).name'}, 'direction_iterator';
               'period', {'loop_variables.periods_longsOnly'}, 'period_iterator';
               }; 

parameters.instancesDim = 1; % for calculating SEM
parameters.ylimits = [0 4];

% Inputs
% average
parameters.loop_list.things_to_load.average.dir = {[parameters.dir_exper 'behavior\gait analysis\normalized with std only\stride segmentations from own x depressions\all mice\data\'],'body_part', '\', 'velocity_direction', '\'};
parameters.loop_list.things_to_load.average.filename = {'average_','period_iterator','_', 'period', '.mat'};
parameters.loop_list.things_to_load.average.variable = {'average'}; 
parameters.loop_list.things_to_load.average.level = 'period';
% std dev
parameters.loop_list.things_to_load.std_dev.dir = {[parameters.dir_exper 'behavior\gait analysis\normalized with std only\stride segmentations from own x depressions\all mice\data\'],'body_part', '\', 'velocity_direction', '\'};
parameters.loop_list.things_to_load.std_dev.filename = {'std_dev_','period_iterator','_', 'period', '.mat'};
parameters.loop_list.things_to_load.std_dev.variable = {'std_dev'}; 
parameters.loop_list.things_to_load.std_dev.level = 'period';
% data_all (to get the number of mice used for SEM)
parameters.loop_list.things_to_load.concatenated_data.dir = {[parameters.dir_exper 'behavior\gait analysis\normalized with std only\stride segmentations from own x depressions\all mice\data\'],'body_part', '\', 'velocity_direction', '\'};
parameters.loop_list.things_to_load.concatenated_data.filename = {'data_all_','period_iterator','_', 'period', '.mat'};
parameters.loop_list.things_to_load.concatenated_data.variable = {'data_all'}; 
parameters.loop_list.things_to_load.concatenated_data.level = 'period';

% Outputs
% figure
parameters.loop_list.things_to_save.fig.dir = {[parameters.dir_exper 'behavior\gait analysis\normalized with std only\stride segmentations from own x depressions\all mice\average figures\'],'body_part', '\', 'velocity_direction', '\'};
parameters.loop_list.things_to_save.fig.filename = {'overlay_','period_iterator','_', 'period'};
parameters.loop_list.things_to_save.fig.variable = {'fig'}; 
parameters.loop_list.things_to_save.fig.level = 'period';
parameters.loop_list.things_to_save.fig.saveas_type = 'svg';
% SEM
parameters.loop_list.things_to_save.SEM.dir = {[parameters.dir_exper 'behavior\gait analysis\normalized with std only\stride segmentations from own x depressions\all mice\data\'],'body_part', '\', 'velocity_direction', '\'};
parameters.loop_list.things_to_save.SEM.filename = {'SEM_','period_iterator','_', 'period', '.mat'};
parameters.loop_list.things_to_save.SEM.variable = {'SEM'}; 
parameters.loop_list.things_to_save.SEM.level = 'period';

RunAnalysis({@PlotMiceStrideAverages}, parameters);
close all;

%% Own x: average motorized within animals 
if isfield(parameters, 'loop_list')
parameters = rmfield(parameters,'loop_list');
end

% Is so you can use a single loop for calculations. 
parameters.loop_list.iterators = {
               'paw', {'loop_variables.paws'}, 'paw_iterator';
                'mouse', {'loop_variables.mice_all(:).name'}, 'mouse_iterator';
               'period', {'loop_variables.periods_longsOnly'}, 'period_iterator';      
               }; 


% Don't include spontaneous in this average
parameters.evaluation_instructions = {{'period_iterator = parameters.values{strcmp(parameters.keywords, "period_iterator")};'...
                                     'if  period_iterator == 1;' ...
                                     'data_evaluated = [];'...
                                     'else;'...
                                     'data_evaluated = parameters.data;'...
                                     'end'}};
parameters.concatDim = 1;
parameters.concatenation_level = 'period';
parameters.averageDim = 1;

% Inputs
% each mouse
parameters.loop_list.things_to_load.data.dir = {[parameters.dir_exper 'behavior\gait analysis\normalized with std only\stride segmentations from own x depressions\concatenated periods\'],'paw', '\total_magnitude\', 'mouse', '\'};
parameters.loop_list.things_to_load.data.filename = {'stride_segmentations_average', '.mat'};
parameters.loop_list.things_to_load.data.variable = {'stride_segmentations_average{', 'period_iterator', ', 1}'}; 
parameters.loop_list.things_to_load.data.level = 'period';

% Outputs
% concatenated data
parameters.loop_list.things_to_save.concatenated_data.dir = {[parameters.dir_exper 'behavior\gait analysis\normalized with std only\stride segmentations from own x depressions\average of all motor speeds\'],'mouse', '\'};
parameters.loop_list.things_to_save.concatenated_data.filename = {'data_all_', 'paw', '.mat'};
parameters.loop_list.things_to_save.concatenated_data.variable = {'data_all'}; 
parameters.loop_list.things_to_save.concatenated_data.level = 'mouse';
% average
parameters.loop_list.things_to_save.average.dir = {[parameters.dir_exper 'behavior\gait analysis\normalized with std only\stride segmentations from own x depressions\average of all motor speeds\'],'mouse', '\'};
parameters.loop_list.things_to_save.average.filename = {'average_','paw', '.mat'};
parameters.loop_list.things_to_save.average.variable = {'average'}; 
parameters.loop_list.things_to_save.average.level = 'mouse';


parameters.loop_list.things_to_rename = {{'data_evaluated', 'data'}
                                            {'concatenated_data', 'data'}};

RunAnalysis({@EvaluateOnData, @ConcatenateData, @AverageData}, parameters);

%% Own x: Average single motor period across mice
if isfield(parameters, 'loop_list')
parameters = rmfield(parameters,'loop_list');
end

% Is so you can use a single loop for calculations. 
parameters.loop_list.iterators = {
               'paw', {'loop_variables.paws'}, 'paw_iterator';
               'mouse', {'loop_variables.mice_all(:).name'}, 'mouse_iterator';
               }; 

parameters.concatDim = 1;
parameters.concatenation_level = 'mouse';
parameters.averageDim = 1;

% Inputs
% each mouse
parameters.loop_list.things_to_load.data.dir = {[parameters.dir_exper 'behavior\gait analysis\normalized with std only\stride segmentations from own x depressions\average of all motor speeds\'],'mouse', '\'};
parameters.loop_list.things_to_load.data.filename = {'average_','paw', '.mat'};
parameters.loop_list.things_to_load.data.variable = {'average'}; 
parameters.loop_list.things_to_load.data.level = 'mouse';

% Outputs
% concatenated data
parameters.loop_list.things_to_save.concatenated_data.dir = {[parameters.dir_exper 'behavior\gait analysis\normalized with std only\stride segmentations from own x depressions\average of all motor speeds\all mice\']};
parameters.loop_list.things_to_save.concatenated_data.filename = {'data_all_', 'paw', '.mat'};
parameters.loop_list.things_to_save.concatenated_data.variable = {'data_all'}; 
parameters.loop_list.things_to_save.concatenated_data.level = 'paw';
% average
parameters.loop_list.things_to_save.average.dir = {[parameters.dir_exper 'behavior\gait analysis\normalized with std only\stride segmentations from own x depressions\average of all motor speeds\all mice\']};
parameters.loop_list.things_to_save.average.filename = {'average_','paw', '.mat'};
parameters.loop_list.things_to_save.average.variable = {'average'}; 
parameters.loop_list.things_to_save.average.level = 'paw';
% std_dev
parameters.loop_list.things_to_save.std_dev.dir = {[parameters.dir_exper 'behavior\gait analysis\normalized with std only\stride segmentations from own x depressions\average of all motor speeds\all mice\']};
parameters.loop_list.things_to_save.std_dev.filename = {'std_dev_', 'paw', '.mat'};
parameters.loop_list.things_to_save.std_dev.variable = {'std_dev'}; 
parameters.loop_list.things_to_save.std_dev.level = 'paw';

parameters.loop_list.things_to_rename = {{'concatenated_data', 'data'}};

RunAnalysis({ @ConcatenateData, @AverageData}, parameters);

%% Own x: Plot single motorized average across mice
if isfield(parameters, 'loop_list')
parameters = rmfield(parameters,'loop_list');
end

% Is so you can use a single loop for calculations. 
parameters.loop_list.iterators = {
               'paw', {'loop_variables.paws'}, 'paw_iterator';
               }; 
parameters.errorType = 'SEM';
parameters.instancesDim = 1; % for calculating SEM
parameters.ylimits = [0 2];
parameters.colors = [ 1, 0, 0;  % spontaneous
                     0, 0, 0;  % 1600
                     0, 0, 0.586;  % 2000
                     0, 0.234, 1;  % 2400
                     0, 0.625, 1  % 2800
                     ]; 
parameters.lineTypes = {'-';  % spontaneous
                       '-';  % 1600
                    '-.'; % 2000
                    ':'; % 2400
                    '--'; % 2800
                     }; 

parameters.do_allPeriods = false; 

% Inputs
% concatenated data
parameters.loop_list.things_to_load.concatenated_data.dir = {[parameters.dir_exper 'behavior\gait analysis\normalized with std only\stride segmentations from own x depressions\average of all motor speeds\all mice\']};
parameters.loop_list.things_to_load.concatenated_data.filename = {'data_all_', 'paw', '.mat'};
parameters.loop_list.things_to_load.concatenated_data.variable = {'data_all'}; 
parameters.loop_list.things_to_load.concatenated_data.level = 'paw';
% average
parameters.loop_list.things_to_load.average.dir = {[parameters.dir_exper 'behavior\gait analysis\normalized with std only\stride segmentations from own x depressions\average of all motor speeds\all mice\']};
parameters.loop_list.things_to_load.average.filename = {'average_','paw', '.mat'};
parameters.loop_list.things_to_load.average.variable = {'average'}; 
parameters.loop_list.things_to_load.average.level = 'paw';
% std_dev
parameters.loop_list.things_to_load.std_dev.dir = {[parameters.dir_exper 'behavior\gait analysis\normalized with std only\stride segmentations from own x depressions\average of all motor speeds\all mice\']};
parameters.loop_list.things_to_load.std_dev.filename = {'std_dev_', 'paw', '.mat'};
parameters.loop_list.things_to_load.std_dev.variable = {'std_dev'}; 
parameters.loop_list.things_to_load.std_dev.level = 'paw';

% Outputs
% figure
parameters.loop_list.things_to_save.fig.dir = {[parameters.dir_exper 'behavior\gait analysis\normalized with std only\stride segmentations from own x depressions\average of all motor speeds\all mice\']};
parameters.loop_list.things_to_save.fig.filename = {'overlay_','paw', ['_' parameters.errorType]};
parameters.loop_list.things_to_save.fig.variable = {'fig'}; 
parameters.loop_list.things_to_save.fig.level = 'paw';
parameters.loop_list.things_to_save.fig.saveas_type = 'svg';
% SEM
parameters.loop_list.things_to_save.SEM.dir = {[parameters.dir_exper 'behavior\gait analysis\normalized with std only\stride segmentations from own x depressions\average of all motor speeds\all mice\']};
parameters.loop_list.things_to_save.SEM.filename = {'SEM_','paw', '.mat'};
parameters.loop_list.things_to_save.SEM.variable = {'SEM'}; 
parameters.loop_list.things_to_save.SEM.level = 'paw';

RunAnalysis({@PlotMiceStrideAverages}, parameters);
close all;

%% Own x: plot one motorized and one spontaneous 
if isfield(parameters, 'loop_list')
parameters = rmfield(parameters,'loop_list');
end

% Is so you can use a single loop for calculations. 
parameters.loop_list.iterators = {
               'paw', {'loop_variables.paws'}, 'paw_iterator';
               }; 

parameters.errorType = 'SEM';
parameters.instancesDim = 1; % for calculating SEM
parameters.ylimits = [0 5];
parameters.colors = [1 0 0; % spontaneous
                    0 0 1]; % motorized
% Inputs
% Motorized
% concatenated data
parameters.loop_list.things_to_load.concatenated_data_motor.dir = {[parameters.dir_exper 'behavior\gait analysis\normalized with std only\stride segmentations from own x depressions\average of all motor speeds\all mice\']};
parameters.loop_list.things_to_load.concatenated_data_motor.filename = {'data_all_', 'paw', '.mat'};
parameters.loop_list.things_to_load.concatenated_data_motor.variable = {'data_all'}; 
parameters.loop_list.things_to_load.concatenated_data_motor.level = 'paw';
% average
parameters.loop_list.things_to_load.average_motor.dir = {[parameters.dir_exper 'behavior\gait analysis\normalized with std only\stride segmentations from own x depressions\average of all motor speeds\all mice\']};
parameters.loop_list.things_to_load.average_motor.filename = {'average_','paw', '.mat'};
parameters.loop_list.things_to_load.average_motor.variable = {'average'}; 
parameters.loop_list.things_to_load.average_motor.level = 'paw';
% std_dev
parameters.loop_list.things_to_load.std_dev_motor.dir = {[parameters.dir_exper 'behavior\gait analysis\normalized with std only\stride segmentations from own x depressions\average of all motor speeds\all mice\']};
parameters.loop_list.things_to_load.std_dev_motor.filename = {'std_dev_', 'paw', '.mat'};
parameters.loop_list.things_to_load.std_dev_motor.variable = {'std_dev'}; 
parameters.loop_list.things_to_load.std_dev_motor.level = 'paw';

% Spontaneous
% average
parameters.loop_list.things_to_load.average_spon.dir = {[parameters.dir_exper 'behavior\gait analysis\normalized with std only\stride segmentations from own x depressions\all mice\data\'],'paw', '\total_magnitude\'};
parameters.loop_list.things_to_load.average_spon.filename = {'average_1_walkLong_spon.mat'};
parameters.loop_list.things_to_load.average_spon.variable = {'average'}; 
parameters.loop_list.things_to_load.average_spon.level = 'paw';
% std dev
parameters.loop_list.things_to_load.std_dev_spon.dir = {[parameters.dir_exper 'behavior\gait analysis\normalized with std only\stride segmentations from own x depressions\all mice\data\'],'paw', '\total_magnitude\'};
parameters.loop_list.things_to_load.std_dev_spon.filename = {'std_dev_1_walkLong_spon.mat'};
parameters.loop_list.things_to_load.std_dev_spon.variable = {'std_dev'}; 
parameters.loop_list.things_to_load.std_dev_spon.level = 'paw';
% data_all (to get the number of mice used for SEM)
parameters.loop_list.things_to_load.concatenated_data_spon.dir = {[parameters.dir_exper 'behavior\gait analysis\normalized with std only\stride segmentations from own x depressions\all mice\data\'],'paw', '\total_magnitude\'};
parameters.loop_list.things_to_load.concatenated_data_spon.filename = {'data_all_1_walkLong_spon.mat'};
parameters.loop_list.things_to_load.concatenated_data_spon.variable = {'data_all'}; 
parameters.loop_list.things_to_load.concatenated_data_spon.level = 'paw';

% Outputs
parameters.loop_list.things_to_save.fig.dir = {[parameters.dir_exper 'behavior\gait analysis\normalized with std only\stride segmentations from own x depressions\average of all motor speeds\all mice\']};
parameters.loop_list.things_to_save.fig.filename = {'oneMotor_oneSpon_','paw'};
parameters.loop_list.things_to_save.fig.variable = {'fig'}; 
parameters.loop_list.things_to_save.fig.level = 'paw';
parameters.loop_list.things_to_save.fig.saveas_type = 'svg';

RunAnalysis({@PlotMiceStrideOneMotorAndSpon}, parameters)

close all;

%% run phase difference detections with phdiffmeasure -- spontaneous
% % use fillmissing Matlab build-in function first ('movmean', 10)
% 
% % Compare HL and tail x to FL x. (tail x is larger on average than tail y... or at least I didn't segment based on y depressions) 
% 
% if isfield(parameters, 'loop_list')
% parameters = rmfield(parameters,'loop_list');
% end
% 
% % Is so you can use a single loop for calculations. 
% parameters.loop_list.iterators = {
%                'paw', {'loop_variables.paws_sublist'}, 'paw_iterator'; % iterate through HL and tail only
%                'mouse', {'loop_variables.mice_all(:).name'}, 'mouse_iterator';
%                }; 
% 
% parameters.isLongWalk = true;
% parameters.fillMissing_window = 10; % The width of the window to run 'movmean' over in the 'fillmissing' step
% parameters.minimumLength = 60; % The minimum length of the timeseries to try to calculate phase on (3 seconds)
% 
% % Inputs 
% % reference (FL x)
% parameters.loop_list.things_to_load.reference.dir = {[parameters.dir_exper 'behavior\body\normalized with std only\concatenated velocity\FL\x\spontaneous\'], 'mouse', '\'};
% parameters.loop_list.things_to_load.reference.filename = {'concatenated_velocity_longPeriods_walk.mat'};
% parameters.loop_list.things_to_load.reference.variable = {'velocity_all'}; 
% parameters.loop_list.things_to_load.reference.level = 'mouse';
% % compare to (other paw x)
% parameters.loop_list.things_to_load.data.dir = {[parameters.dir_exper 'behavior\body\normalized with std only\concatenated velocity\'], 'paw', '\x\spontaneous\', 'mouse', '\'};
% parameters.loop_list.things_to_load.data.filename = {'concatenated_velocity_longPeriods_walk.mat'};
% parameters.loop_list.things_to_load.data.variable = {'velocity_all'}; 
% parameters.loop_list.things_to_load.data.level = 'mouse';
% 
% % Outputs
% parameters.loop_list.things_to_save.phase_differences.dir = {[parameters.dir_exper 'behavior\gait analysis\normalized with std only\phase difference\'], 'paw', '\', 'mouse', '\'};
% parameters.loop_list.things_to_save.phase_differences.filename = {'phaseDifference_longWalk_spontaneous.mat'};
% parameters.loop_list.things_to_save.phase_differences.variable = {'phase_difference'}; 
% parameters.loop_list.things_to_save.phase_differences.level = 'mouse';
% 
% RunAnalysis({@FindPhaseDifference}, parameters);
% 
% 
% %% run phase difference detections with phdiffmeasure -- motorized
% % use fillmissing Matlab build-in function first ('movmean', 10)
% 
% % Compare HL and tail x to FL x. (tail x is larger on average than tail y... or at least I didn't segment based on y depressions) 
% 
% if isfield(parameters, 'loop_list')
% parameters = rmfield(parameters,'loop_list');
% end
% 
% % Is so you can use a single loop for calculations. 
% parameters.loop_list.iterators = {
%                'paw', {'loop_variables.paws_sublist'}, 'paw_iterator'; % iterate through HL and tail only
%                'mouse', {'loop_variables.mice_all(:).name'}, 'mouse_iterator';
%                'motorSpeed', {'loop_variables.motorSpeeds'}, 'motorSpeed_iterator'
%                }; 
% 
% parameters.isLongWalk = true;
% parameters.fillMissing_window = 10; % The width of the window to run 'movmean' over in the 'fillmissing' step
% parameters.minimumLength = 60; % The minimum length of the timeseries to try to calculate phase on (3 seconds)
% 
% % Inputs 
% % reference (FL x)
% parameters.loop_list.things_to_load.reference.dir = {[parameters.dir_exper 'behavior\body\normalized with std only\concatenated velocity\FL\x\motorized\'], 'mouse', '\'};
% parameters.loop_list.things_to_load.reference.filename = {'concatenated_velocity_longPeriods_walk_', 'motorSpeed', '.mat'};
% parameters.loop_list.things_to_load.reference.variable = {'velocity_all'}; 
% parameters.loop_list.things_to_load.reference.level = 'motorSpeed';
% % compare to (other paw x)
% parameters.loop_list.things_to_load.data.dir = {[parameters.dir_exper 'behavior\body\normalized with std only\concatenated velocity\'], 'paw', '\x\motorized\', 'mouse', '\'};
% parameters.loop_list.things_to_load.data.filename = {'concatenated_velocity_longPeriods_walk_', 'motorSpeed', '.mat'};
% parameters.loop_list.things_to_load.data.variable = {'velocity_all'}; 
% parameters.loop_list.things_to_load.data.level = 'motorSpeed';
% 
% % Outputs
% parameters.loop_list.things_to_save.phase_differences.dir = {[parameters.dir_exper 'behavior\gait analysis\normalized with std only\phase difference\'], 'paw', '\', 'mouse', '\'};
% parameters.loop_list.things_to_save.phase_differences.filename = {'phaseDifference_longWalk_motorized', 'motorSpeed', '.mat'};
% parameters.loop_list.things_to_save.phase_differences.variable = {'phase_difference'}; 
% parameters.loop_list.things_to_save.phase_differences.level = 'motorSpeed';
% 
% RunAnalysis({@FindPhaseDifference}, parameters);
% 
% 
% %% Phases: Concatenate
% if isfield(parameters, 'loop_list')
% parameters = rmfield(parameters,'loop_list');
% end
% 
% % Is so you can use a single loop for calculations. 
% parameters.loop_list.iterators = {
%                'paw', {'loop_variables.paws_sublist'}, 'paw_iterator';
%                'mouse', {'loop_variables.mice_all(:).name'}, 'mouse_iterator';
%                'type_tag', {'loop_variables.type_tags'}, 'type_tag_iterator'
%                 };
% 
% parameters.concatDim = 1;
% parameters.concatenation_level = 'type_tag';
% parameters.concatenate_across_cells = true;
% 
% % Inputs
% % data
% parameters.loop_list.things_to_load.data.dir = {[parameters.dir_exper 'behavior\gait analysis\normalized with std only\phase difference\'],'paw', '\', 'mouse', '\'};
% parameters.loop_list.things_to_load.data.filename = {'phaseDifference_', 'type_tag', '.mat'};
% parameters.loop_list.things_to_load.data.variable = {'phase_difference'}; 
% parameters.loop_list.things_to_load.data.level = 'type_tag';
% 
% % Outputs
% parameters.loop_list.things_to_save.concatenated_data.dir = {[parameters.dir_exper 'behavior\gait analysis\normalized with std only\phase difference\'],'paw', '\', 'mouse', '\'};
% parameters.loop_list.things_to_save.concatenated_data.filename = {'phaseDifference_concatenated.mat'};
% parameters.loop_list.things_to_save.concatenated_data.variable = {'phase_difference'}; 
% parameters.loop_list.things_to_save.concatenated_data.level = 'mouse';
% 
% RunAnalysis({@ConcatenateData}, parameters);
% 
%% Get lengths of the long velocity segments: Spontaneous
% so when you average the phase differences, you can weight by the length
% of the walking segment

if isfield(parameters, 'loop_list')
parameters = rmfield(parameters,'loop_list');
end

% Is so you can use a single loop for calculations. 
parameters.loop_list.iterators = {
               'paw', {'loop_variables.paws_sublist'}, 'paw_iterator'; % iterate through HL and tail only
               'mouse', {'loop_variables.mice_all(:).name'}, 'mouse_iterator';
               }; 

parameters.evaluation_instructions = {{'data_evaluated  = {cellfun("size", parameters.data, 1)};'}};
% Inputs 
parameters.loop_list.things_to_load.data.dir = {[parameters.dir_exper 'behavior\body\normalized with std only\concatenated velocity\'], 'paw', '\x\spontaneous\', 'mouse', '\'};
parameters.loop_list.things_to_load.data.filename = {'concatenated_velocity_longPeriods_walk.mat'};
parameters.loop_list.things_to_load.data.variable = {'velocity_all'}; 
parameters.loop_list.things_to_load.data.level = 'mouse';

% Outputs
parameters.loop_list.things_to_save.data_evaluated.dir = {[parameters.dir_exper 'behavior\gait analysis\normalized with std only\phase difference\'], 'paw', '\', 'mouse', '\'};
parameters.loop_list.things_to_save.data_evaluated.filename = {'istanceLengths_longWalk_spontaneous.mat'};
parameters.loop_list.things_to_save.data_evaluated.variable = {'lengths'}; 
parameters.loop_list.things_to_save.data_evaluated.level = 'mouse';

RunAnalysis({@EvaluateOnData}, parameters);

%% Get lengths of the long velocity segments: Motorized
% so when you average the phase differences, you can weight by the length
% of the walking segment

if isfield(parameters, 'loop_list')
parameters = rmfield(parameters,'loop_list');
end

% Is so you can use a single loop for calculations. 
parameters.loop_list.iterators = {
               'paw', {'loop_variables.paws_sublist'}, 'paw_iterator'; % iterate through HL and tail only
               'mouse', {'loop_variables.mice_all(:).name'}, 'mouse_iterator';
               'motorSpeed', {'loop_variables.motorSpeeds'}, 'motorSpeed_iterator'
               }; 

parameters.evaluation_instructions = {{'data_evaluated  = {cellfun("size", parameters.data, 1)};'}};
% Inputs 
parameters.loop_list.things_to_load.data.dir = {[parameters.dir_exper 'behavior\body\normalized with std only\concatenated velocity\'], 'paw', '\x\motorized\', 'mouse', '\'};
parameters.loop_list.things_to_load.data.filename = {'concatenated_velocity_longPeriods_walk_', 'motorSpeed', '.mat'};
parameters.loop_list.things_to_load.data.variable = {'velocity_all'}; 
parameters.loop_list.things_to_load.data.level = 'motorSpeed';

% Outputs
parameters.loop_list.things_to_save.data_evaluated.dir = {[parameters.dir_exper 'behavior\gait analysis\normalized with std only\phase difference\'], 'paw', '\', 'mouse', '\'};
parameters.loop_list.things_to_save.data_evaluated.filename = {'istanceLengths_longWalk_motorized', 'motorSpeed', '.mat'};
parameters.loop_list.things_to_save.data_evaluated.variable = {'lengths'}; 
parameters.loop_list.things_to_save.data_evaluated.level = 'motorSpeed';

RunAnalysis({@EvaluateOnData}, parameters);

%% Lengths: Concatenate
if isfield(parameters, 'loop_list')
parameters = rmfield(parameters,'loop_list');
end

% Is so you can use a single loop for calculations. 
parameters.loop_list.iterators = {
               'paw', {'loop_variables.paws_sublist'}, 'paw_iterator';
               'mouse', {'loop_variables.mice_all(:).name'}, 'mouse_iterator';
               'type_tag', {'loop_variables.type_tags'}, 'type_tag_iterator'
                };

parameters.concatDim = 1;
parameters.concatenation_level = 'type_tag';
parameters.concatenate_across_cells = true;

% Inputs
% data
parameters.loop_list.things_to_load.data.dir = {[parameters.dir_exper 'behavior\gait analysis\normalized with std only\phase difference\'],'paw', '\', 'mouse', '\'};
parameters.loop_list.things_to_load.data.filename = {'istanceLengths_', 'type_tag', '.mat'};
parameters.loop_list.things_to_load.data.variable = {'lengths'}; 
parameters.loop_list.things_to_load.data.level = 'type_tag';

% Outputs
parameters.loop_list.things_to_save.concatenated_data.dir = {[parameters.dir_exper 'behavior\gait analysis\normalized with std only\phase difference\'],'paw', '\', 'mouse', '\'};
parameters.loop_list.things_to_save.concatenated_data.filename = {'instanceLengths_concatenated.mat'};
parameters.loop_list.things_to_save.concatenated_data.variable = {'lengths_all'}; 
parameters.loop_list.things_to_save.concatenated_data.level = 'mouse';

RunAnalysis({@ConcatenateData}, parameters);


%% Phases: Average within mouse using circular statistics
% if isfield(parameters, 'loop_list')
% parameters = rmfield(parameters,'loop_list');
% end
% 
% % Is so you can use a single loop for calculations. 
% parameters.loop_list.iterators = {
%                'paw', {'loop_variables.paws_sublist'}, 'paw_iterator';
%                'mouse', {'loop_variables.mice_all(:).name'}, 'mouse_iterator';
%                'period', {'loop_variables.periods_longsOnly'}, 'period_iterator';
%                }; 
% 
% parameters.averageDim = 1;
% parameters.useWeights = true;
% 
% % Inputs
% % data
% parameters.loop_list.things_to_load.data.dir = {[parameters.dir_exper 'behavior\gait analysis\normalized with std only\phase difference\'],'paw', '\', 'mouse', '\'};
% parameters.loop_list.things_to_load.data.filename = {'phaseDifference_concatenated.mat'};
% parameters.loop_list.things_to_load.data.variable = {'phase_difference{', 'period_iterator', '}'}; 
% parameters.loop_list.things_to_load.data.level = 'mouse';
% % weights
% parameters.loop_list.things_to_load.weights.dir = {[parameters.dir_exper 'behavior\gait analysis\normalized with std only\phase difference\'],'paw', '\', 'mouse', '\'};
% parameters.loop_list.things_to_load.weights.filename = {'instanceLengths_concatenated.mat'};
% parameters.loop_list.things_to_load.weights.variable = {'lengths_all{', 'period_iterator', '}'}; 
% parameters.loop_list.things_to_load.weights.level = 'mouse';
% 
% % Outputs
% parameters.loop_list.things_to_save.average.dir = {[parameters.dir_exper 'behavior\gait analysis\normalized with std only\phase difference\'],'paw', '\', 'mouse', '\'};
% parameters.loop_list.things_to_save.average.filename = {'within_mouse_average.mat'};
% parameters.loop_list.things_to_save.average.variable = {'average{', 'period_iterator', ', 1}'}; 
% parameters.loop_list.things_to_save.average.level = 'mouse';
% 
% parameters.loop_list.things_to_save.std_dev.dir = {[parameters.dir_exper 'behavior\gait analysis\normalized with std only\phase difference\'],'paw', '\', 'mouse', '\'};
% parameters.loop_list.things_to_save.std_dev.filename = {'within_mouse_std_dev.mat'};
% parameters.loop_list.things_to_save.std_dev.variable = {'std_dev{', 'period_iterator', ', 1}'}; 
% parameters.loop_list.things_to_save.std_dev.level = 'mouse';
% 
% RunAnalysis({@AverageCircularData}, parameters);
% 
% parameters.useWeights = false;

%% Find stride durations (from stride velocity segmentations)
if isfield(parameters, 'loop_list')
parameters = rmfield(parameters,'loop_list');
end

% Don't need to do different segmentation methods or velocity directions (just do x of each paw)
parameters.loop_list.iterators = {
               'paw', {'loop_variables.paws'}, 'paw_iterator';
               'mouse', {'loop_variables.mice_all(:).name'}, 'mouse_iterator';
               'period', {'loop_variables.periods_longsOnly'}, 'period_iterator';
               }; 

parameters.evaluation_instructions = {{'data_evaluated = cellfun(@numel, parameters.data)./parameters.fps;'}};

% Inputs
parameters.loop_list.things_to_load.data.dir = {[parameters.dir_exper 'behavior\gait analysis\normalized with std only\stride segmentations\from depressions\concatenated periods\'],'paw', '\', 'mouse', '\'};
parameters.loop_list.things_to_load.data.filename = {'stride_segmentations_together', '.mat'};
parameters.loop_list.things_to_load.data.variable = {'stride_segmentations_together{', 'period_iterator', ', 1}'}; 
parameters.loop_list.things_to_load.data.level = 'mouse';

% Outputs
parameters.loop_list.things_to_save.data_evaluated.dir = {[parameters.dir_exper 'behavior\gait analysis\normalized with std only\stride durations\'],'paw', '\', 'mouse', '\'};
parameters.loop_list.things_to_save.data_evaluated.filename = {'stride_durations', '.mat'};
parameters.loop_list.things_to_save.data_evaluated.variable = {'stride_durations{', 'period_iterator', ', 1}'}; 
parameters.loop_list.things_to_save.data_evaluated.level = 'mouse';

RunAnalysis({@EvaluateOnData}, parameters);

%% Stride durations: plot histograms
if isfield(parameters, 'loop_list')
parameters = rmfield(parameters,'loop_list');
end

parameters.loop_list.iterators = {
               'paw', {'loop_variables.paws'}, 'paw_iterator';
               'mouse', {'loop_variables.mice_all(:).name'}, 'mouse_iterator';
               'period', {'loop_variables.periods_longsOnly'}, 'period_iterator';
               }; 

parameters.bin_edges = 0:0.25:5;
parameters.title_string = 'stride durations';

% Inputs
parameters.loop_list.things_to_load.data.dir = {[parameters.dir_exper 'behavior\gait analysis\normalized with std only\stride durations\'],'paw', '\', 'mouse', '\'};
parameters.loop_list.things_to_load.data.filename = {'stride_durations', '.mat'};
parameters.loop_list.things_to_load.data.variable = {'stride_durations{', 'period_iterator', ', 1}'}; 
parameters.loop_list.things_to_load.data.level = 'mouse';

% Outputs
parameters.loop_list.things_to_save.fig.dir = {[parameters.dir_exper 'behavior\gait analysis\normalized with std only\stride durations\'],'paw', '\', 'mouse', '\'};
parameters.loop_list.things_to_save.fig.filename = {'stride_durations_', 'period', '_', 'period_iterator', '.fig'};
parameters.loop_list.things_to_save.fig.variable = {'fig'}; 
parameters.loop_list.things_to_save.fig.level = 'period';

RunAnalysis({@StrideDurationHistograms}, parameters);

close all;

%% *** Stride durations: fancy plots might do the mean within and across mice in Prism Graphpad ***


%% Segment positiosn by behavior period (long walk periods)-- Spontaneous

% (similar to paw_velocity_pipeline_code.m)
if isfield(parameters, 'loop_list')
parameters = rmfield(parameters,'loop_list');
end

% Iterators
parameters.loop_list.iterators = {'mouse', {'loop_variables.mice_all(:).name'}, 'mouse_iterator'; 
               'day', {'loop_variables.mice_all(', 'mouse_iterator', ').days(:).name'}, 'day_iterator';
                   'stack', {'loop_variables.mice_all(',  'mouse_iterator', ').days(', 'day_iterator', ').spontaneous'}, 'stack_iterator';
                    'body_part', {'loop_variables.body_parts_all_position(:).name'}, 'body_part_iterator';
               'velocity_direction', {'loop_variables.body_parts_all_position(', 'body_part_iterator', ').directions(:).name'}, 'direction_iterator';
                };
% Skip any files that don't exist (spontaneous or problem files)
parameters.load_abort_flag = false; 

% Dimension of different time range pairs.
parameters.rangePairs = 1; 

% 
parameters.segmentDim = 1;
parameters.concatDim = 2;

% Are the segmentations the same length? (If false, will put outputs into
% cell array)
parameters.uniformSegments = false;

% Input values. 
% Extracted timeseries.
parameters.loop_list.things_to_load.timeseries.dir = {[parameters.dir_exper 'behavior\body\normalized with std only\paw position normalized\'], 'mouse', '\', 'day', '\'};
parameters.loop_list.things_to_load.timeseries.filename= {'position', 'stack', '.mat'};
parameters.loop_list.things_to_load.timeseries.variable= {'position.', 'body_part', '.', 'velocity_direction'}; 
parameters.loop_list.things_to_load.timeseries.level = 'stack';
% Time ranges
parameters.loop_list.things_to_load.time_ranges.dir = {[parameters.dir_exper 'behavior\spontaneous\segmented behavior periods\'], 'mouse', '\', 'day', '\'};
parameters.loop_list.things_to_load.time_ranges.filename= {'long_periods_', 'stack', '.mat'};
parameters.loop_list.things_to_load.time_ranges.variable= {'long_periods.walk'}; 
parameters.loop_list.things_to_load.time_ranges.level = 'stack';

% Output Values
parameters.loop_list.things_to_save.segmented_timeseries.dir = {[parameters.dir_exper 'behavior\body\normalized with std only\segmented positions\'], 'body_part', '\', 'velocity_direction', '\', 'mouse', '\', 'day', '\'};
parameters.loop_list.things_to_save.segmented_timeseries.filename= {'segmented_timeseries_longWalk_spontaneous', '_', 'stack', '.mat'};
parameters.loop_list.things_to_save.segmented_timeseries.variable= {'segmented_timeseries'}; 
parameters.loop_list.things_to_save.segmented_timeseries.level = 'velocity_direction';

RunAnalysis({@SegmentTimeseriesData}, parameters);

%% Segment positiosn by behavior period (long walk periods)-- Motorized

% (similar to paw_velocity_pipeline_code.m)
if isfield(parameters, 'loop_list')
parameters = rmfield(parameters,'loop_list');
end

% Iterators
parameters.loop_list.iterators = {'mouse', {'loop_variables.mice_all(:).name'}, 'mouse_iterator'; 
               'day', {'loop_variables.mice_all(', 'mouse_iterator', ').days(:).name'}, 'day_iterator';
               'stack', {'loop_variables.mice_all(',  'mouse_iterator', ').days(', 'day_iterator', ').stacks'}, 'stack_iterator';
               'body_part', {'loop_variables.body_parts_all_position(:).name'}, 'body_part_iterator';
               'velocity_direction', {'loop_variables.body_parts_all_position(', 'body_part_iterator', ').directions(:).name'}, 'direction_iterator';
               'motorSpeed', {'loop_variables.motorSpeeds'}, 'motorSpeed_iterator';
               };

% Skip any files that don't exist (spontaneous or problem files)
parameters.load_abort_flag = true; 

% Dimension of different time range pairs.
parameters.rangePairs = 1; 

% 
parameters.segmentDim = 1;
parameters.concatDim = 2;

% Are the segmentations the same length? (If false, will put outputs into
% cell array)
parameters.uniformSegments = false;

% Input values. 
% Extracted timeseries.
parameters.loop_list.things_to_load.timeseries.dir = {[parameters.dir_exper 'behavior\body\normalized with std only\paw position normalized\'], 'mouse', '\', 'day', '\'};
parameters.loop_list.things_to_load.timeseries.filename= {'position', 'stack', '.mat'};
parameters.loop_list.things_to_load.timeseries.variable= {'position.', 'body_part', '.', 'velocity_direction'}; 
parameters.loop_list.things_to_load.timeseries.level = 'stack';
% Time ranges
parameters.loop_list.things_to_load.time_ranges.dir = {[parameters.dir_exper 'behavior\motorized\period instances\'], 'mouse', '\', 'day', '\'};
parameters.loop_list.things_to_load.time_ranges.filename= {'long_periods_', 'stack', '.mat'};
parameters.loop_list.things_to_load.time_ranges.variable= {'long_periods.walk_', 'motorSpeed'}; 
parameters.loop_list.things_to_load.time_ranges.level = 'stack';

% Output Values
parameters.loop_list.things_to_save.segmented_timeseries.dir = {[parameters.dir_exper 'behavior\body\normalized with std only\segmented positions\'], 'body_part', '\', 'velocity_direction', '\', 'mouse', '\', 'day', '\'};
parameters.loop_list.things_to_save.segmented_timeseries.filename= {'segmented_timeseries_longWalk_motorized', 'motorSpeed', '_', 'stack', '.mat'};
parameters.loop_list.things_to_save.segmented_timeseries.variable= {'segmented_timeseries'}; 
parameters.loop_list.things_to_save.segmented_timeseries.level = 'motorSpeed';

RunAnalysis({@SegmentTimeseriesData}, parameters);

%% Concatenate positions (long walk periods)
% use type_tag to avoid separating spon & motorized
if isfield(parameters, 'loop_list')
parameters = rmfield(parameters,'loop_list');
end

% Iterators
parameters.loop_list.iterators = {
                    'body_part', {'loop_variables.body_parts_all_position(:).name'}, 'body_part_iterator';
                    'velocity_direction', {'loop_variables.body_parts_all_position(', 'body_part_iterator', ').directions(:).name'}, 'direction_iterator';
                   'type_tag', {'loop_variables.type_tags'}, 'type_tag_iterator';
                   'mouse', {'loop_variables.mice_all(:).name'}, 'mouse_iterator'; 
                   'day', {'loop_variables.mice_all(', 'mouse_iterator', ').days(:).name'}, 'day_iterator';
                   'stack', {'getfield(loop_variables, {1}, "mice_all", {',  'mouse_iterator', '}, "days", {', 'day_iterator', '}, ', 'loop_variables.conditions_stack_locations_long{', 'type_tag_iterator', '})'}, 'stack_iterator'; 
                   };

parameters.concatDim = 1;
parameters.concatenate_across_cells = true;

% Inputs 
parameters.loop_list.things_to_load.data.dir = {[parameters.dir_exper 'behavior\body\normalized with std only\segmented positions\'], 'body_part', '\', 'velocity_direction', '\', 'mouse', '\', 'day', '\'};
parameters.loop_list.things_to_load.data.filename= {'segmented_timeseries_', 'type_tag', '_', 'stack', '.mat'};
parameters.loop_list.things_to_load.data.variable= {'segmented_timeseries'}; 
parameters.loop_list.things_to_load.data.level = 'stack';

% Outputs
parameters.loop_list.things_to_save.concatenated_data.dir = {[parameters.dir_exper 'behavior\body\normalized with std only\concatenated positions\'], 'body_part', '\', 'velocity_direction', '\', 'mouse', '\'};
parameters.loop_list.things_to_save.concatenated_data.filename = {'concatenated_positions_', 'type_tag', '.mat'};
parameters.loop_list.things_to_save.concatenated_data.variable = {'position_all'}; 
parameters.loop_list.things_to_save.concatenated_data.level = 'mouse';

RunAnalysis({@ConcatenateData}, parameters);

%% Segment positions by stride -- own x-direction 
% To find stride lengths in space/distance 
if isfield(parameters, 'loop_list')
parameters = rmfield(parameters,'loop_list');
end

% Is so you can use a single loop for calculations. 
parameters.loop_list.iterators = {
                'body_part', {'loop_variables.body_parts_all_position(:).name'}, 'body_part_iterator';
                    'velocity_direction', {'loop_variables.body_parts_all_position(', 'body_part_iterator', ').directions(:).name'}, 'direction_iterator';
               'mouse', {'loop_variables.mice_all(:).name'}, 'mouse_iterator';
               'type_tag', {'loop_variables.type_tags'}, 'type_tag_iterator';
                };

parameters.segmentDim = 1;
parameters.concatDim = 2;
parameters.instancesAsCells = true; 
parameters.uniformSegments = false;
parameters.add_extra_cell_layer = false;

% Inputs
% timeseries
parameters.loop_list.things_to_load.timeseries.dir = {[parameters.dir_exper 'behavior\body\normalized with std only\concatenated positions\'], 'body_part', '\', 'velocity_direction', '\', 'mouse', '\'};
parameters.loop_list.things_to_load.timeseries.filename = {'concatenated_positions_', 'type_tag', '.mat'};
parameters.loop_list.things_to_load.timeseries.variable = {'position_all'}; 
parameters.loop_list.things_to_load.timeseries.level = 'type_tag';
%time ranges
parameters.loop_list.things_to_load.time_ranges.dir = {[parameters.dir_exper 'behavior\gait analysis\normalized with std only\x peaks\all periods\'], 'body_part', '\', 'mouse', '\'};
parameters.loop_list.things_to_load.time_ranges.filename = {'x_peaks_', 'type_tag', '.mat'};
parameters.loop_list.things_to_load.time_ranges.variable = {'x_peaks.depression_ranges{1}'}; 
parameters.loop_list.things_to_load.time_ranges.level = 'type_tag';

% Outputs
% segmented timseries
parameters.loop_list.things_to_save.segmented_timeseries.dir = {[parameters.dir_exper 'behavior\gait analysis\normalized with std only\position strides from own x depressions\all periods\'],'body_part', '\', 'velocity_direction', '\' 'mouse', '\'};
parameters.loop_list.things_to_save.segmented_timeseries.filename = {'stride_segmentations.mat'};
parameters.loop_list.things_to_save.segmented_timeseries.variable = {'stride_segmentations_depression{', 'type_tag_iterator', ', 1}'}; 
parameters.loop_list.things_to_save.segmented_timeseries.level = 'mouse';

RunAnalysis({@StrideSegmentationLooper}, parameters);

%% Segment positions by stride -- FL x (might not do?)

%% plot position strides
% plot x vs y 
if isfield(parameters, 'loop_list')
parameters = rmfield(parameters,'loop_list');
end

% Is so you can use a single loop for calculations. 
parameters.loop_list.iterators = {
               'paw', {'loop_variables.paws'}, 'paw_iterator';
               'mouse', {'loop_variables.mice_all(:).name'}, 'mouse_iterator';
               'type_tag', {'loop_variables.type_tags'}, 'type_tag_iterator';
                };
parameters.resampleLength = 10;
parameters.closeFigures = true;
% Inputs
% x data
parameters.loop_list.things_to_load.x.dir = {[parameters.dir_exper 'behavior\gait analysis\normalized with std only\position strides from own x depressions\all periods\'],'paw', '\x\' 'mouse', '\'};
parameters.loop_list.things_to_load.x.filename = {'stride_segmentations.mat'};
parameters.loop_list.things_to_load.x.variable = {'stride_segmentations_depression{', 'type_tag_iterator', ', 1}'}; 
parameters.loop_list.things_to_load.x.level = 'mouse';
% y data
parameters.loop_list.things_to_load.y.dir = {[parameters.dir_exper 'behavior\gait analysis\normalized with std only\position strides from own x depressions\all periods\'],'paw', '\y\' 'mouse', '\'};
parameters.loop_list.things_to_load.y.filename = {'stride_segmentations.mat'};
parameters.loop_list.things_to_load.y.variable = {'stride_segmentations_depression{', 'type_tag_iterator', ', 1}'}; 
parameters.loop_list.things_to_load.y.level = 'mouse';

% Outputs
parameters.loop_list.things_to_save.x_together.dir = {[parameters.dir_exper 'behavior\gait analysis\normalized with std only\position strides from own x depressions\figures\'],'paw', '\y\' 'mouse', '\'};
parameters.loop_list.things_to_save.x_together.filename = {'stride_segmentations_x_together.mat'};
parameters.loop_list.things_to_save.x_together.variable = {'stride_segmentations_x_together{', 'type_tag_iterator', ', 1}'}; 
parameters.loop_list.things_to_save.x_together.level = 'mouse';

parameters.loop_list.things_to_save.y_together.dir = {[parameters.dir_exper 'behavior\gait analysis\normalized with std only\position strides from own x depressions\figures\'],'paw', '\y\' 'mouse', '\'};
parameters.loop_list.things_to_save.y_together.filename = {'stride_segmentations_y_together.mat'};
parameters.loop_list.things_to_save.y_together.variable = {'stride_segmentations_y_together{', 'type_tag_iterator', ', 1}'}; 
parameters.loop_list.things_to_save.y_together.level = 'mouse';

parameters.loop_list.things_to_save.x_resampled.dir = {[parameters.dir_exper 'behavior\gait analysis\normalized with std only\position strides from own x depressions\figures\'],'paw', '\y\' 'mouse', '\'};
parameters.loop_list.things_to_save.x_resampled.filename = {'stride_segmentations_x_resampled.mat'};
parameters.loop_list.things_to_save.x_resampled.variable = {'stride_segmentations_x_resampled{', 'type_tag_iterator', ', 1}'}; 
parameters.loop_list.things_to_save.x_resampled.level = 'mouse';

parameters.loop_list.things_to_save.y_resampled.dir = {[parameters.dir_exper 'behavior\gait analysis\normalized with std only\position strides from own x depressions\figures\'],'paw', '\y\' 'mouse', '\'};
parameters.loop_list.things_to_save.y_resampled.filename = {'stride_segmentations_y_resampled.mat'};
parameters.loop_list.things_to_save.y_resampled.variable = {'stride_segmentations_y_resampled{', 'type_tag_iterator', ', 1}'}; 
parameters.loop_list.things_to_save.y_resampled.level = 'mouse';

parameters.loop_list.things_to_save.average.dir = {[parameters.dir_exper 'behavior\gait analysis\normalized with std only\position strides from own x depressions\figures\'],'paw', '\y\' 'mouse', '\'};
parameters.loop_list.things_to_save.average.filename = {'stride_segmentations_average.mat'};
parameters.loop_list.things_to_save.average.variable = {'stride_segmentations_average{', 'type_tag_iterator', ', 1}'}; 
parameters.loop_list.things_to_save.average.level = 'mouse';

parameters.loop_list.things_to_save.std_dev.dir = {[parameters.dir_exper 'behavior\gait analysis\normalized with std only\position strides from own x depressions\figures\'],'paw', '\y\' 'mouse', '\'};
parameters.loop_list.things_to_save.std_dev.filename = {'stride_segmentations_std_dev.mat'};
parameters.loop_list.things_to_save.std_dev.variable = {'stride_segmentations_std_dev{', 'type_tag_iterator', ', 1}'}; 
parameters.loop_list.things_to_save.std_dev.level = 'mouse';

% parameters.loop_list.things_to_save.fig_positions_together.dir = {[parameters.dir_exper 'behavior\gait analysis\normalized with std only\position strides from own x depressions\figures\'],'paw', '\y\' 'mouse', '\'};
% parameters.loop_list.things_to_save.fig_positions_together.filename = {'stride_segmentations_together.fig'};
% parameters.loop_list.things_to_save.fig_positions_together.variable = {'fig_positions_together'}; 
% parameters.loop_list.things_to_save.fig_positions_together.level = 'type_tag';

parameters.loop_list.things_to_save.fig_resampled.dir = {[parameters.dir_exper 'behavior\gait analysis\normalized with std only\position strides from own x depressions\figures\'],'paw', '\y\' 'mouse', '\'};
parameters.loop_list.things_to_save.fig_resampled.filename = {'stride_segmentations_resampled.fig'};
parameters.loop_list.things_to_save.fig_resampled.variable = {'fig_resampled'}; 
parameters.loop_list.things_to_save.fig_resampled.level = 'type_tag';

parameters.loop_list.things_to_save.fig_average.dir = {[parameters.dir_exper 'behavior\gait analysis\normalized with std only\position strides from own x depressions\figures\'],'paw', '\y\' 'mouse', '\'};
parameters.loop_list.things_to_save.fig_average.filename = {'stride_segmentations_average.fig'};
parameters.loop_list.things_to_save.fig_average.variable = {'fig_average'}; 
parameters.loop_list.things_to_save.fig_average.level = 'type_tag';

RunAnalysis({@PlotStridePosition}, parameters);

parameters.closeFigures = false;

%% Average paw trajectories across mice 
if isfield(parameters, 'loop_list')
parameters = rmfield(parameters,'loop_list');
end

% Is so you can use a single loop for calculations. 
parameters.loop_list.iterators = {
               'paw', {'loop_variables.paws'}, 'paw_iterator';
               'type_tag', {'loop_variables.type_tags'}, 'type_tag_iterator';
               'mouse', {'loop_variables.mice_all(:).name'}, 'mouse_iterator';
                };
paramters.concatDim = 3;
parameters.averageDim = 3; 

% Inputs 
% average trajectory per mouse
parameters.loop_list.things_to_load.data_organized.dir = {[parameters.dir_exper 'behavior\gait analysis\normalized with std only\position strides from own x depressions\figures\'],'paw', '\y\' 'mouse', '\'};
parameters.loop_list.things_to_load.data_organized.filename = {'stride_segmentations_average.mat'};
parameters.loop_list.things_to_load.data_organized.variable = {'stride_segmentations_average{', 'type_tag_iterator', ', 1}'}; 
parameters.loop_list.things_to_load.data_organized.level = 'mouse';

% Outputs
% average trajectory

% std of trajectory

% concatenated trajectory

%% Get stride lengths (from stride position segmentations)
% load stride_positions_together
if isfield(parameters, 'loop_list')
parameters = rmfield(parameters,'loop_list');
end

% Is so you can use a single loop for calculations. 
parameters.loop_list.iterators = {
               'paw', {'loop_variables.paws'}, 'paw_iterator';
               'mouse', {'loop_variables.mice_all(:).name'}, 'mouse_iterator';
               'type_tag', {'loop_variables.type_tags'}, 'type_tag_iterator';
                };

parameters.evaluation_instructions = {{'data_evaluated = cellfun( @range, parameters.data);'}};
% Input
parameters.loop_list.things_to_load.data.dir = {[parameters.dir_exper 'behavior\gait analysis\normalized with std only\position strides from own x depressions\figures\'],'paw', '\y\' 'mouse', '\'};
parameters.loop_list.things_to_load.data.filename = {'stride_segmentations_x_together.mat'};
parameters.loop_list.things_to_load.data.variable = {'stride_segmentations_x_together{', 'type_tag_iterator', ', 1}'}; 
parameters.loop_list.things_to_load.data.level = 'mouse';

% Outputs
parameters.loop_list.things_to_save.data_evaluated.dir = {[parameters.dir_exper 'behavior\gait analysis\normalized with std only\position strides from own x depressions\figures\'],'paw', '\' 'mouse', '\'};
parameters.loop_list.things_to_save.data_evaluated.filename = {'stride_lengths.mat'};
parameters.loop_list.things_to_save.data_evaluated.variable = {'stride_lengths{', 'type_tag_iterator', ', 1}'}; 
parameters.loop_list.things_to_save.data_evaluated.level = 'mouse';

RunAnalysis({@EvaluateOnData}, parameters);

%% Stride lengths: plot histograms
if isfield(parameters, 'loop_list')
parameters = rmfield(parameters,'loop_list');
end

parameters.loop_list.iterators = {
               'paw', {'loop_variables.paws'}, 'paw_iterator';
               'mouse', {'loop_variables.mice_all(:).name'}, 'mouse_iterator';
               'period', {'loop_variables.periods_longsOnly'}, 'period_iterator';
               }; 

parameters.bin_edges = 0:0.5:10;
parameters.title_string = 'lengths';

% Inputs
parameters.loop_list.things_to_load.data.dir = {[parameters.dir_exper 'behavior\gait analysis\normalized with std only\position strides from own x depressions\figures\'],'paw', '\' 'mouse', '\'};
parameters.loop_list.things_to_load.data.filename = {'stride_lengths', '.mat'};
parameters.loop_list.things_to_load.data.variable = {'stride_lengths{', 'period_iterator', ', 1}'}; 
parameters.loop_list.things_to_load.data.level = 'mouse';

% Outputs
parameters.loop_list.things_to_save.fig.dir = {[parameters.dir_exper 'behavior\gait analysis\normalized with std only\stride lengths\'],'paw', '\', 'mouse', '\'};
parameters.loop_list.things_to_save.fig.filename = {'stride_lengths_', 'period', '_', 'period_iterator', '.fig'};
parameters.loop_list.things_to_save.fig.variable = {'fig'}; 
parameters.loop_list.things_to_save.fig.level = 'period';

RunAnalysis({@StrideDurationHistograms}, parameters);

close all;

%% Stride durations: Organize for Prism 
if isfield(parameters, 'loop_list')
parameters = rmfield(parameters,'loop_list');
end

parameters.loop_list.iterators = {
               'paw', {'loop_variables.paws'}, 'paw_iterator';
               'mouse', {'loop_variables.mice_all(:).name'}, 'mouse_iterator';
               'period', {'loop_variables.periods_longsOnly'}, 'period_iterator';
               }; 

parameters.number_of_mice = 7;
parameters.number_of_periods = 5;
parameters.instancesDim = 1;
parameters.NisMice = true;

% Inputs
parameters.loop_list.things_to_load.data.dir = {[parameters.dir_exper 'behavior\gait analysis\normalized with std only\stride durations\'],'paw', '\' 'mouse', '\'};
parameters.loop_list.things_to_load.data.filename = {'stride_durations', '.mat'};
parameters.loop_list.things_to_load.data.variable = {'stride_durations{', 'period_iterator', ', 1}'}; 
parameters.loop_list.things_to_load.data.level = 'mouse';

% Outputs
parameters.loop_list.things_to_save.data_organized_out.dir = {[parameters.dir_exper 'behavior\gait analysis\normalized with std only\stride durations\organized for Prism\']};
parameters.loop_list.things_to_save.data_organized_out.filename = {'stride_durations_organized_', 'paw', '.mat'};
parameters.loop_list.things_to_save.data_organized_out.variable = {'stride_durations_organized'}; 
parameters.loop_list.things_to_save.data_organized_out.level = 'paw';

RunAnalysis({@organize_gait_data_forPrism}, parameters)

%% Stride lengths: Organize for Prism 
if isfield(parameters, 'loop_list')
parameters = rmfield(parameters,'loop_list');
end

parameters.loop_list.iterators = {
               'paw', {'loop_variables.paws'}, 'paw_iterator';
               'mouse', {'loop_variables.mice_all(:).name'}, 'mouse_iterator';
               'period', {'loop_variables.periods_longsOnly'}, 'period_iterator';
               }; 

parameters.number_of_mice = 7;
parameters.number_of_periods = 5;
parameters.instancesDim = 1;
parameters.NisMice = true;

% Inputs
parameters.loop_list.things_to_load.data.dir = {[parameters.dir_exper 'behavior\gait analysis\normalized with std only\position strides from own x depressions\figures\'],'paw', '\' 'mouse', '\'};
parameters.loop_list.things_to_load.data.filename = {'stride_lengths', '.mat'};
parameters.loop_list.things_to_load.data.variable = {'stride_lengths{', 'period_iterator', ', 1}'}; 
parameters.loop_list.things_to_load.data.level = 'mouse';
% Outputs
parameters.loop_list.things_to_save.data_organized_out.dir = {[parameters.dir_exper 'behavior\gait analysis\normalized with std only\position strides from own x depressions\organized for Prism\']};
parameters.loop_list.things_to_save.data_organized_out.filename = {'stride_lengths_organized_', 'paw', '.mat'};
parameters.loop_list.things_to_save.data_organized_out.variable = {'stride_lengths_organized'}; 
parameters.loop_list.things_to_save.data_organized_out.level = 'paw';

RunAnalysis({@organize_gait_data_forPrism}, parameters)

%% Count number of occurances of long motorized walk periods 
% across speeds, per mouse
if isfield(parameters, 'loop_list')
parameters = rmfield(parameters,'loop_list');
end

parameters.loop_list.iterators = {
               'mouse', {'loop_variables.mice_all(:).name'}, 'mouse_iterator'; % don't include m1100
               'speed', {'loop_variables.motorSpeeds'}, 'speed_iterator';
               }; 

parameters.evaluation_instructions = {{'data = parameters.data;'...
                                        'data_evaluated = size(data,1);'
                                        };
                                       {};
                                       {'data = parameters.data;'...
                                        'data_evaluated = sum(data);';
                                       };
                                       };

parameters.concatDim = 1; 
parameters.concatenation_level = 'speed';

% Inputs
parameters.loop_list.things_to_load.data.dir = {[parameters.dir_exper 'behavior\body\normalized with std only\concatenated velocity\FL\total_magnitude\motorized\'] ,'mouse', '\' };
parameters.loop_list.things_to_load.data.filename = {'concatenated_velocity_longPeriods_walk_', 'speed', '.mat'};
parameters.loop_list.things_to_load.data.variable = {'velocity_all'}; 
parameters.loop_list.things_to_load.data.level = 'speed';

% Outputs 
parameters.loop_list.things_to_save.data_evaluated.dir = {[parameters.dir_exper 'behavior\motorized\long periods counts\within mice\'], 'mouse', '\'};
parameters.loop_list.things_to_save.data_evaluated.filename = {'counts_walk.mat'};
parameters.loop_list.things_to_save.data_evaluated.variable = {'counts'}; 
parameters.loop_list.things_to_save.data_evaluated.level = 'mouse';

parameters.loop_list.things_to_rename = {{'data_evaluated', 'data'}, {'concatenated_data', 'data'}};

RunAnalysis({@EvaluateOnData, @ConcatenateData, @EvaluateOnData}, parameters);

%% walk counts across mice
if isfield(parameters, 'loop_list')
parameters = rmfield(parameters,'loop_list');
end

% Iterations.
parameters.loop_list.iterators = {   
               'mouse', {'loop_variables.mice_all(:).name'}, 'mouse_iterator';
               };

parameters.concatDim = 1; 
parameters.concatenation_level = 'mouse';
parameters.averageDim = 1; 
% Inputs
parameters.loop_list.things_to_load.data.dir = {[parameters.dir_exper 'behavior\motorized\long periods counts\within mice\'], 'mouse', '\'};
parameters.loop_list.things_to_load.data.filename= {'counts_walk.mat'};
parameters.loop_list.things_to_load.data.variable= {'counts'}; 
parameters.loop_list.things_to_load.data.level = 'mouse';

% Outputs
% concatenated across mice
parameters.loop_list.things_to_save.concatenated_data.dir = {[parameters.dir_exper 'behavior\motorized\long periods counts\across mice\']};
parameters.loop_list.things_to_save.concatenated_data.filename= {'counts_allmice_walk.mat'};
parameters.loop_list.things_to_save.concatenated_data.variable= {'counts_all'}; 
parameters.loop_list.things_to_save.concatenated_data.level = 'end';
% average
parameters.loop_list.things_to_save.average.dir = {[parameters.dir_exper 'behavior\motorized\long periods counts\across mice\']};
parameters.loop_list.things_to_save.average.filename= {'counts_averageAcrossMice_walk.mat'};
parameters.loop_list.things_to_save.average.variable= {'counts_averageAcrossMice'}; 
parameters.loop_list.things_to_save.average.level = 'end';
% std dev
parameters.loop_list.things_to_save.std_dev.dir = {[parameters.dir_exper 'behavior\motorized\long periods counts\across mice\']};
parameters.loop_list.things_to_save.std_dev.filename= {'counts_stdAcrossMice_walk.mat'};
parameters.loop_list.things_to_save.std_dev.variable= {'counts_stdAcrossMice'}; 
parameters.loop_list.things_to_save.std_dev.level = 'end';

parameters.loop_list.things_to_rename = {{'concatenated_data', 'data'}};

RunAnalysis({@ConcatenateData, @AverageData}, parameters)




%% motorized rest counts across mice 
% to get the counts, extract the timesereies of 1 body part
% (similar to paw_velocity_pipeline_code.m)
if isfield(parameters, 'loop_list')
parameters = rmfield(parameters,'loop_list');
end

% Iterators
parameters.loop_list.iterators = {'mouse', {'loop_variables.mice_all(:).name'}, 'mouse_iterator'; 
               'day', {'loop_variables.mice_all(', 'mouse_iterator', ').days(:).name'}, 'day_iterator';
                   'stack', {'loop_variables.mice_all(',  'mouse_iterator', ').days(', 'day_iterator', ').stacks'}, 'stack_iterator';
                   'body_part', {'loop_variables.body_parts_all(1).name'}, 'body_part_iterator';
                   'velocity_direction', {'loop_variables.body_parts_all(', 'body_part_iterator', ').directions(1).name'}, 'direction_iterator';
                   };

% Skip any files that don't exist (spontaneous or problem files)
parameters.load_abort_flag = true; 

% Dimension of different time range pairs.
parameters.rangePairs = 1; 

% 
parameters.segmentDim = 1;
parameters.concatDim = 2;

% Are the segmentations the same length? (If false, will put outputs into
% cell array)
parameters.uniformSegments = false;

% Input values. 
% Extracted timeseries.
parameters.loop_list.things_to_load.timeseries.dir = {[parameters.dir_exper 'behavior\body\normalized with std only\paw velocity normalized\'], 'mouse', '\', 'day', '\'};
parameters.loop_list.things_to_load.timeseries.filename= {'velocity', 'stack', '.mat'};
parameters.loop_list.things_to_load.timeseries.variable= {'velocity.', 'body_part', '.', 'velocity_direction'}; 
parameters.loop_list.things_to_load.timeseries.level = 'stack';
% Time ranges
parameters.loop_list.things_to_load.time_ranges.dir = {[parameters.dir_exper 'behavior\motorized\period instances\'], 'mouse', '\', 'day', '\'};
parameters.loop_list.things_to_load.time_ranges.filename= {'long_periods_', 'stack', '.mat'};
parameters.loop_list.things_to_load.time_ranges.variable= {'long_periods.rest'}; 
parameters.loop_list.things_to_load.time_ranges.level = 'stack';

% Output Values
parameters.loop_list.things_to_save.segmented_timeseries.dir = {[parameters.dir_exper 'behavior\body\normalized with std only\segmented velocities\'], 'body_part', '\', 'velocity_direction', '\motorized\', 'mouse', '\', 'day', '\'};
parameters.loop_list.things_to_save.segmented_timeseries.filename= {'segmented_timeseries_longPeriods_rest_', 'stack', '.mat'};
parameters.loop_list.things_to_save.segmented_timeseries.variable= {'segmented_timeseries'}; 
parameters.loop_list.things_to_save.segmented_timeseries.level = 'velocity_direction';

RunAnalysis({@SegmentTimeseriesData}, parameters);


%% concatenate rest velocities
if isfield(parameters, 'loop_list')
parameters = rmfield(parameters,'loop_list');
end

% Iterators
parameters.loop_list.iterators = {
               'body_part', {'loop_variables.body_parts_all(1).name'}, 'body_part_iterator';
                'velocity_direction', {'loop_variables.body_parts_all(', 'body_part_iterator', ').directions(1).name'}, 'velocity_direction_iterator';
               'mouse', {'loop_variables.mice_all(:).name'}, 'mouse_iterator'; 
               %'motorSpeed', {'loop_variables.motorSpeeds'}, 'motorSpeed_iterator';
               'day', {'loop_variables.mice_all(', 'mouse_iterator', ').days(:).name'}, 'day_iterator';
               'stack', {'loop_variables.mice_all(',  'mouse_iterator', ').days(', 'day_iterator', ').stacks'}, 'stack_iterator';
                    };

parameters.concatDim = 1;
%parameters.concatenation_level = 'day';
parameters.concatenate_across_cells = true;

% Input
parameters.loop_list.things_to_load.data.dir = {[parameters.dir_exper 'behavior\body\normalized with std only\segmented velocities\'], 'body_part', '\', 'velocity_direction', '\motorized',  '\' 'mouse', '\', 'day', '\'};
parameters.loop_list.things_to_load.data.filename = {'segmented_timeseries_longPeriods_rest_', 'stack', '.mat'};
parameters.loop_list.things_to_load.data.variable = {'segmented_timeseries'}; 
parameters.loop_list.things_to_load.data.level = 'stack';

% Output
parameters.loop_list.things_to_save.concatenated_data.dir = {[parameters.dir_exper 'behavior\body\normalized with std only\concatenated velocity\'], 'body_part', '\', 'velocity_direction', '\motorized\', 'mouse', '\'};
parameters.loop_list.things_to_save.concatenated_data.filename = {'concatenated_velocity_longPeriods_rest.mat'};
parameters.loop_list.things_to_save.concatenated_data.variable = {'velocity_all'}; 
parameters.loop_list.things_to_save.concatenated_data.level = 'mouse';

RunAnalysis({@ConcatenateData}, parameters);

%% now count the rest numbers
if isfield(parameters, 'loop_list')
parameters = rmfield(parameters,'loop_list');
end

% Iterations.
parameters.loop_list.iterators = {     
               'mouse', {'loop_variables.mice_all(:).name'}, 'mouse_iterator'; % don't include mouse 1100
               };

parameters.evaluation_instructions = {{'data = parameters.data;'...
                                      'data_evaluated = size(data, 1);'
                                       }};
parameters.concatDim = 1; 
parameters.concatenation_level = 'mouse';
parameters.averageDim = 1; 
% Inputs
parameters.loop_list.things_to_load.data.dir = {[parameters.dir_exper 'behavior\body\normalized with std only\concatenated velocity\FL\total_magnitude\motorized\'], 'mouse', '\'};
parameters.loop_list.things_to_load.data.filename= {'concatenated_velocity_longPeriods_rest.mat'};
parameters.loop_list.things_to_load.data.variable= {'velocity_all'}; 
parameters.loop_list.things_to_load.data.level = 'mouse';

% Outputs
% concatenated across mice
parameters.loop_list.things_to_save.concatenated_data.dir = {[parameters.dir_exper 'behavior\motorized\long periods counts\across mice\']};
parameters.loop_list.things_to_save.concatenated_data.filename= {'counts_allmice_rest.mat'};
parameters.loop_list.things_to_save.concatenated_data.variable= {'counts_all'}; 
parameters.loop_list.things_to_save.concatenated_data.level = 'end';
% average
parameters.loop_list.things_to_save.average.dir = {[parameters.dir_exper 'behavior\motorized\long periods counts\across mice\']};
parameters.loop_list.things_to_save.average.filename= {'counts_averageAcrossMice_rest.mat'};
parameters.loop_list.things_to_save.average.variable= {'counts_averageAcrossMice'}; 
parameters.loop_list.things_to_save.average.level = 'end';
% std dev
parameters.loop_list.things_to_save.std_dev.dir = {[parameters.dir_exper 'behavior\motorized\long periods counts\across mice\']};
parameters.loop_list.things_to_save.std_dev.filename= {'counts_stdAcrossMice_rest.mat'};
parameters.loop_list.things_to_save.std_dev.variable= {'counts_stdAcrossMice'}; 
parameters.loop_list.things_to_save.std_dev.level = 'end';

parameters.loop_list.things_to_rename = {{'data_evaluated', 'data' }, {'concatenated_data', 'data'}};

RunAnalysis({@EvaluateOnData, @ConcatenateData, @AverageData}, parameters)
