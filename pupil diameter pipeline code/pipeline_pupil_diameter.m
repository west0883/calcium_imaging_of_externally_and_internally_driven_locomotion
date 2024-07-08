% pipeline_pupil_diameter.m
% Sarah West
% 6/14/22 

% Use "create_mice_all.m" before using this.

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
% If you want to change the list of stacks, use ListStacks function.
% Ex: numberVector=2:12; digitNumber=2;
% Ex cont: stackList=ListStacks(numberVector,digitNumber); 
% Ex cont: mice_all(1).stacks(1)=stackList;

parameters.mice_all = parameters.mice_all;    

% Give the number of digits that should be included in each stack number.
parameters.digitNumber=2; 

% *************************************************************************
% Parameters

% Sampling frequency of eye video, in fps.
parameters.eye_fps = 20;

% Sampling frequency of collected brain data (per channel), in Hz or frames per
% second.
parameters.fps= 20; 

% Number of channels from brain data (need this to calculate correct
% "skip" time length).
parameters.channelNumber = 2;

% Number of frames you recorded from brain and want to keep (don't make chunks longer than this)  
% (after skipped frames are removed)
parameters.frames=6000; 

% Number of initial brain frames to skip, allows for brightness/image
% stabilization of camera. Need this to know how much to skip in the
% behavior.
parameters.skip = 1200; 

% Load names of motorized periods
load([parameters.dir_exper 'periods_nametable.mat']);
periods_motorized = periods;

% Load names of spontaneous periods
load([parameters.dir_exper 'periods_nametable_spontaneous.mat']);
periods_spontaneous = periods(1:6, :);
clear periods; 

% Create a shared motorized & spontaneous list.
periods_bothConditions = [periods_motorized; periods_spontaneous]; 
parameters.periods_bothConditions = periods_bothConditions.condition;
parameters.periods_motorized = periods_motorized;
parameters.periods_spontaneous = periods_spontaneous;

% Loop variables.
parameters.loop_variables.mice_all = parameters.mice_all;
parameters.loop_variables.conditions = {'motorized'; 'spontaneous'};
parameters.loop_variables.conditions_stack_locations = {'stacks'; 'spontaneous'};
parameters.loop_variables.periods_bothConditions = parameters.periods_bothConditions;

% If it exists, load mice_all_no_missing_data.m
if isfile([parameters.dir_exper 'behavior\eye\mice_all_no_missing_data.mat'])
    load([parameters.dir_exper 'behavior\eye\mice_all_no_missing_data.mat']);

    % put into loop_variables
    parameters.loop_variables.mice_all_no_missing_data = mice_all_no_missing_data;

end

%% Import DeepLabCut pupil extraction data. 
% Calls ImportDLCPupilData.m, but that function doesn't really do anything,
% as RunAnalysis can import it in with a load function.

% Always clear loop list first. 
if isfield(parameters, 'loop_list')
parameters = rmfield(parameters,'loop_list');
end

% Iterators
parameters.loop_list.iterators = {
               'mouse', {'loop_variables.mice_all(:).name'}, 'mouse_iterator'; 
               'day', {'loop_variables.mice_all(', 'mouse_iterator', ').days(:).name'}, 'day_iterator';
               'stack_name', { 'dir("Y:\Sarah\Data\Random Motorized Treadmill\', 'day', '\', 'mouse', '\eye\eye*filtered.csv").name'}, 'stack_name_iterator'; 
               };

% Abort analysis if there's no corresponding file.
parameters.load_abort_flag = true; 

% Input values
parameters.loop_list.things_to_load.import_in.dir = {'Y:\Sarah\Data\Random Motorized Treadmill\', 'day', '\', 'mouse', '\eye\'};
parameters.loop_list.things_to_load.import_in.filename= {'stack_name'}; 
parameters.loop_list.things_to_load.import_in.variable= {'trial_in'}; 
parameters.loop_list.things_to_load.import_in.level = 'stack_name';
parameters.loop_list.things_to_load.import_in.load_function = @importdata;

% Output
parameters.loop_list.things_to_save.import_out.dir = {[parameters.dir_exper 'behavior\eye\extracted pupil tracking\'], 'mouse', '\', 'day', '\'};
parameters.loop_list.things_to_save.import_out.filename= {'trial', 'stack_name', '.mat'};
parameters.loop_list.things_to_save.import_out.variable= {'trial'}; 
parameters.loop_list.things_to_save.import_out.level = 'stack_name';

RunAnalysis({@ImportDLCPupilData}, parameters)

%% Search for data that wasn't imported

% Iterators
parameters.loop_list.iterators = {
               'mouse', {'loop_variables.mice_all(:).name'}, 'mouse_iterator'; 
               'day', {'loop_variables.mice_all(', 'mouse_iterator', ').days(:).name'}, 'day_iterator';
               'condition', {'loop_variables.conditions'}, 'condition_iterator';
               'stack', {'getfield(loop_variables, {1}, "mice_all", {',  'mouse_iterator', '}, "days", {', 'day_iterator', '}, ', 'loop_variables.conditions_stack_locations{', 'condition_iterator', '})'}, 'stack_iterator'; 
               };


% Input values
parameters.loop_list.things_to_check.dir = {[parameters.dir_exper 'behavior\eye\extracted pupil tracking\'], 'mouse', '\', 'day', '\'};
parameters.loop_list.things_to_check.filename= {'trialeye', 'stack', '*.mat'};  

% Output
parameters.loop_list.missing_data.dir = {[parameters.dir_exper 'behavior\eye\']};
parameters.loop_list.missing_data.filename= {'missing_data.mat'};

SearchForData(parameters)

%% Make a mice_all that doesn't have the missing eye data in it.
% create_mice_all_no_missing_eyedata.m


%% Fit circles to pupil edges
% Using mice_all_no_missing_data
% Always clear loop list first. 
if isfield(parameters, 'loop_list')
parameters = rmfield(parameters,'loop_list');
end

% Iterators
parameters.loop_list.iterators = {
               'mouse', {'loop_variables.mice_all_no_missing_data(:).name'}, 'mouse_iterator'; 
               'day', {'loop_variables.mice_all_no_missing_data(', 'mouse_iterator', ').days(:).name'}, 'day_iterator';
               'condition', {'loop_variables.conditions'}, 'condition_iterator';
               'stack', {'getfield(loop_variables, {1}, "mice_all_no_missing_data", {',  'mouse_iterator', '}, "days", {', 'day_iterator', '}, ', 'loop_variables.conditions_stack_locations{', 'condition_iterator', '})'}, 'stack_iterator'; 
               };

% Number of points plotted
parameters.numberOfPoints = 8;

% Abort analysis if there's no corresponding file.
parameters.load_abort_flag = true; 

% Columns of position data.
parameters.xPositionColumns = 2:3: parameters.numberOfPoints * 3 + 1;
parameters.yPositionColumns = 3:3: parameters.numberOfPoints * 3 + 1;
parameters.confidenceColumns = 4:3: parameters.numberOfPoints * 3 + 1;

% Confidence thrshold that the tracking confidence of a point needs to be
% over to be included in the diameter calculation.
parameters.confidenceThreshold = 0.3;

% Info about format of the input data. Should always be the same if using
% "importdata.m" on DLC output .csv files in step above. 
parameters.framesDim = 1;
parameters.dataDim = 2;

% Inputs
parameters.loop_list.things_to_load.DLC_data.dir = {[parameters.dir_exper 'behavior\eye\extracted pupil tracking\'], 'mouse', '\', 'day', '\'};
parameters.loop_list.things_to_load.DLC_data.filename= {'trialeye', 'stack', '*.mat'}; % Has a variable middle part to the name
parameters.loop_list.things_to_load.DLC_data.variable= {'trial.data'}; 
parameters.loop_list.things_to_load.DLC_data.level = 'stack';

% Outputs
parameters.loop_list.things_to_save.circle_info.dir = {[parameters.dir_exper 'behavior\eye\pupil diameters\'], 'mouse', '\', 'day', '\'};
parameters.loop_list.things_to_save.circle_info.filename= {'trial', 'stack', '.mat'};
parameters.loop_list.things_to_save.circle_info.variable= {'trial'}; 
parameters.loop_list.things_to_save.circle_info.level = 'stack';

RunAnalysis({@FitCircles}, parameters);

%% Find max pupil diameter per day
% Concatenate the max pupil diameters of each stack in a day, then take
% maximum. (Max is taken after each concatenation, but only the last is
% saved.)
% Use mice_all_no_missing_data
if isfield(parameters, 'loop_list')
parameters = rmfield(parameters,'loop_list');
end

% Iterators   
% Both motorized & spontaneous stacks are concatenated together.
parameters.loop_list.iterators = {
               'mouse', {'loop_variables.mice_all_no_missing_data(:).name'}, 'mouse_iterator'; 
               'day', {'loop_variables.mice_all_no_missing_data(', 'mouse_iterator', ').days(:).name'}, 'day_iterator';
               'stack', {'[loop_variables.mice_all_no_missing_data(',  'mouse_iterator', ').days(', 'day_iterator', ').stacks; loop_variables.mice_all_no_missing_data(',  'mouse_iterator', ').days(', 'day_iterator', ').spontaneous]'}, 'stack_iterator';
                 };

parameters.concatDim = 1;
parameters.concatenation_level = 'stack';
parameters.evaluation_instructions = {{};
                                      {'data_evaluated = max(parameters.data);'}};

% Input
parameters.loop_list.things_to_load.data.dir = {[parameters.dir_exper 'behavior\eye\pupil diameters\'], 'mouse', '\', 'day', '\'};
parameters.loop_list.things_to_load.data.filename = {'trial', 'stack', '.mat'};
parameters.loop_list.things_to_load.data.variable = {'trial.max_diameter'}; 
parameters.loop_list.things_to_load.data.level = 'stack';

% Outputs
parameters.loop_list.things_to_save.data_evaluated.dir = {[parameters.dir_exper 'behavior\eye\pupil diameters\'], 'mouse', '\', 'day', '\'};
parameters.loop_list.things_to_save.data_evaluated.filename = {'day_max_diameter.mat'};
parameters.loop_list.things_to_save.data_evaluated.variable = {'max_diameter'}; 
parameters.loop_list.things_to_save.data_evaluated.level = 'day';

parameters.loop_list.things_to_rename = {{'concatenated_data', 'data'}};

RunAnalysis({@ConcatenateData, @EvaluateOnData}, parameters);

parameters = rmfield(parameters, 'concatenation_level');

%% Normalize pupil diameters by proportion of max per day
% Always clear loop list first. 
% Use mice_all_nomissing_data
if isfield(parameters, 'loop_list')
parameters = rmfield(parameters,'loop_list');
end

% Iterators   
% Both motorized & spontaneous stacks are concatenated together.
parameters.loop_list.iterators = {
               'mouse', {'loop_variables.mice_all_no_missing_data(:).name'}, 'mouse_iterator'; 
               'day', {'loop_variables.mice_all_no_missing_data(', 'mouse_iterator', ').days(:).name'}, 'day_iterator';
               'stack', {'[loop_variables.mice_all_no_missing_data(',  'mouse_iterator', ').days(', 'day_iterator', ').stacks; loop_variables.mice_all_no_missing_data(',  'mouse_iterator', ').days(', 'day_iterator', ').spontaneous]'}, 'stack_iterator';
                 };

parameters.evaluation_instructions = {{'data_evaluated = parameters.data./parameters.max_diameter;'}};

% Inputs
% max diameter per day
parameters.loop_list.things_to_load.max_diameter.dir = {[parameters.dir_exper 'behavior\eye\pupil diameters\'], 'mouse', '\', 'day', '\'};
parameters.loop_list.things_to_load.max_diameter.filename = {'day_max_diameter.mat'};
parameters.loop_list.things_to_load.max_diameter.variable = {'max_diameter'}; 
parameters.loop_list.things_to_load.max_diameter.level = 'day';
% timeseries of pupil diameters
parameters.loop_list.things_to_load.data.dir = {[parameters.dir_exper 'behavior\eye\pupil diameters\'], 'mouse', '\', 'day', '\'};
parameters.loop_list.things_to_load.data.filename= {'trial', 'stack', '.mat'};
parameters.loop_list.things_to_load.data.variable= {'trial.diameter'}; 
parameters.loop_list.things_to_load.data.level = 'stack';

% Ouputs
parameters.loop_list.things_to_save.data_evaluated.dir = {[parameters.dir_exper 'behavior\eye\pupil diameters normalized\'], 'mouse', '\', 'day', '\'};
parameters.loop_list.things_to_save.data_evaluated.filename = {'diameters', 'stack', '.mat'};
parameters.loop_list.things_to_save.data_evaluated.variable = {'diameters'}; 
parameters.loop_list.things_to_save.data_evaluated.level = 'stack';

RunAnalysis({@EvaluateOnData}, parameters);

%% Pad short stacks with NaNs.
% Sometimes the behavior cameras were 50-100 frames short, but we still
% want the data that DID get collected. Without this, the segmentation
% steps will throw errors. 

% Not running with RunAnalysis because you don't have to load each of them
% to check their length.

% Folder/filenames you're checking.
parameters.filename_forcheck =  {[parameters.dir_exper 'behavior\eye\pupil diameters normalized\'], 'mouse', '\', 'day', '\diameters', 'stack', '.mat'};

parameters.loop_list.iterators = {
               'mouse', {'loop_variables.mice_all_no_missing_data(:).name'}, 'mouse_iterator'; 
               'day', {'loop_variables.mice_all_no_missing_data(', 'mouse_iterator', ').days(:).name'}, 'day_iterator';
               'condition', {'loop_variables.conditions'}, 'condition_iterator';
               'stack', {'getfield(loop_variables, {1}, "mice_all_no_missing_data", {',  'mouse_iterator', '}, "days", {', 'day_iterator', '}, ', 'loop_variables.conditions_stack_locations{', 'condition_iterator', '})'}, 'stack_iterator'; 
               };

looping_output_list = LoopGenerator(parameters.loop_list, parameters.loop_variables); 

% For each element of looping_output_list, 
for itemi = 1:size(looping_output_list,1)

    % Get keywords, like in RunAnalysis
    parameters.keywords = [parameters.loop_list.iterators(:,1); parameters.loop_list.iterators(:,3)];
    
    % Get values, like in RunAnalysis
    parameters.values = cell(size(parameters.keywords));
    for i = 1: numel(parameters.keywords)
        parameters.values{i} = looping_output_list(itemi).(cell2mat(parameters.keywords(i)));
    end

    MessageToUser('Checking ', parameters);

    % Get the filename 
    filestring = CreateStrings(parameters.filename_forcheck, parameters.keywords, parameters.values);

    % Usin matfile objects, check size of the file
    mat_object = matfile(filestring);

    % If less than frames
    if size(mat_object.diameters, 1) < parameters.frames

        % Load.
        load(filestring, 'diameters');

        % Pad 
        short_number = parameters.frames - size(diameters,1);
        diameters = [diameters; NaN(short_number, 1)]; 

        % Tell user.
        disp(['Stack short by ' num2str(short_number) ' frames.']);

        % Save
        save(filestring, 'diameters');
    end
end

%% If a stack of pupil diameters is missing, a vector of NaNs is created.
% This is so the instances stay properly aligned with fluorescence data
% Checks in \eye\pupil diameters\, puts into \eye\pupil diameters normalized

missing_eye_data.m
%% Motirized: Segment by behavior
% Always clear loop list first. 
if isfield(parameters, 'loop_list')
parameters = rmfield(parameters,'loop_list');
end

% Iterators
parameters.loop_list.iterators = {'mouse', {'loop_variables.mice_all(:).name'}, 'mouse_iterator'; 
               'day', {'loop_variables.mice_all(', 'mouse_iterator', ').days(:).name'}, 'day_iterator';
                   'stack', {'loop_variables.mice_all(',  'mouse_iterator', ').days(', 'day_iterator', ').stacks'}, 'stack_iterator'};
parameters.loop_variables.periods_nametable = periods_motorized; 

% Skip any files that don't exist (spontaneous or problem files)
parameters.load_abort_flag = true; 

% Dimension of different time range pairs.
parameters.rangePairs = 1; 

% 
parameters.segmentDim = 1;
parameters.concatDim = 2;

% Input values. 
% Extracted timeseries.
parameters.loop_list.things_to_load.timeseries.dir = {[parameters.dir_exper 'behavior\eye\pupil diameters normalized\'], 'mouse', '\', 'day', '\'};
parameters.loop_list.things_to_load.timeseries.filename= {'diameters', 'stack', '.mat'};
parameters.loop_list.things_to_load.timeseries.variable= {'diameters'}; 
parameters.loop_list.things_to_load.timeseries.level = 'stack';
% Time ranges
parameters.loop_list.things_to_load.time_ranges.dir = {[parameters.dir_exper 'behavior\motorized\period instances table format\'], 'mouse', '\', 'day', '\'};
parameters.loop_list.things_to_load.time_ranges.filename= {'all_periods_', 'stack', '.mat'};
parameters.loop_list.things_to_load.time_ranges.variable= {'all_periods.time_ranges'}; 
parameters.loop_list.things_to_load.time_ranges.level = 'stack';

% Output Values
parameters.loop_list.things_to_save.segmented_timeseries.dir = {[parameters.dir_exper 'behavior\eye\segmented eye pupil diameters\motorized\'], 'mouse', '\', 'day', '\'};
parameters.loop_list.things_to_save.segmented_timeseries.filename= {'segmented_timeseries_', 'stack', '.mat'};
parameters.loop_list.things_to_save.segmented_timeseries.variable= {'segmented_timeseries'}; 
parameters.loop_list.things_to_save.segmented_timeseries.level = 'stack';

RunAnalysis({@SegmentTimeseriesData}, parameters);

%% Spontaneous: Segment by behavior
% Always clear loop list first. 
if isfield(parameters, 'loop_list')
parameters = rmfield(parameters,'loop_list');
end

% Iterators
parameters.loop_list.iterators = {'mouse', {'loop_variables.mice_all(:).name'}, 'mouse_iterator'; 
               'day', {'loop_variables.mice_all(', 'mouse_iterator', ').days(:).name'}, 'day_iterator';
                   'stack', {'loop_variables.mice_all(',  'mouse_iterator', ').days(', 'day_iterator', ').spontaneous'}, 'stack_iterator';
                   'period', {'loop_variables.periods_spontaneous{:}'}, 'period_iterator'};

parameters.loop_variables.periods_spontaneous = periods_spontaneous.condition; 

% Skip any files that don't exist (spontaneous or problem files)
parameters.load_abort_flag = true; 

% Dimension of different time range pairs.
parameters.rangePairs = 1; 

% 
parameters.segmentDim = 1;
parameters.concatDim = 2;

% Input values. 
% Extracted timeseries.
parameters.loop_list.things_to_load.timeseries.dir = {[parameters.dir_exper 'behavior\eye\pupil diameters normalized\'], 'mouse', '\', 'day', '\'};
parameters.loop_list.things_to_load.timeseries.filename= {'diameters', 'stack', '.mat'};
parameters.loop_list.things_to_load.timeseries.variable= {'diameters'}; 
parameters.loop_list.things_to_load.timeseries.level = 'stack';
% Time ranges
parameters.loop_list.things_to_load.time_ranges.dir = {[parameters.dir_exper 'behavior\spontaneous\segmented behavior periods\'], 'mouse', '\', 'day', '\'};
parameters.loop_list.things_to_load.time_ranges.filename= {'behavior_periods_', 'stack', '.mat'};
parameters.loop_list.things_to_load.time_ranges.variable= {'behavior_periods.', 'period'}; 
parameters.loop_list.things_to_load.time_ranges.level = 'stack';

% Output Values
% (Convert to cell format to be compatible with motorized in below code)
parameters.loop_list.things_to_save.segmented_timeseries.dir = {[parameters.dir_exper 'behavior\eye\segmented eye pupil diameters\spontaneous\'], 'mouse', '\', 'day', '\'};
parameters.loop_list.things_to_save.segmented_timeseries.filename= {'segmented_timeseries_', 'stack', '.mat'};
parameters.loop_list.things_to_save.segmented_timeseries.variable= {'segmented_timeseries{', 'period_iterator',',1}'}; 
parameters.loop_list.things_to_save.segmented_timeseries.level = 'stack';

RunAnalysis({@SegmentTimeseriesData}, parameters);

%% Look for stacks when number of instances don't match those of fluorescence.
% Don't need to load, so don't use RunAnalysis.

parameters.loop_list.iterators = {
               'mouse', {'loop_variables.mice_all(:).name'}, 'mouse_iterator'; 
               'day', {'loop_variables.mice_all(', 'mouse_iterator', ').days(:).name'}, 'day_iterator';
               'condition', {'loop_variables.conditions'}, 'condition_iterator';
               'stack', {'getfield(loop_variables, {1}, "mice_all", {',  'mouse_iterator', '}, "days", {', 'day_iterator', '}, ', 'loop_variables.conditions_stack_locations{', 'condition_iterator', '})'}, 'stack_iterator'; 
               };
parameters.checkingDim = 2;
parameters.check_againstDim = 3; 
parameters.use_substructure = false;
% Input values
parameters.loop_list.things_to_check.dir = {[parameters.dir_exper 'behavior\eye\segmented eye pupil diameters\'], 'condition', '\', 'mouse', '\', 'day', '\'};
parameters.loop_list.things_to_check.filename= {'segmented_timeseries_', 'stack', '.mat'};  
parameters.loop_list.things_to_check.variable = {'segmented_timeseries'};

parameters.loop_list.check_against.dir = {[parameters.dir_exper 'fluorescence analysis\segmented timeseries\'], 'condition', '\', 'mouse', '\', 'day', '\'};
parameters.loop_list.check_against.filename= {'segmented_timeseries_', 'stack', '.mat'};  
parameters.loop_list.check_against.variable = {'segmented_timeseries'};

% Output
parameters.loop_list.mismatched_data.dir = {[parameters.dir_exper 'behavior\eye\']};
parameters.loop_list.mismatched_data.filename= {'mismatched_data.mat'};

CheckSizes2(parameters);

%% Notes for removal.
% for these, it's okay to remove ending values to match fluorescence

% load mismatched data 
load('Y:\Sarah\Analysis\Experiments\Random Motorized Treadmill\behavior\eye\mismatched_data.mat')

% For each item in mismatched data, 
for itemi = 2:size(mismatched_data, 1)
   
    % load the pupil data
    load(['Y:\Sarah\Analysis\Experiments\Random Motorized Treadmill\behavior\eye\segmented eye pupil diameters\' mismatched_data{itemi, 3} '\' mismatched_data{itemi, 1} '\' mismatched_data{itemi, 2} '\segmented_timeseries_' mismatched_data{itemi, 4} '.mat'])
    pupil = segmented_timeseries;
    
    % load the flourescence data 
    load(['Y:\Sarah\Analysis\Experiments\Random Motorized Treadmill\fluorescence analysis\segmented timeseries\' mismatched_data{itemi, 3} '\' mismatched_data{itemi, 1} '\' mismatched_data{itemi, 2} '\segmented_timeseries_' mismatched_data{itemi, 4} '.mat'])
    fluor = segmented_timeseries;

    % find number of instances you should have
    % if empty, make it set = 0 (because 3rd dim will be 1)
    if  ~isempty(fluor{mismatched_data{itemi, 5}})
        correct_num = size(fluor{mismatched_data{itemi, 5}}, 3);
    else
        correct_num = 0;
    end
    pupil_new = pupil;
    pupil_new{mismatched_data{itemi, 5}} = pupil{mismatched_data{itemi, 5}}(:, 1:correct_num);

    % Save over previous 
    segmented_timeseries = pupil_new;
    save(['Y:\Sarah\Analysis\Experiments\Random Motorized Treadmill\behavior\eye\segmented eye pupil diameters\' mismatched_data{itemi, 3} '\' mismatched_data{itemi, 1} '\' mismatched_data{itemi, 2} '\segmented_timeseries_' mismatched_data{itemi, 4} '.mat'], 'segmented_timeseries')

end 
% %From first round: 
% '1087'	'011222'	'motorized'	'11'	180
% % --> get rid of the last instance, the fluoresence was short by ~20 frames
% load('Y:\Sarah\Analysis\Experiments\Random Motorized Treadmill\behavior\eye\segmented eye pupil diameters\motorized\1087\011222\segmented_timeseries_11.mat')
% segmented_timeseries{180} = segmented_timeseries{180}(:, 1:end-1);
% save('Y:\Sarah\Analysis\Experiments\Random Motorized Treadmill\behavior\eye\segmented eye pupil diameters\motorized\1087\011222\segmented_timeseries_11.mat', 'segmented_timeseries')
% 
% %'1088'	'010522'	'motorized'	'11'	180 
% % --> same as above
% load('Y:\Sarah\Analysis\Experiments\Random Motorized Treadmill\behavior\eye\segmented eye pupil diameters\motorized\1088\010522\segmented_timeseries_11.mat')
% segmented_timeseries{180} = segmented_timeseries{180}(:, 1:end-1);
% save('Y:\Sarah\Analysis\Experiments\Random Motorized Treadmill\behavior\eye\segmented eye pupil diameters\motorized\1088\010522\segmented_timeseries_11.mat', 'segmented_timeseries')
% 
% % '1088'	'011322'	'motorized'	'06'	180
% % short by 2 frames; remove last 2 instances. 
% load('Y:\Sarah\Analysis\Experiments\Random Motorized Treadmill\behavior\eye\segmented eye pupil diameters\motorized\1088\011322\segmented_timeseries_06.mat')
% segmented_timeseries{180} = segmented_timeseries{180}(:, 1:end-2);
% save('Y:\Sarah\Analysis\Experiments\Random Motorized Treadmill\behavior\eye\segmented eye pupil diameters\motorized\1088\011322\segmented_timeseries_06.mat', 'segmented_timeseries')


%% Concatenate within behavior (spon & motorized independently) 
% Always clear loop list first. 
if isfield(parameters, 'loop_list')
parameters = rmfield(parameters,'loop_list');
end

% Iterators
parameters.loop_list.iterators = {
               'condition', {'loop_variables.conditions'}, 'condition_iterator';
               'mouse', {'loop_variables.mice_all(:).name'}, 'mouse_iterator'; 
               'day', {'loop_variables.mice_all(', 'mouse_iterator', ').days(:).name'}, 'day_iterator';
               'stack', {'getfield(loop_variables, {1}, "mice_all", {',  'mouse_iterator', '}, "days", {', 'day_iterator', '}, ', 'loop_variables.conditions_stack_locations{', 'condition_iterator', '})'}, 'stack_iterator'; 
               };

% Dimension to concatenate the timeseries across.
parameters.concatDim = 2; 
parameters.concatenate_across_cells = false; 

% Clear any reshaping instructions 
if isfield(parameters, 'reshapeDims')
    parameters = rmfield(parameters,'reshapeDims');
end

% Input Values
parameters.loop_list.things_to_load.data.dir = {[parameters.dir_exper 'behavior\eye\segmented eye pupil diameters\'],'condition', '\' 'mouse', '\', 'day', '\'};
parameters.loop_list.things_to_load.data.filename= {'segmented_timeseries_', 'stack', '.mat'};
parameters.loop_list.things_to_load.data.variable= {'segmented_timeseries'}; 
parameters.loop_list.things_to_load.data.level = 'stack';

% Output values
parameters.loop_list.things_to_save.concatenated_data.dir = {[parameters.dir_exper 'behavior\eye\concatenated diameters\'], 'condition', '\', 'mouse', '\'};
parameters.loop_list.things_to_save.concatenated_data.filename= {'concatenated_diameters_all_periods.mat'};
parameters.loop_list.things_to_save.concatenated_data.variable= {'diameters'}; 
parameters.loop_list.things_to_save.concatenated_data.level = 'mouse';

RunAnalysis({@ConcatenateData}, parameters);

%% Concatenate spon & motorized into same cell array.
if isfield(parameters, 'loop_list')
parameters = rmfield(parameters,'loop_list');
end

% Is so you can use a single loop for calculations. 
parameters.loop_list.iterators = {'mouse', {'loop_variables.mice_all(:).name'}, 'mouse_iterator'; 
               'condition', 'loop_variables.conditions', 'condition_iterator';
                };

% Tell it to concatenate across cells, not within cells. 
parameters.concatenate_across_cells = true; 
parameters.concatDim = 1;

% Input Values 
parameters.loop_list.things_to_load.data.dir = {[parameters.dir_exper 'behavior\eye\concatenated diameters\'], 'condition', '\', 'mouse', '\'};
parameters.loop_list.things_to_load.data.filename= {'concatenated_diameters_all_periods.mat'};
parameters.loop_list.things_to_load.data.variable= {'diameters'}; 
parameters.loop_list.things_to_load.data.level = 'condition';

% Output values
parameters.loop_list.things_to_save.concatenated_data.dir = {[parameters.dir_exper 'behavior\eye\concatenated diameters\both conditions\'], 'mouse', '\'};
parameters.loop_list.things_to_save.concatenated_data.filename= {'concatenated_diameters_all_periods.mat'};
parameters.loop_list.things_to_save.concatenated_data.variable= {'diameters_all'}; 
parameters.loop_list.things_to_save.concatenated_data.level = 'mouse';

RunAnalysis({@ConcatenateData}, parameters);

parameters.concatenate_across_cells = false;

%% Roll diameter timeseries
% Always clear loop list first. 
if isfield(parameters, 'loop_list')
parameters = rmfield(parameters,'loop_list');
end

% Iterators
parameters.loop_list.iterators = {'mouse', {'loop_variables.mice_all(:).name'}, 'mouse_iterator'; 
               'period', {'loop_variables.periods'}, 'period_iterator';            
               };
parameters.loop_variables.mice_all = parameters.mice_all;
parameters.loop_variables.periods = periods_bothConditions.condition;

% Dimension to roll across (time dimension). Will automatically add new
% data to the last + 1 dimension. 
parameters.rollDim = 1; 

% Window and step sizes (in frames)
parameters.windowSize = 20;
parameters.stepSize = 5; 

% Input 
parameters.loop_list.things_to_load.data.dir = {[parameters.dir_exper 'behavior\eye\concatenated diameters\both conditions\'], 'mouse', '\'};
parameters.loop_list.things_to_load.data.filename= {'concatenated_diameters_all_periods.mat'};
parameters.loop_list.things_to_load.data.variable= {'diameters_all{', 'period_iterator', '}'}; 
parameters.loop_list.things_to_load.data.level = 'mouse';

% Output
parameters.loop_list.things_to_save.data_rolled.dir = {[parameters.dir_exper 'behavior\eye\rolled concatenated diameters\'], 'mouse', '\'};
parameters.loop_list.things_to_save.data_rolled.filename= {'diameters_rolled.mat'};
parameters.loop_list.things_to_save.data_rolled.variable= {'diameters_rolled{', 'period_iterator', ',1}'}; 
parameters.loop_list.things_to_save.data_rolled.level = 'mouse';

parameters.loop_list.things_to_save.roll_number.dir = {[parameters.dir_exper 'behavior\eye\rolled concatenated diameters\'], 'mouse', '\'};
parameters.loop_list.things_to_save.roll_number.filename= {'diameter_rolled_rollnumber.mat'};
parameters.loop_list.things_to_save.roll_number.variable= {'roll_number{', 'period_iterator', ',1}'}; 
parameters.loop_list.things_to_save.roll_number.level = 'mouse';

RunAnalysis({@RollData}, parameters);

%% Average diameter per roll & instance
% Permute so instances are in last dimension

% Always clear loop list first. 
if isfield(parameters, 'loop_list')
parameters = rmfield(parameters,'loop_list');
end

% Iterations.
parameters.loop_list.iterators = {'mouse', {'loop_variables.mice_all(:).name'}, 'mouse_iterator'; 
               'period', {'loop_variables.periods{:}'}, 'period_iterator'};

parameters.loop_variables.mice_all = parameters.mice_all;
parameters.loop_variables.periods = periods_bothConditions.condition;

% Permute data so instances are in last dimension 
parameters.DimOrder = [1, 3, 2];

% Dimension to average across 
parameters.averageDim  = 1; 

% Load & put in the "true" roll number there's supposed to be.
load([parameters.dir_exper 'behavior\eye\rolled concatenated diameters\1087\diameter_rolled_rollnumber.mat'], 'roll_number'); 
parameters.roll_number = roll_number;
clear roll_number;

% Evaluation instructions (put instances in last dimension)
parameters.evaluation_instructions = {{}; {};{'if size(parameters.data,1) ~= parameters.roll_number{', 'period_iterator', '};'...
                                        'data_evaluated = transpose(parameters.data);'...
                                        'else;'...
                                         'data_evaluated = parameters.data;'...
                                         'end'}};


parameters.loop_list.things_to_load.data.dir = {[parameters.dir_exper 'behavior\eye\rolled concatenated diameters\'], 'mouse', '\'};
parameters.loop_list.things_to_load.data.filename= {'diameters_rolled.mat'};
parameters.loop_list.things_to_load.data.variable= {'diameters_rolled{', 'period_iterator', '}'}; 
parameters.loop_list.things_to_load.data.level = 'mouse';

parameters.loop_list.things_to_save.data_evaluated.dir = {[parameters.dir_exper 'behavior\eye\rolled concatenated diameters\'], 'mouse', '\'};
parameters.loop_list.things_to_save.data_evaluated.filename= {'diameter_averaged_by_instance.mat'};
parameters.loop_list.things_to_save.data_evaluated.variable= {'diameter_averaged_by_instance{', 'period_iterator',',1}'}; 
parameters.loop_list.things_to_save.data_evaluated.level = 'mouse';

parameters.loop_list.things_to_rename = {{'data_permuted', 'data'}; 
                                         { 'average', 'data'}};
RunAnalysis({@PermuteData, @AverageData, @EvaluateOnData}, parameters);

%% Figure out how much of the pupil data is missing per mouse, per behavior
% Always clear loop list first. 
if isfield(parameters, 'loop_list')
parameters = rmfield(parameters,'loop_list');
end

% Iterations.
parameters.loop_list.iterators = {'mouse', {'loop_variables.mice_all(:).name'}, 'mouse_iterator'; 
               'period', {'loop_variables.periods{:}'}, 'period_iterator'};

parameters.loop_variables.mice_all = parameters.mice_all;
parameters.loop_variables.periods = periods_bothConditions.condition;

parameters.evaluation_instructions = {{'data_evaluated = sum(isnan(parameters.data),"all")/numel(parameters.data);'}}; 
                                        
% Input 
parameters.loop_list.things_to_load.data.dir = {[parameters.dir_exper 'behavior\eye\rolled concatenated diameters\'], 'mouse', '\'};
parameters.loop_list.things_to_load.data.filename= {'diameter_averaged_by_instance.mat'};
parameters.loop_list.things_to_load.data.variable= {'diameter_averaged_by_instance{', 'period_iterator',',1}'}; 
parameters.loop_list.things_to_load.data.level = 'mouse';

% Output 
parameters.loop_list.things_to_save.data_evaluated.dir = {[parameters.dir_exper 'behavior\eye\rolled concatenated diameters\'], 'mouse', '\'};
parameters.loop_list.things_to_save.data_evaluated.filename= {'diameter_isnan_ratio.mat'};
parameters.loop_list.things_to_save.data_evaluated.variable= {'diameter_isnan_ratio{', 'period_iterator',',1}'}; 
parameters.loop_list.things_to_save.data_evaluated.level = 'mouse';

RunAnalysis({@EvaluateOnData}, parameters);


%% Plot histograms of the ratio of missing data.
% Always clear loop list first. 
if isfield(parameters, 'loop_list')
parameters = rmfield(parameters,'loop_list');
end

% Iterations.
parameters.loop_list.iterators = {'mouse', {'loop_variables.mice_all(:).name'}, 'mouse_iterator'};

parameters.loop_variables.mice_all = parameters.mice_all;

% Input 
parameters.loop_list.things_to_load.data.dir = {[parameters.dir_exper 'behavior\eye\rolled concatenated diameters\'], 'mouse', '\'};
parameters.loop_list.things_to_load.data.filename= {'diameter_isnan_ratio.mat'};
parameters.loop_list.things_to_load.data.variable= {'diameter_isnan_ratio'}; 
parameters.loop_list.things_to_load.data.level = 'mouse';

% Output
parameters.loop_list.things_to_save.fig.dir = {[parameters.dir_exper 'behavior\eye\rolled concatenated diameters\'], 'mouse', '\'};
parameters.loop_list.things_to_save.fig.filename= {'diameter_isnan_ratio_histogram.fig'};
parameters.loop_list.things_to_save.fig.variable= {'fig'}; 
parameters.loop_list.things_to_save.fig.level = 'mouse';

RunAnalysis({@HistogramOfRatios}, parameters);

%% Pupil diameter for fluorescence PLSR
% Instead of rolling the diameters, reshape so each time point is its own instance
% (Is for fluorescence PLSR)

% Put into same cell array, to match other formatting
% Always clear loop list first. 
if isfield(parameters, 'loop_list')
parameters = rmfield(parameters,'loop_list');
end

% Iterators
parameters.loop_list.iterators = {'mouse', {'loop_variables.mice_all(:).name'}, 'mouse_iterator'; 
               'period', {'loop_variables.periods_bothConditions'}, 'period_iterator';            
               };

% one column, timepoints of each instance are together 
parameters.evaluation_instructions = {{'data_evaluated = transpose(reshape(parameters.data, [], 1));'}};
% Input 
parameters.loop_list.things_to_load.data.dir = {[parameters.dir_exper 'behavior\eye\concatenated diameters\both conditions\'], 'mouse', '\'};
parameters.loop_list.things_to_load.data.filename= {'concatenated_diameters_all_periods.mat'};
parameters.loop_list.things_to_load.data.variable= {'diameters_all{', 'period_iterator', '}'}; 
parameters.loop_list.things_to_load.data.level = 'mouse';

% Output
parameters.loop_list.things_to_save.data_evaluated.dir = {[parameters.dir_exper 'behavior\eye\diameter for fluorescence PLSR\'], 'mouse', '\'};
parameters.loop_list.things_to_save.data_evaluated.filename= {'diameter_forFluorescence.mat'};
parameters.loop_list.things_to_save.data_evaluated.variable= {'diameter_forFluorescence{', 'period_iterator', ', 1}'}; 
parameters.loop_list.things_to_save.data_evaluated.level = 'mouse';

RunAnalysis({@EvaluateOnData}, parameters);

%% Pupil diameter for fluorescence PLSR 
% Shorten the warning peridos

% Put into same cell array, to match other formatting
% Always clear loop list first. 
if isfield(parameters, 'loop_list')
parameters = rmfield(parameters,'loop_list');
end

% Iterators
parameters.loop_list.iterators = {'mouse', {'loop_variables.mice_all(:).name'}, 'mouse_iterator'; 
               'period', {'loop_variables.periods_bothConditions'}, 'period_iterator';            
               };

parameters.shorten_dimensions = '41:end, :';

% one column, timepoints of each instance are together 
parameters.evaluation_instructions = {{}; {'data_evaluated = transpose(reshape(parameters.data, [], 1));'}};
% Input
% indices to shorten.
parameters.loop_list.things_to_load.indices_to_shorten.dir = {[parameters.dir_exper 'PLSR Warning Periods\']};
parameters.loop_list.things_to_load.indices_to_shorten.filename= {'indices_to_shorten.mat'};
parameters.loop_list.things_to_load.indices_to_shorten.variable= {'indices_to_shorten_original_index'}; 
parameters.loop_list.things_to_load.indices_to_shorten.level = 'start';
% data
parameters.loop_list.things_to_load.data.dir = {[parameters.dir_exper 'behavior\eye\concatenated diameters\both conditions\'], 'mouse', '\'};
parameters.loop_list.things_to_load.data.filename= {'concatenated_diameters_all_periods.mat'};
parameters.loop_list.things_to_load.data.variable= {'diameters_all{', 'period_iterator', '}'}; 
parameters.loop_list.things_to_load.data.level = 'mouse';

% Output
parameters.loop_list.things_to_save.data_evaluated.dir = {[parameters.dir_exper 'behavior\eye\diameter for fluorescence PLSR Warning Periods\'], 'mouse', '\'};
parameters.loop_list.things_to_save.data_evaluated.filename= {'diameter_forFluorescence.mat'};
parameters.loop_list.things_to_save.data_evaluated.variable= {'diameter_forFluorescence{', 'period_iterator', ', 1}'}; 
parameters.loop_list.things_to_save.data_evaluated.level = 'mouse';

parameters.loop_list.things_to_rename = {{'data_shortened', 'data'}};

RunAnalysis({@ShortenFluorescenceWarningPeriods, @EvaluateOnData}, parameters);