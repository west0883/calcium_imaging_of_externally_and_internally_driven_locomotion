% pipeline_fluorescence_analysis.m

%% Initial Setup  
% Put all needed paramters in a structure called "parameters", which you
% can then easily feed into your functions. 
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

% Include stacks from a "spontaneous" field of mice_all?
parameters.use_spontaneous_also = true;

% Other parameters
parameters.digitNumber = 2;
parameters.yDim = 256;
parameters.xDim = 256;
number_of_sources = 32; 
parameters.number_of_sources = number_of_sources;

% Lower triangle only 
parameters.indices = find(tril(ones(number_of_sources), -1));

% Load names of motorized periods
load([parameters.dir_exper 'periods_nametable.mat']);
periods_motorized = periods;

% Load names of spontaneous periods
load([parameters.dir_exper 'periods_nametable_spontaneous.mat']);
periods_spontaneous = periods(1:6, :);
clear periods; 

% Create a shared motorized & spontaneous list.
periods_bothConditions = [periods_motorized; periods_spontaneous]; 

parameters.loop_variables.mice_all = parameters.mice_all;
parameters.loop_variables.transformations = {'not transformed'; 'Fisher transformed'};
parameters.loop_variables.conditions = {'motorized'; 'spontaneous'};
parameters.loop_variables.conditions_stack_locations = {'stacks'; 'spontaneous'};
parameters.loop_variables.periods_bothConditions = periods_bothConditions.condition;

%% Run fluorescence extraction. 
% Always clear loop list first. 
if isfield(parameters, 'loop_list')
parameters = rmfield(parameters,'loop_list');
end

% Iterators
parameters.loop_list.iterators = {'mouse', {'loop_variables.mice_all(:).name'}, 'mouse_iterator'; 
               'day', {'loop_variables.mice_all(', 'mouse_iterator', ').days(:).name'}, 'day_iterator';
               'stack', {'[loop_variables.mice_all(',  'mouse_iterator', ').days(', 'day_iterator', ').stacks; loop_variables.mice_all(',  'mouse_iterator', ').days(', 'day_iterator', ').spontaneous]'}, 'stack_iterator'};

% Dimension different sources are in
parameters.sourcesDim = 2; 

% If the mean timeseries should be weighted by the weights of pixels in the sources (default is uniform mask)
parameters.weightedMean = true; 

% Input values
% Source masks
parameters.loop_list.things_to_load.sources.dir = {[parameters.dir_exper 'spatial segmentation\500 SVD components\manual assignments\'], 'mouse', '\'};
parameters.loop_list.things_to_load.sources.filename= {'sources_reordered_masked.mat'};
parameters.loop_list.things_to_load.sources.variable= {'sources_masked'};
parameters.loop_list.things_to_load.sources.level = 'mouse';

% Preprocessed fluorescence data videos
parameters.loop_list.things_to_load.data.dir = {[parameters.dir_exper 'preprocessing\fully preprocessed stacks\'], 'mouse', '\', 'day', '\'};
parameters.loop_list.things_to_load.data.filename= {'data', 'stack', '.mat'};
parameters.loop_list.things_to_load.data.variable= {'data'}; 
parameters.loop_list.things_to_load.data.level = 'stack';

% Output values. 
parameters.loop_list.things_to_save.timeseries.dir = {[parameters.dir_exper 'fluorescence analysis\timeseries\'], 'mouse', '\', 'day', '\'};
parameters.loop_list.things_to_save.timeseries.filename= {'timeseries', 'stack', '.mat'};
parameters.loop_list.things_to_save.timeseries.variable= {'timeseries'}; 
parameters.loop_list.things_to_save.timeseries.level = 'stack';

% Run 
RunAnalysis({@ExtractFluorescenceTimeseries}, parameters);

%% Motorized: Segment fluorescence by behavior period
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
parameters.concatDim = 3;

% Input values. 
% Extracted timeseries.
parameters.loop_list.things_to_load.timeseries.dir = {[parameters.dir_exper 'fluorescence analysis\timeseries\'], 'mouse', '\', 'day', '\'};
parameters.loop_list.things_to_load.timeseries.filename= {'timeseries', 'stack', '.mat'};
parameters.loop_list.things_to_load.timeseries.variable= {'timeseries'}; 
parameters.loop_list.things_to_load.timeseries.level = 'stack';
% Time ranges
parameters.loop_list.things_to_load.time_ranges.dir = {[parameters.dir_exper 'behavior\motorized\period instances table format\'], 'mouse', '\', 'day', '\'};
parameters.loop_list.things_to_load.time_ranges.filename= {'all_periods_', 'stack', '.mat'};
parameters.loop_list.things_to_load.time_ranges.variable= {'all_periods.time_ranges'}; 
parameters.loop_list.things_to_load.time_ranges.level = 'stack';

% Output Values
parameters.loop_list.things_to_save.segmented_timeseries.dir = {[parameters.dir_exper 'fluorescence analysis\segmented timeseries\motorized\'], 'mouse', '\', 'day', '\'};
parameters.loop_list.things_to_save.segmented_timeseries.filename= {'segmented_timeseries_', 'stack', '.mat'};
parameters.loop_list.things_to_save.segmented_timeseries.variable= {'segmented_timeseries'}; 
parameters.loop_list.things_to_save.segmented_timeseries.level = 'stack';

RunAnalysis({@SegmentTimeseriesData}, parameters);

%% SPONTANEOUS-- Segment fluorescence by behavior period
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
parameters.concatDim = 3;

% Input values. 
% Extracted timeseries.
parameters.loop_list.things_to_load.timeseries.dir = {[parameters.dir_exper 'fluorescence analysis\timeseries\'], 'mouse', '\', 'day', '\'};
parameters.loop_list.things_to_load.timeseries.filename= {'timeseries', 'stack', '.mat'};
parameters.loop_list.things_to_load.timeseries.variable= {'timeseries'}; 
parameters.loop_list.things_to_load.timeseries.level = 'stack';
% Time ranges
parameters.loop_list.things_to_load.time_ranges.dir = {[parameters.dir_exper 'behavior\spontaneous\segmented behavior periods\'], 'mouse', '\', 'day', '\'};
parameters.loop_list.things_to_load.time_ranges.filename= {'behavior_periods_', 'stack', '.mat'};
parameters.loop_list.things_to_load.time_ranges.variable= {'behavior_periods.', 'period'}; 
parameters.loop_list.things_to_load.time_ranges.level = 'stack';

% Output Values
% (Convert to cell format to be compatible with motorized in below code)
parameters.loop_list.things_to_save.segmented_timeseries.dir = {[parameters.dir_exper 'fluorescence analysis\segmented timeseries\spontaneous\'], 'mouse', '\', 'day', '\'};
parameters.loop_list.things_to_save.segmented_timeseries.filename= {'segmented_timeseries_', 'stack', '.mat'};
parameters.loop_list.things_to_save.segmented_timeseries.variable= {'segmented_timeseries{', 'period_iterator',',1}'}; 
parameters.loop_list.things_to_save.segmented_timeseries.level = 'stack';

RunAnalysis({@SegmentTimeseriesData}, parameters);

%% Concatenate fluorescence by behavior per mouse 
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
parameters.concatDim = 3; 
parameters.concatenate_across_cells = false; 

% Clear any reshaping instructions 
if isfield(parameters, 'reshapeDims')
    parameters = rmfield(parameters,'reshapeDims');
end

% Input Values
parameters.loop_list.things_to_load.data.dir = {[parameters.dir_exper 'fluorescence analysis\segmented timeseries\'],'condition', '\' 'mouse', '\', 'day', '\'};
parameters.loop_list.things_to_load.data.filename= {'segmented_timeseries_', 'stack', '.mat'};
parameters.loop_list.things_to_load.data.variable= {'segmented_timeseries'}; 
parameters.loop_list.things_to_load.data.level = 'stack';

% Output values
parameters.loop_list.things_to_save.concatenated_data.dir = {[parameters.dir_exper 'fluorescence analysis\concatenated timeseries\'], 'condition', '\', 'mouse', '\'};
parameters.loop_list.things_to_save.concatenated_data.filename= {'concatenated_timeseries_all_periods.mat'};
parameters.loop_list.things_to_save.concatenated_data.variable= {'timeseries'}; 
parameters.loop_list.things_to_save.concatenated_data.level = 'mouse';

RunAnalysis({@ConcatenateData}, parameters);

%% Concatenate motorized & spontaneous together. 
% Is so you can use a single loop for calculations. 
parameters.loop_list.iterators = {'mouse', {'loop_variables.mice_all(:).name'}, 'mouse_iterator'; 
               'condition', 'loop_variables.conditions', 'condition_iterator';
                };

% Tell it to concatenate across cells, not within cells. 
parameters.concatenate_across_cells = true; 
parameters.concatDim = 1;
parameters.concatenation_level = 'condition';

% Input Values (use a trick to concatenate just the 2 conditions)
parameters.loop_list.things_to_load.data.dir = {[parameters.dir_exper 'fluorescence analysis\concatenated timeseries\'], 'condition', '\', 'mouse', '\'};
parameters.loop_list.things_to_load.data.filename= {'concatenated_timeseries_all_periods.mat'};
parameters.loop_list.things_to_load.data.variable= {'timeseries'}; 
parameters.loop_list.things_to_load.data.level = 'condition';

% Output values
parameters.loop_list.things_to_save.concatenated_data.dir = {[parameters.dir_exper 'fluorescence analysis\concatenated timeseries both conditions\'], 'mouse', '\'};
parameters.loop_list.things_to_save.concatenated_data.filename= {'concatenated_timeseries_all_periods.mat'};
parameters.loop_list.things_to_save.concatenated_data.variable= {'timeseries_all'}; 
parameters.loop_list.things_to_save.concatenated_data.level = 'mouse';

RunAnalysis({@ConcatenateData}, parameters);


%% [FROM HERE DOWN YOU CAN COMBINE MOTORIZED & SPONTANOUS SECTIONS]
% Because they're concatenated.

%% Take average of fluorescence by behavior 
% Always clear loop list first. 
if isfield(parameters, 'loop_list')
parameters = rmfield(parameters,'loop_list');
end

% Iterators
parameters.loop_list.iterators = {
                'mouse', {'loop_variables.mice_all(:).name'}, 'mouse_iterator'; 
                'period', {'loop_variables.periods'}, 'period_iterator';              
               };

parameters.loop_variables.periods = periods_bothConditions.condition; 
parameters.loop_variables.mice_all = parameters.mice_all;

% Dimension to average across
parameters.averageDim = 3; 

% Input 
parameters.loop_list.things_to_load.data.dir = {[parameters.dir_exper 'fluorescence analysis\concatenated timeseries both conditions\'], 'mouse', '\'};
parameters.loop_list.things_to_load.data.filename= {'concatenated_timeseries_all_periods.mat'};
parameters.loop_list.things_to_load.data.variable= {'timeseries_all{', 'period_iterator', '}'}; 
parameters.loop_list.things_to_load.data.level = 'mouse';

% Output
parameters.loop_list.things_to_save.average.dir = {[parameters.dir_exper 'fluorescence analysis\concatenated timeseries both conditions\'], 'mouse', '\'};
parameters.loop_list.things_to_save.average.filename= {'average_timeseries_all_periods_mean.mat'};
parameters.loop_list.things_to_save.average.variable= {'average{', 'period_iterator', ',1}'}; 
parameters.loop_list.things_to_save.average.level = 'mouse';

parameters.loop_list.things_to_save.std_dev.dir = {[parameters.dir_exper 'fluorescence analysis\concatenated timeseries both conditions\'], 'mouse', '\'};
parameters.loop_list.things_to_save.std_dev.filename= {'average_timeseries_all_periods_std.mat'};
parameters.loop_list.things_to_save.std_dev.variable= {'std_dev{', 'period_iterator', ',1}'}; 
parameters.loop_list.things_to_save.std_dev.level = 'mouse';

RunAnalysis({@AverageData}, parameters);

%% Make a "true" roll number vector to refer to, based on periods tables. 
% Window and step sizes (in frames)
parameters.windowSize = 20;
parameters.stepSize = 5; 
parameters.duration = periods_bothConditions.duration; 

windowSize = cell(size(parameters.duration));
stepSize = cell(size(parameters.duration));
windowSize(:) = {parameters.windowSize};
stepSize(:) = {parameters.stepSize}; 

roll_number = cellfun(@CountRolls, parameters.duration, windowSize, stepSize, 'UniformOutput',false);

save([parameters.dir_exper 'roll_number.mat'], 'roll_number');

clear roll_number;

%% Roll data 
% Always clear loop list first. 
if isfield(parameters, 'loop_list')
parameters = rmfield(parameters,'loop_list');
end

% Iterators
parameters.loop_list.iterators = {'mouse', {'loop_variables.mice_all(:).name'}, 'mouse_iterator'; 
               'period', {'loop_variables.periods'}, 'period_iterator';            
               };

parameters.loop_variables.periods = periods_bothConditions.condition; 

% Dimension to roll across (time dimension). Will automatically add new
% data to the last + 1 dimension. 
parameters.rollDim = 1; 

% Window and step sizes (in frames)
parameters.windowSize = 20;
parameters.stepSize = 5; 

% Input 
parameters.loop_list.things_to_load.data.dir = {[parameters.dir_exper 'fluorescence analysis\concatenated timeseries both conditions\'], 'mouse', '\'};
parameters.loop_list.things_to_load.data.filename= {'concatenated_timeseries_all_periods.mat'};
parameters.loop_list.things_to_load.data.variable= {'timeseries_all{', 'period_iterator', '}'}; 
parameters.loop_list.things_to_load.data.level = 'mouse';

% Output
parameters.loop_list.things_to_save.data_rolled.dir = {[parameters.dir_exper 'fluorescence analysis\rolled timeseries\'], 'mouse', '\'};
parameters.loop_list.things_to_save.data_rolled.filename= {'timeseries_rolled.mat'};
parameters.loop_list.things_to_save.data_rolled.variable= {'timeseries_rolled{', 'period_iterator', ',1}'}; 
parameters.loop_list.things_to_save.data_rolled.level = 'mouse';

parameters.loop_list.things_to_save.roll_number.dir = {[parameters.dir_exper 'fluorescence analysis\rolled timeseries\'], 'mouse', '\'};
parameters.loop_list.things_to_save.roll_number.filename= {'roll_number.mat'};
parameters.loop_list.things_to_save.roll_number.variable= {'roll_number{', 'period_iterator', ', 1}'}; 
parameters.loop_list.things_to_save.roll_number.level = 'mouse';

RunAnalysis({@RollData}, parameters);

%% Timeseries for fluorescence PLSR
% Instead of rolling the timeseries, reshape so each time point is its own instance
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

% permute so it's: time x instance x source 
% one column, timepoints of each instance are together 
parameters.evaluation_instructions = {{'holder = permute(parameters.data, [1 3 2]);' ...  
                                       'data_evaluated = squeeze(reshape(holder, [], 1, parameters.number_of_sources));'}};

% Input 
parameters.loop_list.things_to_load.data.dir = {[parameters.dir_exper 'fluorescence analysis\concatenated timeseries both conditions\'], 'mouse', '\'};
parameters.loop_list.things_to_load.data.filename= {'concatenated_timeseries_all_periods.mat'};
parameters.loop_list.things_to_load.data.variable= {'timeseries_all{', 'period_iterator', '}'}; 
parameters.loop_list.things_to_load.data.level = 'mouse';

% Output
parameters.loop_list.things_to_save.data_evaluated.dir = {[parameters.dir_exper 'fluorescence analysis\fluorescence for fluorescence PLSR\'], 'mouse', '\'};
parameters.loop_list.things_to_save.data_evaluated.filename= {'forFluorescence.mat'};
parameters.loop_list.things_to_save.data_evaluated.variable= {'forFluorescence{', 'period_iterator', ', 1}'}; 
parameters.loop_list.things_to_save.data_evaluated.level = 'mouse';

RunAnalysis({@EvaluateOnData}, parameters);

%% Timeseries for WARNING PERIOD fluorescence (trim first)
% Put into same cell array, to match other formatting
% Always clear loop list first. 
if isfield(parameters, 'loop_list')
parameters = rmfield(parameters,'loop_list');
end

% Iterators
parameters.loop_list.iterators = {'mouse', {'loop_variables.mice_all(:).name'}, 'mouse_iterator'; 
               'period', {'loop_variables.periods_bothConditions'}, 'period_iterator';            
               };

% shortening instructions 
parameters.shorten_dimensions = '41:end, :, :';

% data starts as time x source x instance
parameters.evaluation_instructions = {{};
                                     % permute so it's: time x instance x source 
                                     % one column, timepoints of each instance are together 
                                     {'holder = permute(parameters.data, [1 3 2]);' ...  
                                     'data_evaluated = squeeze(reshape(holder, [], 1, parameters.number_of_sources));'}};
% Input
% indices to shorten.
parameters.loop_list.things_to_load.indices_to_shorten.dir = {[parameters.dir_exper 'PLSR Warning Periods\']};
parameters.loop_list.things_to_load.indices_to_shorten.filename= {'indices_to_shorten.mat'};
parameters.loop_list.things_to_load.indices_to_shorten.variable= {'indices_to_shorten_original_index'}; 
parameters.loop_list.things_to_load.indices_to_shorten.level = 'start';
% data
parameters.loop_list.things_to_load.data.dir = {[parameters.dir_exper 'fluorescence analysis\concatenated timeseries both conditions\'], 'mouse', '\'};
parameters.loop_list.things_to_load.data.filename= {'concatenated_timeseries_all_periods.mat'};
parameters.loop_list.things_to_load.data.variable= {'timeseries_all{', 'period_iterator', '}'}; 
parameters.loop_list.things_to_load.data.level = 'mouse';

% Output
parameters.loop_list.things_to_save.data_evaluated.dir = {[parameters.dir_exper 'fluorescence analysis\fluorescence for fluorescence PLSR Warning Periods\'], 'mouse', '\'};
parameters.loop_list.things_to_save.data_evaluated.filename= {'forFluorescence.mat'};
parameters.loop_list.things_to_save.data_evaluated.variable= {'forFluorescence{', 'period_iterator', ', 1}'}; 
parameters.loop_list.things_to_save.data_evaluated.level = 'mouse';

parameters.loop_list.things_to_rename = {{'data_shortened', 'data'}};

RunAnalysis({@ShortenFluorescenceWarningPeriods, @EvaluateOnData}, parameters);

%% Correlate data
% Separating these into smaller files because correlating takes a long time
% & I want the progress to be saved periodically.

% Always clear loop list first. 
if isfield(parameters, 'loop_list')
parameters = rmfield(parameters,'loop_list');
end

% Iterators
parameters.loop_list.iterators = {'mouse', {'loop_variables.mice_all(:).name'}, 'mouse_iterator'; 
               'period', {'loop_variables.periods'}, 'period_iterator';            
               };

parameters.loop_variables.periods = periods_bothConditions.condition; 

% Dimension to correlate across (dimensions where different sources are). 
parameters.sourceDim = 2; 

% Time dimension (the dimension of the timeseries that will be correlated)
parameters.timeDim = 1; 

% Input 
parameters.loop_list.things_to_load.data.dir = {[parameters.dir_exper 'fluorescence analysis\rolled timeseries\'], 'mouse', '\'};
parameters.loop_list.things_to_load.data.filename= {'timeseries_rolled.mat'};
parameters.loop_list.things_to_load.data.variable= {'timeseries_rolled{', 'period_iterator', ',1}'}; 
parameters.loop_list.things_to_load.data.level = 'mouse';

% Output
parameters.loop_list.things_to_save.correlation.dir = {[parameters.dir_exper 'fluorescence analysis\correlations\not transformed\'], 'mouse', '\instances\'};
parameters.loop_list.things_to_save.correlation.filename= {'correlations.mat'};
parameters.loop_list.things_to_save.correlation.variable= {'correlations{', 'period_iterator', '}'}; 
parameters.loop_list.things_to_save.correlation.level = 'mouse';

RunAnalysis({@CorrelateTimeseriesData}, parameters);

%% Run Fisher z - transformation 
% Save as separate files so you match the structure of the correlation
% matrices for looping in the next steps. 

% Always clear loop list first. 
if isfield(parameters, 'loop_list')
parameters = rmfield(parameters,'loop_list');
end

% Iterators
parameters.loop_list.iterators = {'mouse', {'loop_variables.mice_all(:).name'}, 'mouse_iterator'; 
               'period', {'loop_variables.periods'}, 'period_iterator';            
               };

parameters.loop_variables.periods = periods_bothConditions.condition; 

% Input
parameters.loop_list.things_to_load.data.dir = {[parameters.dir_exper 'fluorescence analysis\correlations\not transformed\'], 'mouse', '\instances\'};
parameters.loop_list.things_to_load.data.filename= {'correlations.mat'};
parameters.loop_list.things_to_load.data.variable= {'correlations{', 'period_iterator', '}'}; 
parameters.loop_list.things_to_load.data.level = 'mouse';

% Output 
parameters.loop_list.things_to_save.data_transformed.dir = {[parameters.dir_exper 'fluorescence analysis\correlations\Fisher transformed\'], 'mouse', '\instances\'};
parameters.loop_list.things_to_save.data_transformed.filename= {'correlations.mat'};
parameters.loop_list.things_to_save.data_transformed.variable= {'correlations{', 'period_iterator', '}'}; 
parameters.loop_list.things_to_save.data_transformed.level = 'mouse';

RunAnalysis({@FisherTransform}, parameters);

%% Calculate contralateral-ipsalateral averages 
% Calculate on the Fisher transformed & on not transformed
% Always clear loop list first. 
if isfield(parameters, 'loop_list')
parameters = rmfield(parameters,'loop_list');
end

% Iterators
parameters.loop_list.iterators = {
               'transformation', 'loop_variables.transformations', 'transformation_iterator';
               'mouse', {'loop_variables.mice_all(:).name'}, 'mouse_iterator'; 
               'period', {'loop_variables.periods'}, 'period_iterator';            
               };

parameters.loop_variables.periods = periods_bothConditions.condition; 

% Input
parameters.loop_list.things_to_load.data.dir = {[parameters.dir_exper 'fluorescence analysis\correlations\'], 'transformation', '\', 'mouse', '\instances\'};
parameters.loop_list.things_to_load.data.filename= {'correlations.mat'};
parameters.loop_list.things_to_load.data.variable= {'correlations{', 'period_iterator', '}'}; 
parameters.loop_list.things_to_load.data.level = 'mouse';

% Output
parameters.loop_list.things_to_save.data_ipsacontra.dir = {[parameters.dir_exper 'fluorescence analysis\correlations\'] 'transformation','\', 'mouse', '\instances\'};
parameters.loop_list.things_to_save.data_ipsacontra.filename= {'correlations_IpsaContra.mat'};
parameters.loop_list.things_to_save.data_ipsacontra.variable= {'correlations_{', 'period_iterator', '}'}; 
parameters.loop_list.things_to_save.data_ipsacontra.level = 'mouse';

RunAnalysis({@IpsaContraAverage}, parameters);

%% 
% From here on, can run everything with a "transform" iterator -- "not
% transformed" or "Fisher transformed". (Because I want to see how the
% analyses look with & without transformation.)


%% Save reshaped data (2D + roll dim)
% You end up using this more than once, so might as well save it.
% Always clear loop list first. 
% Also permute so instances are in the last dimension. 
% Keep only lower triangle of matrix.

if isfield(parameters, 'loop_list')
parameters = rmfield(parameters,'loop_list');
end

% Iterators
parameters.loop_list.iterators = {
                'transformation', {'loop_variables.transformations'}, 'transformation_iterator';
               'mouse', {'loop_variables.mice_all(:).name'}, 'mouse_iterator'; 
               'period', {'loop_variables.periods'}, 'period_iterator';            
               };

parameters.loop_variables.periods = periods_bothConditions.condition; 

% Lower triangle only.
parameters.indices = find(tril(ones(number_of_sources), -1));

% Load & put in the "true" roll number there's supposed to be.
load([parameters.dir_exper 'roll_number.mat'], 'roll_number'); 
parameters.roll_number = roll_number;
clear roll_number;

% Variable/data you want to reshape. 
parameters.toReshape = {'parameters.data'}; 

% Dimensions for reshaping, before removing data & before cnocatenation.
% Turning it into 2 dims + roll dim. 
parameters.reshapeDims = {'{size(parameters.data, 1) * size(parameters.data,2), [], parameters.roll_number{', 'period_iterator', '},}'};

% Permute data instructions/dimensions. Puts instances in last dimension. 
parameters.DimOrder = [1, 3, 2]; 

% Load & put in roll number.

% Evaluation instructions.
parameters.evaluation_instructions = {{}, {},{'if ~isempty(parameters.data);'...
          'data_evaluated = parameters.data(parameters.indices, :,:);' ...
          'else;'...
          'data_evaluated = [];'...
          'end'
           }};

% Input 
parameters.loop_list.things_to_load.data.dir = {[parameters.dir_exper 'fluorescence analysis\correlations\'], 'transformation', '\', 'mouse', '\instances\'};
parameters.loop_list.things_to_load.data.filename= {'correlations.mat'};
parameters.loop_list.things_to_load.data.variable= {'correlations{', 'period_iterator', '}'}; 
parameters.loop_list.things_to_load.data.level = 'mouse';

% Output
parameters.loop_list.things_to_save.data_evaluated.dir = {[parameters.dir_exper 'fluorescence analysis\correlations\'], 'transformation', '\', 'mouse', '\instances reshaped\'};
parameters.loop_list.things_to_save.data_evaluated.filename= {'values.mat'};
parameters.loop_list.things_to_save.data_evaluated.variable= {'values{', 'period_iterator', ', 1}'}; 
parameters.loop_list.things_to_save.data_evaluated.level = 'mouse';

parameters.loop_list.things_to_rename = {{'data_reshaped', 'data'};
                                         {'data_permuted', 'data'}}; 

RunAnalysis({@ReshapeData, @PermuteData, @EvaluateOnData}, parameters); 

%% Reshape the ipsa-contra averaged data.

if isfield(parameters, 'loop_list')
parameters = rmfield(parameters,'loop_list');
end

% Iterators
parameters.loop_list.iterators = {
                'transformation', {'loop_variables.transformations'}, 'transformation_iterator';
               'mouse', {'loop_variables.mice_all(:).name'}, 'mouse_iterator'; 
               'period', {'loop_variables.periods'}, 'period_iterator';            
               };

parameters.loop_variables.periods = periods_bothConditions.condition; 

% Lower triangle only.
parameters.indices = find(tril(ones(number_of_sources), -1));

% Load & put in the "true" roll number there's supposed to be.
load([parameters.dir_exper 'roll_number.mat'], 'roll_number'); 
parameters.roll_number = roll_number;
clear roll_number;

% Variable/data you want to reshape. 
parameters.toReshape = {'parameters.data'}; 

% Dimensions for reshaping, before removing data & before cnocatenation.
% Turning it into 2 dims + roll dim. 
parameters.reshapeDims = {'{size(parameters.data, 1) * size(parameters.data,2), [], parameters.roll_number{', 'period_iterator', '},}'};

% Permute data instructions/dimensions. Puts instances in last dimension. 
parameters.DimOrder = [1, 3, 2]; 

% Load & put in roll number.

% Evaluation instructions.
parameters.evaluation_instructions = {{}, {},{'if ~isempty(parameters.data);'...
          'data_evaluated = parameters.data(parameters.indices, :,:);' ...
          'else;'...
          'data_evaluated = [];'...
          'end'
           }};

% Input 
parameters.loop_list.things_to_load.data.dir = {[parameters.dir_exper 'fluorescence analysis\correlations\'], 'transformation', '\', 'mouse', '\instances\'};
parameters.loop_list.things_to_load.data.filename= {'correlations_IpsaContra.mat'};
parameters.loop_list.things_to_load.data.variable= {'correlations_{', 'period_iterator', '}'}; 
parameters.loop_list.things_to_load.data.level = 'mouse';

% Output
parameters.loop_list.things_to_save.data_evaluated.dir = {[parameters.dir_exper 'fluorescence analysis\correlations\'], 'transformation', '\', 'mouse', '\instances reshaped\'};
parameters.loop_list.things_to_save.data_evaluated.filename= {'values_IpsaContra.mat'};
parameters.loop_list.things_to_save.data_evaluated.variable= {'values{', 'period_iterator', ', 1}'}; 
parameters.loop_list.things_to_save.data_evaluated.level = 'mouse';

parameters.loop_list.things_to_rename = {{'data_reshaped', 'data'};
                                         {'data_permuted', 'data'}}; 

RunAnalysis({@ReshapeData, @PermuteData, @EvaluateOnData}, parameters); 
