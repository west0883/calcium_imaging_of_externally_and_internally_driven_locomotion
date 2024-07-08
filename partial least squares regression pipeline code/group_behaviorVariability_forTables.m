
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
parameters.mice_all = parameters.mice_all;   %([1:3 5:end]);

% Other parameters
parameters.digitNumber = 2;
parameters.yDim = 256;
parameters.xDim = 256;
parameters.number_of_sources = 32; 
parameters.indices = 1:16;

% Load the motorized/spontaneous list of periods, to fit in with
% correlations pipeline 

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

% Load periods_nametable_PLSR.m, if it exists yet. (Otherwise is created in
% first step).
if isfile([parameters.dir_exper 'PLSR\periods_nametable_forPLSR_specials.mat'])
    load([parameters.dir_exper 'PLSR\periods_nametable_forPLSR_specials.mat']);
    parameters.periods = periods;

    % Also load the indices to remove
    load([parameters.dir_exper 'PLSR\indices_to_remove_Specials.mat']);
    parameters.indices_to_remove = indices_to_remove;

    % Load lists of response categories
    load([parameters.dir_exper 'PLSR\response_categories.mat']);
    parameters.loop_variables.response_categories = categories;
    parameters.categories = categories;

    clear periods indices_to_remove categories;

end

% Load comparisons for first level continuous, if it exists yet.
if isfile([parameters.dir_exper 'PLSR\comparisons_special_continuous.mat'])
    load([parameters.dir_exper 'PLSR\comparisons_special_continuous.mat']);
    parameters.comparisons_continuous = comparisons([1:5 7:end]);   
    parameters.loop_variables.comparisons_continuous = parameters.comparisons_continuous; 
    clear comparisons;
end

% Load comparisons for first level categorical, if it exists yet.
if isfile([parameters.dir_exper 'PLSR\comparisons_special_categorical.mat'])
    load([parameters.dir_exper 'PLSR\comparisons_special_categorical.mat']);
    parameters.comparisons_categorical = comparisons; %([9 10]);
    parameters.loop_variables.comparisons_categorical = parameters.comparisons_categorical;
    clear comparisons;
end

% Load list of variables to subtract from level 1 categoricals, if it
% exists yet.
if isfile([parameters.dir_exper 'PLSR\variablesToSubtract_level1_categorical.mat'])
    load([parameters.dir_exper 'PLSR\variablesToSubtract_level1_categorical.mat']);
    parameters.variablesToSubtract = variablesToSubtract;
    clear variablesToSubtract;
end
 
% Names of all continuous variables.
parameters.continuous_variable_names = {'speed', 'accel', 'duration', 'pupil_diameter', 'tail', 'nose', 'FL', 'HL', 'x'};

% Put relevant variables into loop_variables.
parameters.loop_variables.mice_all = parameters.mice_all;
parameters.loop_variables.periods = parameters.periods.condition; 
parameters.loop_variables.periods_bothConditions = parameters.periods_bothConditions.condition; 
parameters.loop_variables.categories.type = parameters.categories.type;
parameters.loop_variables.output_types = {'BETA'};% {'Cov', 'BETA'}; % For continuous variables, I want to also see the actual regressors

parameters.loop_variables.data_types = {'', 'fluorescence '};
parameters.average_and_std_together = false;

if isfield(parameters, 'loop_list')
parameters = rmfield(parameters,'loop_list');
end

%% fluorescence 
if isfield(parameters, 'loop_list')
    parameters = rmfield(parameters,'loop_list');
end
% Iterators
parameters.loop_list.iterators = {
               'comparison', {'loop_variables.comparisons_continuous(:).name'}, 'comparison_iterator';
               % 'data_type', {'loop_variables.data_types'}, 'data_type_iterator';
               }; 

parameters.evaluation_instructions = {{
                                       'data = transpose(parameters.data);'
                                        'data_evaluated = [data parameters.average parameters.std_dev] .* 100 ;'
                                      }};
parameters.concatDim = 1; 
% parameters.concatenate_across = 'comparison';
% Inputs 
% all the numbers
parameters.loop_list.things_to_load.data.dir = {[parameters.dir_exper 'PLSR fluorescence Specials Inverted Intercepts\results\total percent variance behavior\continuous'], '\', 'comparison','\'};
parameters.loop_list.things_to_load.data.filename= {'all_percent_var.mat'};
parameters.loop_list.things_to_load.data.variable= {'all_percent_var'}; 
parameters.loop_list.things_to_load.data.level = 'comparison';
% average
parameters.loop_list.things_to_load.average.dir =  {[parameters.dir_exper 'PLSR fluorescence Specials Inverted Intercepts\results\total percent variance behavior\continuous'], '\', 'comparison','\'};
parameters.loop_list.things_to_load.average.filename= {'average_percent_var.mat'};
parameters.loop_list.things_to_load.average.variable= {'average_percent_var'}; 
parameters.loop_list.things_to_load.average.level = 'comparison';
% std dev
parameters.loop_list.things_to_load.std_dev.dir =  {[parameters.dir_exper 'PLSR fluorescence Specials Inverted Intercepts\results\total percent variance behavior\continuous'], '\', 'comparison','\'};
parameters.loop_list.things_to_load.std_dev.filename= {'std_dev_percent_var.mat'};
parameters.loop_list.things_to_load.std_dev.variable= {'std_dev_percent_var'}; 
parameters.loop_list.things_to_load.std_dev.level = 'comparison';

% Outputs
parameters.loop_list.things_to_save.concatenated_data.dir =  {[parameters.dir_exper 'PLSR fluorescence Specials Inverted Intercepts\results\total percent variance behavior\continuous'], '\'};
parameters.loop_list.things_to_save.concatenated_data.filename= {'organized.mat'};
parameters.loop_list.things_to_save.concatenated_data.variable= {'organized'}; 
parameters.loop_list.things_to_save.concatenated_data.level = 'end';

parameters.loop_list.things_to_rename = {{'data_evaluated', 'data'}}; 

RunAnalysis({@EvaluateOnData, @ConcatenateData}, parameters);

%% reorder to match comparison order of table

% old position is the index, new position is value
old_order = [2; % walk
             1; % rest
              4;% start
              5 ;% stop 
             6  ;% finished stop
             3 ; % pre walk
             7 ; % accel
             9 ; % decel 
              8 ;% post accel
              10 ;% post decel
             ];
[~, new_order] = sort(old_order);

load([parameters.dir_exper 'PLSR fluorescence Specials Inverted Intercepts\results\total percent variance behavior\continuous\organized.mat'])
organized2 = organized(new_order, :);
save([parameters.dir_exper 'PLSR fluorescence Specials Inverted Intercepts\results\total percent variance behavior\continuous\organized2.mat'], 'organized2');

%% correlations
if isfield(parameters, 'loop_list')
    parameters = rmfield(parameters,'loop_list');
end
% Iterators
parameters.loop_list.iterators = {
               'comparison', {'loop_variables.comparisons_continuous(:).name'}, 'comparison_iterator';
               % 'data_type', {'loop_variables.data_types'}, 'data_type_iterator';
               }; 

parameters.evaluation_instructions = {{
                                       'data = transpose(parameters.data);'
                                        'data_evaluated = [data parameters.average parameters.std_dev] .* 100 ;'
                                      }};
parameters.concatDim = 1; 
% parameters.concatenate_across = 'comparison';
% Inputs 
% all the numbers
parameters.loop_list.things_to_load.data.dir = {[parameters.dir_exper 'PLSR Specials Inverted Intercepts\results\total percent variance behavior\continuous'], '\', 'comparison','\'};
parameters.loop_list.things_to_load.data.filename= {'all_percent_var.mat'};
parameters.loop_list.things_to_load.data.variable= {'all_percent_var'}; 
parameters.loop_list.things_to_load.data.level = 'comparison';
% average
parameters.loop_list.things_to_load.average.dir =  {[parameters.dir_exper 'PLSR Specials Inverted Intercepts\results\total percent variance behavior\continuous'], '\', 'comparison','\'};
parameters.loop_list.things_to_load.average.filename= {'average_percent_var.mat'};
parameters.loop_list.things_to_load.average.variable= {'average_percent_var'}; 
parameters.loop_list.things_to_load.average.level = 'comparison';
% std dev
parameters.loop_list.things_to_load.std_dev.dir =  {[parameters.dir_exper 'PLSR Specials Inverted Intercepts\results\total percent variance behavior\continuous'], '\', 'comparison','\'};
parameters.loop_list.things_to_load.std_dev.filename= {'std_dev_percent_var.mat'};
parameters.loop_list.things_to_load.std_dev.variable= {'std_dev_percent_var'}; 
parameters.loop_list.things_to_load.std_dev.level = 'comparison';

% Outputs
parameters.loop_list.things_to_save.concatenated_data.dir =  {[parameters.dir_exper 'PLSR Specials Inverted Intercepts\results\total percent variance behavior\continuous'], '\'};
parameters.loop_list.things_to_save.concatenated_data.filename= {'organized.mat'};
parameters.loop_list.things_to_save.concatenated_data.variable= {'organized'}; 
parameters.loop_list.things_to_save.concatenated_data.level = 'end';

parameters.loop_list.things_to_rename = {{'data_evaluated', 'data'}}; 

RunAnalysis({@EvaluateOnData, @ConcatenateData}, parameters);

%% reorder to match comparison order of table

% old position is the index, new position is value
old_order = [2; % walk
             1; % rest
              4;% start
              5 ;% stop 
             6  ;% finished stop
             3 ; % pre walk
             7 ; % accel
             9 ; % decel 
              8 ;% post accel
              10 ;% post decel
             ];
[~, new_order] = sort(old_order);

load([parameters.dir_exper 'PLSR Specials Inverted Intercepts\results\total percent variance behavior\continuous\organized.mat'])
organized2 = organized(new_order, :);
save([parameters.dir_exper 'PLSR Specials Inverted Intercepts\results\total percent variance behavior\continuous\organized2.mat'], 'organized2');