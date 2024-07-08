% pipeline_PLSR_variability.m
% Sarah West
% 10/3/23

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
parameters.mice_all = parameters.mice_all;

% Other parameters
parameters.digitNumber = 2;
parameters.yDim = 256;
parameters.xDim = 256;
parameters.number_of_sources = 32; 
parameters.indices = find(tril(ones(parameters.number_of_sources), -1));

% Load periods_nametable_PLSR.m, if it exists yet. (Otherwise is created in
% first step).
if isfile([parameters.dir_exper 'PLSR\periods_nametable_forPLSR_specials_correlations.mat'])
    load([parameters.dir_exper 'PLSR\periods_nametable_forPLSR_specials_correlations.mat']);
    parameters.periods = periods;

    % Also load the indices to remove
    load([parameters.dir_exper 'PLSR\indices_to_remove_Specials.mat']);
    parameters.indices_to_remove = indices_to_remove;

    clear periods indices_to_remove categories;

end

% Load comparisons for first level continuous
load([parameters.dir_exper 'PLSR\comparisons_continuous_inverted_intercepts.mat']);
parameters.comparisons_continuous = comparisons;
clear comparisons;


% make variablesToUse2, without the "_vector" in it

for i = 1:size(parameters.comparisons_continuous, 2)

    old = parameters.comparisons_continuous(i).variablesToUse;

    old_str = repmat({'_vector'}, 1, numel(old));
    new_str = repmat({''}, 1, numel(old));

    new = cellfun(@strrep, old, old_str, new_str, 'UniformOutput',false);
    
    parameters.comparisons_continuous(i).variablesToUse2 = new; 

end 
clear old old_str new new_str i


% Names of all continuous variables.
parameters.continuous_variable_names = {'speed', 'accel', 'duration', 'pupil_diameter', 'tail', 'nose', 'FL', 'HL', 'x'};

dir_stems(1).name = 'PLSR fluorescence Inverted Intercepts'; 
dir_stems(1).comparisons = parameters.comparisons_continuous;
dir_stems(2).name = 'PLSR Inverted Intercepts'; 
dir_stems(2).comparisons = parameters.comparisons_continuous;


% change dir_stems here
parameters.dir_stems = dir_stems;
%parameters.dir_stems(1).comparisons = parameters.dir_stems(1).comparisons(5:end);
%parameters.dir_stems(1).comparisons.variablesToUse2 = parameters.dir_stems(1).comparisons.variablesToUse2(5:end);

% Put relevant variables into loop_variables.
parameters.loop_variables.mice_all = parameters.mice_all;
parameters.loop_variables.continuous_variables = parameters.continuous_variable_names;
parameters.loop_variables.comparisons_continuous = parameters.comparisons_continuous; 
parameters.loop_variables.singleActions = {'singleVariable', 'omittedVariable'};
parameters.loop_variables.dir_stems = parameters.dir_stems;
parameters.average_and_std_together = false;

%% fluorescence

if isfield(parameters, 'loop_list')
parameters = rmfield(parameters,'loop_list');
end

% Iterators
parameters.loop_list.iterators = {
               'mouse', {'loop_variables.mice_all(:).name'}, 'mouse_iterator'; 
               'comparison', {'loop_variables.comparisons_continuous(:).name'}, 'comparison_iterator'     
               };

parameters.comparisons_continuous = parameters.comparisons_continuous;

% Run single and omitted variables?
parameters.runSingleVariables = true;
parameters.omitSingleVariables = true;
parameters.shuffleOthers = false; 
parameters.shuffleNumber = 10;

% PLSR parameters
parameters.findBestNComponents = true;
parameters.ncomponents_max_default = 8; 
parameters.contiguous_partitions = true; 
parameters.kFolds = 10;
parameters.MonteCarloReps = 10;
parameters.comparison_type = 'continuous';
parameters.stratify = false;
parameters.permutationGeneration = false;

% Inputs 
% dataset for PLSR (from DatasetPrep)
parameters.loop_list.things_to_load.dataset.dir = {[parameters.dir_exper 'PLSR fluorescence Inverted Intercepts\variable prep\datasets\level 1 continuous\'], 'comparison', '\' 'mouse', '\'};
parameters.loop_list.things_to_load.dataset.filename= {'PLSR_dataset_info.mat'};
parameters.loop_list.things_to_load.dataset.variable= {'dataset_info'}; 
parameters.loop_list.things_to_load.dataset.level = 'comparison';

% Outputs
parameters.loop_list.things_to_save.all_variability.dir = {[parameters.dir_exper 'PLSR fluorescence Inverted Intercepts\results\variability by parameter shuffles\'], 'comparison', '\' 'mouse', '\'};
parameters.loop_list.things_to_save.all_variability.filename= {'all_variability.mat'};
parameters.loop_list.things_to_save.all_variability.variable= {'all_variability'}; 
parameters.loop_list.things_to_save.all_variability.level = 'comparison';

RunAnalysis({@VariabilityByVariable}, parameters);


%% correlations

if isfield(parameters, 'loop_list')
parameters = rmfield(parameters,'loop_list');
end

% Iterators
parameters.loop_list.iterators = {
               'mouse', {'loop_variables.mice_all(:).name'}, 'mouse_iterator'; 
               'comparison', {'loop_variables.comparisons_continuous(:).name'}, 'comparison_iterator'     
               };

parameters.comparisons_continuous = parameters.comparisons_continuous;

%parameters.loop_variables.dir_stems = {'PLSR Inverted Intercepts', 'PLSR fluorescence Inverted Intercepts'};
% Run single and omitted variables?
parameters.runSingleVariables = true;
parameters.omitSingleVariables = true;
parameters.shuffleOthers = false; 
parameters.shuffleNumber = 10;

% PLSR parameters
parameters.findBestNComponents = true;
parameters.ncomponents_max_default = 8; 
parameters.contiguous_partitions = true; 
parameters.kFolds = 10;
parameters.MonteCarloReps = 10;
parameters.comparison_type = 'continuous';
parameters.stratify = false;
parameters.permutationGeneration = false;

% Inputs 
% dataset for PLSR (from DatasetPrep)
parameters.loop_list.things_to_load.dataset.dir = {[parameters.dir_exper 'PLSR Inverted Intercepts\variable prep\datasets\level 1 continuous\'], 'comparison', '\' 'mouse', '\'};
parameters.loop_list.things_to_load.dataset.filename= {'PLSR_dataset_info.mat'};
parameters.loop_list.things_to_load.dataset.variable= {'dataset_info'}; 
parameters.loop_list.things_to_load.dataset.level = 'comparison';

% Outputs
parameters.loop_list.things_to_save.all_variability.dir = {[parameters.dir_exper 'PLSR Inverted Intercepts\results\variability by parameter shuffles\'], 'comparison', '\' 'mouse', '\'};
parameters.loop_list.things_to_save.all_variability.filename= {'all_variability.mat'};
parameters.loop_list.things_to_save.all_variability.variable= {'all_variability'}; 
parameters.loop_list.things_to_save.all_variability.level = 'comparison';

RunAnalysis({@VariabilityByVariable}, parameters);

%% Sums of variability 
if isfield(parameters, 'loop_list')
parameters = rmfield(parameters,'loop_list');
end

% Iterators
parameters.loop_list.iterators = {
               'dir_stem', {'loop_variables.dir_stems(:).name'}, 'dir_stem_iterator'; 
               'mouse', {'loop_variables.mice_all(:).name'}, 'mouse_iterator'; 
               'comparison', {'loop_variables.dir_stems(', 'dir_stem_iterator', ').comparisons(:).name'}, 'comparison_iterator';     
               'singleAction', {'loop_variables.singleActions'}, 'singleAction_iterator';
               'continuous_variable', {'loop_variables.dir_stems(', 'dir_stem_iterator', ').comparisons(', 'comparison_iterator', ').variablesToUse2'}, 'continuous_variable_iterator';
               };

% for brain data, take the sum of the second row
parameters.evaluation_instructions = {{'data = parameters.data;'...
                                       'data_evaluated = sum(data(2, :));'
                                       }};

% Inputs
parameters.loop_list.things_to_load.data.dir = {[parameters.dir_exper], 'dir_stem',  '\results\variability by parameter shuffles\', 'comparison', '\' 'mouse', '\'};
parameters.loop_list.things_to_load.data.filename= {'all_variability.mat'};
parameters.loop_list.things_to_load.data.variable= {'all_variability.', 'continuous_variable', '.', 'singleAction', '.PCTVAR'}; 
parameters.loop_list.things_to_load.data.level = 'comparison';

% Outputs
parameters.loop_list.things_to_save.data_evaluated.dir = {[parameters.dir_exper], 'dir_stem',  '\results\variability by parameter shuffles\', 'comparison', '\' 'mouse', '\'};
parameters.loop_list.things_to_save.data_evaluated.filename= {'singleAction', '_', 'continuous_variable', '.mat'};
parameters.loop_list.things_to_save.data_evaluated.variable= {'PCTVAR'}; 
parameters.loop_list.things_to_save.data_evaluated.level = 'continuous_variable';

RunAnalysis({@EvaluateOnData}, parameters);


%% Average across mice 
% skip mouse 1100
if isfield(parameters, 'loop_list')
parameters = rmfield(parameters,'loop_list');
end

% Iterators
parameters.loop_list.iterators = {
              'dir_stem', {'loop_variables.dir_stems(:).name'}, 'dir_stem_iterator'; 
               'comparison', {'loop_variables.dir_stems(', 'dir_stem_iterator', ').comparisons(:).name'}, 'comparison_iterator';     
               'singleAction', {'loop_variables.singleActions'}, 'singleAction_iterator';
               'continuous_variable', {'loop_variables.dir_stems(', 'dir_stem_iterator', ').comparisons(', 'comparison_iterator', ').variablesToUse2'}, 'continuous_variable_iterator';
               'mouse', {'loop_variables.mice_all(:).name'}, 'mouse_iterator'; 
               };

parameters.concatDim = 1;
parameters.averageDim = 1;
                        

% Inputs
parameters.loop_list.things_to_load.data.dir = {[parameters.dir_exper], 'dir_stem',  '\results\variability by parameter shuffles\', 'comparison', '\' 'mouse', '\'};
parameters.loop_list.things_to_load.data.filename= {'singleAction', '_', 'continuous_variable', '.mat'};
parameters.loop_list.things_to_load.data.variable= {'PCTVAR'}; 
parameters.loop_list.things_to_load.data.level = 'mouse';

% Outputs
% all concatenated 
parameters.loop_list.things_to_save.concatenated_data.dir = {[parameters.dir_exper], 'dir_stem',  '\results\variability by parameter shuffles\', 'comparison', '\mean across mice\'};
parameters.loop_list.things_to_save.concatenated_data.filename= {'singleAction', '_', 'continuous_variable', '_all.mat'};
parameters.loop_list.things_to_save.concatenated_data.variable= {'all_PCTVAR'}; 
parameters.loop_list.things_to_save.concatenated_data.level = 'continuous_variable';
% average
parameters.loop_list.things_to_save.average.dir = {[parameters.dir_exper], 'dir_stem',  '\results\variability by parameter shuffles\', 'comparison', '\mean across mice\'};
parameters.loop_list.things_to_save.average.filename= {'singleAction', '_', 'continuous_variable', '_average.mat'};
parameters.loop_list.things_to_save.average.variable= {'average_PCTVAR'}; 
parameters.loop_list.things_to_save.average.level = 'continuous_variable';
% std_dev
parameters.loop_list.things_to_save.std_dev.dir = {[parameters.dir_exper], 'dir_stem',  '\results\variability by parameter shuffles\', 'comparison', '\mean across mice\'};
parameters.loop_list.things_to_save.std_dev.filename= {'singleAction', '_', 'continuous_variable', '_std.mat'};
parameters.loop_list.things_to_save.std_dev.variable= {'std_PCTVAR'}; 
parameters.loop_list.things_to_save.std_dev.level = 'continuous_variable';

parameters.loop_list.things_to_rename = {{'data_out', 'data'}
                                         {'concatenated_data', 'data'}};

RunAnalysis({@MouseNotToUseVariability, @ConcatenateData, @AverageData}, parameters);


%% For omitted variables, get difference from total variance explained
% within mice

if isfield(parameters, 'loop_list')
parameters = rmfield(parameters,'loop_list');
end

% Iterators
parameters.loop_list.iterators = {
               'dir_stem', {'loop_variables.dir_stems(:).name'}, 'dir_stem_iterator'; 
               'mouse', {'loop_variables.mice_all(:).name'}, 'mouse_iterator'; 
               'comparison', {'loop_variables.dir_stems(', 'dir_stem_iterator', ').comparisons(:).name'}, 'comparison_iterator';     
               'continuous_variable', {'loop_variables.dir_stems(', 'dir_stem_iterator', ').comparisons(', 'comparison_iterator', ').variablesToUse2'}, 'continuous_variable_iterator';
               };

% Sum total variability of all components 
parameters.evaluation_instructions = {{'total = parameters.total;'...
                                       'data_evaluated = sum(total(2,:));'
                                        }
% Subtract individual variability from total
                                       {'data = parameters.data;'...
                                        'total = parameters.total;'...
                                         'data_evaluated = total - data;'}
                                      }; 


% Input
% summed percents for variable 
parameters.loop_list.things_to_load.data.dir = {[parameters.dir_exper], 'dir_stem',  '\results\variability by parameter shuffles\', 'comparison', '\' 'mouse', '\'};
parameters.loop_list.things_to_load.data.filename= {'omittedVariable_', 'continuous_variable', '.mat'};
parameters.loop_list.things_to_load.data.variable= {'PCTVAR'}; 
parameters.loop_list.things_to_load.data.level = 'continuous_variable';
% Total variability
parameters.loop_list.things_to_load.total.dir = {[parameters.dir_exper], 'dir_stem',  '\results\level 1 continuous\', 'comparison', '\' 'mouse', '\'};
parameters.loop_list.things_to_load.total.filename= {'PLSR_results.mat'};
parameters.loop_list.things_to_load.total.variable= {'PLSR_results.PCTVAR'}; 
parameters.loop_list.things_to_load.total.level = 'comparison';

% Output
parameters.loop_list.things_to_save.data_evaluated.dir = {[parameters.dir_exper], 'dir_stem',  '\results\variability by parameter shuffles\', 'comparison', '\' 'mouse', '\'};
parameters.loop_list.things_to_save.data_evaluated.filename= {'omittedVariable_difference_', 'continuous_variable', '.mat'};
parameters.loop_list.things_to_save.data_evaluated.variable= {'PCTVAR'}; 
parameters.loop_list.things_to_save.data_evaluated.level = 'continuous_variable';

parameters.loop_list.things_to_rename = {{'data_evaluated', 'total'};
                                          };

RunAnalysis({@EvaluateOnData, @EvaluateOnData}, parameters);

%% Average omitted variable difference across mice 

if isfield(parameters, 'loop_list')
parameters = rmfield(parameters,'loop_list');
end

% Iterators
parameters.loop_list.iterators = {
              'dir_stem', {'loop_variables.dir_stems(:).name'}, 'dir_stem_iterator'; 
               'comparison', {'loop_variables.dir_stems(', 'dir_stem_iterator', ').comparisons(:).name'}, 'comparison_iterator';     
               'continuous_variable', {'loop_variables.dir_stems(', 'dir_stem_iterator', ').comparisons(', 'comparison_iterator', ').variablesToUse2'}, 'continuous_variable_iterator';
               'mouse', {'loop_variables.mice_all(:).name'}, 'mouse_iterator'; 
               };

parameters.concatDim = 1;
parameters.averageDim = 1;
                        

% Inputs
parameters.loop_list.things_to_load.data.dir = {[parameters.dir_exper], 'dir_stem',  '\results\variability by parameter shuffles\', 'comparison', '\' 'mouse', '\'};
parameters.loop_list.things_to_load.data.filename= {'omittedVariable_difference_', 'continuous_variable', '.mat'};
parameters.loop_list.things_to_load.data.variable= {'PCTVAR'}; 
parameters.loop_list.things_to_load.data.level = 'mouse';

% Outputs
% all concatenated 
parameters.loop_list.things_to_save.concatenated_data.dir = {[parameters.dir_exper], 'dir_stem',  '\results\variability by parameter shuffles\', 'comparison', '\mean across mice\'};
parameters.loop_list.things_to_save.concatenated_data.filename= {'omittedVariable_difference_', 'continuous_variable', '_all.mat'};
parameters.loop_list.things_to_save.concatenated_data.variable= {'allPCTVAR'}; 
parameters.loop_list.things_to_save.concatenated_data.level = 'continuous_variable';
% average
parameters.loop_list.things_to_save.average.dir = {[parameters.dir_exper], 'dir_stem',  '\results\variability by parameter shuffles\', 'comparison', '\mean across mice\'};
parameters.loop_list.things_to_save.average.filename= {'omittedVariable_difference_', 'continuous_variable', '_average.mat'};
parameters.loop_list.things_to_save.average.variable= {'averagePCTVAR'}; 
parameters.loop_list.things_to_save.average.level = 'continuous_variable';
% std_dev
parameters.loop_list.things_to_save.std_dev.dir = {[parameters.dir_exper], 'dir_stem',  '\results\variability by parameter shuffles\', 'comparison', '\mean across mice\'};
parameters.loop_list.things_to_save.std_dev.filename= {'omittedVariable_difference_', 'continuous_variable', '_std.mat'};
parameters.loop_list.things_to_save.std_dev.variable= {'stdPCTVAR'}; 
parameters.loop_list.things_to_save.std_dev.level = 'continuous_variable';

parameters.loop_list.things_to_rename = {{'data_out', 'data'}
                                         {'concatenated_data', 'data'}};

RunAnalysis({@MouseNotToUseVariability, @ConcatenateData, @AverageData}, parameters);
