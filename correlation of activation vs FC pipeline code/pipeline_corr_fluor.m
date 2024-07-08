% pipeline_corr_fluor.m
% Sarah West
% 10/24/23
% Runs correlations between results of PLSR on fluorescece and PLSR on correlaions
% per node. Is to validate if average node correlation is dependent on node
% fluorescence.

% Updated 10/24/23 for the "specials" and "other pairs" comparisons

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

% normal vs Warning periods, comparison types
parameters.loop_variables.categories = {'normal', 'warningPeriods'};
parameters.loop_variables.comparison_types = {'categorical', 'continuous'};

% Load comparisons

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

% Load comparisons for first level continuous, specials
load([parameters.dir_exper 'PLSR\comparisons_special_continuous.mat']);
parameters.comparisons_continuous_Specials = comparisons;
clear comparisons;

% other pairs
load([parameters.dir_exper 'PLSR Specials Other Pairs\comparisons_special_otherPairs_continuous.mat']);
parameters.comparisons_continuous_OtherPairs = comparisons;
clear comparisons

% Load comparisons for first level categorical, special
load([parameters.dir_exper 'PLSR\comparisons_special_categorical.mat']);
parameters.comparisons_categorical_Specials = comparisons;
clear comparisons;

% other pairs
load([parameters.dir_exper 'PLSR Specials Other Pairs\comparisons_special_otherPairs_categorical.mat']);
parameters.comparisons_categorical_OtherPairs = comparisons;
clear comparisons;


dir_stems(1).name = 'PLSR fluorescence Specials Inverted Intercepts'; 
dir_stems(1).comparisons = parameters.comparisons_categorical_Specials([1:27 30:end]);
dir_stems(2).name = 'PLSR Specials Inverted Intercepts'; 
dir_stems(2).comparisons = parameters.comparisons_categorical_Specials([1:27 30:end]);
% dir_stems(3).name = 'PLSR fluorescence Specials Other Pairs Inverted Intercepts'; 
% dir_stems(3).comparisons = parameters.comparisons_categorical_OtherPairs(1:8);
% dir_stems(4).name = 'PLSR Specials Other Pairs Inverted Intercepts'; 
% dir_stems(4).comparisons = parameters.comparisons_categorical_OtherPairs(1:8);


dir_stem_types(1).name = 'correlations'; 
dir_stem_types(1).sub_names = {'PLSR Specials Inverted Intercepts', 'PLSR Specials Other Pairs Inverted Intercepts'};
dir_stem_types(2).name = 'fluorescence';
dir_stem_types(2).sub_names = {'PLSR fluorescence Specials Inverted Intercepts', 'PLSR fluorescence Specials Other Pairs Inverted Intercepts'};
% Put relevant variables into loop_variables.
parameters.loop_variables.dir_stems = dir_stems;
parameters.loop_variables.dir_stem_types = dir_stem_types;
parameters.loop_variables.mice_all = parameters.mice_all;
parameters.loop_variables.data_type = {'correlations'; 'fluorescence'};
parameters.average_and_std_together = false;



% Renumbering. Go down to only 16 nodes (one hemisphere only). 
load('Y:\Sarah\Analysis\Experiments\Random Motorized Treadmill\node_renumbering.mat');
parameters.node_renumbering = node_renumbering(2:2:end, :)/2;

%% *** Start with categorical only ***

%% reshape level 1 correlation results to not have double representation (16 instead of 32 values, matches fluorescence)
% also transpose so dimensions match fluorescence

if isfield(parameters, 'loop_list')
parameters = rmfield(parameters,'loop_list');
end

% get only the correlation dir_stems
parameters.dir_stems = dir_stems([2]); % 4]); 
parameters.loop_variables.dir_stems = parameters.dir_stems;

% Iterators
parameters.loop_list.iterators = {
               'dir_stem', {'loop_variables.dir_stems(:).name'}, 'dir_stem_iterator'; 
               'mouse', {'loop_variables.mice_all(:).name'}, 'mouse_iterator'; 
               'comparison', {'loop_variables.dir_stems(', 'dir_stem_iterator', ').comparisons(:).name'}, 'comparison_iterator';   
               };

parameters.evaluation_instructions = {{ 'data = parameters.data;'...
                                       'data_evaluated = transpose(data(:, 1:2:end));'}};

% Inputs 
parameters.loop_list.things_to_load.data.dir = {[parameters.dir_exper], 'dir_stem', '\results\level 1 categorical\', 'comparison', '\', 'mouse', '\'};
parameters.loop_list.things_to_load.data.filename= {'average_by_nodes_BETA.mat'};
parameters.loop_list.things_to_load.data.variable= {'average_by_nodes'}; 
parameters.loop_list.things_to_load.data.level = 'comparison';

% Outputs
parameters.loop_list.things_to_save.data_evaluated.dir = {[parameters.dir_exper 'corr fluor\data reshaped\'], 'dir_stem', '\', 'comparison', '\', 'mouse', '\'};
parameters.loop_list.things_to_save.data_evaluated.filename= {'average_by_nodes.mat'};
parameters.loop_list.things_to_save.data_evaluated.variable= {'average_by_nodes'}; 
parameters.loop_list.things_to_save.data_evaluated.level = 'comparison';

RunAnalysis({@EvaluateOnData}, parameters);


%% keep only odd nodes for fluorescence
% multiply by sigmas to account for variability per comparison

if isfield(parameters, 'loop_list')
parameters = rmfield(parameters,'loop_list');
end

% get only the fluorescence dir_stems
parameters.dir_stems = dir_stems([1]); %  3]); 
parameters.loop_variables.dir_stems = parameters.dir_stems;

% Iterators
parameters.loop_list.iterators = {
               'dir_stem', {'loop_variables.dir_stems(:).name'}, 'dir_stem_iterator'; 
               'mouse', {'loop_variables.mice_all(:).name'}, 'mouse_iterator'; 
               'comparison', {'loop_variables.dir_stems(', 'dir_stem_iterator', ').comparisons(:).name'}, 'comparison_iterator';   
               };

parameters.evaluation_instructions = {{ 'data = parameters.data;'...
                                        'sigmas = parameters.sigmas;'...
                                       'holder = data .* sigmas;'...
                                       'data_evaluated = transpose(holder(2, 1:2:end));'}}; % remove intercepts

% Inputs 
% data
parameters.loop_list.things_to_load.data.dir = {[parameters.dir_exper], 'dir_stem', '\results\level 1 categorical\', 'comparison', '\', 'mouse', '\'};
parameters.loop_list.things_to_load.data.filename= {'PLSR_results.mat'};
parameters.loop_list.things_to_load.data.variable= {'PLSR_results.BETA'}; 
parameters.loop_list.things_to_load.data.level = 'comparison';
% sigmas
parameters.loop_list.things_to_load.sigmas.dir = {[parameters.dir_exper], 'dir_stem', '\variable prep\datasets\level 1 categorical\', 'comparison', '\',  'mouse', '\'};
parameters.loop_list.things_to_load.sigmas.filename= {'sigmas.mat'};
parameters.loop_list.things_to_load.sigmas.variable= {'sigmas'}; 
parameters.loop_list.things_to_load.sigmas.level = 'comparison';

% Outputs
parameters.loop_list.things_to_save.data_evaluated.dir = {[parameters.dir_exper 'corr fluor\data reshaped\'], 'dir_stem', '\', 'comparison', '\', 'mouse', '\'};
parameters.loop_list.things_to_save.data_evaluated.filename= {'average_by_nodes.mat'}; % make this the file name to keep it consistent with above
parameters.loop_list.things_to_save.data_evaluated.variable= {'average_by_nodes'}; 
parameters.loop_list.things_to_save.data_evaluated.level = 'comparison';

RunAnalysis({@EvaluateOnData}, parameters);

%% Concatenate across comparisons within mice 
% (categorical only)

if isfield(parameters, 'loop_list')
parameters = rmfield(parameters,'loop_list');
end

% get all dir_stems
parameters.dir_stems = dir_stems; 
parameters.loop_variables.dir_stems = parameters.dir_stems;

% Iterators
parameters.loop_list.iterators = {
           'dir_stem', {'loop_variables.dir_stems(:).name'}, 'dir_stem_iterator'; 
           'mouse', {'loop_variables.mice_all(:).name'}, 'mouse_iterator'; 
           'comparison', {'loop_variables.dir_stems(', 'dir_stem_iterator', ').comparisons(:).name'}, 'comparison_iterator';   
   };

parameters.concatenation_level = 'comparison'; 
parameters.concatDim = 2;

% Inputs
parameters.loop_list.things_to_load.data.dir = {[parameters.dir_exper 'corr fluor\data reshaped\'], 'dir_stem', '\', 'comparison', '\', 'mouse', '\'};
parameters.loop_list.things_to_load.data.filename= {'average_by_nodes.mat'};
parameters.loop_list.things_to_load.data.variable= {'average_by_nodes'}; 
parameters.loop_list.things_to_load.data.level = 'comparison';

% Output
parameters.loop_list.things_to_save.concatenated_data.dir = {[parameters.dir_exper 'corr fluor\concatenated within mice across comparisons\'], 'dir_stem', '\categorical\', 'mouse', '\'}; 
parameters.loop_list.things_to_save.concatenated_data.filename= {'data_concatenated.mat'};
parameters.loop_list.things_to_save.concatenated_data.variable= {'data_concatenated'}; 
parameters.loop_list.things_to_save.concatenated_data.level = 'mouse';

RunAnalysis({@ConcatenateData}, parameters);
   

%% Concatenate across dir_stem within mice
% (categorical only) 

if isfield(parameters, 'loop_list')
parameters = rmfield(parameters,'loop_list');
end

% Iterators
parameters.loop_list.iterators = {
           'mouse', {'loop_variables.mice_all(:).name'}, 'mouse_iterator'; 
           'dir_stem_type', {'loop_variables.dir_stem_types(:).name'}, 'dir_stem_type_iterator';       
           'dir_stem', {'loop_variables.dir_stem_types(', 'dir_stem_type_iterator', ').sub_names{:}'}, 'dir_stem_iterator';   
   };

parameters.concatenation_level = 'dir_stem'; 
parameters.concatDim = 2;

% Inputs
parameters.loop_list.things_to_load.data.dir = {[parameters.dir_exper 'corr fluor\concatenated within mice across comparisons\'], 'dir_stem', '\categorical\', 'mouse', '\'}; 
parameters.loop_list.things_to_load.data.filename= {'data_concatenated.mat'};
parameters.loop_list.things_to_load.data.variable= {'data_concatenated'}; 
parameters.loop_list.things_to_load.data.level = 'dir_stem';

% Outputs
parameters.loop_list.things_to_save.concatenated_data.dir = {[parameters.dir_exper 'corr fluor\concatenated within mice across dir stems\categorical\'], 'dir_stem_type', '\', 'mouse', '\'}; 
parameters.loop_list.things_to_save.concatenated_data.filename= {'data_concatenated.mat'};
parameters.loop_list.things_to_save.concatenated_data.variable= {'data_concatenated'}; 
parameters.loop_list.things_to_save.concatenated_data.level = 'dir_stem_type';

RunAnalysis({@ConcatenateData}, parameters)

%% Concatenate across mice 
% (categorical only) 
if isfield(parameters, 'loop_list')
parameters = rmfield(parameters,'loop_list');
end

% Iterators
parameters.loop_list.iterators = {
           'dir_stem_type', {'loop_variables.dir_stem_types(:).name'}, 'dir_stem_type_iterator';  
           'mouse', {'loop_variables.mice_all(:).name'}, 'mouse_iterator'; 
   };
parameters.concatenation_level = 'mouse';
parameters.concatDim = 3; 

% Inputs 
parameters.loop_list.things_to_load.data.dir = {[parameters.dir_exper 'corr fluor\concatenated within mice across dir stems\categorical\'], 'dir_stem_type', '\', 'mouse', '\'}; 
parameters.loop_list.things_to_load.data.filename= {'data_concatenated.mat'};
parameters.loop_list.things_to_load.data.variable= {'data_concatenated'}; 
parameters.loop_list.things_to_load.data.level = 'mouse';

% Outputs
parameters.loop_list.things_to_save.concatenated_data.dir = {[parameters.dir_exper 'corr fluor\concatenated across mice\categorical\'], 'dir_stem_type', '\'}; 
parameters.loop_list.things_to_save.concatenated_data.filename= {'data_concatenated.mat'};
parameters.loop_list.things_to_save.concatenated_data.variable= {'data_concatenated'}; 
parameters.loop_list.things_to_save.concatenated_data.level = 'dir_stem_type';

RunAnalysis({@ConcatenateData}, parameters)

%% Permute correlations and fluorescence
% so they match expected format for CorrFluor.m
if isfield(parameters, 'loop_list')
parameters = rmfield(parameters,'loop_list');
end

% Iterators
parameters.loop_list.iterators = {
           'dir_stem_type', {'loop_variables.dir_stem_types(:).name'}, 'dir_stem_type_iterator';  
   };

parameters.evaluation_instructions = {{'data = parameters.data;'...
                                        'data_evaluated = permute(data, [3 1 2]);'}};
% Inputs
parameters.loop_list.things_to_load.data.dir = {[parameters.dir_exper 'corr fluor\concatenated across mice\categorical\'], 'dir_stem_type', '\'}; 
parameters.loop_list.things_to_load.data.filename= {'data_concatenated.mat'};
parameters.loop_list.things_to_load.data.variable= {'data_concatenated'}; 
parameters.loop_list.things_to_load.data.level = 'dir_stem_type';

% Outputs
parameters.loop_list.things_to_save.data_evaluated.dir = {[parameters.dir_exper 'corr fluor\concatenated across mice\categorical\'], 'dir_stem_type', '\'}; 
parameters.loop_list.things_to_save.data_evaluated.filename= {'data_concatenated_permuted.mat'};
parameters.loop_list.things_to_save.data_evaluated.variable= {'data_concatenated_permuted'}; 
parameters.loop_list.things_to_save.data_evaluated.level = 'dir_stem_type';

RunAnalysis({@EvaluateOnData}, parameters);

%% Run correlations per mouse
% (categorical
% Correlate within mice
% Average across mice

if isfield(parameters, 'loop_list')
parameters = rmfield(parameters,'loop_list');
end

% Iterators
parameters.loop_list.iterators = 'none';

parameters.mouseDim = 1;
parameters.nodeDim = 2;
parameters.comparison_type = 'categorical';

% Inputs 
% Correlations
parameters.loop_list.things_to_load.correlations.dir = {[parameters.dir_exper 'corr fluor\concatenated across mice\categorical\correlations\']}; 
parameters.loop_list.things_to_load.correlations.filename= {'data_concatenated_permuted.mat'};
parameters.loop_list.things_to_load.correlations.variable= {'data_concatenated_permuted'}; 
parameters.loop_list.things_to_load.correlations.level = 'start';

% Fluorescence
parameters.loop_list.things_to_load.fluorescence.dir = {[parameters.dir_exper 'corr fluor\concatenated across mice\categorical\fluorescence\']}; 
parameters.loop_list.things_to_load.fluorescence.filename= {'data_concatenated_permuted.mat'};
parameters.loop_list.things_to_load.fluorescence.variable= {'data_concatenated_permuted'}; 
parameters.loop_list.things_to_load.fluorescence.level = 'start';

% Outputs
% correlations per mouse p values 
parameters.loop_list.things_to_save.corrs_per_mouse.dir = {[parameters.dir_exper 'corr fluor\results\categorical\']}; 
parameters.loop_list.things_to_save.corrs_per_mouse.filename= {'corrs_per_mouse.mat'};
parameters.loop_list.things_to_save.corrs_per_mouse.variable= {'corrs_per_mouse'}; 
parameters.loop_list.things_to_save.corrs_per_mouse.level = 'end';

parameters.loop_list.things_to_save.corrs_per_mouse_pvalues.dir = {[parameters.dir_exper 'corr fluor\results\categorical\']}; 
parameters.loop_list.things_to_save.corrs_per_mouse_pvalues.filename= {'corrs_per_mouse_pvalues.mat'};
parameters.loop_list.things_to_save.corrs_per_mouse_pvalues.variable= {'corrs_per_mouse_pvalues'}; 
parameters.loop_list.things_to_save.corrs_per_mouse_pvalues.level = 'end';

% average correlations across mice 
% parameters.loop_list.things_to_save.corrs_across_mice.dir = {[parameters.dir_exper 'corr fluor\results\categorical\']}; 
% parameters.loop_list.things_to_save.corrs_across_mice.filename= {'corrs_across_mice.mat'};
% parameters.loop_list.things_to_save.corrs_across_mice.variable= {'corrs_across_mice'}; 
% parameters.loop_list.things_to_save.corrs_across_mice.level = 'end';

RunAnalysis({@CorrFluor}, parameters); 

%% create null distributions
% shuffle correlations relative to fluorescence 

if isfield(parameters, 'loop_list')
parameters = rmfield(parameters,'loop_list');
end

% Iterators
parameters.loop_list.iterators = 'none';

parameters.mouseDim = 1;
parameters.nodeDim = 2;
parameters.comparison_type = 'categorical';

parameters.shuffleDim = 4; 
parameters.shuffle_along_dimension = 3;
parameters.shuffleNumber = 5000;
parameters.saveShuffled = true;
parameters.mix_distributions = false;
parameters.reorder_one_distribution = true;

% Inputs 
% data1: Fluorescence
parameters.loop_list.things_to_load.data1.dir = {[parameters.dir_exper 'corr fluor\concatenated across mice\categorical\fluorescence\']}; 
parameters.loop_list.things_to_load.data1.filename= {'data_concatenated_permuted.mat'};
parameters.loop_list.things_to_load.data1.variable= {'data_concatenated_permuted'}; 
parameters.loop_list.things_to_load.data1.level = 'start';
% data2: Correlations
parameters.loop_list.things_to_load.data2.dir = {[parameters.dir_exper 'corr fluor\concatenated across mice\categorical\correlations\']}; 
parameters.loop_list.things_to_load.data2.filename= {'data_concatenated_permuted.mat'};
parameters.loop_list.things_to_load.data2.variable= {'data_concatenated_permuted'}; 
parameters.loop_list.things_to_load.data2.level = 'start';

% Outputs
parameters.loop_list.things_to_save.null_distributions.dir = {[parameters.dir_exper 'corr fluor\concatenated across mice\categorical\null distributions\']}; 
parameters.loop_list.things_to_save.null_distributions.filename= {'null_distributions_new_values.mat'};
parameters.loop_list.things_to_save.null_distributions.variable= {'new_values'}; 
parameters.loop_list.things_to_save.null_distributions.level = 'end';

RunAnalysis({@NullDistribution}, parameters);

%% Correlate the new null distributions 
% Correlate within mice

if isfield(parameters, 'loop_list')
parameters = rmfield(parameters,'loop_list');
end

% Iterators
parameters.loop_list.iterators = {
                                'shuffle', {'loop_variables.shuffles'}, 'shuffle_iterator';
                                 }; 

parameters.loop_variables.shuffles = 1:5000; 
parameters.mouseDim = 1;
parameters.nodeDim = 2;
parameters.comparison_type = 'categorical';

% Inputs 
% Correlations
parameters.loop_list.things_to_load.correlations.dir = {[parameters.dir_exper 'corr fluor\concatenated across mice\categorical\null distributions\']}; 
parameters.loop_list.things_to_load.correlations.filename= {'null_distributions_new_values.mat'};
parameters.loop_list.things_to_load.correlations.variable= {'new_values{2}(:,:,:,', 'shuffle_iterator', ')'}; 
parameters.loop_list.things_to_load.correlations.level = 'start';

% Fluorescence
parameters.loop_list.things_to_load.fluorescence.dir = {[parameters.dir_exper 'corr fluor\concatenated across mice\categorical\null distributions\']}; 
parameters.loop_list.things_to_load.fluorescence.filename= {'null_distributions_new_values.mat'};
parameters.loop_list.things_to_load.fluorescence.variable= {'new_values{1}(:,:,:,', 'shuffle_iterator', ')'}; 
parameters.loop_list.things_to_load.fluorescence.level = 'start';

% Outputs
% correlations per mouse p values 
parameters.loop_list.things_to_save.corrs_per_mouse.dir = {[parameters.dir_exper 'corr fluor\results\categorical\']}; 
parameters.loop_list.things_to_save.corrs_per_mouse.filename= {'null_distributions_corrs_per_mouse.mat'};
parameters.loop_list.things_to_save.corrs_per_mouse.variable= {'null_distributions_corrs_per_mouse{', 'shuffle_iterator', ',1}'}; 
parameters.loop_list.things_to_save.corrs_per_mouse.level = 'end';

% average correlations across mice 
% parameters.loop_list.things_to_save.corrs_across_mice.dir = {[parameters.dir_exper 'corr fluor\results\categorical\']}; 
% parameters.loop_list.things_to_save.corrs_across_mice.filename= {'corrs_across_mice.mat'};
% parameters.loop_list.things_to_save.corrs_across_mice.variable= {'corrs_across_mice'}; 
% parameters.loop_list.things_to_save.corrs_across_mice.level = 'end';

RunAnalysis({@CorrFluor}, parameters); 


%% Take average of correlations without mouse 1100: test values 

% load corrs_per_mouse
load('Y:\Sarah\Analysis\Experiments\Random Motorized Treadmill\corr fluor\results\categorical\corrs_per_mouse.mat');

% find average & std across mice
corrs_average = mean(corrs_per_mouse([1:3 5:7], :), 1, 'omitnan');
corrs_std_dev = std(corrs_per_mouse([1:3 5:7], :), [], 1, 'omitnan');

% save average & std
save('Y:\Sarah\Analysis\Experiments\Random Motorized Treadmill\corr fluor\results\categorical\corrs_average.mat', 'corrs_average');
save('Y:\Sarah\Analysis\Experiments\Random Motorized Treadmill\corr fluor\results\categorical\corrs_std_dev.mat', 'corrs_std_dev');

clear corrs_per_mouse corrs_average corrs_std_dev;

%% Take average of correlations without mouse 1100: null distributions

% load null_distributions_corrs_per_mouse
load('Y:\Sarah\Analysis\Experiments\Random Motorized Treadmill\corr fluor\results\categorical\null_distributions_corrs_per_mouse.mat');

% put all null distributions into 1 matrix.
null_distributions = cat(3, null_distributions_corrs_per_mouse{:});

% find average of null distributions across mice
null_distributions_average = mean(null_distributions([1:3 5:7], :, :), 1, 'omitnan');

% save average & std
save('Y:\Sarah\Analysis\Experiments\Random Motorized Treadmill\corr fluor\results\categorical\null_distribution_average.mat', 'null_distributions_average');

clear null_distributions null_distributions_average null_distributions_corrs_per_mouse;

%% Calculate significance with null distribution

if isfield(parameters, 'loop_list')
parameters = rmfield(parameters,'loop_list');
end

% Iterators
parameters.loop_list.iterators = 'none';

parameters.useNormalDistribution = true;
parameters.alphaValue = 0.05;
parameters.useFDR = true;
parameters.shufflesDim = 3;

% Inputs
% null distribution
parameters.loop_list.things_to_load.null_distribution.dir = {[parameters.dir_exper 'corr fluor\results\categorical\']}; 
parameters.loop_list.things_to_load.null_distribution.filename= {'null_distribution_average.mat'};
parameters.loop_list.things_to_load.null_distribution.variable= {'null_distributions_average'}; 
parameters.loop_list.things_to_load.null_distribution.level = 'start';
% test values
parameters.loop_list.things_to_load.test_values.dir = {[parameters.dir_exper 'corr fluor\results\categorical\']}; 
parameters.loop_list.things_to_load.test_values.filename= {'corrs_average.mat'};
parameters.loop_list.things_to_load.test_values.variable= {'corrs_average'}; 
parameters.loop_list.things_to_load.test_values.level = 'start';

% Output
parameters.loop_list.things_to_save.significance.dir = {[parameters.dir_exper 'corr fluor\results\categorical\']}; 
parameters.loop_list.things_to_save.significance.filename= {'significance.mat'};
parameters.loop_list.things_to_save.significance.variable= {'significance'}; 
parameters.loop_list.things_to_save.significance.level = 'end';

RunAnalysis({@SignificanceCalculation}, parameters);

%% Renumber

% test values --average
load('Y:\Sarah\Analysis\Experiments\Random Motorized Treadmill\corr fluor\results\categorical\corrs_average.mat');
corrs_average_renumbered= ArrangeNewNumbering(corrs_average, parameters.node_renumbering, false, 2, [], []);
save('Y:\Sarah\Analysis\Experiments\Random Motorized Treadmill\corr fluor\results\categorical\corrs_averaged_renumbered.mat', 'corrs_average_renumbered');

% test values --std
load('Y:\Sarah\Analysis\Experiments\Random Motorized Treadmill\corr fluor\results\categorical\corrs_std_dev.mat');
corrs_std_dev_renumbered= ArrangeNewNumbering(corrs_std_dev, parameters.node_renumbering, false, 2, [], []);
save('Y:\Sarah\Analysis\Experiments\Random Motorized Treadmill\corr fluor\results\categorical\corrs_std_dev_renumbered.mat', 'corrs_std_dev_renumbered');

% significance (for completeness)
load('Y:\Sarah\Analysis\Experiments\Random Motorized Treadmill\corr fluor\results\categorical\significance.mat');
significance_renumbered = ArrangeNewNumbering(significance.all, parameters.node_renumbering, false, 2, [], []);
save('Y:\Sarah\Analysis\Experiments\Random Motorized Treadmill\corr fluor\results\categorical\significance_renumbered.mat', 'significance_renumbered');


%% plot 

load('Y:\Sarah\Analysis\Experiments\Random Motorized Treadmill\corr fluor\results\categorical\corrs_averaged_renumbered.mat')
load('Y:\Sarah\Analysis\Experiments\Random Motorized Treadmill\corr fluor\results\categorical\corrs_std_dev_renumbered.mat')
fig = figure; bar(corrs_average_renumbered);
hold on;
% calculate SEM
SEM = corrs_std_dev_renumbered./sqrt(6); 
errorbar(corrs_average_renumbered, SEM, 'LineStyle', 'none');

% x ticks 
xticks(1:16)
xticklabels([1:2:32]);

title('correlation between fluorescence and FC, with SEM')
ylabel('correlation');
xlabel('node')

saveas(fig, 'Y:\Sarah\Analysis\Experiments\Random Motorized Treadmill\corr fluor\results\categorical\plot', 'fig' );
saveas(fig, 'Y:\Sarah\Analysis\Experiments\Random Motorized Treadmill\corr fluor\results\categorical\plot', 'svg' );

%% *** Repeat for continuous ***

%% reshape level 1 correlation results to not have double representation (16 instead of 32 values, matches fluorescence)
% For each category (couldn't put in interators because of folder names)
for categoryi = 1:numel(parameters.loop_variables.categories)

    category = parameters.loop_variables.categories{categoryi};

    if isfield(parameters, 'loop_list')
    parameters = rmfield(parameters,'loop_list');
    end
    
    % Iterators
    parameters.loop_list.iterators = {
                   'comparison', {'loop_variables.comparisons_continuous.' category '(:).name'}, 'comparison_iterator' ;    
                   };

    parameters.evaluation_instructions = {'data_evaluated = permute(parameters.data(1:2:end, :, :), [2 1 3]);'};

    % Inputs 
    if strcmp(category, 'normal')
        parameters.loop_list.things_to_load.data.dir = {[parameters.dir_exper 'PLSR\results\level 2 continuous\Ipsa Contra\'], 'comparison', '\'};
    else
        parameters.loop_list.things_to_load.data.dir = {[parameters.dir_exper 'PLSR Warning Periods\results\level 2 continuous\'], 'comparison', '\'};
    end
    parameters.loop_list.things_to_load.data.filename= {'average_by_nodes_Cov.mat'};
    parameters.loop_list.things_to_load.data.variable= {'average_by_nodes'}; 
    parameters.loop_list.things_to_load.data.level = 'comparison';

    % Outputs
    parameters.loop_list.things_to_save.data_evaluated.dir = {[parameters.dir_exper 'corr fluor\only relevant corrs\continuous\' category '\'], 'comparison', '\'};
    parameters.loop_list.things_to_save.data_evaluated.filename= {'average_by_nodes.mat'};
    parameters.loop_list.things_to_save.data_evaluated.variable= {'average_by_nodes'}; 
    parameters.loop_list.things_to_save.data_evaluated.level = 'comparison';

    RunAnalysis({@EvaluateOnData}, parameters);
end 

%% Reshape Continuous data
% For each category (couldn't put in interators because of folder names)
for categoryi = 1:numel(parameters.loop_variables.categories)

    category = parameters.loop_variables.categories{categoryi};

    parameters.this_comparison_set = parameters.comparisons_continuous.(category);

    for data_typei = 1:numel(parameters.loop_variables.data_type)
        data_type = parameters.loop_variables.data_type{data_typei}; 

        if isfield(parameters, 'loop_list')
        parameters = rmfield(parameters,'loop_list');
        end
        
        % Iterators
        parameters.loop_list.iterators = {
                       'comparison', {'loop_variables.comparisons_continuous.', category, '(:).name'}, 'comparison_iterator' ;    
                       };
        
        parameters.comparison_type = 'continuous'; 
        parameters.variablesDimIn = 2;
    
        % Inputs 
        if strcmp(data_type, 'fluors')
            % fluorescence normal
            if strcmp(category, 'normal') 
                parameters.loop_list.things_to_load.data.dir = {[parameters.dir_exper 'PLSR fluorescence\variable prep\datasets\level 2 continuous\'], 'comparison', '\'};
            % fluorescence warning periods
            else 
                parameters.loop_list.things_to_load.data.dir = {[parameters.dir_exper 'PLSR fluorescence Warning Periods\variable prep\datasets\level 2 continuous\'], 'comparison', '\'};
            end
    
            parameters.loop_list.things_to_load.data.filename= {'PLSR_dataset_info_Cov.mat'};
            parameters.loop_list.things_to_load.data.variable= {'dataset_info.responseVariables'}; 
            parameters.loop_list.things_to_load.data.level = 'comparison';
        else
            % correlations 
            parameters.loop_list.things_to_load.data.dir = {[parameters.dir_exper 'corr fluor\only relevant corrs\continuous\' category '\'], 'comparison', '\'};
            parameters.loop_list.things_to_load.data.filename= {'average_by_nodes.mat'};
            parameters.loop_list.things_to_load.data.variable= {'average_by_nodes'}; 
            parameters.loop_list.things_to_load.data.level = 'comparison';
        end
        
        % Outputs
        parameters.loop_list.things_to_save.data_reshaped.dir = {[parameters.dir_exper 'corr fluor\reshaped continuous\' category '\'], 'comparison', '\'};
        parameters.loop_list.things_to_save.data_reshaped.filename= {data_type, '_reshaped.mat'};
        parameters.loop_list.things_to_save.data_reshaped.variable= {'data'}; 
        parameters.loop_list.things_to_save.data_reshaped.level = 'comparison';
       
        RunAnalysis({@ReshapeContinuousData}, parameters); 

    end 
end

%% Pad missing mouse
if isfield(parameters, 'loop_list')
parameters = rmfield(parameters,'loop_list');
end

% Iterators
parameters.loop_list.iterators = {
               'category', {'loop_variables.categories'}, 'category_iterator';
               'comparison', {'loop_variables.comparisons_continuous.', 'category', '(:).name'}, 'comparison_iterator' ;    
               };

% number of mice 
parameters.number_of_mice = 7; 

% placement where mouse 1100 is missing
parameters.placement = 4; 

% dimension difference mice are in 
parameters.mouseDim = 1; 

% Inputs
% Average correlations
% 7 x 16 x 4
parameters.loop_list.things_to_load.correlations.dir = {[parameters.dir_exper 'corr fluor\reshaped continuous\'], 'category', '\', 'comparison', '\'};
parameters.loop_list.things_to_load.correlations.filename= {'corrs_reshaped.mat'};
parameters.loop_list.things_to_load.correlations.variable= {'data'}; 
parameters.loop_list.things_to_load.correlations.level = 'comparison';

% Fluorescence
% 7 x 16 x 4
parameters.loop_list.things_to_load.fluorescence.dir = {[parameters.dir_exper 'corr fluor\reshaped continuous\'], 'category', '\', 'comparison', '\'};
parameters.loop_list.things_to_load.fluorescence.filename= {'fluors_reshaped.mat'};
parameters.loop_list.things_to_load.fluorescence.variable= {'data'}; 
parameters.loop_list.things_to_load.fluorescence.level = 'comparison';

% Outputs 
parameters.loop_list.things_to_save.corrs_padded.dir = {[parameters.dir_exper 'corr fluor\data padded\continuous\']}; 
parameters.loop_list.things_to_save.corrs_padded.filename= {'corrs_padded_', 'comparison', '.mat'};
parameters.loop_list.things_to_save.corrs_padded.variable= {'corrs'}; 
parameters.loop_list.things_to_save.corrs_padded.level = 'comparison';

parameters.loop_list.things_to_save.fluors_padded.dir = {[parameters.dir_exper 'corr fluor\data padded\continuous\']}; 
parameters.loop_list.things_to_save.fluors_padded.filename= {'fluors_padded_', 'comparison', '.mat'};
parameters.loop_list.things_to_save.fluors_padded.variable= {'fluors'}; 
parameters.loop_list.things_to_save.fluors_padded.level = 'comparison';

RunAnalysis({@PadMissingMouse}, parameters);  

%% Concatenate across comparisons
if isfield(parameters, 'loop_list')
parameters = rmfield(parameters,'loop_list');
end

% Iterators
parameters.loop_list.iterators = {
   'data_type', {'loop_variables.data_type'}, 'data_type_iterator';
   'comparison', {'loop_variables.comparisons_continuous_both(:).name'}, 'comparison_iterator' ; 
   };

parameters.concatenation_level = 'comparison'; 
parameters.concatDim = 4;

% Inputs
parameters.loop_list.things_to_load.data.dir = {[parameters.dir_exper 'corr fluor\data padded\continuous\']}; 
parameters.loop_list.things_to_load.data.filename= {'data_type', '_padded_', 'comparison', '.mat'};
parameters.loop_list.things_to_load.data.variable= {'data_type'}; 
parameters.loop_list.things_to_load.data.level = 'comparison';

% Output
parameters.loop_list.things_to_save.concatenated_data.dir = {[parameters.dir_exper 'corr fluor\data concatenated across comparisons\continuous\']}; 
parameters.loop_list.things_to_save.concatenated_data.filename= {'data_type', '_concatenated.mat'};
parameters.loop_list.things_to_save.concatenated_data.variable= {'data_type'}; 
parameters.loop_list.things_to_save.concatenated_data.level = 'data_type';

RunAnalysis({@ConcatenateData}, parameters);

%% Run correlations 
% Correlate within mice
% Average across mice

if isfield(parameters, 'loop_list')
parameters = rmfield(parameters,'loop_list');
end

% Iterators
parameters.loop_list.iterators = 'none';

parameters.comparison_type = 'continuous'; 
parameters.mouseDim = 1;
parameters.observationDim = 4;
parameters.variableDim = 3; 
parameters.nodeDim = 2;
               
% Inputs 
% Correlations
parameters.loop_list.things_to_load.correlations.dir = {[parameters.dir_exper 'corr fluor\data concatenated across comparisons\continuous\']}; 
parameters.loop_list.things_to_load.correlations.filename= {'corrs_concatenated.mat'};
parameters.loop_list.things_to_load.correlations.variable= {'corrs'}; 
parameters.loop_list.things_to_load.correlations.level = 'corrs';

% Fluorescence
parameters.loop_list.things_to_load.fluorescence.dir = {[parameters.dir_exper 'corr fluor\data concatenated across comparisons\continuous\']}; 
parameters.loop_list.things_to_load.fluorescence.filename= {'fluors_concatenated.mat'};
parameters.loop_list.things_to_load.fluorescence.variable= {'fluors'}; 
parameters.loop_list.things_to_load.fluorescence.level = 'fluors';

% Outputs
% correlations per mouse p values 
parameters.loop_list.things_to_save.corrs_per_mouse.dir = {[parameters.dir_exper 'corr fluor\results\continuous\']}; 
parameters.loop_list.things_to_save.corrs_per_mouse.filename= {'corrs_per_mouse.mat'};
parameters.loop_list.things_to_save.corrs_per_mouse.variable= {'corrs_per_mouse'}; 
parameters.loop_list.things_to_save.corrs_per_mouse.level = 'end';

parameters.loop_list.things_to_save.corrs_per_mouse_pvalues.dir = {[parameters.dir_exper 'corr fluor\results\continuous\']}; 
parameters.loop_list.things_to_save.corrs_per_mouse_pvalues.filename= {'corrs_per_mouse_pvalues.mat'};
parameters.loop_list.things_to_save.corrs_per_mouse_pvalues.variable= {'corrs_per_mouse_pvalues'}; 
parameters.loop_list.things_to_save.corrs_per_mouse_pvalues.level = 'end';

% average correlations across mice 
parameters.loop_list.things_to_save.corrs_across_mice.dir = {[parameters.dir_exper 'corr fluor\results\continuous\']}; 
parameters.loop_list.things_to_save.corrs_across_mice.filename= {'corrs_across_mice.mat'};
parameters.loop_list.things_to_save.corrs_across_mice.variable= {'corrs_across_mice'}; 
parameters.loop_list.things_to_save.corrs_across_mice.level = 'end';

RunAnalysis({@CorrFluor}, parameters); 

%% Renumber continuous
load('Y:\Sarah\Analysis\Experiments\Random Motorized Treadmill\corr fluor\results\continuous\corrs_across_mice.mat');
corrs_across_mice = ArrangeNewNumbering(corrs_across_mice, parameters.node_renumbering, false, 1, []);
save('Y:\Sarah\Analysis\Experiments\Random Motorized Treadmill\corr fluor\results\continuous\corrs_across_mice_renumbered.mat', 'corrs_across_mice');

%% Continuous significance

load('Y:\Sarah\Analysis\Experiments\Random Motorized Treadmill\corr fluor\results\continuous\corrs_per_mouse.mat')
corrs_per_mouse_fisher = atanh(corrs_per_mouse);

hs = NaN(16, 4);
ps = NaN(16, 4);

for variablei = 1:4

    [h, p, ci, stats1] = ttest(corrs_per_mouse_fisher(:,:, variablei), [], 'Alpha', 0.05/16/5);
    hs(:, variablei) = h';
    ps(:, variablei) = p';
    stats(variablei) = stats1;
end
[hs_fdr_cont, crit_p_cont ] = fdr_bh(ps);
% renumber 
hs_fdr_cont = ArrangeNewNumbering(hs_fdr_cont, parameters.node_renumbering, false, 1, []);

%% 

fig = figure;
hold on;

% 0 axis
plot([0 16], [0 0]);
labels = {'origin'};

% categorical
load('Y:\Sarah\Analysis\Experiments\Random Motorized Treadmill\corr fluor\results\categorical\corrs_across_mice_renumbered.mat');

% not significant
indices = find(hs_fdr_cat== 0);
marker_shape = 'o';
marker_color = 'k';
if ~isempty(indices)
    plot(indices, corrs_across_mice(indices), marker_shape, 'MarkerEdgeColor', 'k', "MarkerFaceColor", marker_color, 'MarkerSize', 12 );
    labels = [labels {'categorical'}];
end 
% significant
indices = find(hs_fdr_cat == 1);
if ~isempty(indices)
    marker_shape = 'hexagram';
    plot(indices, corrs_across_mice(indices), marker_shape, 'MarkerEdgeColor', 'k', "MarkerFaceColor", marker_color, 'MarkerSize', 14 );
    labels = [labels {'categorical significant'}]; 
end 


load('Y:\Sarah\Analysis\Experiments\Random Motorized Treadmill\corr fluor\results\continuous\corrs_across_mice_renumbered.mat')

hold on;
% speed
for variablei = 1:4
    switch variablei
        case 1
            marker_color = 'r';
        case 2
            marker_color = 'y';
        case 3
            marker_color = 'b';
        case 4
            marker_color = 'g';
    end 
   
    % not significant
    indices = find(hs_fdr_cont(:, variablei) == 0);
    marker_shape = 'o';
    if ~isempty(indices)
        plot(indices, corrs_across_mice(indices, variablei),marker_shape, 'MarkerEdgeColor', 'k', "MarkerFaceColor", marker_color );
        labels = [labels parameters.continuous_variable_names(variablei)];
    end 
    % significant
    indices = find(hs_fdr_cont(:, variablei) == 1);
    if ~isempty(indices)
        marker_shape = 'hexagram';
        plot(indices, corrs_across_mice(indices, variablei),marker_shape, 'MarkerEdgeColor', 'k', "MarkerFaceColor", marker_color );
        labels = [labels {[parameters.continuous_variable_names{variablei} ' significant']}]; 
    end 
end


xlim([0 17]);
ylim([-0.6 0.6]);

legend(labels);

xticks([2:2:16]);
xticklabels({'3&4', '7&8', '11&12', '15&16', '19&20', '23&24', '27&28', '31&32'});

ylabel('correlation'); xlabel('node');
title('relationship between fluorescence and correlation');

%% across nodes
% categorical, averaged
load('Y:\Sarah\Analysis\Experiments\Random Motorized Treadmill\corr fluor\results\categorical\corrs_across_mice.mat');
a = reshape(corrs_across_mice, [], 1);
[hs1,ps1] = ttest(atanh(a));

% categorical, all
load('Y:\Sarah\Analysis\Experiments\Random Motorized Treadmill\corr fluor\results\categorical\corrs_per_mouse.mat');
a = reshape(corrs_per_mouse, [], 1);
[hs2,ps2] = ttest(atanh(a));

% anovan
group_mouse = repmat([1:7]', 1, 16);
group_node = repmat([1:16], 7, 1);
group_mouse = reshape(group_mouse, [], 1);
group_node = reshape(group_node, [], 1);

[p, tbl, stats] = anovan(atanh(a), {group_mouse, group_node});
 figure; c = multcompare(stats, 'Dimension', 2);
 figure; c = multcompare(stats, 'Dimension', 1);
