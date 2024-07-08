% pipeline_means_across_mice.m 
% Sarah West
% 11/13/23


%% Initial Setup  
% Put all needed paramters in a structure called "parameters", which you

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
if isfile([parameters.dir_exper 'PLSR\periods_nametable_forPLSR_specials.mat'])
    load([parameters.dir_exper 'PLSR\periods_nametable_forPLSR_specials.mat']);
    parameters.periods = periods;

    % Also load the indices to remove
    load([parameters.dir_exper 'PLSR\indices_to_remove_Specials.mat']);
    parameters.indices_to_remove = indices_to_remove;

    clear periods indices_to_remove;

end

% Load comparisons for first level categorical, if it exists yet.
if isfile([parameters.dir_exper 'PLSR\comparisons_special_categorical.mat'])
    load([parameters.dir_exper 'PLSR\comparisons_special_categorical.mat']);
    parameters.comparisons_categorical = comparisons(1:31);
    parameters.loop_variables.comparisons_categorical = parameters.comparisons_categorical;
    clear comparisons;
end

% Load comparisons for first level continuous, if it exists yet.
if isfile([parameters.dir_exper 'PLSR\comparisons_special_continuous.mat'])
    load([parameters.dir_exper 'PLSR\comparisons_special_continuous.mat']);
    parameters.comparisons_continuous = comparisons;
    parameters.loop_variables.comparisons_continuous = parameters.comparisons_continuous;
    clear comparisons;
end

parameters.loop_variables.mice_all = parameters.mice_all;
parameters.loop_variables.data_types = {'correlations', 'fluorescence'};
parameters.loop_variables.data_types2 = {'correlations', 'fluorescence DFF', 'fluorescence'};
parameters.loop_variables.change_types = {'mean', 'beta'};
parameters.loop_variables.comparison_types = {'categorical', 'continuous'};

% For plots
parameters.region_labels = {'M2', 'M1', 'S1', 'LP', 'MV', 'Rs'};
parameters.region_label_locations_range = [3.5 8.5 13.5 18.5 24.5 30.5] ;
parameters.grid_width_major = 2.7;
parameters.grid_width_minor = 0.25;
parameters.node_label_text = true;

% Names of all continuous variables.
parameters.continuous_variable_names = {'speed', 'accel', 'duration', 'pupil_diameter', 'tail', 'nose', 'FL', 'HL', 'x'};

%% find mean changes, after PSLR vertical correction
% use the datasets from dataset prep
parameters.comparison_type = 'categorical';

for data_typei = 1:numel(parameters.loop_variables.data_types)

    data_type = parameters.loop_variables.data_types{data_typei};
    if strcmp(data_type, 'correlations')
        input_string_name = [];
    else 
        input_string_name = ['fluorescence '];
    end 

    if isfield(parameters, 'loop_list')
    parameters = rmfield(parameters,'loop_list');
    end
    
    % Iterators
    parameters.loop_list.iterators = {
                   'mouse', {'loop_variables.mice_all(:).name'}, 'mouse_iterator'; 
                   'comparison', {'loop_variables.comparisons_categorical(:).name'}, 'comparison_iterator'     
                   };
    
    % Input 
    % dataset
    parameters.loop_list.things_to_load.dataset.dir = {[parameters.dir_exper 'PLSR ' input_string_name 'Specials Inverted Intercepts\variable prep\datasets\level 1 categorical\'], 'comparison', '\', 'mouse', '\'};
    parameters.loop_list.things_to_load.dataset.filename= {'PLSR_dataset_info.mat'};
    parameters.loop_list.things_to_load.dataset.variable= {'dataset_info'}; 
    parameters.loop_list.things_to_load.dataset.level = 'comparison';
    
    % Output 
    % mean changes
    parameters.loop_list.things_to_save.average_change.dir = {[parameters.dir_exper 'mean checks\mean changes\' data_type '\'], 'comparison', '\', 'mouse', '\'};
    parameters.loop_list.things_to_save.average_change.filename= {'mean_change.mat'};
    parameters.loop_list.things_to_save.average_change.variable= {'mean_change'}; 
    parameters.loop_list.things_to_save.average_change.level = 'comparison';
    % std of each set [set1; set2]
    parameters.loop_list.things_to_save.std_change.dir = {[parameters.dir_exper 'mean checks\mean changes\' data_type '\'], 'comparison', '\', 'mouse', '\'};
    parameters.loop_list.things_to_save.std_change.filename= {'std_of_each.mat'};
    parameters.loop_list.things_to_save.std_change.variable= {'std_of_each'}; 
    parameters.loop_list.things_to_save.std_change.level = 'comparison';

    % save each set separately to find null distributions later
    % set 1 
    parameters.loop_list.things_to_save.set1.dir = {[parameters.dir_exper 'mean checks\mean changes\' data_type '\'], 'comparison', '\', 'mouse', '\'};
    parameters.loop_list.things_to_save.set1.filename= {'set1.mat'};
    parameters.loop_list.things_to_save.set1.variable= {'set1'}; 
    parameters.loop_list.things_to_save.set1.level = 'comparison';
    % set 2 
    parameters.loop_list.things_to_save.set2.dir = {[parameters.dir_exper 'mean checks\mean changes\' data_type '\'], 'comparison', '\', 'mouse', '\'};
    parameters.loop_list.things_to_save.set2.filename= {'set2.mat'};
    parameters.loop_list.things_to_save.set2.variable= {'set2'}; 
    parameters.loop_list.things_to_save.set2.level = 'comparison';

    RunAnalysis({@FindMeanChanges}, parameters)
end 
parameters = rmfield(parameters,'comparison_type');

%% Find   for fluorescence mean changes
if isfield(parameters, 'loop_list')
parameters = rmfield(parameters,'loop_list');
end

% Iterators
parameters.loop_list.iterators = {
               'mouse', {'loop_variables.mice_all(:).name'}, 'mouse_iterator'; 
               'comparison', {'loop_variables.comparisons_categorical(:).name'}, 'comparison_iterator'     
               };
% Inputs
% data
parameters.loop_list.things_to_load.data.dir = {[parameters.dir_exper 'mean checks\mean changes\fluorescence\'], 'comparison', '\', 'mouse', '\'};
parameters.loop_list.things_to_load.data.filename= {'mean_change.mat'};
parameters.loop_list.things_to_load.data.variable= {'mean_change'}; 
parameters.loop_list.things_to_load.data.level = 'comparison';
% fluorescence means
parameters.loop_list.things_to_load.fluorescence_mean.dir = {[parameters.dir_exper '\preprocessing\stack means\']};
parameters.loop_list.things_to_load.fluorescence_mean.filename = {'IC_means_acrossMice_homologousTogether_hemo_corrected.mat'};
parameters.loop_list.things_to_load.fluorescence_mean.variable = {'source_mean'};
parameters.loop_list.things_to_load.fluorescence_mean.level = 'start';

% Ouputs
parameters.loop_list.things_to_save.DFF.dir = {[parameters.dir_exper 'mean checks\mean changes\fluorescence DFF\'], 'comparison', '\', 'mouse', '\'};
parameters.loop_list.things_to_save.DFF.filename= {'mean_change.mat'};
parameters.loop_list.things_to_save.DFF.variable= {'mean_change'}; 
parameters.loop_list.things_to_save.DFF.level = 'comparison';

RunAnalysis({@DFF}, parameters);

%% Find change per mouse from PLSR, sigma adjusted 

for comparison_typei = 1:numel(parameters.loop_variables.comparison_types)
    comparison_type = parameters.loop_variables.comparison_types{comparison_typei};
    parameters.comparison_type = comparison_type;
    for data_typei = 1:numel(parameters.loop_variables.data_types)
    
        data_type = parameters.loop_variables.data_types{data_typei};
        if strcmp(data_type, 'correlations')
            input_string_name = [];
        else 
            input_string_name = ['fluorescence '];
        end 
        if isfield(parameters, 'loop_list')
        parameters = rmfield(parameters,'loop_list');
        end
        
        % Iterators
        parameters.loop_list.iterators = {
                       'mouse', {'loop_variables.mice_all(:).name'}, 'mouse_iterator'; 
                       'comparison', {['loop_variables.comparisons_' comparison_type '(:).name']}, 'comparison_iterator'     
                       };
        
        % Inputs
        % dataset (for sigmas)
        parameters.loop_list.things_to_load.dataset.dir = {[parameters.dir_exper 'PLSR ' input_string_name 'Specials Inverted Intercepts\variable prep\datasets\level 1 ' comparison_type '\'], 'comparison', '\', 'mouse', '\'};
        parameters.loop_list.things_to_load.dataset.filename= {'PLSR_dataset_info.mat'};
        parameters.loop_list.things_to_load.dataset.variable= {'dataset_info'}; 
        parameters.loop_list.things_to_load.dataset.level = 'comparison';    
        % PSLR Betas
        parameters.loop_list.things_to_load.betas.dir = {[parameters.dir_exper 'PLSR ' input_string_name 'Specials Inverted Intercepts\results\level 1 ' comparison_type '\'], 'comparison', '\', 'mouse', '\'};
        parameters.loop_list.things_to_load.betas.filename= {'PLSR_results.mat'};
        parameters.loop_list.things_to_load.betas.variable= {'PLSR_results.BETA'}; 
        parameters.loop_list.things_to_load.betas.level = 'comparison';
        
        % Outputs
        parameters.loop_list.things_to_save.beta_change.dir = {[parameters.dir_exper 'mean checks\beta changes\' comparison_type, '\' data_type '\'], 'comparison', '\', 'mouse', '\'};
        parameters.loop_list.things_to_save.beta_change.filename= {'beta_change.mat'};
        parameters.loop_list.things_to_save.beta_change.variable= {'beta_change'}; 
        parameters.loop_list.things_to_save.beta_change.level = 'comparison';
    
        RunAnalysis({@FindBetaChanges}, parameters)
    end
end 
parameters = rmfield(parameters,'comparison_type');

%% Find DFF for fluorescence Betas
for comparison_typei = 2 %1:numel(parameters.loop_variables.comparison_types)
    
    comparison_type = parameters.loop_variables.comparison_types{comparison_typei};
    parameters.comparison_type = comparison_type;
    
    if isfield(parameters, 'loop_list')
    parameters = rmfield(parameters,'loop_list');
    end
    
    % Iterators
    parameters.loop_list.iterators = {
                   'mouse', {'loop_variables.mice_all(:).name'}, 'mouse_iterator'; 
                   'comparison', {['loop_variables.comparisons_' comparison_type '(:).name']}, 'comparison_iterator'     
                   };
    % Inputs
    % data
    parameters.loop_list.things_to_load.data.dir = {[parameters.dir_exper 'mean checks\beta changes\' comparison_type, '\fluorescence\'], 'comparison', '\', 'mouse', '\'};
    parameters.loop_list.things_to_load.data.filename= {'beta_change.mat'};
    parameters.loop_list.things_to_load.data.variable= {'beta_change'}; 
    parameters.loop_list.things_to_load.data.level = 'comparison';
    % fluorescence means
    parameters.loop_list.things_to_load.fluorescence_mean.dir = {[parameters.dir_exper '\preprocessing\stack means\']};
    parameters.loop_list.things_to_load.fluorescence_mean.filename = {'IC_means_acrossMice_homologousTogether_hemo_corrected.mat'};
    parameters.loop_list.things_to_load.fluorescence_mean.variable = {'source_mean'};
    parameters.loop_list.things_to_load.fluorescence_mean.level = 'start';
    
    % Ouputs
    parameters.loop_list.things_to_save.DFF.dir = {[parameters.dir_exper 'mean checks\beta changes\' comparison_type '\fluorescence DFF\'], 'comparison', '\', 'mouse', '\'};
    parameters.loop_list.things_to_save.DFF.filename= {'beta_change.mat'};
    parameters.loop_list.things_to_save.DFF.variable= {'beta_change'}; 
    parameters.loop_list.things_to_save.DFF.level = 'comparison';
    
    RunAnalysis({@DFF}, parameters);
end
parameters = rmfield(parameters,'comparison_type');

%% means of changes across mice: categorical
% do just beta changes (ran mean changes previously)
if isfield(parameters, 'loop_list')
parameters = rmfield(parameters,'loop_list');
end

% Iterators
parameters.loop_list.iterators = {
               'data_type', 'loop_variables.data_types2(3)', 'data_type_iterator';
               'change_type', 'loop_variables.change_types(2)', 'change_type_iterator';
               'comparison', {'loop_variables.comparisons_categorical(:).name'}, 'comparison_iterator';
               'mouse', {'loop_variables.mice_all(:).name'}, 'mouse_iterator'; 
               };
parameters.comparison_type = 'categorical';
parameters.concatDim = 1;
parameters.concatenation_level = 'mouse';
parameters.averageDim = 1;

% Inputs
parameters.loop_list.things_to_load.data.dir = {[parameters.dir_exper 'mean checks\'], 'change_type', ' changes\categorical\', 'data_type', '\', 'comparison', '\', 'mouse', '\'};
parameters.loop_list.things_to_load.data.filename= {'change_type', '_change.mat'};
parameters.loop_list.things_to_load.data.variable= {'change_type', '_change'}; 
parameters.loop_list.things_to_load.data.level = 'mouse';

% Outputs
% all together
parameters.loop_list.things_to_save.concatenated_data.dir = {[parameters.dir_exper 'mean checks\'], 'change_type', ' changes\categorical\', 'data_type', '\', 'comparison', '\mean across mice\'};
parameters.loop_list.things_to_save.concatenated_data.filename= {'changes_all.mat'};
parameters.loop_list.things_to_save.concatenated_data.variable= {'changes_all'}; 
parameters.loop_list.things_to_save.concatenated_data.level = 'comparison';
% average
parameters.loop_list.things_to_save.average.dir = {[parameters.dir_exper 'mean checks\'], 'change_type', ' changes\categorical\', 'data_type', '\', 'comparison', '\mean across mice\'};
parameters.loop_list.things_to_save.average.filename= {'changes_average.mat'};
parameters.loop_list.things_to_save.average.variable= {'changes_average'}; 
parameters.loop_list.things_to_save.average.level = 'comparison';
% std_dev
parameters.loop_list.things_to_save.std_dev.dir = {[parameters.dir_exper 'mean checks\'], 'change_type', ' changes\categorical\', 'data_type', '\', 'comparison', '\mean across mice\'};
parameters.loop_list.things_to_save.std_dev.filename= {'changes_std_dev.mat'};
parameters.loop_list.things_to_save.std_dev.variable= {'changes_std_dev'}; 
parameters.loop_list.things_to_save.std_dev.level = 'comparison';

parameters.loop_list.things_to_rename = {{};
                                         {'concatenated_data', 'data'}};

RunAnalysis({@UseThisMouse, @ConcatenateData, @AverageData}, parameters)

parameters = rmfield(parameters,'comparison_type');

%% means of changes across mice: continuous 
if isfield(parameters, 'loop_list')
parameters = rmfield(parameters,'loop_list');
end

% Iterators
parameters.loop_list.iterators = {
               'data_type', 'loop_variables.data_types2', 'data_type_iterator';
               'comparison', {'loop_variables.comparisons_continuous(:).name'}, 'comparison_iterator';
               'mouse', {'loop_variables.mice_all(:).name'}, 'mouse_iterator'; 
               };
parameters.comparison_type = 'continuous';
parameters.concatDim = 3;
parameters.concatenation_level = 'mouse';
parameters.averageDim = 3;

% Inputs
parameters.loop_list.things_to_load.data.dir = {[parameters.dir_exper 'mean checks\beta changes\continuous\'], 'data_type', '\', 'comparison', '\', 'mouse', '\'};
parameters.loop_list.things_to_load.data.filename= {'beta_change.mat'};
parameters.loop_list.things_to_load.data.variable= {'beta_change'}; 
parameters.loop_list.things_to_load.data.level = 'mouse';

% Outputs
% all together
parameters.loop_list.things_to_save.concatenated_data.dir = {[parameters.dir_exper 'mean checks\beta changes\continuous\'], 'data_type', '\', 'comparison', '\mean across mice\'};
parameters.loop_list.things_to_save.concatenated_data.filename= {'changes_all.mat'};
parameters.loop_list.things_to_save.concatenated_data.variable= {'changes_all'}; 
parameters.loop_list.things_to_save.concatenated_data.level = 'comparison';
% average
parameters.loop_list.things_to_save.average.dir = {[parameters.dir_exper 'mean checks\beta changes\continuous\'],'data_type', '\', 'comparison', '\mean across mice\'};
parameters.loop_list.things_to_save.average.filename= {'changes_average.mat'};
parameters.loop_list.things_to_save.average.variable= {'changes_average'}; 
parameters.loop_list.things_to_save.average.level = 'comparison';
% std_dev
parameters.loop_list.things_to_save.std_dev.dir = {[parameters.dir_exper 'mean checks\beta changes\continuous\'],'data_type', '\', 'comparison', '\mean across mice\'};
parameters.loop_list.things_to_save.std_dev.filename= {'changes_std_dev.mat'};
parameters.loop_list.things_to_save.std_dev.variable= {'changes_std_dev'}; 
parameters.loop_list.things_to_save.std_dev.level = 'comparison';

parameters.loop_list.things_to_rename = {{};
                                         {'concatenated_data', 'data'}};

RunAnalysis({@UseThisMouse, @ConcatenateData, @AverageData}, parameters)
parameters = rmfield(parameters,'comparison_type');

%% Continuous fluorescence data: Average fluorescence accross hemispheres for some hemispheres
if isfield(parameters, 'loop_list')
parameters = rmfield(parameters,'loop_list');
end

% Iterators
parameters.loop_list.iterators = {
               'comparison', {'loop_variables.comparisons_continuous(:).name'}, 'comparison_iterator';     
               };
parameters.output_type = 'BETA';                      
parameters.variables_to_average_across_hemispheres = {'speed_vector', 'accel_vector', 'duration_vector', 'pupil_diameter_vector', 'tail_vector', 'nose_vector'};
parameters.nodeDim = 2;
parameters.variableDim = 1; 

% Input 
parameters.loop_list.things_to_load.data.dir = {[parameters.dir_exper 'mean checks\beta changes\continuous\fluorescence DFF\'], 'comparison', '\mean across mice\'};
parameters.loop_list.things_to_load.data.filename= {'changes_average.mat'};
parameters.loop_list.things_to_load.data.variable= {'changes_average'}; 
parameters.loop_list.things_to_load.data.level = 'comparison';
 
% Output
% Keep same name so it can be used for iterations below
parameters.loop_list.things_to_save.output.dir = {[parameters.dir_exper 'mean checks\beta changes\continuous\fluorescence DFF\'], '\', 'comparison', '\mean across mice\'};
parameters.loop_list.things_to_save.output.filename= {'changes_average.mat'};
parameters.loop_list.things_to_save.output.variable= {'changes_average'}; 
parameters.loop_list.things_to_save.output.level = 'comparison';

RunAnalysis({@FluorescenceContinuousAverageHemispheres_Inverted_Modified}, parameters)

%% Correlations: Calculate ipsa-contra averaging on random permutations
% (this wasn't done when you run random permutations on the inverted PLSR;
% the correlations were suffled independently)


%% Average null distributions by each mouse's sigmas 

for comparison_typei = 1:numel(parameters.loop_variables.comparison_types)
    comparison_type = parameters.loop_variables.comparison_types{comparison_typei};
    parameters.comparison_type = comparison_type;
    for data_typei = 1:numel(parameters.loop_variables.data_types)
    
        data_type = parameters.loop_variables.data_types{data_typei};
        if strcmp(data_type, 'correlations')
            input_string_name = [];
        else 
            input_string_name = ['fluorescence '];
        end 
        if isfield(parameters, 'loop_list')
        parameters = rmfield(parameters,'loop_list');
        end
        
        % Iterators
        parameters.loop_list.iterators = {
                       'mouse', {'loop_variables.mice_all(:).name'}, 'mouse_iterator'; 
                       'comparison', {['loop_variables.comparisons_' comparison_type '(:).name']}, 'comparison_iterator'     
                       };

        parameters.onPermutations = true;
        
        % Inputs
        % dataset (for sigmas)
        parameters.loop_list.things_to_load.dataset.dir = {[parameters.dir_exper 'PLSR ' input_string_name 'Specials Inverted Intercepts\variable prep\datasets\level 1 ' comparison_type '\'], 'comparison', '\', 'mouse', '\'};
        parameters.loop_list.things_to_load.dataset.filename= {'PLSR_dataset_info.mat'};
        parameters.loop_list.things_to_load.dataset.variable= {'dataset_info'}; 
        parameters.loop_list.things_to_load.dataset.level = 'comparison';    
        % PSLR Betas
        parameters.loop_list.things_to_load.betas.dir = {[parameters.dir_exper 'PLSR ' input_string_name 'Specials Inverted Intercepts\results\level 1 ' comparison_type '\'], 'comparison', '\', 'mouse', '\'};
        parameters.loop_list.things_to_load.betas.filename= {'PLSR_BETAs_randomPermutations.mat'};
        parameters.loop_list.things_to_load.betas.variable= {'BETAs_randomPermutations'}; 
        parameters.loop_list.things_to_load.betas.level = 'comparison';
        
        % Outputs
        parameters.loop_list.things_to_save.beta_change.dir = {[parameters.dir_exper 'mean checks\beta changes\' comparison_type, '\' data_type '\'], 'comparison', '\', 'mouse', '\'};
        parameters.loop_list.things_to_save.beta_change.filename= {'beta_change_randomPermutations.mat'};
        parameters.loop_list.things_to_save.beta_change.variable= {'beta_change'}; 
        parameters.loop_list.things_to_save.beta_change.level = 'comparison';
    
        RunAnalysis({@FindBetaChanges}, parameters)
    end
end 
parameters = rmfield(parameters,'comparison_type');

%% Fluorescence: find DFF of null distributions 
for comparison_typei = 1:numel(parameters.loop_variables.comparison_types)
    
    comparison_type = parameters.loop_variables.comparison_types{comparison_typei};
    parameters.comparison_type = comparison_type;
    
    if isfield(parameters, 'loop_list')
    parameters = rmfield(parameters,'loop_list');
    end
    
    % Iterators
    parameters.loop_list.iterators = {
                   'mouse', {'loop_variables.mice_all(:).name'}, 'mouse_iterator'; 
                   'comparison', {['loop_variables.comparisons_' comparison_type '(:).name']}, 'comparison_iterator'     
                   };
    % Inputs
    % data
    parameters.loop_list.things_to_load.data.dir = {[parameters.dir_exper 'mean checks\beta changes\' comparison_type, '\fluorescence\'], 'comparison', '\', 'mouse', '\'};
    parameters.loop_list.things_to_load.data.filename= {'beta_change_randomPermutations.mat'};
    parameters.loop_list.things_to_load.data.variable= {'beta_change'}; 
    parameters.loop_list.things_to_load.data.level = 'comparison';
    % fluorescence means
    parameters.loop_list.things_to_load.fluorescence_mean.dir = {[parameters.dir_exper '\preprocessing\stack means\']};
    parameters.loop_list.things_to_load.fluorescence_mean.filename = {'IC_means_acrossMice_homologousTogether_hemo_corrected.mat'};
    parameters.loop_list.things_to_load.fluorescence_mean.variable = {'source_mean'};
    parameters.loop_list.things_to_load.fluorescence_mean.level = 'start';
    
    % Ouputs
    parameters.loop_list.things_to_save.DFF.dir = {[parameters.dir_exper 'mean checks\beta changes\' comparison_type '\fluorescence DFF\'], 'comparison', '\', 'mouse', '\'};
    parameters.loop_list.things_to_save.DFF.filename= {'beta_change_randomPermutations.mat'};
    parameters.loop_list.things_to_save.DFF.variable= {'beta_change'}; 
    parameters.loop_list.things_to_save.DFF.level = 'comparison';
    
    RunAnalysis({@DFF}, parameters);
end
parameters = rmfield(parameters,'comparison_type');

%% Average null distributions across mice

if isfield(parameters, 'loop_list')
parameters = rmfield(parameters,'loop_list');
end

% Iterators
parameters.loop_list.iterators = {
               'data_type', {'loop_variables.data_types2'}, 'data_type_iterator';
               'comparison_type', {'loop_variables.comparison_types'}, 'comparison_type_iterator';
               'comparison', {'loop_variables.comparisons_', 'comparison_type', '(:).name'}, 'comparison_iterator';
               'mouse', {'loop_variables.mice_all(:).name'}, 'mouse_iterator'; 
               };

parameters.concatDim = 4;
parameters.concatenation_level = 'mouse';
parameters.averageDim = 4;

% Inputs
parameters.loop_list.things_to_load.data.dir = {[parameters.dir_exper 'mean checks\beta changes\'], 'comparison_type', '\', 'data_type', '\', 'comparison', '\', 'mouse', '\'};
parameters.loop_list.things_to_load.data.filename= {'beta_change_randomPermutations.mat'};
parameters.loop_list.things_to_load.data.variable= {'beta_change'}; 
parameters.loop_list.things_to_load.data.level = 'comparison';

% Outputs
% all together
parameters.loop_list.things_to_save.concatenated_data.dir = {[parameters.dir_exper 'mean checks\beta changes\'], 'comparison_type', '\', 'data_type', '\', 'comparison', '\mean across mice\'};
parameters.loop_list.things_to_save.concatenated_data.filename= {'randomPermutations_all.mat'};
parameters.loop_list.things_to_save.concatenated_data.variable= {'randomPermutations_all'}; 
parameters.loop_list.things_to_save.concatenated_data.level = 'comparison';
% average
parameters.loop_list.things_to_save.average.dir = {[parameters.dir_exper 'mean checks\beta changes\'], 'comparison_type', '\','data_type', '\', 'comparison', '\mean across mice\'};
parameters.loop_list.things_to_save.average.filename= {'randomPermutations_average.mat'};
parameters.loop_list.things_to_save.average.variable= {'randomPermutations_average'}; 
parameters.loop_list.things_to_save.average.level = 'comparison';

parameters.loop_list.things_to_rename = {{};
                                         {'concatenated_data', 'data'}};

RunAnalysis({@UseThisMouse, @ConcatenateData, @AverageData}, parameters)

%% Continuous fluorescence data: Average null distributions of fluorescence accross hemispheres for some hemispheres
if isfield(parameters, 'loop_list')
parameters = rmfield(parameters,'loop_list');
end

% Iterators
parameters.loop_list.iterators = {
               'comparison', {'loop_variables.comparisons_continuous(:).name'}, 'comparison_iterator';     
               'mouse', {'loop_variables.mice_all(:).name'}, 'mouse_iterator';
               };
parameters.output_type = 'BETA';                      
parameters.variables_to_average_across_hemispheres = {'speed_vector', 'accel_vector', 'duration_vector', 'pupil_diameter_vector', 'tail_vector', 'nose_vector'};
parameters.nodeDim = 2;
parameters.variableDim = 1; 

% Input 
parameters.loop_list.things_to_load.data.dir = {[parameters.dir_exper 'mean checks\beta changes\continuous\fluorescence DFF\'], 'comparison', '\mean across mice\'};
parameters.loop_list.things_to_load.data.filename= {'randomPermutations_average.mat'};
parameters.loop_list.things_to_load.data.variable= {'randomPermutations_average'}; 
parameters.loop_list.things_to_load.data.level = 'comparison';
 
% Output
% Make the same name so you can use it in loop iterations below
parameters.loop_list.things_to_save.output.dir = {[parameters.dir_exper 'mean checks\beta changes\continuous\fluorescence DFF\'], '\', 'comparison', '\mean across mice\'};
parameters.loop_list.things_to_save.output.filename= {'randomPermutations_average.mat'};
parameters.loop_list.things_to_save.output.variable= {'randomPermutations_average'}; 
parameters.loop_list.things_to_save.output.level = 'comparison';

RunAnalysis({@FluorescenceContinuousAverageHemispheres_Inverted_Modified}, parameters)

%%  Check significance 
% Always clear loop list first. 

for comparison_typei = 1:numel(parameters.loop_variables.comparison_types) 

    comparison_type = parameters.loop_variables.comparison_types{comparison_typei};

    % If fluorescence, use _acrossHemispheres version 

    for data_typei = 1:numel(parameters.loop_variables.data_types2)
        data_type = parameters.loop_variables.data_types2{data_typei};

        if isfield(parameters, 'loop_list')
        parameters = rmfield(parameters,'loop_list');
        end
    
        if strcmp(comparison_type, 'categorical')
            parameters.shufflesDim = 2;
        else 
            parameters.shufflesDim = 3; 
        end
    
        % Iterators
        parameters.loop_list.iterators = {
                       'comparison', {'loop_variables.comparisons_', comparison_type, '(:).name'}, 'comparison_iterator';
                       }; 
       
        parameters.find_significance = true;
        
        % The statistical alpha value
        parameters.alphaValue = 0.05;  %0.05/496; %0.001; %/numel(parameters.comparisons_continuous);
        
        % Say that you do want to use bootstrapping.
        parameters.useBootstrapping = false;
        
        % If you want to fit a normal distribution before t-test (default = true)
        parameters.useNormalDistribution = true; 
        
        % Use false discovery rate correction?
        parameters.useFDR = true;
        % Inputs:
        %Test values (will grab only the intercepts with EvaluateOnData)
        parameters.loop_list.things_to_load.test_values.dir = {[parameters.dir_exper 'mean checks\beta changes\' comparison_type '\' data_type, '\'], 'comparison', '\mean across mice\'};
        if strcmp(data_type, 'fluorescence DFF') & strcmp(comparison_type, 'continuous')
            parameters.loop_list.things_to_load.test_values.filename= {'changes_average_acrossHemispheres.mat'};
        else
            parameters.loop_list.things_to_load.test_values.filename= {'changes_average.mat'};
        end
        parameters.loop_list.things_to_load.test_values.variable= {'changes_average'}; 
        parameters.loop_list.things_to_load.test_values.level = 'comparison';
        % Null distribution
        parameters.loop_list.things_to_load.null_distribution.dir = {[parameters.dir_exper 'mean checks\beta changes\' comparison_type '\', data_type '\'], 'comparison', '\mean across mice\'};
        if strcmp(data_type, 'fluorescence DFF') & strcmp(comparison_type, 'continuous')
            parameters.loop_list.things_to_load.null_distribution.filename= {'randomPermutations_average_acrossHemispheres.mat'};
        else
            parameters.loop_list.things_to_load.null_distribution.filename= {'randomPermutations_average.mat'};
        end 
        parameters.loop_list.things_to_load.null_distribution.variable= {'randomPermutations_average'}; 
        parameters.loop_list.things_to_load.null_distribution.level = 'comparison';
        
        % Outputs
        parameters.loop_list.things_to_save.significance.dir = {[parameters.dir_exper 'mean checks\beta changes\' comparison_type '\', data_type, '\'], 'comparison', '\mean across mice\'};
        parameters.loop_list.things_to_save.significance.filename= {'significance.mat'};
        parameters.loop_list.things_to_save.significance.variable= {'significance'}; 
        parameters.loop_list.things_to_save.significance.level = 'comparison';
        
        RunAnalysis({@SignificanceCalculation}, parameters);
    end 
end    
parameters.useFDR = false;

%% Correlations: average for dot plots

for comparison_typei =1:numel(parameters.loop_variables.comparison_types) 

    comparison_type = parameters.loop_variables.comparison_types{comparison_typei};

    % Dimension to average across AFTER data has gone through AverageByNode
    % code.
    if strcmp(comparison_type, 'categorical')
        parameters.averageDim = 2;
    else 
        parameters.averageDim = 3; 
    end

    if isfield(parameters, 'loop_list')
        parameters = rmfield(parameters,'loop_list');
    end
    % Iterators
    parameters.loop_list.iterators = {
                   
                   'comparison', {['loop_variables.comparisons_' comparison_type '(:).name']}, 'comparison_iterator' };
    
    parameters.fromPLSR = true; 
    
    % Parameters for AverageByNode code.
    parameters.isVector = true;
    parameters.corrsDim = 2;
    
    % Multiply by average sigmas?
    parameters.multiply_by_average_sigma = false;
    
    % Input
    parameters.loop_list.things_to_load.data.dir = {[parameters.dir_exper 'mean checks\beta changes\' comparison_type '\'], '\correlations\', 'comparison', '\mean across mice\'};
    parameters.loop_list.things_to_load.data.filename = {'changes_all.mat'};
    parameters.loop_list.things_to_load.data.variable = {'changes_all'};
    parameters.loop_list.things_to_load.data.level = 'comparison';
    
    % Output 
    % each mouse, as a matrix
    parameters.loop_list.things_to_save.node_averages.dir = {[parameters.dir_exper 'mean checks\beta changes\' comparison_type '\'], '\correlations\', 'comparison', '\mean across mice\'};
    parameters.loop_list.things_to_save.node_averages.filename = {'average_by_nodes_all.mat'};
    parameters.loop_list.things_to_save.node_averages.variable = {'average_by_nodes'};
    parameters.loop_list.things_to_save.node_averages.level = 'comparison';
    % Across mice.
    parameters.loop_list.things_to_save.average.dir =  {[parameters.dir_exper 'mean checks\beta changes\' comparison_type '\'], '\correlations\', 'comparison', '\mean across mice\'};
    parameters.loop_list.things_to_save.average.filename = {'average_by_nodes_across_mice.mat'};
    parameters.loop_list.things_to_save.average.variable = {'average_by_nodes'};
    parameters.loop_list.things_to_save.average.level = 'comparison';
    
    parameters.loop_list.things_to_rename = {;
                                              {'node_averages', 'data'}}; 
    
    RunAnalysis({@AverageByNode_Inverted, @AverageData}, parameters);

end

%% Correlations: Permute random permutations for next steps
if isfield(parameters, 'loop_list')
parameters = rmfield(parameters,'loop_list');
end

% Iterators
parameters.loop_list.iterators = {
               'comparison_type', {'loop_variables.comparison_types'}, 'comparison_type_iterator';
               'comparison', {'loop_variables.comparisons_', 'comparison_type', '(:).name'}, 'comparison_iterator';
               };

parameters.evaluation_instructions = {{'data = parameters.data;'...
                                        'data_evaluated = permute(data, [1 2 4 3]);'}};

% Input
parameters.loop_list.things_to_load.data.dir = {[parameters.dir_exper 'mean checks\beta changes\'],  'comparison_type', '\correlations\', 'comparison', '\mean across mice\'};
parameters.loop_list.things_to_load.data.filename = {'randomPermutations_all.mat'};
parameters.loop_list.things_to_load.data.variable = {'randomPermutations_all'};
parameters.loop_list.things_to_load.data.level = 'comparison';

% output 
parameters.loop_list.things_to_save.data_evaluated.dir = {[parameters.dir_exper 'mean checks\beta changes\'],  'comparison_type', '\correlations\', 'comparison', '\mean across mice\'};
parameters.loop_list.things_to_save.data_evaluated.filename = {'randomPermutations_all_permuted.mat'};
parameters.loop_list.things_to_save.data_evaluated.variable = {'randomPermutations_all'};
parameters.loop_list.things_to_save.data_evaluated.level = 'comparison';

RunAnalysis({@EvaluateOnData}, parameters);

%% Correlations: average null distribtions for dot plots

parameters.onPermutations = true; 
parameters.fromPLSR = false;
for comparison_typei = 1:numel(parameters.loop_variables.comparison_types) 

    comparison_type = parameters.loop_variables.comparison_types{comparison_typei};

    % Dimension to average across AFTER data has gone through AverageByNode
    % code.

    parameters.averageDim = 2; 

    if isfield(parameters, 'loop_list')
        parameters = rmfield(parameters,'loop_list');
    end
    % Iterators
    parameters.loop_list.iterators = {
                   
                   'comparison', {['loop_variables.comparisons_' comparison_type '(:).name']}, 'comparison_iterator' };
    
    parameters.fromPLSR = true; 
    
    % Parameters for AverageByNode code.
    parameters.isVector = true;
    parameters.corrsDim = 2;
    
    % Multiply by average sigmas?
    parameters.multiply_by_average_sigma = false;
    
    % Input
    parameters.loop_list.things_to_load.data.dir = {[parameters.dir_exper 'mean checks\beta changes\' comparison_type '\'], '\correlations\', 'comparison', '\mean across mice\'};
    parameters.loop_list.things_to_load.data.filename = {'randomPermutations_all_permuted.mat'};
    parameters.loop_list.things_to_load.data.variable = {'randomPermutations_all'};
    parameters.loop_list.things_to_load.data.level = 'comparison';
    
    % Output 
    % each mouse, as a matrix
    parameters.loop_list.things_to_save.node_averages.dir = {[parameters.dir_exper 'mean checks\beta changes\' comparison_type '\'], '\correlations\', 'comparison', '\mean across mice\'};
    parameters.loop_list.things_to_save.node_averages.filename = {'average_by_nodes_randomPermutations_all.mat'};
    parameters.loop_list.things_to_save.node_averages.variable = {'average_by_nodes'};
    parameters.loop_list.things_to_save.node_averages.level = 'comparison';
    % Across mice.
    parameters.loop_list.things_to_save.average.dir =  {[parameters.dir_exper 'mean checks\beta changes\' comparison_type '\'], '\correlations\', 'comparison', '\mean across mice\'};
    parameters.loop_list.things_to_save.average.filename = {'average_by_nodes_randomPermutations_across_mice.mat'};
    parameters.loop_list.things_to_save.average.variable = {'average_by_nodes'};
    parameters.loop_list.things_to_save.average.level = 'comparison';
    
    parameters.loop_list.things_to_rename = {
                                             {'node_averages', 'data'}}; 
    
    RunAnalysis({@AverageByNode_Inverted_Modified, @AverageData}, parameters);

end
parameters.onPermutations = false;

%% Reshape continuous average by nodes correlations for significance &  dot plots
if isfield(parameters, 'loop_list')
parameters = rmfield(parameters,'loop_list');
end

% Iterators
parameters.loop_list.iterators = {
               'comparison', {'loop_variables.comparisons_', comparison_type, '(:).name'}, 'comparison_iterator';
               }; 
parameters.evaluation_instructions = {{  'data = parameters.data;'...
                                         'data_evaluated = reshape(data, size(data,1) * size(data,2), 1);'
                                        }};

% Inputs
parameters.loop_list.things_to_load.data.dir = {[parameters.dir_exper 'mean checks\beta changes\continuous\correlations\'], 'comparison', '\mean across mice\'};
parameters.loop_list.things_to_load.data.filename= {'average_by_nodes_across_mice.mat'};
parameters.loop_list.things_to_load.data.variable= {'average_by_nodes'}; 
parameters.loop_list.things_to_load.data.level = 'comparison';

% Outputs
parameters.loop_list.things_to_save.data_evaluated.dir = {[parameters.dir_exper 'mean checks\beta changes\continuous\correlations\'], 'comparison', '\mean across mice\'};
parameters.loop_list.things_to_save.data_evaluated.filename= {'average_by_nodes_across_mice_reshaped.mat'};
parameters.loop_list.things_to_save.data_evaluated.variable= {'average_by_nodes_reshaped'}; 
parameters.loop_list.things_to_save.data_evaluated.level = 'comparison';

RunAnalysis({@EvaluateOnData}, parameters)

%% Correlations average by nodes: check significance for dot plots
% Always clear loop list first. 

if isfield(parameters, 'loop_list')
parameters = rmfield(parameters,'loop_list');
end

for comparison_typei = 1:numel(parameters.loop_variables.comparison_types) 

    comparison_type = parameters.loop_variables.comparison_types{comparison_typei};
    parameters.shufflesDim = 2;

    % Iterators
    parameters.loop_list.iterators = {
                   'comparison', {'loop_variables.comparisons_', comparison_type, '(:).name'}, 'comparison_iterator';
                   }; 
   
    parameters.find_significance = true;
    
    % The statistical alpha value
    parameters.alphaValue = 0.05;  %0.05/496; %0.001; %/numel(parameters.comparisons_continuous);
    
    % Say that you do want to use bootstrapping.
    parameters.useBootstrapping = false;
    
    % If you want to fit a normal distribution before t-test (default = true)
    parameters.useNormalDistribution = true; 
    
    % Use false discovery rate correction?
    parameters.useFDR = true;
    % Inputs:
    %Test values (will grab only the intercepts with EvaluateOnData)
    parameters.loop_list.things_to_load.test_values.dir = {[parameters.dir_exper 'mean checks\beta changes\' comparison_type '\'],'correlations\', 'comparison', '\mean across mice\'};
    if strcmp(comparison_type, 'categorical')
    parameters.loop_list.things_to_load.test_values.filename= {'average_by_nodes_across_mice.mat'};
    parameters.loop_list.things_to_load.test_values.variable= {'average_by_nodes'}; 
    else
     parameters.loop_list.things_to_load.test_values.filename= {'average_by_nodes_across_mice_reshaped.mat'};
    parameters.loop_list.things_to_load.test_values.variable= {'average_by_nodes_reshaped'}; 
    end 
    parameters.loop_list.things_to_load.test_values.level = 'comparison';
    % Null distribution
    parameters.loop_list.things_to_load.null_distribution.dir = {[parameters.dir_exper 'mean checks\beta changes\' comparison_type '\'], 'correlations\', 'comparison', '\mean across mice\'};
    parameters.loop_list.things_to_load.null_distribution.filename= {'average_by_nodes_randomPermutations_across_mice.mat'};
    parameters.loop_list.things_to_load.null_distribution.variable= {'average_by_nodes'}; 
    parameters.loop_list.things_to_load.null_distribution.level = 'comparison';
    
    % Outputs
    parameters.loop_list.things_to_save.significance.dir = {[parameters.dir_exper 'mean checks\beta changes\' comparison_type '\'],'correlations\', 'comparison', '\mean across mice\'};
    parameters.loop_list.things_to_save.significance.filename= {'average_by_nodes_significance.mat'};
    parameters.loop_list.things_to_save.significance.variable= {'significance'}; 
    parameters.loop_list.things_to_save.significance.level = 'comparison';
    
    RunAnalysis({@SignificanceCalculation}, parameters);
end    
parameters.useFDR = false;

