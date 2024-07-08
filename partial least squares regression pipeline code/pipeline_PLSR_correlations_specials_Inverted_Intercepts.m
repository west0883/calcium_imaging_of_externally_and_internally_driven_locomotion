% pipeline_PLSR_correlations_specials_Inverted_Intercepts.m
% Sarah West
% 6/13/23

% Sets up variables for & runs partial least squares regression analysis
% for Random Motorized Treadmill experiements---correlations for non-warning 
% periods. 

% This is the edits made after peer review from Nature Communications. The
% original pipeline script was pipeline_PLSR_multilevel.m, which is now in
% the "old pipeline scripts_before peer review" subfolder.

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
    parameters.comparisons_categorical = comparisons; % don't run the motorized vs spontaneous rest periods;
    parameters.loop_variables.comparisons_categorical = parameters.comparisons_categorical;
    clear comparisons;
end

% Names of all continuous variables.
parameters.continuous_variable_names = {'speed', 'accel', 'duration', 'pupil_diameter', 'tail', 'nose', 'FL', 'HL', 'x'};

% Intercept values for each continuous variable
% (All are 0, except pupil diameter which is 0.7)
parameters.intercepts = [0 0 0 0.7 0 0 0 0 0];

% Put relevant variables into loop_variables.
parameters.loop_variables.mice_all = parameters.mice_all;
parameters.loop_variables.periods = parameters.periods.condition; 
parameters.loop_variables.conditions = {'motorized'; 'spontaneous'};
parameters.loop_variables.conditions_stack_locations = {'stacks'; 'spontaneous'};
parameters.loop_variables.variable_type = {'response variables', 'correlations'};
parameters.loop_variables.categories.type = parameters.categories.type;
parameters.loop_variables.comparison_types = {'categorical', 'continuous'};
parameters.loop_variables.output_types =  {'BETA'}; %{'Cov', % For continuous variables, I want to also see the actual regressors

parameters.average_and_std_together = false;

%% Put in response variables, no vertical concatenation

if isfield(parameters, 'loop_list')
parameters = rmfield(parameters,'loop_list');
end

% Iterators
parameters.loop_list.iterators = {'mouse', {'loop_variables.mice_all(:).name'}, 'mouse_iterator'};

% Variables to replicate
parameters.response_variable_names = {'motorized_vs_spon_dummyvars_vector', 'type_dummyvars_vector', 'transition_or_not_dummyvars_vector', 'speed_vector', 'accel_vector', 'duration_vector', 'pupil_diameter_vector', 'tail_vector', 'nose_vector', 'FL_vector', 'HL_vector', 'x_vector'};
parameters.variables_static = {'motorized_vs_spon_dummyvars_vector', 'type_dummyvars_vector', 'transition_or_not_dummyvars_vector', 'duration_vector'};
parameters.motorized_variables_static = {'speed_vector', 'accel_vector'}; % These are the ones that are static in motorized, not static in spontaneous
% Additional variables -- pupil, tail, nose, FL, HL; always present & loaded in
parameters.additional_variables = parameters.response_variable_names(7:end);
% Original order of spontaneous (for velocity & accel indexing)
parameters.spontaneous_periods_order = {'rest', 'walk', 'prewalk', 'startwalk', 'stopwalk', 'postwalk'};

parameters.concatenate_vertically = false;

% Input
% Correlations (for instances count)
parameters.loop_list.things_to_load.data.dir = {[parameters.dir_exper 'PLSR Specials\variable prep\correlations\'], 'mouse', '\'};
parameters.loop_list.things_to_load.data.filename= {'values.mat'};
parameters.loop_list.things_to_load.data.variable= {'values'}; 
parameters.loop_list.things_to_load.data.level = 'mouse';

% Spontaneous velocity
parameters.loop_list.things_to_load.speed_vector.dir = {[parameters.dir_exper 'behavior\spontaneous\rolled concatenated velocity\'], 'mouse', '\'};
parameters.loop_list.things_to_load.speed_vector.filename= {'velocity_averaged_by_instance.mat'};
parameters.loop_list.things_to_load.speed_vector.variable= {'velocity_averaged_by_instance'}; 
parameters.loop_list.things_to_load.speed_vector.level = 'mouse';

% Spontaneous accel.
parameters.loop_list.things_to_load.accel_vector.dir = {[parameters.dir_exper 'behavior\spontaneous\rolled concatenated velocity\'], 'mouse', '\'};
parameters.loop_list.things_to_load.accel_vector.filename= {'accel_averaged_by_instance.mat'};
parameters.loop_list.things_to_load.accel_vector.variable= {'accel_averaged_by_instance'}; 
parameters.loop_list.things_to_load.accel_vector.level = 'mouse';

% Pupil diameter
parameters.loop_list.things_to_load.pupil_diameter_vector.dir = {[parameters.dir_exper 'behavior\eye\rolled concatenated diameters\'], 'mouse', '\'};
parameters.loop_list.things_to_load.pupil_diameter_vector.filename= {'diameter_averaged_by_instance.mat'};
parameters.loop_list.things_to_load.pupil_diameter_vector.variable= {'diameter_averaged_by_instance'}; 
parameters.loop_list.things_to_load.pupil_diameter_vector.level = 'mouse';

% Tail 
parameters.loop_list.things_to_load.tail_vector.dir = {[parameters.dir_exper 'behavior\body\normalized with std only\value per roll velocity\tail\total_magnitude\'], 'mouse', '\'};
parameters.loop_list.things_to_load.tail_vector.filename= {'velocity_averaged_by_instance.mat'};
parameters.loop_list.things_to_load.tail_vector.variable= {'velocity_averaged_by_instance'}; 
parameters.loop_list.things_to_load.tail_vector.level = 'mouse';

% Nose 
parameters.loop_list.things_to_load.nose_vector.dir = {[parameters.dir_exper 'behavior\body\normalized with std only\value per roll velocity\nose\total_magnitude\'], 'mouse', '\'};
parameters.loop_list.things_to_load.nose_vector.filename= {'velocity_averaged_by_instance.mat'};
parameters.loop_list.things_to_load.nose_vector.variable= {'velocity_averaged_by_instance'}; 
parameters.loop_list.things_to_load.nose_vector.level = 'mouse';

% FL 
parameters.loop_list.things_to_load.FL_vector.dir = {[parameters.dir_exper 'behavior\body\normalized with std only\value per roll velocity\FL\total_magnitude\'], 'mouse', '\'};
parameters.loop_list.things_to_load.FL_vector.filename= {'velocity_averaged_by_instance.mat'};
parameters.loop_list.things_to_load.FL_vector.variable= {'velocity_averaged_by_instance'}; 
parameters.loop_list.things_to_load.FL_vector.level = 'mouse';

% HL 
parameters.loop_list.things_to_load.HL_vector.dir = {[parameters.dir_exper 'behavior\body\normalized with std only\value per roll velocity\HL\total_magnitude\'], 'mouse', '\'};
parameters.loop_list.things_to_load.HL_vector.filename= {'velocity_averaged_by_instance.mat'};
parameters.loop_list.things_to_load.HL_vector.variable= {'velocity_averaged_by_instance'}; 
parameters.loop_list.things_to_load.HL_vector.level = 'mouse';

% x 
parameters.loop_list.things_to_load.x_vector.dir = {[parameters.dir_exper 'behavior\body\normalized with std only\value per roll velocity\FL\x\'], 'mouse', '\'};
parameters.loop_list.things_to_load.x_vector.filename= {'velocity_averaged_by_instance.mat'};
parameters.loop_list.things_to_load.x_vector.variable= {'velocity_averaged_by_instance'}; 
parameters.loop_list.things_to_load.x_vector.level = 'mouse';

% rest & walk duration 
parameters.loop_list.things_to_load.duration_place.dir = {[parameters.dir_exper 'behavior\duration place concatenated\both conditions\'], 'mouse', '\'};
parameters.loop_list.things_to_load.duration_place.filename= {'all_duration_place.mat'};
parameters.loop_list.things_to_load.duration_place.variable= {'all_duration_place'}; 
parameters.loop_list.things_to_load.duration_place.level = 'mouse';

% Output 
parameters.loop_list.things_to_save.response_variables.dir = {[parameters.dir_exper 'PLSR Specials Inverted Intercepts\variable prep\response variables\'], 'mouse', '\'};
parameters.loop_list.things_to_save.response_variables.filename= {'response_variables_table.mat'};
parameters.loop_list.things_to_save.response_variables.variable= {'response_variables'}; 
parameters.loop_list.things_to_save.response_variables.level = 'mouse';

RunAnalysis({@PopulateResponseVariables}, parameters);

%% Prepare datasets per continuous comparison. 
if isfield(parameters, 'loop_list')
parameters = rmfield(parameters,'loop_list');
end

% Iterators
parameters.loop_list.iterators = {
               'mouse', {'loop_variables.mice_all(:).name'}, 'mouse_iterator'; 
               'comparison', {'loop_variables.comparisons_continuous(:).name'}, 'comparison_iterator'     
               };

% Specify which comparisons should be used for this dataset prep. 
parameters.this_comparison_set = parameters.comparisons_continuous;
parameters.comparison_type = 'continuous';

% Remove outliers from explanatory variables.
parameters.removeOutliers = true;

% Flag for whether or not missing data (NaNs) should be imputed.
parameters.imputeMissing = true; 

% Number of PLSR components that should be used for imputing missing data.
% Using just 85% instead of 90% usually cuts number of components needed by
% half.
parameters.imputation_components_variance_explained = 85; % in percents
parameters.imputation_max_components = 9; 

% Input 
parameters.loop_list.things_to_load.response.dir = {[parameters.dir_exper 'PLSR Specials Inverted Intercepts\variable prep\response variables\'], 'mouse', '\'};
parameters.loop_list.things_to_load.response.filename= {'response_variables_table.mat'};
parameters.loop_list.things_to_load.response.variable= {'response_variables'}; 
parameters.loop_list.things_to_load.response.level = 'mouse';

parameters.loop_list.things_to_load.explanatory.dir = {[parameters.dir_exper 'PLSR Specials\variable prep\correlations\'], 'mouse', '\'};
parameters.loop_list.things_to_load.explanatory.filename= {'values.mat'};
parameters.loop_list.things_to_load.explanatory.variable= {'values'}; 
parameters.loop_list.things_to_load.explanatory.level = 'mouse';

% Output
parameters.loop_list.things_to_save.dataset.dir = {[parameters.dir_exper 'PLSR Specials Inverted Intercepts\variable prep\datasets\level 1 continuous\'], 'comparison', '\' 'mouse', '\'};
parameters.loop_list.things_to_save.dataset.filename= {'PLSR_dataset_info.mat'};
parameters.loop_list.things_to_save.dataset.variable= {'dataset_info'}; 
parameters.loop_list.things_to_save.dataset.level = 'comparison';

RunAnalysis({@DatasetPrep}, parameters);

parameters.removeOutliers = false;
parameters.imputeMissing = false;

%% Find the average ratio of NaNs in pupil diameter, tail, nose, FL, HL, x, across continuous comparisons per mouse.
if isfield(parameters, 'loop_list')
parameters = rmfield(parameters,'loop_list');
end

% Is so you can use a single loop for calculations. 
parameters.loop_list.iterators = {
                'mouse', {'loop_variables.mice_all(:).name'}, 'mouse_iterator'; 
                'comparison', {'loop_variables.comparisons_continuous(:).name'}, 'comparison_iterator' };

parameters.evaluation_instructions = {{'data_evaluated = parameters.data.responseVariables(end-5:end);'}}; % DLC data is always the last 5 variables entered.
parameters.concatDim = 1;
parameters.concatenation_level = 'comparison';
parameters.averageDim = 1;
parameters.average_and_std_together = true;

% Input Values 
parameters.loop_list.things_to_load.data.dir = {[parameters.dir_exper 'PLSR Specials Inverted Intercepts\variable prep\datasets\level 1 continuous\'], 'comparison', '\' 'mouse', '\'};
parameters.loop_list.things_to_load.data.filename= {'PLSR_dataset_info.mat'};
parameters.loop_list.things_to_load.data.variable= {'dataset_info.NaN_ratios'}; 
parameters.loop_list.things_to_load.data.level = 'comparison';

% Output values
parameters.loop_list.things_to_save.average.dir = {[parameters.dir_exper 'PLSR Specials Inverted Intercepts\variable prep\datasets\level 1 continuous\NaN ratios\'], 'mouse', '\'};
parameters.loop_list.things_to_save.average.filename= {'average_missing_pupil_data_ratios.mat'};
parameters.loop_list.things_to_save.average.variable= {'average_missing_pupil_data_ratios'}; 
parameters.loop_list.things_to_save.average.level = 'mouse';

parameters.loop_list.things_to_rename = {{'data_evaluated', 'data'}; 
                                         { 'concatenated_data', 'data'}};

RunAnalysis({@EvaluateOnData,@ConcatenateData, @AverageData}, parameters);

%% PLSR Level 1, continuous: optimize components.
% Don't run any permutations yet.
% Always clear loop list first. 
if isfield(parameters, 'loop_list')
parameters = rmfield(parameters,'loop_list');
end

% Iterators
parameters.loop_list.iterators = {
               'mouse', {'loop_variables.mice_all(:).name'}, 'mouse_iterator'; 
               'comparison', {'loop_variables.comparisons_continuous(:).name'}, 'comparison_iterator' };

% Parameters for calculating best number of components. If
% "findBestNComponents" = false, just run the ncomponents_max
parameters.findBestNComponents = true;
parameters.ncomponents_max = 9; 
parameters.contiguous_partitions = true; 
parameters.kFolds = 10;
parameters.MonteCarloReps = 10;
parameters.comparison_type = 'continuous';
parameters.stratify = false;

% Do you want permutations?
parameters.permutationGeneration = false;

% Input 
parameters.loop_list.things_to_load.dataset.dir = {[parameters.dir_exper 'PLSR Specials Inverted Intercepts\variable prep\datasets\level 1 continuous\'], 'comparison', '\' 'mouse', '\'};
parameters.loop_list.things_to_load.dataset.filename= {'PLSR_dataset_info.mat'};
parameters.loop_list.things_to_load.dataset.variable= {'dataset_info'}; 
parameters.loop_list.things_to_load.dataset.level = 'comparison';

% Output
parameters.loop_list.things_to_save.results.dir = {[parameters.dir_exper 'PLSR Specials Inverted Intercepts\results\level 1 continuous\'], 'comparison','\', 'mouse', '\'};
parameters.loop_list.things_to_save.results.filename= {'PLSR_results.mat'};
parameters.loop_list.things_to_save.results.variable= {'PLSR_results'}; 
parameters.loop_list.things_to_save.results.level = 'comparison';

RunAnalysis({@PLSR_forRunAnalysis_Inverted}, parameters);  

parameters.findBestNComponents = false;


%% Remove continuous variables effects from each behavior type. 
% only fill indices that you used in the given continuous comparison (leave
% other cells blank) 

% Always clear loop list first. 
if isfield(parameters, 'loop_list')
parameters = rmfield(parameters,'loop_list');
end

% Iterators
parameters.loop_list.iterators = {
               'mouse', {'loop_variables.mice_all(:).name'}, 'mouse_iterator'; 
               'comparison', {'loop_variables.comparisons_continuous(:).name'}, 'comparison_iterator' };

parameters.removeOutliers = false; 
parameters.imputeMissing = false; 

parameters.useIntercepts = true;

% Input 
% The variables from the comparison
parameters.loop_list.things_to_load.dataset.dir = {[parameters.dir_exper 'PLSR Specials Inverted Intercepts\variable prep\datasets\level 1 continuous\'], 'comparison', '\' 'mouse', '\'};
parameters.loop_list.things_to_load.dataset.filename= {'PLSR_dataset_info.mat'};
parameters.loop_list.things_to_load.dataset.variable= {'dataset_info'}; 
parameters.loop_list.things_to_load.dataset.level = 'comparison';
% The results from the continuous regression (for the Betas)
parameters.loop_list.things_to_load.PLSR_results.dir = {[parameters.dir_exper 'PLSR Specials Inverted Intercepts\results\level 1 continuous\'], 'comparison', '\' 'mouse', '\'};
parameters.loop_list.things_to_load.PLSR_results.filename= {'PLSR_results.mat'};
parameters.loop_list.things_to_load.PLSR_results.variable= {'PLSR_results'}; 
parameters.loop_list.things_to_load.PLSR_results.level = 'comparison';
% Old correlation values (for the number of instances)
parameters.loop_list.things_to_load.values_old.dir = {[parameters.dir_exper 'PLSR Specials\variable prep\correlations\'], 'mouse', '\'};
parameters.loop_list.things_to_load.values_old.filename= {'values.mat'};
parameters.loop_list.things_to_load.values_old.variable= {'values'}; 
parameters.loop_list.things_to_load.values_old.level = 'mouse';

% Output
parameters.loop_list.things_to_save.values_new.dir = {[parameters.dir_exper 'PLSR Specials Inverted Intercepts\variable prep\correlations\'], 'mouse', '\'};
parameters.loop_list.things_to_save.values_new.filename= {'correlations_continuousSubtracted.mat'};
parameters.loop_list.things_to_save.values_new.variable= {'correlations'}; 
parameters.loop_list.things_to_save.values_new.level = 'mouse';

RunAnalysis({@ResidualsFromContinuous_Inverted}, parameters); 



%% Level 1 categorical -- Prepare datasets, continuous subtracted.
if isfield(parameters, 'loop_list')
parameters = rmfield(parameters,'loop_list');
end

% Iterators
parameters.loop_list.iterators = {
               'mouse', {'loop_variables.mice_all(:).name'}, 'mouse_iterator'; 
               'comparison', {'loop_variables.comparisons_categorical(:).name'}, 'comparison_iterator'     
               };

% Specify which comparisons should be used for this dataset prep. 
parameters.this_comparison_set = parameters.comparisons_categorical;
parameters.comparison_type = 'categorical';

% Flag for whether or not missing data (NaNs) should be imputed. (Don't
% need it for these comparisons, already did it at residual level in previous step.)
parameters.removeOutliers = false;
parameters.imputeMissing = false; 
parameters.imputation_components_variance_explained = 85; % in percents
parameters.imputation_max_components = 10; 

% Input 
% Don't need any outliers-removed responses with categorical
parameters.loop_list.things_to_load.response.dir = {[parameters.dir_exper 'PLSR Specials Inverted Intercepts\variable prep\response variables\'], 'mouse', '\'};
parameters.loop_list.things_to_load.response.filename= {'response_variables_table.mat'};
parameters.loop_list.things_to_load.response.variable= {'response_variables'}; 
parameters.loop_list.things_to_load.response.level = 'mouse';

parameters.loop_list.things_to_load.explanatory.dir = {[parameters.dir_exper 'PLSR Specials Inverted Intercepts\variable prep\correlations\'], 'mouse', '\'};
parameters.loop_list.things_to_load.explanatory.filename= {'correlations_continuousSubtracted.mat'};
parameters.loop_list.things_to_load.explanatory.variable= {'correlations'}; 
parameters.loop_list.things_to_load.explanatory.level = 'mouse';

% Output
parameters.loop_list.things_to_save.dataset.dir = {[parameters.dir_exper 'PLSR Specials Inverted Intercepts\variable prep\datasets\level 1 categorical\'], 'comparison', '\' 'mouse', '\'};
parameters.loop_list.things_to_save.dataset.filename= {'PLSR_dataset_info.mat'};
parameters.loop_list.things_to_save.dataset.variable= {'dataset_info'}; 
parameters.loop_list.things_to_save.dataset.level = 'comparison';

RunAnalysis({@DatasetPrep}, parameters);

%% Level 1 categorical -- optimize number of components
% Will look at the outputs from 10 calculated components.

% Always clear loop list first. 
if isfield(parameters, 'loop_list')
parameters = rmfield(parameters,'loop_list');
end

% Iterators
parameters.loop_list.iterators = {
               'mouse', {'loop_variables.mice_all(:).name'}, 'mouse_iterator'; 
               'comparison', {'loop_variables.comparisons_categorical(:).name'}, 'comparison_iterator' };

% Parameters for calculating best number of components. If
% "findBestNComponents" = false, just run the ncomponents_max
parameters.findBestNComponents = true;
parameters.ncomponents_max = 2;  
parameters.contiguous_partitions = true; 
parameters.kFolds = 10;
parameters.MonteCarloReps = 10;
parameters.comparison_type = 'categorical';
parameters.stratify = true;

% Do you want permutations?
parameters.permutationGeneration = false;

% Input 
parameters.loop_list.things_to_load.dataset.dir = {[parameters.dir_exper 'PLSR Specials Inverted Intercepts\variable prep\datasets\level 1 categorical\'], 'comparison', '\' 'mouse', '\'};
parameters.loop_list.things_to_load.dataset.filename= {'PLSR_dataset_info.mat'};
parameters.loop_list.things_to_load.dataset.variable= {'dataset_info'}; 
parameters.loop_list.things_to_load.dataset.level = 'comparison';

% Output
parameters.loop_list.things_to_save.results.dir = {[parameters.dir_exper 'PLSR Specials Inverted Intercepts\results\level 1 categorical\'], 'comparison', '\' 'mouse', '\'};
parameters.loop_list.things_to_save.results.filename= {'PLSR_results.mat'};
parameters.loop_list.things_to_save.results.variable= {'PLSR_results'}; 
parameters.loop_list.things_to_save.results.level = 'comparison';

RunAnalysis({@PLSR_forRunAnalysis_Inverted}, parameters);  

parameters.findBestNComponents = false;





%% SIGNIFICANCE STUFF 

%% Level 1 continuous -- run random permutations.
% With best number of components.
% Always clear loop list first. 
do = true;
if do 
    if isfield(parameters, 'loop_list')
    parameters = rmfield(parameters,'loop_list');
    end
    
    % Iterators
    parameters.loop_list.iterators = {
                   'mouse', {'loop_variables.mice_all(:).name'}, 'mouse_iterator'; 
                   'comparison', {'loop_variables.comparisons_continuous(:).name'}, 'comparison_iterator' };
    
    % Do you want permutations?
    parameters.permutationGeneration = true;
    parameters.useBootstrapping = false;
    parameters.n_permutations = 1000;
    parameters.stratify = false;
    parameters.comparison_type = 'continuous';
    
    % Input 
    % dataset
    parameters.loop_list.things_to_load.dataset.dir = {[parameters.dir_exper 'PLSR Specials Inverted Intercepts\variable prep\datasets\level 1 continuous\'], 'comparison', '\' 'mouse', '\'};
    parameters.loop_list.things_to_load.dataset.filename= {'PLSR_dataset_info.mat'};
    parameters.loop_list.things_to_load.dataset.variable= {'dataset_info'}; 
    parameters.loop_list.things_to_load.dataset.level = 'comparison';
    % optimized number of components to use.
    parameters.loop_list.things_to_load.ncomponents_max.dir = {[parameters.dir_exper 'PLSR Specials Inverted Intercepts\results\level 1 continuous\'], 'comparison', '\', 'mouse', '\'};
    parameters.loop_list.things_to_load.ncomponents_max.filename= {'PLSR_results.mat'};
    parameters.loop_list.things_to_load.ncomponents_max.variable= {'PLSR_results.ncomponents_used'}; 
    parameters.loop_list.things_to_load.ncomponents_max.level = 'comparison';

    % Output
    % Covs
    parameters.loop_list.things_to_save.Covs_randomPermutations.dir = {[parameters.dir_exper 'PLSR Specials Inverted Intercepts\results\level 1 continuous\'], 'comparison', '\' 'mouse', '\'};
    parameters.loop_list.things_to_save.Covs_randomPermutations.filename= {'PLSR_Covs_randomPermutations.mat'};
    parameters.loop_list.things_to_save.Covs_randomPermutations.variable= {'Covs_randomPermutations'}; 
    parameters.loop_list.things_to_save.Covs_randomPermutations.level = 'comparison';
    % Betas
    parameters.loop_list.things_to_save.BETAs_randomPermutations.dir = {[parameters.dir_exper 'PLSR Specials Inverted Intercepts\results\level 1 continuous\'], 'comparison', '\' 'mouse', '\'};
    parameters.loop_list.things_to_save.BETAs_randomPermutations.filename= {'PLSR_BETAs_randomPermutations.mat'};
    parameters.loop_list.things_to_save.BETAs_randomPermutations.variable= {'BETAs_randomPermutations'}; 
    parameters.loop_list.things_to_save.BETAs_randomPermutations.level = 'comparison';

    RunAnalysis({@PLSR_forRunAnalysis_Inverted}, parameters);  
    
    parameters.permutationGeneration = false;
end 


%% Level 1 categorical -- run random permutations.
do = true;
if do 
    % Always clear loop list first. 
    if isfield(parameters, 'loop_list')
    parameters = rmfield(parameters,'loop_list');
    end
    
    % Iterators
    parameters.loop_list.iterators = {
                   'mouse', {'loop_variables.mice_all(:).name'}, 'mouse_iterator'; 
                   'comparison', {'loop_variables.comparisons_categorical(:).name'}, 'comparison_iterator' }; 
    
    % Do you want permutations?
    parameters.permutationGeneration = true;
    parameters.useBootstrapping = false;
    parameters.n_permutations = 1000;
    parameters.stratify = true;
    parameters.comparison_type = 'categorical';
    
    % Input 
    parameters.loop_list.things_to_load.dataset.dir = {[parameters.dir_exper 'PLSR Specials Inverted Intercepts\variable prep\datasets\level 1 categorical\'], 'comparison', '\' 'mouse', '\'};
    parameters.loop_list.things_to_load.dataset.filename= {'PLSR_dataset_info.mat'};
    parameters.loop_list.things_to_load.dataset.variable= {'dataset_info'}; 
    parameters.loop_list.things_to_load.dataset.level = 'comparison';
    % optimized number of components to use.
    parameters.loop_list.things_to_load.ncomponents_max.dir = {[parameters.dir_exper 'PLSR Specials Inverted Intercepts\results\level 1 categorical\'], 'comparison', '\', 'mouse', '\'};
    parameters.loop_list.things_to_load.ncomponents_max.filename= {'PLSR_results.mat'};
    parameters.loop_list.things_to_load.ncomponents_max.variable= {'PLSR_results.ncomponents_used'}; 
    parameters.loop_list.things_to_load.ncomponents_max.level = 'comparison';
    
    % Output
    % Covs
    parameters.loop_list.things_to_save.Covs_randomPermutations.dir = {[parameters.dir_exper 'PLSR Specials Inverted Intercepts\results\level 1 categorical\'], 'comparison', '\' 'mouse', '\'};
    parameters.loop_list.things_to_save.Covs_randomPermutations.filename= {'PLSR_Covs_randomPermutations.mat'};
    parameters.loop_list.things_to_save.Covs_randomPermutations.variable= {'Covs_randomPermutations'}; 
    parameters.loop_list.things_to_save.Covs_randomPermutations.level = 'comparison';
    % Betas
    parameters.loop_list.things_to_save.BETAs_randomPermutations.dir = {[parameters.dir_exper 'PLSR Specials Inverted Intercepts\results\level 1 categorical\'], 'comparison', '\' 'mouse', '\'};
    parameters.loop_list.things_to_save.BETAs_randomPermutations.filename= {'PLSR_BETAs_randomPermutations.mat'};
    parameters.loop_list.things_to_save.BETAs_randomPermutations.variable= {'BETAs_randomPermutations'}; 
    parameters.loop_list.things_to_save.BETAs_randomPermutations.level = 'comparison';
    
    RunAnalysis({@PLSR_forRunAnalysis_Inverted}, parameters);  
    
    parameters.permutationGeneration = false;
end



%% Brain: Average % variance of correlations explained per comparison
if isfield(parameters, 'loop_list')
parameters = rmfield(parameters,'loop_list');
end
        
% Iterators
parameters.loop_list.iterators = {
               'comparison_type', {'loop_variables.comparison_types'}, 'comparison_type_iterator';
               'comparison', {'loop_variables.comparisons_', 'comparison_type', '(:).name'}, 'comparison_iterator';
                % Don't use m1100
               'mouse', {'loop_variables.mice_all([1:3 5:7]).name'}, 'mouse_iterator'; 
               }; 

parameters.evaluation_instructions = {{'data = parameters.data;'...
                                        'data_evaluated = sum(data(2,:));'
                                      }};
parameters.concatDim = 1; 
parameters.averageDim = 1;

% Inputs
% PLSR percent variance
parameters.loop_list.things_to_load.data.dir = {[parameters.dir_exper 'PLSR Specials Inverted Intercepts\results\level 1 '], 'comparison_type', '\', 'comparison','\', 'mouse', '\'};
parameters.loop_list.things_to_load.data.filename= {'PLSR_results.mat'};
parameters.loop_list.things_to_load.data.variable= {'PLSR_results.PCTVAR'}; 
parameters.loop_list.things_to_load.data.level = 'mouse';

% Outputs
% concatenated data
parameters.loop_list.things_to_save.concatenated_data.dir = {[parameters.dir_exper 'PLSR Specials Inverted Intercepts\results\total percent variance brain\'],  'comparison_type', '\' 'comparison','\'};
parameters.loop_list.things_to_save.concatenated_data.filename= {'all_percent_var.mat'};
parameters.loop_list.things_to_save.concatenated_data.variable= {'all_percent_var'}; 
parameters.loop_list.things_to_save.concatenated_data.level = 'comparison';
% average
parameters.loop_list.things_to_save.average.dir = {[parameters.dir_exper 'PLSR Specials Inverted Intercepts\results\total percent variance brain\'],  'comparison_type', '\' 'comparison','\'};
parameters.loop_list.things_to_save.average.filename= {'average_percent_var.mat'};
parameters.loop_list.things_to_save.average.variable= {'average_percent_var'}; 
parameters.loop_list.things_to_save.average.level = 'comparison';
% std_dev
parameters.loop_list.things_to_save.std_dev.dir = {[parameters.dir_exper 'PLSR Specials Inverted Intercepts\results\total percent variance brain\'],  'comparison_type', '\' 'comparison','\'};
parameters.loop_list.things_to_save.std_dev.filename= {'std_dev_percent_var.mat'};
parameters.loop_list.things_to_save.std_dev.variable= {'std_dev_percent_var'}; 
parameters.loop_list.things_to_save.std_dev.level = 'comparison';

parameters.loop_list.things_to_rename = {{'data_evaluated', 'data'};
                                         {'concatenated_data','data'}
                                        };
RunAnalysis({@EvaluateOnData, @ConcatenateData, @AverageData}, parameters);


%% Brain: Get % variance of total for categorical 
if isfield(parameters, 'loop_list')
parameters = rmfield(parameters,'loop_list');
end
parameters.loop_variables.sub_comparisons = {'rest', 'start', 'walk', 'stop'};
% Iterators
parameters.loop_list.iterators = {
               'comparison', {'loop_variables.sub_comparisons'}, 'comparison_iterator';         
               }; 

parameters.evaluation_instructions = {{'continuous = parameters.continuous_data;'...
                                        'categorical = parameters.categorical_data;'...
                                        'data_evaluated = (1 - continuous) .* categorical .* 100;'
                                      }};

% Inputs
% continuous data 
parameters.loop_list.things_to_load.continuous_data.dir = {[parameters.dir_exper 'PLSR Specials Inverted Intercepts\results\total percent variance brain\continuous\all_'], 'comparison','_continuousVars', '\'};
parameters.loop_list.things_to_load.continuous_data.filename= {'all_percent_var.mat'};
parameters.loop_list.things_to_load.continuous_data.variable= {'all_percent_var'}; 
parameters.loop_list.things_to_load.continuous_data.level = 'comparison';
% categorical data
parameters.loop_list.things_to_load.categorical_data.dir = {[parameters.dir_exper 'PLSR Specials Inverted Intercepts\results\total percent variance brain\categorical\'], 'comparison','_motorizedvsspon_categorical', '\'};
parameters.loop_list.things_to_load.categorical_data.filename= {'all_percent_var.mat'};
parameters.loop_list.things_to_load.categorical_data.variable= {'all_percent_var'}; 
parameters.loop_list.things_to_load.categorical_data.level = 'comparison';
% Outputs
parameters.loop_list.things_to_save.data_evaluated.dir = {[parameters.dir_exper 'PLSR Specials Inverted Intercepts\results\total percent variance brain\categorical\'], 'comparison','_motorizedvsspon_categorical', '\'};
parameters.loop_list.things_to_save.data_evaluated.filename= {'all_percent_var_ofTotal.mat'};
parameters.loop_list.things_to_save.data_evaluated.variable= {'all_percent_var'}; 
parameters.loop_list.things_to_save.data_evaluated.level = 'comparison';

RunAnalysis({@EvaluateOnData}, parameters);

%% Behavior: Average % variance of fluorescence explained per comparison
if isfield(parameters, 'loop_list')
parameters = rmfield(parameters,'loop_list');
end
        
% Iterators
parameters.loop_list.iterators = {
               'comparison_type', {'loop_variables.comparison_types'}, 'comparison_type_iterator';
               'comparison', {'loop_variables.comparisons_', 'comparison_type', '(:).name'}, 'comparison_iterator';
                % Don't use m1100
               'mouse', {'loop_variables.mice_all([1:3 5:7]).name'}, 'mouse_iterator'; 
               }; 

parameters.evaluation_instructions = {{'data = parameters.data;'...
                                        'data_evaluated = sum(data(1,:));'
                                      }};
parameters.concatDim = 1; 
parameters.averageDim = 1;

% Inputs
% PLSR percent variance
parameters.loop_list.things_to_load.data.dir = {[parameters.dir_exper 'PLSR Specials Inverted Intercepts\results\level 1 '], 'comparison_type', '\', 'comparison','\', 'mouse', '\'};
parameters.loop_list.things_to_load.data.filename= {'PLSR_results.mat'};
parameters.loop_list.things_to_load.data.variable= {'PLSR_results.PCTVAR'}; 
parameters.loop_list.things_to_load.data.level = 'mouse';

% Outputs
% concatenated data
parameters.loop_list.things_to_save.concatenated_data.dir = {[parameters.dir_exper 'PLSR Specials Inverted Intercepts\results\total percent variance behavior\'],  'comparison_type', '\' 'comparison','\'};
parameters.loop_list.things_to_save.concatenated_data.filename= {'all_percent_var.mat'};
parameters.loop_list.things_to_save.concatenated_data.variable= {'all_percent_var'}; 
parameters.loop_list.things_to_save.concatenated_data.level = 'comparison';
% average
parameters.loop_list.things_to_save.average.dir = {[parameters.dir_exper 'PLSR Specials Inverted Intercepts\results\total percent variance behavior\'],  'comparison_type', '\' 'comparison','\'};
parameters.loop_list.things_to_save.average.filename= {'average_percent_var.mat'};
parameters.loop_list.things_to_save.average.variable= {'average_percent_var'}; 
parameters.loop_list.things_to_save.average.level = 'comparison';
% std_dev
parameters.loop_list.things_to_save.std_dev.dir = {[parameters.dir_exper 'PLSR Specials Inverted Intercepts\results\total percent variance behavior\'],  'comparison_type', '\' 'comparison','\'};
parameters.loop_list.things_to_save.std_dev.filename= {'std_dev_percent_var.mat'};
parameters.loop_list.things_to_save.std_dev.variable= {'std_dev_percent_var'}; 
parameters.loop_list.things_to_save.std_dev.level = 'comparison';

parameters.loop_list.things_to_rename = {{'data_evaluated', 'data'};
                                         {'concatenated_data','data'}
                                        };
RunAnalysis({@EvaluateOnData, @ConcatenateData, @AverageData}, parameters);
