%% pipeline_ICs_MotorizedTreadmill.m
% Sarah West
% 9/1/21

% After preprocessing data and SVD compressing all data from a mouse into one compression
% with the supercomputers, use this code to calculate and clean the ICs of each animal. 

% **Use and run create_days_all.m before using this***

%% Initial Setup  

clear all; 

% ***********************************
% Directories

% Create the experiment name. This is used to name the output folder. 
parameters.experiment_name='Random Motorized Treadmill';

% Create the input directory of the SVD compressed datasets for each mouse
parameters.dir_dataset=['Y:\Sarah\Analysis\Experiments\' parameters.experiment_name '\spatial segmentation\500 SVD components\SVD compressions\'];

% Establish the format of the file names of compressed data. Each piece
% needs to be a separate entry in a cell array. Put the string 'mouse', 'day',
% or 'stack number' where the mouse, day, or stack number will be. If you 
% concatenated this as a sigle string, it should create a file name, with the 
% correct mouse/day/stack name inserted accordingly. 
parameters.compressed_data_name={parameters.dir_dataset, 'm', 'mouse number', '_SVD_compressed.mat'}; 

% Output directory name bases
parameters.dir_base='Y:\Sarah\Analysis\Experiments\';
parameters.dir_exper=[parameters.dir_base parameters.experiment_name '\']; 

% Was the data masked in preprocessing? (Masking that removed pixels, so
% masks have to be loaded in to reshape the data into proper images).
% If it was masked, the flag = 1; If not masked, the flag = 0.
parameters.masked_flag=1; 

% Directory of where the masks are saved. If not masked, can leave this
% empty. If masked, use the cell array format as above. 
parameters.dir_input_mask=[parameters.dir_exper 'preprocessing\masks\'];
parameters.mask_filename={'masks_m', 'mouse number', '.mat' };
parameters.mask_variable = {'indices_of_mask'};

% Compressed component that represents [spatial dimension]. Will be 'S' or
% 'V', with the number of (masked) pixels as one of the dimensions. Put in 
% as a character with quotes.
parameters.spatial_component='V'; 

% (DON'T EDIT). Load the "mice_all" variable you've created with "create_mice_all.m"
load([parameters.dir_exper 'mice_all.mat']);

% ****Change here if there are specific mice, days, and/or stacks you want to work with**** 
parameters.mice_all=mice_all;

parameters.mice_all=parameters.mice_all;

% ****************************************
% ***Parameters.*** 

% Set image dimensions.
parameters.yDim=256;
parameters.xDim=256;

% The number of ICs you want to calculate. 
parameters.num_sources=100;

% Number of digits you want in your stack number names (default is 3).
parameters.digitNumber = 2; 

% Use split domain ICs?
parameters.splitDomains = true;

% Always clear loop list first. 
if isfield(parameters, 'loop_list')
parameters = rmfield(parameters,'loop_list');
end


%% Calculate ICs
% Calculates ICs from SVD compressed data. Assumes one compressed dataset
% per mouse.

% output directory & filename
parameters.dir_out = {[parameters.dir_exper 'spatial segmentation\500 SVD components\raw ICs\'], 'mouse number', '\'};
parameters.output_filename = {['sources' num2str(parameters.num_sources) '.mat']};

% Use a gpu for this calculation? (t/f)
parameters.use_gpu = true;

% (DON'T EDIT). Run code. 
calculate_ICs(parameters); 


%% Plot raw ICs

% Input directory 
parameters.dir_input_base = {[parameters.dir_exper 'spatial segmentation\500 SVD components\raw ICs\'], 'mouse number', '\'};
parameters.input_filename = {['sources' num2str(parameters.num_sources) '.mat']};
parameters.input_variable = {'sources'};

% output directory
parameters.dir_output_base = {[parameters.dir_exper 'spatial segmentation\500 SVD components\raw ICs\'], 'mouse number', '\'};
parameters.output_filename = {['sources' num2str(parameters.num_sources)]};

% Run code
plot_rawICs(parameters); 

%% Regularize ICs 
% For cleaning the ICs
% Applies a (raw) threshold to the ICs.

parameters.amplitude_threshold = 3.5; %2.0; % 3.5

% Minimim size in pixels of an IC.
parameters.minPixels = 150; % 150

% Indicate if you want to z-score your ICs before regularizing (true/false)
%parameters.zscore_flag = false;

parameters.large_component_conditional_zscore_flag = true;
parameters.maxPixels = 5000; 
parameters.large_component_conditional_zscore_thresh = 1; %2; % 1

parameters.small_component_conditional_zscore_flag = true;
parameters.small_component_conditional_zscore_thresh =2.5; %1; % 2.5

% Input directory 
parameters.dir_input_base = {[parameters.dir_exper 'spatial segmentation\500 SVD components\raw ICs\'], 'mouse number', '\'};
parameters.input_filename = {['sources' num2str(parameters.num_sources) '.mat']};
parameters.input_variable = {'sources'};

% output directory
parameters.dir_output_base = {[parameters.dir_exper 'spatial segmentation\500 SVD components\regularized ICs ' ... 
                             num2str(parameters.minPixels) ' amp ' num2str(parameters.amplitude_threshold) ' two conditionals small ' num2str(parameters.small_component_conditional_zscore_thresh) ...
                              ' large ' num2str(num2str(parameters.large_component_conditional_zscore_thresh)) '\'], 'mouse number', '\'};
parameters.output_filename = {['sources' num2str(parameters.num_sources) '.mat']};
parameters.output_variable = {'sources'};

% (DON'T EDIT). Run code. 
regularize_ICs(parameters);


%% ***PRE ADDBACK***: Remove IC artifacts (interactive)
% 
% parameters.amplitude_threshold = 3.5; % 3.5
% parameters.minPixels = 150; % 150
% parameters.large_component_conditional_zscore_flag = true;
% parameters.maxPixels = 5000; 
% parameters.large_component_conditional_zscore_thresh = 1;
% parameters.small_component_conditional_zscore_flag = true;
% parameters.small_component_conditional_zscore_thresh = 2.5;
% 
% % Use second (right) monitor for figures, if available?
% parameters.second_monitor = false;
% 
% % Always clear loop list first. 
% if isfield(parameters, 'loop_list')
% parameters = rmfield(parameters,'loop_list');
% end
% 
% % Dimension different sources are in.
% parameters.sourcesDim = 3;
% 
% % Dimension the pixels dimension, different sources dimension of original data.
% parameters.originalSourcesPixelsDim = 2; 
% parameters.originalSourcesDim = 1; 
% 
% % Loop variables
% parameters.loop_list.iterators = {'mouse', {'loop_variables.mice_all(:).name'}, 'mouse_iterator';
%                                   'source', {'1:70'}, 'source_iterator'};
% parameters.loop_variables.mice_all = parameters.mice_all;
% 
% parameters.remove_entire_sources = true;
% 
% 
% % Input values
% parameters.loop_list.things_to_load.sources.dir = {[parameters.dir_exper 'spatial segmentation\500 SVD components\regularized ICs ' ... 
%                              num2str(parameters.minPixels) ' amp ' num2str(parameters.amplitude_threshold) ' two conditionals small ' num2str(parameters.small_component_conditional_zscore_thresh) ...
%                               ' large ' num2str(num2str(parameters.large_component_conditional_zscore_thresh)) '\'], 'mouse', '\'};
% parameters.loop_list.things_to_load.sources.filename= {['sources' num2str(parameters.num_sources) '.mat']};
% parameters.loop_list.things_to_load.sources.variable= {'sources'};
% parameters.loop_list.things_to_load.sources.level = 'mouse';
% 
% parameters.loop_list.things_to_load.indices_of_mask.dir = {[parameters.dir_exper 'preprocessing\masks\']};
% parameters.loop_list.things_to_load.indices_of_mask.filename= {'masks_m', 'mouse', '.mat'};
% parameters.loop_list.things_to_load.indices_of_mask.variable= {'indices_of_mask'}; 
% parameters.loop_list.things_to_load.indices_of_mask.level = 'mouse';
% 
% parameters.loop_list.things_to_load.reference_image.dir = {[parameters.dir_exper 'preprocessing\representative images\'], 'mouse', '\'};
% parameters.loop_list.things_to_load.reference_image.filename= {'reference_image.mat'};
% parameters.loop_list.things_to_load.reference_image.variable= {'reference_image'};
% parameters.loop_list.things_to_load.reference_image.level = 'mouse';
% 
% % [Right now, code assumes raw sources are in same file, I think it still
% % works if it's each source separate.]
% parameters.loop_list.things_to_load.original_sources.dir = {[parameters.dir_exper 'spatial segmentation\500 SVD components\raw ICs\'], 'mouse', '\'};
% parameters.loop_list.things_to_load.original_sources.filename= {['sources' num2str(parameters.num_sources) '.mat']};
% parameters.loop_list.things_to_load.original_sources.variable= {'sources'};
% parameters.loop_list.things_to_load.original_sources.level = 'mouse';
% 
% % (for any existing artifact removals)
% parameters.loop_list.things_to_load.sources_artifacts_removed_old.dir = {[parameters.dir_exper 'spatial segmentation\500 SVD components\artifacts removed conditional thresholding\pre addback\'], 'mouse', '\'};
% parameters.loop_list.things_to_load.sources_artifacts_removed_old.filename = {'sources.mat'};
% parameters.loop_list.things_to_load.sources_artifacts_removed_old.variable= {'sources'};
% parameters.loop_list.things_to_load.sources_artifacts_removed_old.level = 'mouse';
% 
% % Output values
% parameters.loop_list.things_to_save.sources_artifacts_removed.dir = {[parameters.dir_exper 'spatial segmentation\500 SVD components\artifacts removed conditional thresholding\pre addback\'], 'mouse', '\'};
% parameters.loop_list.things_to_save.sources_artifacts_removed.filename = {'sources.mat'};
% parameters.loop_list.things_to_save.sources_artifacts_removed.variable= {'sources'};
% parameters.loop_list.things_to_save.sources_artifacts_removed.level = 'source';
% 
% RunAnalysis({@RemoveArtifacts}, parameters);

%% Plot IC overlays together

% figure; 
% for i = 1:size(parameters.mice_all(:))
%     mouse = parameters.mice_all(i).name;
%     load(['Y:\Sarah\Analysis\Experiments\Random Motorized Treadmill\spatial segmentation\500 SVD components\artifacts removed conditional thresholding\' mouse '\sources.mat']);
%     subplot(2, 3, i); imagesc(sources.overlay); axis square; title(mouse); colorbar;
%     xticks([]); yticks([]);
% end 

%% Add back some ICs 
% *****RUN ONLY ONCE **** 
% raw IC segments that you think are real but didn't survive the
% regularization process. Wil run as a script, saves to end of list of
% regularizing ICs so it doesn't disrupt any previously saved artifact
% removals.
% Puts ouput in "post addback staging" folder to keep you from overwriting
% any artifact removal you want to do to the new sources. Have to copy the
% output into a new folder manually to call it in the next steps.

%add_back_ICs_RandomMotorizedTreadmill.m

%% ***POST ADDBACK***: Remove IC artifacts (interactive)
% ONLY GET RID OF THE OBVIOUS ARTIFACTS -- NO BLOOD VESSELS YET (to make
% figures later)
% parameters.amplitude_threshold = 3.5; % 3.5
% parameters.minPixels = 150; % 150
% parameters.large_component_conditional_zscore_flag = true;
% parameters.maxPixels = 5000; 
% parameters.large_component_conditional_zscore_thresh = 1;
% parameters.small_component_conditional_zscore_flag = true;
% parameters.small_component_conditional_zscore_thresh = 2.5;
% 
% % Use second (right) monitor for figures, if available?
% parameters.second_monitor = false;
% 
% % Always clear loop list first. 
% if isfield(parameters, 'loop_list')
% parameters = rmfield(parameters,'loop_list');
% end
% 
% % Dimension different sources are in.
% parameters.sourcesDim = 3;
% 
% parameters.remove_entire_sources = true;
% parameters.use_darkness_threshold = false;
% parameters.draw_artifact_masks = true;
% 
% % Dimension the pixels dimension, different sources dimension of original data.
% parameters.originalSourcesPixelsDim = 2; 
% parameters.originalSourcesDim = 1; 
% 
% % Loop variables % delete 16 & 45 (post-add back & removed)
% parameters.loop_list.iterators = {'mouse', {'loop_variables.mice_all(:).name'}, 'mouse_iterator';
%                                   'source', {'45:70'}, 'source_iterator'};
% parameters.loop_variables.mice_all = parameters.mice_all;
% 
% % Input values
% parameters.loop_list.things_to_load.sources.dir = {[parameters.dir_exper 'spatial segmentation\500 SVD components\regularized ICs ' ... 
%                              num2str(parameters.minPixels) ' amp ' num2str(parameters.amplitude_threshold) ' two conditionals small ' num2str(parameters.small_component_conditional_zscore_thresh) ...
%                               ' large ' num2str(num2str(parameters.large_component_conditional_zscore_thresh)) '\post addback\'], 'mouse', '\'};
% parameters.loop_list.things_to_load.sources.filename= {['sources' num2str(parameters.num_sources) '.mat']};
% parameters.loop_list.things_to_load.sources.variable= {'sources'};
% parameters.loop_list.things_to_load.sources.level = 'mouse';
% 
% parameters.loop_list.things_to_load.indices_of_mask.dir = {[parameters.dir_exper 'preprocessing\masks\']};
% parameters.loop_list.things_to_load.indices_of_mask.filename= {'masks_m', 'mouse', '.mat'};
% parameters.loop_list.things_to_load.indices_of_mask.variable= {'indices_of_mask'}; 
% parameters.loop_list.things_to_load.indices_of_mask.level = 'mouse';
% 
% parameters.loop_list.things_to_load.reference_image.dir = {[parameters.dir_exper 'preprocessing\representative images\'], 'mouse', '\'};
% parameters.loop_list.things_to_load.reference_image.filename= {'reference_image.mat'};
% parameters.loop_list.things_to_load.reference_image.variable= {'reference_image'};
% parameters.loop_list.things_to_load.reference_image.level = 'mouse';
% 
% parameters.loop_list.things_to_load.original_sources.dir = {[parameters.dir_exper 'spatial segmentation\500 SVD components\raw ICs\'], 'mouse', '\'};
% parameters.loop_list.things_to_load.original_sources.filename= {['sources' num2str(parameters.num_sources) '.mat']};
% parameters.loop_list.things_to_load.original_sources.variable= {'sources'};
% parameters.loop_list.things_to_load.original_sources.level = 'mouse';
% 
% % (for any existing artifact removals)
% parameters.loop_list.things_to_load.sources_artifacts_removed_old.dir = {[parameters.dir_exper 'spatial segmentation\500 SVD components\artifacts removed conditional thresholding\post addback without high fine tuning\'], 'mouse', '\'};
% parameters.loop_list.things_to_load.sources_artifacts_removed_old.filename = {'sources.mat'};
% parameters.loop_list.things_to_load.sources_artifacts_removed_old.variable= {'sources'};
% parameters.loop_list.things_to_load.sources_artifacts_removed_old.level = 'mouse';
% 
% % Output values
% parameters.loop_list.things_to_save.sources_artifacts_removed.dir = {[parameters.dir_exper 'spatial segmentation\500 SVD components\artifacts removed conditional thresholding\post addback without high fine tuning\'], 'mouse', '\'};
% parameters.loop_list.things_to_save.sources_artifacts_removed.filename = {'sources.mat'};
% parameters.loop_list.things_to_save.sources_artifacts_removed.variable= {'sources'};
% parameters.loop_list.things_to_save.sources_artifacts_removed.level = 'source';
% 
% RunAnalysis({@RemoveArtifacts}, parameters);


%% POST ADD BACK: Remove Artifacts INCLUDING BLOOD VESSELS

% Use second (right) monitor for figures, if available?
parameters.second_monitor = false;

% Always clear loop list first. 
if isfield(parameters, 'loop_list')
parameters = rmfield(parameters,'loop_list');
end

parameters.amplitude_threshold = 3.5; % 3.5
parameters.minPixels = 150; % 150
parameters.large_component_conditional_zscore_flag = true;
parameters.maxPixels = 5000; 
parameters.large_component_conditional_zscore_thresh = 1;
parameters.small_component_conditional_zscore_flag = true;
parameters.small_component_conditional_zscore_thresh = 2.5;

% Dimension different sources are in.
parameters.sourcesDim = 3;

parameters.remove_entire_sources = false;
parameters.use_darkness_threshold = true;
parameters.draw_artifact_masks = true;

% Dimension the pixels dimension, different sources dimension of original data.
parameters.originalSourcesPixelsDim = 2; 
parameters.originalSourcesDim = 1; 

% Loop variables
parameters.loop_list.iterators = {'mouse', {'loop_variables.mice_all(:).name'}, 'mouse_iterator';
                                  'source', {'32:70'}, 'source_iterator'};
parameters.loop_variables.mice_all = parameters.mice_all;

% Input values
parameters.loop_list.things_to_load.sources.dir = {[parameters.dir_exper 'spatial segmentation\500 SVD components\regularized ICs ' ... 
                             num2str(parameters.minPixels) ' amp ' num2str(parameters.amplitude_threshold) ' two conditionals small ' num2str(parameters.small_component_conditional_zscore_thresh) ...
                              ' large ' num2str(num2str(parameters.large_component_conditional_zscore_thresh)) '\post addback\'], 'mouse', '\'};
parameters.loop_list.things_to_load.sources.filename= {['sources' num2str(parameters.num_sources) '.mat']};
parameters.loop_list.things_to_load.sources.variable= {'sources'};
parameters.loop_list.things_to_load.sources.level = 'mouse';

parameters.loop_list.things_to_load.indices_of_mask.dir = {[parameters.dir_exper 'preprocessing\masks\']};
parameters.loop_list.things_to_load.indices_of_mask.filename= {'masks_m', 'mouse', '.mat'};
parameters.loop_list.things_to_load.indices_of_mask.variable= {'indices_of_mask'}; 
parameters.loop_list.things_to_load.indices_of_mask.level = 'mouse';

parameters.loop_list.things_to_load.reference_image.dir = {[parameters.dir_exper 'preprocessing\representative images\'], 'mouse', '\'};
parameters.loop_list.things_to_load.reference_image.filename= {'reference_image.mat'};
parameters.loop_list.things_to_load.reference_image.variable= {'reference_image'};
parameters.loop_list.things_to_load.reference_image.level = 'mouse';

parameters.loop_list.things_to_load.original_sources.dir = {[parameters.dir_exper 'spatial segmentation\500 SVD components\raw ICs\'], 'mouse', '\'};
parameters.loop_list.things_to_load.original_sources.filename= {['sources' num2str(parameters.num_sources) '.mat']};
parameters.loop_list.things_to_load.original_sources.variable= {'sources'};
parameters.loop_list.things_to_load.original_sources.level = 'mouse';

% (for any existing artifact removals)
parameters.loop_list.things_to_load.sources_artifacts_removed_old.dir = {[parameters.dir_exper 'spatial segmentation\500 SVD components\artifacts removed conditional thresholding\post addback WITH high fine tuning\'], 'mouse', '\'};
parameters.loop_list.things_to_load.sources_artifacts_removed_old.filename = {'sources.mat'};
parameters.loop_list.things_to_load.sources_artifacts_removed_old.variable= {'sources'};
parameters.loop_list.things_to_load.sources_artifacts_removed_old.level = 'mouse';

% Output values
parameters.loop_list.things_to_save.sources_artifacts_removed.dir = {[parameters.dir_exper 'spatial segmentation\500 SVD components\artifacts removed conditional thresholding\post addback WITH high fine tuning\'], 'mouse', '\'};
parameters.loop_list.things_to_save.sources_artifacts_removed.filename = {'sources.mat'};
parameters.loop_list.things_to_save.sources_artifacts_removed.variable= {'sources'};
parameters.loop_list.things_to_save.sources_artifacts_removed.level = 'source';

RunAnalysis({@RemoveArtifacts}, parameters);

%% Plot most recent cleaned overlays 
figure; 
for i = 1:size(mice_all(:))
    mouse = mice_all(i).name;
    try
    load(['Y:\Sarah\Analysis\Experiments\Random Motorized Treadmill\spatial segmentation\500 SVD components\artifacts removed conditional thresholding\post addback WITH high fine tuning\' mouse '\sources.mat']);
    subplot(2, 3, i); imagesc(sources.overlay); axis square; title(mouse); colorbar;
    xticks([]); yticks([]);
    catch
        continue
    end
end 


%% Return any deleted ICs you deleted by accident.
% % Loop variables
% parameters.loop_list.iterators = {'mouse', {'loop_variables.mice_all(:).name'}, 'mouse_iterator'};
% parameters.loop_variables.mice_all = parameters.mice_all;
% 
% % Inputs 
% % Sources
% parameters.loop_list.things_to_load.sources_artifacts_removed.dir = {[parameters.dir_exper 'spatial segmentation\500 SVD components\artifacts_removed\'], 'mouse', '\'};
% parameters.loop_list.things_to_load.sources_artifacts_removed.filename = {'sources.mat'};
% parameters.loop_list.things_to_load.sources_artifacts_removed.variable= {'sources'};
% parameters.loop_list.things_to_load.sources_artifacts_removed.level = 'mouse';
% 
% % Masks
% parameters.loop_list.things_to_load.indices_of_mask.dir = {[parameters.dir_exper 'masks/']};
% parameters.loop_list.things_to_load.indices_of_mask.filename= {'masks_m', 'mouse', '.mat'};
% parameters.loop_list.things_to_load.indices_of_mask.variable= {'indices_of_mask'}; 
% parameters.loop_list.things_to_load.indices_of_mask.level = 'mouse';


% Outputs (editing the same file)
% parameters.loop_list.things_to_load.sources_artifacts_removed.dir = {[parameters.dir_exper 'spatial segmentation\500 SVD components\artifacts_removed\'], 'mouse', '\'};
% parameters.loop_list.things_to_load.sources_artifacts_removed.filename = {'sources.mat'};
% parameters.loop_list.things_to_load.sources_artifacts_removed.variable= {'sources'};
% parameters.loop_list.things_to_load.sources_artifacts_removed.level = 'mouse';
% 
% RunAnalysis({@RestoreSources}, parameters);


%% Plot resulting cleaned ICs
% Always clear loop list first. 
% if isfield(parameters, 'loop_list')
% parameters = rmfield(parameters,'loop_list');
% end
% 
% % Dimension different sources are in.
% parameters.sourcesDim = 3;
% 
% % Loop variables
% parameters.loop_list.iterators = {'mouse', {'loop_variables.mice_all(:).name'}, 'mouse_iterator'};
% parameters.loop_variables.mice_all = parameters.mice_all;
% 
% % Input variables
% parameters.loop_list.things_to_load.sources_overlay.dir = {[parameters.dir_exper 'spatial segmentation\artifacts_removed\'], 'mouse', '\'};
% parameters.loop_list.things_to_load.sources_overlay.filename = {'sources.mat'};
% parameters.loop_list.things_to_load.sources_overlay.variable= {'sources.overlay'};
% parameters.loop_list.things_to_load.sources_overlay.level = 'mouse';
% 
% % Output variables
% parameters.loop_list.things_to_save.sources_overlay.dir = {[parameters.dir_exper 'spatial segmentation\artifacts_removed\'], 'mouse', '\'};
% parameters.loop_list.things_to_save.sources_overlay.filename = {'sources_overlay.fig'};
% parameters.loop_list.things_to_save.sources_overlay.level = 'mouse';
% 
% parameters.loop_list.things_to_save.sources_artifacts_removed.dir = {[parameters.dir_exper 'spatial segmentation\artifacts_removed\'], 'mouse', '\'};
% parameters.loop_list.things_to_save.sources_artifacts_removed.filename = {'sources_overlay.fig'};
% parameters.loop_list.things_to_save.sources_artifacts_removed.level = 'mouse';
% 
% % For now, assume everything was saved as structures from RemoveArtifacts.m
% RunAnalysis({@PlotSourceOverlays}, parameters);

%% Align atlases 
% (Get from the LocaNMF preprocessing pipeline folder you made)

%% Make an atlas colorscheme (make the "values" of each region more sensible)

% % Always clear loop list first. 
% if isfield(parameters, 'loop_list')
% parameters = rmfield(parameters,'loop_list');
% end
% 
% % Load UN-aligned atlas. 
% parameters.loop_list.things_to_load.atlas.dir = {[parameters.dir_exper 'spatial segmentation\aligned atlases\']};
% parameters.loop_list.things_to_load.atlas.filename= {'atlas.mat'};
% parameters.loop_list.things_to_load.atlas.variable= {'atlas'}; 
% parameters.loop_list.things_to_load.atlas.level = 'mouse';
% 
% % Save newly-ordered atlas. 
% parameters.loop_list.things_to_save.atlas_reordered.dir = {[parameters.dir_exper 'spatial segmentation\aligned atlases\']};
% parameters.loop_list.things_to_save.atlas_reordered.filename= {'atlas_reordered.mat'};
% parameters.loop_list.things_to_save.atlas_reordered.variable= {'atlas'}; 
% parameters.loop_list.things_to_save.atlas_reordered.level = 'mouse';


%% Find best atlas locations. 

% % Always clear loop list first. 
% if isfield(parameters, 'loop_list')
% parameters = rmfield(parameters,'loop_list');
% end
% 
% % Dimension different sources are in.
% parameters.sourcesDim = 3;
% 
% % Loop variables; iterate through mice, sources
% parameters.loop_list.iterators = {'mouse', {'loop_variables.mice_all(:).name'}, 'mouse_iterator'};
% parameters.loop_variables.mice_all = parameters.mice_all;
% 
% % Load the "new region order" list, if necessary. Is the list of regions
% % sorted by their center-of-mass from fron to back, left side only. Right
% % side is the negative of these values. (Same for all mice)
% parameters.loop_list.things_to_load.region_order.dir = {[parameters.dir_exper 'spatial segmentation\aligned atlases\']};
% parameters.loop_list.things_to_load.region_order.filename= {'all_regions_sorted.mat'};
% parameters.loop_list.things_to_load.region_order.variable= {'all_regions'}; 
% parameters.loop_list.things_to_load.region_order.level = 'mouse';
% 
% % Load aligned atlas. 
% parameters.loop_list.things_to_load.atlas.dir = {[parameters.dir_exper 'spatial segmentation\aligned atlases\'], 'mouse', '\'};
% parameters.loop_list.things_to_load.atlas.filename= {'atlas.mat'};
% parameters.loop_list.things_to_load.atlas.variable= {'atlas'}; 
% parameters.loop_list.things_to_load.atlas.level = 'mouse';
% 
% % Load corresponding atlas region names
% parameters.loop_list.things_to_load.region_names.dir = {[parameters.dir_exper 'spatial segmentation\aligned atlases\'], 'mouse', '\'};
% parameters.loop_list.things_to_load.region_names.filename= {'atlas.mat'};
% parameters.loop_list.things_to_load.region_names.variable= {'areanames'}; 
% parameters.loop_list.things_to_load.region_names.level = 'mouse';
% 
% % Load mouse's mask, if necessary.
% parameters.loop_list.things_to_load.indices_of_mask.dir = {[parameters.dir_exper 'preprocessing\masks\']};
% parameters.loop_list.things_to_load.indices_of_mask.filename= {'masks_m', 'mouse', '.mat'};
% parameters.loop_list.things_to_load.indices_of_mask.variable= {'indices_of_mask'}; 
% parameters.loop_list.things_to_load.indices_of_mask.level = 'mouse';
% 
% % Load sources with artifacts removed
% parameters.loop_list.things_to_load.sources_artifacts_removed.dir = {[parameters.dir_exper 'spatial segmentation\500 SVD components\artifacts_removed\'], 'mouse', '\'};
% parameters.loop_list.things_to_load.sources_artifacts_removed.filename = {'sources.mat'};
% parameters.loop_list.things_to_load.sources_artifacts_removed.variable= {'sources'};
% parameters.loop_list.things_to_load.sources_artifacts_removed.level = 'mouse';
% 
% % Save metrics 
% parameters.loop_list.things_to_save.metrics.dir = {[parameters.dir_exper 'spatial segmentation\500 SVD components\best region fit\'], 'mouse', '\'};
% parameters.loop_list.things_to_save.metrics.filename= {'best_region_fit_metrics.mat'};
% parameters.loop_list.things_to_save.metrics.variable= {'metrics'};
% parameters.loop_list.things_to_save.metrics.level = 'mouse';
% 
% % Save color - coded atlas. 
% parameters.loop_list.things_to_save.atlas_color_coded.dir = {[parameters.dir_exper 'spatial segmentation\aligned atlases\'], 'mouse', '\'};
% parameters.loop_list.things_to_save.atlas_color_coded.filename= {'atlas_color_coded.mat'};
% parameters.loop_list.things_to_save.atlas_color_coded.variable= {'atlas'};
% parameters.loop_list.things_to_save.atlas_color_coded.level = 'mouse';
% 
% % Save figure of sources over atlas that meet a single best region
% parameters.loop_list.things_to_save.figure_best_fit.dir = {[parameters.dir_exper 'spatial segmentation\500 SVD components\best region fit\'], 'mouse', '\'};
% parameters.loop_list.things_to_save.figure_best_fit.filename= {'best_region_fit.fig'};
% parameters.loop_list.things_to_save.figure_best_fit.variable= {};
% parameters.loop_list.things_to_save.figure_best_fit.level = 'mouse';
% 
% % Save figure of plotted metrics matrices. 
% parameters.loop_list.things_to_save.figure_metrics.dir = {[parameters.dir_exper 'spatial segmentation\500 SVD components\best region fit\'], 'mouse', '\'};
% parameters.loop_list.things_to_save.figure_metrics.filename= {'best_region_fit_metrics.fig'};
% parameters.loop_list.things_to_save.figure_metrics.variable= {};
% parameters.loop_list.things_to_save.figure_metrics.level = 'mouse';
% 
% % Run code
% RunAnalysis({@FindAtlasRegions}, parameters);

%% Final manual assignment 
manual_region_assignments.m;

%% Plot manual assignments next to the color-coded atlas.

% % Always clear loop list first. 
% if isfield(parameters, 'loop_list')
% parameters = rmfield(parameters,'loop_list');
% end
% 
% % Loop variables; iterate through mice, sources
% parameters.loop_list.iterators = {'mouse', {'loop_variables.mice_all(:).name'}, 'mouse_iterator'};
% parameters.loop_variables.mice_all = parameters.mice_all;
% 
% % Load list of manual assignments.
% parameters.loop_list.things_to_load.manual_assignments.dir = {[parameters.dir_exper 'spatial segmentation\500 SVD components\manual assignments\'], 'mouse', '\'};
% parameters.loop_list.things_to_load.manual_assignments.filename = {'region_assignments.mat'};
% parameters.loop_list.things_to_load.manual_assignments.variable= {'region_assignments'};
% parameters.loop_list.things_to_load.manual_assignments.level = 'mouse';
% 
% % Load sources with artifacts removed
% parameters.loop_list.things_to_load.sources_artifacts_removed.dir = {[parameters.dir_exper 'spatial segmentation\500 SVD components\artifacts_removed\'], 'mouse', '\'};
% parameters.loop_list.things_to_load.sources_artifacts_removed.filename = {'sources.mat'};
% parameters.loop_list.things_to_load.sources_artifacts_removed.variable= {'sources'};
% parameters.loop_list.things_to_load.sources_artifacts_removed.level = 'mouse';
% 
% % Load aligned atlas. 
% parameters.loop_list.things_to_load.atlas.dir = {[parameters.dir_exper 'spatial segmentation\aligned atlases\'], 'mouse', '\'};
% parameters.loop_list.things_to_load.atlas.filename= {'atlas_color_coded.mat'};
% parameters.loop_list.things_to_load.atlas.variable= {'atlas'}; 
% parameters.loop_list.things_to_load.atlas.level = 'mouse';
% 
% % Save figure of sources over atlas that meet a single best region
% parameters.loop_list.things_to_save.figure_region_assignments.dir = {[parameters.dir_exper 'spatial segmentation\500 SVD components\manual assignments\'], 'mouse', '\'};
% parameters.loop_list.things_to_save.figure_region_assignments.filename= {'region_assignments.fig'};
% parameters.loop_list.things_to_save.figure_region_assignments.variable= {};
% parameters.loop_list.things_to_save.figure_region_assignments.level = 'mouse';
% 
% RunAnalysis({@PlotRegionAssignments}, parameters);

%% Reorder sources.
% Using the manual assignments, reorder the source colormaps.
% (For now assumes sources are listed in 3D); 

% Always clear loop list first. 
if isfield(parameters, 'loop_list')
parameters = rmfield(parameters,'loop_list');
end

% Loop variables; iterate through mice, sources
parameters.loop_list.iterators = {'mouse', {'loop_variables.mice_all(:).name'}, 'mouse_iterator'};
parameters.loop_variables.mice_all = parameters.mice_all;

% Say if you should add together multiple sources that belong to same new
% sources/nodes (otherwise just overwrites previous ones for now)
parameters.add_sources = true; 

% Say how many new sources/nodes there should be.
parameters.num_new_sources = 32;

% Load list of manual assignments.
parameters.loop_list.things_to_load.assigned_region_order.dir = {[parameters.dir_exper 'spatial segmentation\500 SVD components\manual assignments\'], 'mouse', '\'};
parameters.loop_list.things_to_load.assigned_region_order.filename = {'region_assignments.mat'};
parameters.loop_list.things_to_load.assigned_region_order.variable= {'region_assignments'};
parameters.loop_list.things_to_load.assigned_region_order.level = 'mouse';

% Load sources with artifacts removed
parameters.loop_list.things_to_load.sources.dir = {[parameters.dir_exper 'spatial segmentation\500 SVD components\artifacts removed conditional thresholding\post addback WITH high fine tuning\'], 'mouse', '\'};
parameters.loop_list.things_to_load.sources.filename = {'sources.mat'};
parameters.loop_list.things_to_load.sources.variable= {'sources.sources'};
parameters.loop_list.things_to_load.sources.level = 'mouse';

% Save reordered sources. 
parameters.loop_list.things_to_save.sources_reordered.dir = {[parameters.dir_exper 'spatial segmentation\500 SVD components\manual assignments\'], 'mouse', '\'};
parameters.loop_list.things_to_save.sources_reordered.filename = {'sources_reordered.mat'};
parameters.loop_list.things_to_save.sources_reordered.variable= {'sources'};
parameters.loop_list.things_to_save.sources_reordered.level = 'mouse';

RunAnalysis({@ReorderSources}, parameters);

%% Make a reordered overlay. 

fig = figure; 
fig.WindowState = 'maximized';
for i = 1:size(mice_all, 2)
    mouse = mice_all(i).name;
    
    load(['Y:\Sarah\Analysis\Experiments\Random Motorized Treadmill\spatial segmentation\500 SVD components\manual assignments\' mouse '\sources_reordered.mat']);
    
    overlay = zeros(parameters.yDim, parameters.xDim);
    sources_thresholded = sources > 0; 
    for ici = 1:size(sources_thresholded,3)
       indices = find(sources_thresholded(:,:, ici)); 
       overlay(indices) = ici; 
    end 
    subplot(3, 3, i); imagesc(overlay); axis square; title(mouse); colorbar; colormap([1 1 1; parula(size(sources,3))]);
    xticks([]); yticks([]);
    % save data plotted 
    data_plotted =  overlay;
    save(['Y:\Sarah\Analysis\Experiments\Random Motorized Treadmill\spatial segmentation\500 SVD components\manual assignments\m' mouse 'sources_reordered_overlay.mat'], 'data_plotted');
    
   
end 
 savefig(['Y:\Sarah\Analysis\Experiments\Random Motorized Treadmill\spatial segmentation\500 SVD components\manual assignments\sources_reordered_allmice.fig'])
%% Make a reordered overlay, WITHOUT fine tuning. 

fig = figure; 
fig.WindowState = 'maximized';
for i = 1:size(mice_all, 2)
    mouse = mice_all(i).name;
    
    load(['Y:\Sarah\Analysis\Experiments\Random Motorized Treadmill\spatial segmentation\500 SVD components\manual assignments\' mouse '\sources_reordered.mat']);
    
    overlay = zeros(parameters.yDim, parameters.xDim);
    sources_thresholded = sources > 0; 
    for ici = 1:size(sources_thresholded,3)
       indices = find(sources_thresholded(:,:, ici)); 
       overlay(indices) = ici; 
    end 
    subplot(3, 3, i); imagesc(overlay); axis square; title(mouse); colorbar; colormap([1 1 1; parula(size(sources,3))]);
    xticks([]); yticks([]);
    savefig(['Y:\Sarah\Analysis\Experiments\Random Motorized Treadmill\spatial segmentation\500 SVD components\manual assignments\sources_reordered_allmice_withoutFineTuning.fig'])
end 


%% Apply masks to sources.
% Always clear loop list first. 
if isfield(parameters, 'loop_list')
parameters = rmfield(parameters,'loop_list');
end

% Loop variables; iterate through mice, sources
parameters.loop_list.iterators = {'mouse', {'loop_variables.mice_all(:).name'}, 'mouse_iterator'};
parameters.loop_variables.mice_all = parameters.mice_all;

% Input
parameters.loop_list.things_to_load.indices_of_mask.dir = {[parameters.dir_exper 'preprocessing\masks\']};
parameters.loop_list.things_to_load.indices_of_mask.filename= {'masks_m', 'mouse', '.mat'};
parameters.loop_list.things_to_load.indices_of_mask.variable= {'indices_of_mask'}; 
parameters.loop_list.things_to_load.indices_of_mask.level = 'mouse';

parameters.loop_list.things_to_load.data.dir = {[parameters.dir_exper 'spatial segmentation\500 SVD components\manual assignments\'], 'mouse', '\'};
parameters.loop_list.things_to_load.data.filename = {'sources_reordered.mat'};
parameters.loop_list.things_to_load.data.variable= {'sources'};
parameters.loop_list.things_to_load.data.level = 'mouse';

% Output
parameters.loop_list.things_to_save.data_masked.dir = {[parameters.dir_exper 'spatial segmentation\500 SVD components\manual assignments\'], 'mouse', '\'};
parameters.loop_list.things_to_save.data_masked.filename = {'sources_reordered_masked.mat'};
parameters.loop_list.things_to_save.data_masked.variable= {'sources_masked'};
parameters.loop_list.things_to_save.data_masked.level = 'mouse';

RunAnalysis({@ApplyMasks}, parameters);