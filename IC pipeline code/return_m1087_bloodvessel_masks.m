% return_m1087_bloodvessel_masks.m

% Sarah West
% 5/5/22


% Compare artifact removed overlay numbers from:
% Y:\Sarah\Analysis\Experiments\...
% Random Motorized Treadmill\spatial segmentation\500 SVD components\
% artifacts_removed\1087\sources.mat 
% to:
% Y:\Sarah\Analysis\Experiments\Random Motorized Treadmill\spatial
% segmentation\500 SVD components\artifacts removed conditional
% thresholding\post addback without high fine tuning\1087\sources.mat


% Save to "post addback WITH high fine tuning\1087" folder 
% Column 1 = old. Column 2 = new
IC_numbers = [
    13, 14
    19, 20
    29, 33
    12, 13
    7 , 8
    6, 7
    23, 25
    25, 27
    3, 4
    27, 29
    22, 24
    5, 6
    4, 5
    2, 2
    1, 1
    11, 12
    8, 9
    9, 10
    10, 11
    14, 15
    20, 21
    24, 26
    28, 30
    16, 17
    17, 18
    21, 23
    18, 19
    15, 16
    26, 28
];

%%
% Find location/indices of "old" IC numbers with removed ICs taken into
% account. 
% Load artifact removed sources from "new"
 
 load('Y:\Sarah\Analysis\Experiments\Random Motorized Treadmill\spatial segmentation\500 SVD components\artifacts removed conditional thresholding\post addback without high fine tuning\1087\sources.mat')
 indices_new = find(~cellfun(@ischar, sources.indices_to_remove)); 
 
 % do same thing but for "old"
 load('Y:\Sarah\Analysis\Experiments\Random Motorized Treadmill\spatial segmentation\500 SVD components\artifacts_removed\1087\sources.mat')

 indices_old = find(~cellfun(@ischar, sources.indices_to_remove)); 


 % With these indices, pull out appropriate artifact masks and artifact
 % indices
 artifact_masks_new = cell(1,1);
 indices_to_remove_new= cell(1,1);
 for i = 1:size(IC_numbers,1)

        old_ic = IC_numbers(i,1); 
        new_ic = IC_numbers(i,2);
        artifact_masks_new{indices_new(new_ic),1} = sources.artifact_masks{indices_old(old_ic)};
        indices_to_remove_new{indices_new(new_ic),1} = sources.indices_to_remove{indices_old(old_ic)}; 
 
 end 

 %% Load artifacts removed conditional thresholding WITHOUT high fine tuning.
 load('Y:\Sarah\Analysis\Experiments\Random Motorized Treadmill\spatial segmentation\500 SVD components\artifacts removed conditional thresholding\post addback without high fine tuning\1087\sources.mat')

% Concatenate new indices and artifact masks in. 
for i = 1:numel(indices_to_remove_new)
    
    
    % Concatenate
    try
        if ~isempty(indices_to_remove_new{i})
            sources.indices_to_remove{i} =  [sources.indices_to_remove{i}; indices_to_remove_new{i}]; 
            sources.artifact_masks{i} =  cat(3, sources.artifact_masks{i}, artifact_masks_new{i}); 
           
        end
    catch

        disp('here');
    end 
end 

% Save as artifacts removed conditional thresholding WITH high fine tuning.
save('Y:\Sarah\Analysis\Experiments\Random Motorized Treadmill\spatial segmentation\500 SVD components\artifacts removed conditional thresholding\post addback WITH high fine tuning\1087\sources.mat', 'sources');