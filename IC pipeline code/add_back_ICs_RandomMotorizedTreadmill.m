% add_back_ICs_RandomMotorizedTreadmill.m
% Sarah West
% 4/29/22

% A script to add back ICs that you think are real but were taken out by
% regularize_ICs.

% Wil run as a script, saves to end of list of artifact_removed sources
% so it doesn't disrupt any previously saved artifact
% removals.

% All mice used:
% 500 SVD components, 100 ICs.
% parameters.amplitude_threshold = 3.5; 
% parameters.minPixels = 150;
% parameters.maxPixels = 5000; 
% parameters.large_component_conditional_zscore_thresh = 1;
% parameters.small_component_conditional_zscore_thresh = 2.5; 

%%
% Mouse 1087
% Notes:
% Right rostral medial M2 (in same IC as left rostral medial M2, artifacts 
% removed IC #20). 
% Left caudal M2 (in same IC as right caudal M2, artifacts removed IC # 4)
% Right medial parietal (to match artifacts removed IC #1, haven't found
% where it might be yet).
% More of left lateral M1(or S1?) (to add to artifacts removed IC #22)

% load raw ICs & masks to manipulate
add_back = cell(1);
original_IC_numbers = cell(1);
raw_ICs = cell(1);
load('Y:\Sarah\Analysis\Experiments\Random Motorized Treadmill\spatial segmentation\500 SVD components\raw ICs\1087\sources100.mat', 'sources');
load('Y:\Sarah\Analysis\Experiments\Random Motorized Treadmill\preprocessing\masks\masks_m1087.mat');
% Fill masks, convert to absolte value
sources = abs(FillMasks(sources', indices_of_mask, 256, 256));

% Right medial rostral M2
holder = NaN(256);
holder1 = NaN(256);
source = sources(:, :,19);
holder1(1:90, 160:256) = source(1:90, 160:256); % Limit area
inds = find(holder1 > 2.0); % Threshold
holder(inds) = source(inds);
figure; imagesc(holder); colorbar;
add_back{1} = holder;
original_IC_numbers{1} = 19;
raw_ICs{1} = source;

% Left caudal M2
holder = NaN(256);
holder1 = NaN(256);
source = sources(:, :,3);
holder1(88:130, 90:125) = source(88:130, 90:125); % Limit area
inds = find(holder1 > 2.0); % Threshold
holder(inds) = source(inds);
figure; imagesc(holder); colorbar;
add_back{2} = holder;
original_IC_numbers{2} = 3;
raw_ICs{2} = source;

% right medial parietal
holder = NaN(256);
holder1 = NaN(256);
source = sources(:, :,1);
holder1(115:185, 134:165) = source(115:185, 134:165); % Limit area
inds = find(holder1 > 1.5); % Threshold
holder(inds) = source(inds);
figure; imagesc(holder); colorbar; caxis([0 4])
add_back{3} = holder;
original_IC_numbers{3} = 1;
raw_ICs{3} = source;

% more of left lateral M1 (ish)
holder = NaN(256);
holder1 = NaN(256);
source =  nansum(sources(:, :,[13 27 34:36]),3)/5;
holder1(88:126,28:55) = source(88:126,28:55); % Limit area
inds = find(holder1 > 1.5); % Threshold
holder(inds) = source(inds);
figure; imagesc(holder); colorbar; caxis([0 4]);
add_back{4} = holder;
original_IC_numbers{4} = 13;
raw_ICs{4} = source;

add_back_all{1} = add_back;
original_IC_numbers_all{1} = original_IC_numbers;
raw_ICs_all{1} = raw_ICs;

close all;
%% 
% Mouse 1088
% Notes: 
% (done)Right caudal M2 (to go with left caudal M2, artifacts removed IC 31)
% (done) See if there's more of artifacts removed IC 4 (might be more of it in the
% medial direction, to match similar ICs in other mice). 
% (done) Also see if corresponding ICs on the left (13 & 26) have more area between them.
% (done)Find more of IC 12 (cut off by blood vessels) 
% (done) Find more of IC 38 (cut off by blood vessels)
% (done) Fnd more of IC 2 
% (done) See if there's a matching on left for IC 25 (would let me confidently
% count those as retrosplenial & throw out 27 & 28, which might be gunk)
% See if there are corresponding ICs in other mice to 36 & 37 

% load raw ICs & masks to manipulate
add_back = cell(1);
original_IC_numbers = cell(1);
raw_ICs = cell(1);
load('Y:\Sarah\Analysis\Experiments\Random Motorized Treadmill\spatial segmentation\500 SVD components\raw ICs\1088\sources100.mat', 'sources');
load('Y:\Sarah\Analysis\Experiments\Random Motorized Treadmill\preprocessing\masks\masks_m1088.mat');
% Fill masks, convert to absolte value
sources = abs(FillMasks(sources', indices_of_mask, 256, 256));

% Right caudal M2
holder = NaN(256);
holder1 = NaN(256);
source =  sources(:, :,40);
holder1(56:144, 131:161) = source(56:144, 131:161); % Limit area
inds = find(holder1 > 1); % Threshold
holder(inds) = source(inds);
figure; imagesc(holder); colorbar; caxis([0 4]);
add_back{1} = holder;
original_IC_numbers{1} = 40;
raw_ICs{1} = source;

% left anterior-ish M1 (between ICs 13 & 26)
holder = NaN(256);
holder1 = NaN(256);
source =  sources(:, :,58);
holder1(85:130,37:83) = source(85:130,37:83); % Limit area
inds = find(holder1 > 1.5); % Threshold
holder(inds) = source(inds);
figure; imagesc(holder); colorbar; caxis([0 5]);
add_back{2} = holder;
original_IC_numbers{2} = 58;
raw_ICs{2} = source;

% More of 12
holder = NaN(256);
holder1 = NaN(256);
source =  sources(:, :,13);
%holder1() = source(85:130,37:83); % Limit area
inds = find(source > 3.5); % Threshold
holder(inds) = source(inds);
figure; imagesc(holder); colorbar; caxis([0 5]);
add_back{3} = holder;
original_IC_numbers{3} = 13;
raw_ICs{3} = source;

% More of 38
holder = NaN(256);
holder1 = NaN(256);
source =  sources(:, :,68);
holder1(110:155, 1:42) = source(110:155, 1:42); % Limit area
inds = find(holder1 > 3.5); % Threshold
holder(inds) = source(inds);
figure; imagesc(holder); colorbar; caxis([0 5]);
add_back{4} = holder;
original_IC_numbers{4} = 68;
raw_ICs{4} = source;

% More of IC 2 (right caudal M2)
holder = NaN(256);
holder1 = NaN(256);
source =  sources(:, :,1);
holder1(120:190, 138:175) = source(120:190, 138:175); % Limit area
inds = find(holder1 > 2.5); % Threshold
holder(inds) = source(inds);
figure; imagesc(holder); colorbar; caxis([0 5]);
add_back{5} = holder;
original_IC_numbers{5} = 1;
raw_ICs{5} = source;

% Left retrosplenial
holder = NaN(256);
holder1 = NaN(256);
source =  sources(:, :,32);
holder1(160:240, 70:120) = source(160:240, 70:120); % Limit area
inds = find(holder1 > 2.5); % Threshold
holder(inds) = source(inds);
figure; imagesc(holder); colorbar; caxis([0 5]);
add_back{6} = holder;
original_IC_numbers{6} = 32;
raw_ICs{6} = source;

add_back_all{2} = add_back;
original_IC_numbers_all{2} = original_IC_numbers;
raw_ICs_all{2} = raw_ICs; 
close all;
%%
% Mouse 1096
% Notes:
% (all ICs labeled from pre-artifacts removed)
% (done) Right lateral M2, is part of IC 16 
% IC 25 is part of IC 24 (medial rostral M2s, I think), so should be
% combined with IC22
% (done) IC 26 has a left corresponding regions that was small & thresholded out
% (done, might be artifact) IC 30 (a left mid parietal) might have a corresponding domain on rightin same IC
% (done) C 32 (right medial rostral M2) has a left domain that might work better
% than IC 22 --> Actually, is IC 33, but was thresholded very small.
% 34 & 35 are part of 33
% (done) IC 43 has a right side that wasn't saved, is beneath a disctinct part of IC 17. 
% I can't tell if 49 (pre-artifacts removed) is real or not.
% (done) 57 might be real, with corresponding left IC (re-threshold both)
% (done) Right lateral parietal (corresponds to IC 39)
add_back = cell(1);
original_IC_numbers = cell(1);
raw_ICs = cell(1);
load('Y:\Sarah\Analysis\Experiments\Random Motorized Treadmill\spatial segmentation\500 SVD components\raw ICs\1096\sources100.mat', 'sources');
load('Y:\Sarah\Analysis\Experiments\Random Motorized Treadmill\preprocessing\masks\masks_m1096.mat');
% Fill masks, convert to absolte value
sources = abs(FillMasks(sources', indices_of_mask, 256, 256));

% Right lateral rostral M2
holder = NaN(256);
holder1 = NaN(256);
source =  sources(:, :,12);
figure; imagesc(source);
holder1(1:75, 180:256) = source(1:75, 180:256); % Limit area
inds = find(holder1 > 2); % Threshold
holder(inds) = source(inds);
figure; imagesc(holder); colorbar; caxis([0 4]);
add_back{1} = holder;
original_IC_numbers{1} = 12;
raw_ICs{1} = source;

% mid parietal (IC 26 left part)
holder = NaN(256);
holder1 = NaN(256);
source =  sources(:, :, 19);
figure; imagesc(source);
holder1(:, 1:142) = source(:, 1:142); % Limit area
inds = find(holder1 > 3.5); % Threshold
holder(inds) = source(inds);
figure; imagesc(holder); colorbar; caxis([0 4]);
add_back{2} = holder;
original_IC_numbers{2} = 19;
raw_ICs{2} = source;

% right part of 30
holder = NaN(256);
holder1 = NaN(256);
source =  sources(:, :, 23);
figure; imagesc(source);
holder1(150:180, 170:215) = source(150:180, 170:215); % Limit area
inds = find(holder1 > 2.5); % Threshold
holder(inds) = source(inds);
figure; imagesc(holder); colorbar; caxis([0 4]);
add_back{3} = holder;
original_IC_numbers{3} = 23;
raw_ICs{3} = source;

% IC 33, better left medial rostral M2
holder = NaN(256);
holder1 = NaN(256);
source =  sources(:, :, 25);
figure; imagesc(source);
holder1(25:75,40:90) = source(25:75,40:90); % Limit area
inds = find(holder1 > 3.5); % Threshold
holder(inds) = source(inds);
figure; imagesc(holder); colorbar; caxis([0 4]);
add_back{4} = holder;
original_IC_numbers{4} = 25;
raw_ICs{4} = source;

% reight retrosplenial, right part of IC 43 
% holder = NaN(256);
% holder1 = NaN(256);
% source =  sources(:, :, 34);
% figure; imagesc(source);
% holder1(140:215, 136:190) = source(140:215, 136:190); % Limit area
% inds = find(holder1 > 3.5); % Threshold
% holder(inds) = source(inds);
% figure; imagesc(holder); colorbar; caxis([0 4]);
% add_back{5} = holder;
% raw_ICs{5} = source;

% right lateral parietal (corresponds to IC 39)
source =  nansum(sources(:, :, [39 54 59 61]), 3)./4;
figure; imagesc(source);
holder = NaN(256);
holder1 = NaN(256);
holder1(100:150, 1:46) = source(100:150, 1:46); % Limit area
inds = find(holder1 > 2.5); % Threshold
holder(inds) = source(inds);
figure; imagesc(holder); colorbar; caxis([0 4]);
add_back{6} = holder;
original_IC_numbers{6} = source;
raw_ICs{6} = source;

% right part of 57
source =  sources(:, :, 89);
figure; imagesc(source);
holder = NaN(256);
holder1 = NaN(256);
holder1(145:180, 190:235) = source(145:180, 190:235); % Limit area
inds = find(holder1 > 3.5); % Threshold
holder(inds) = source(inds);
figure; imagesc(holder); colorbar; caxis([0 4]);
add_back{7} = holder;
original_IC_numbers{7} = 89;
raw_ICs{7} = source;

% left part of 57
source =  sources(:, :, 89);
figure; imagesc(source);
holder = NaN(256);
holder1 = NaN(256);
holder1(150:190, 20:60) = source(150:190, 20:60); % Limit area
inds = find(holder1 > 3); % Threshold
holder(inds) = source(inds);
figure; imagesc(holder); colorbar; caxis([0 4]);
add_back{8} = holder;
original_IC_numbers{8} = 89;
raw_ICs{8} = source;

add_back_all{3} = add_back;
original_IC_numbers_all{3} = original_IC_numbers; 
raw_ICs_all{3} = raw_ICs; 
close all;
%%

% Mouse 1100
% (you decided to use the same thresholds as the other mice) 

% Notes:
% (all ICs labeled from pre-artifacts removed)
% nodes missing:
% (don't do)Find more of node 1 (left side of raw IC 30? Or raw IC 79, which you originally threw out, or left of raw 38)
% (don't do)2 -- right side of raw IC 21 (maybe + raw IC 18) & don't use 34; or use IC 34 as 2? 
% (done) For 4: right side of 16
% (done) 7 -- left side of raw IC 2, in a pinch
% (done) 9  -- raw IC 12 
% (done) 18 -- IC 1
% (done) 19 -- Right side of raw 11 % might be raw ICs 27 + 28 (rough)
% (done) 25 -- left side of raw IC 17
% (doone) 27 -- (might have been incorporated into IC 12/node 23); divide IC 12
% /node 23/ raw IC 10 based on boundaries suggested by left side of raw IC 29/ node 28?)
% (don't do) 31 & 32 -- might actually be raw IC 4

add_back = cell(1);
original_IC_numbers = cell(1);
raw_ICs = cell(1);
load('Y:\Sarah\Analysis\Experiments\Random Motorized Treadmill\spatial segmentation\500 SVD components\raw ICs\1100_first mask version\sources100.mat', 'sources');
load('Y:\Sarah\Analysis\Experiments\Random Motorized Treadmill\preprocessing\masks\masks_m1100_firstversion.mat');
% Fill masks, convert to absolte value
sources = abs(FillMasks(sources', indices_of_mask, 256, 256));

% more of node 1
% holder = NaN(256);
% holder1 = NaN(256);
% source =  sum(sources(:, :,[30 38]), 3, 'omitnan')./2;
% figure; imagesc(source);
% holder1() = source(); % Limit area
% inds = find(holder1 > 2); % Threshold
% holder(inds) = source(inds);
% figure; imagesc(holder); colorbar; caxis([0 4]);
% add_back{1} = holder;
% original_IC_numbers{1} = NaN;
% raw_ICs{1} = source;

% New node 4
holder = NaN(256);
holder1 = NaN(256);
source = sources(:,:, 16);
figure; imagesc(source);
holder1(1:88, 137:end) = source(1:88, 137:end); % Limit area
inds = find(holder1 > 3); % Threshold
holder(inds) = source(inds);
figure; imagesc(holder); colorbar; caxis([0 4]);
add_back{1} = holder;
original_IC_numbers{1} = 16;
raw_ICs{1} = source;


% Node 7
source = sources(:,:, 2);
figure; imagesc(source);
holder = NaN(256);
holder1 = NaN(256);
holder1(66:116,36:79) = source(66:116,36:79); % Limit area
inds = find(holder1 > 1); % Threshold
holder(inds) = source(inds);
figure; imagesc(holder); colorbar; caxis([0 4]);
add_back{2} = holder;
original_IC_numbers{2} = 2;
raw_ICs{2} = source;

% Node 9
holder = NaN(256);
holder1 = NaN(256);
source = sources(:,:, 12);
figure; imagesc(source);
holder1(30:110,70:126) = source(30:110,70:126); % Limit area
inds = find(holder1 > 2); % Threshold
holder(inds) = source(inds);
figure; imagesc(holder); colorbar; caxis([0 4]);
add_back{3} = holder;
original_IC_numbers{3} = 12;
raw_ICs{3} = source;

% Node 18
source = sources(:,:, 1);
figure; imagesc(source);
holder = NaN(256);
holder1 = NaN(256);
holder1(70:150,125:180) = source(70:150,125:180); % Limit area
inds = find(holder1 > 1.5); % Threshold
holder(inds) = source(inds);
figure; imagesc(holder); colorbar; caxis([0 4]);
add_back{4} = holder;
original_IC_numbers{4} = 1;
raw_ICs{4} = source;

% Node 19
% source =  sum(sources(:, :,[27 28 52 55 57:60 70:73]), 3, 'omitnan')./12;
% figure; imagesc(source);
% holder = NaN(256);
% holder1 = NaN(256);
% holder1(110:155,26:60) = source(110:155,26:60); % Limit area
% inds = find(holder1 > 1.25 & holder1 < 2.5); % Threshold
% holder(inds) = source(inds);
% figure; imagesc(holder); colorbar; caxis([0 4]);
% add_back{5} = holder;
% original_IC_numbers{5} = NaN;
% raw_ICs{5} = source;

% Node 19-- right side of 11
source = sources(:,:, 11);
figure; imagesc(source);
holder = NaN(256);
holder1 = NaN(256);
holder1(120:180, 1:60) = source(120:180, 1:60); % Limit area
inds = find(holder1 > 1.75); % Threshold
holder(inds) = source(inds);
figure; imagesc(holder); colorbar; caxis([0 4]);
add_back{5} = holder;
original_IC_numbers{5} = 11;
raw_ICs{5} = source;

% Node 25
source = sources(:, :,17);
figure; imagesc(source);
holder = NaN(256);
holder1 = NaN(256);
holder1(176:191, 45:72) = source(176:191, 45:72); % Limit area
inds = find((holder1 > 0.3 )); % Threshold
holder(inds) = source(inds);
figure; imagesc(holder); colorbar; caxis([0 4]);
add_back{6} = holder;
original_IC_numbers{6} = 17;
raw_ICs{6} = source;

% Node 27
source = sources(:, :,29);
figure; imagesc(source);
holder = NaN(256);
holder1 = NaN(256);
holder1(173:209, 1:45) = source(173:209, 1:45); % Limit area
inds = find((holder1 > .8 )); % Threshold
holder(inds) = source(inds);
figure; imagesc(holder); colorbar; caxis([0 4]);
add_back{7} = holder;
original_IC_numbers{7} = 29;
raw_ICs{7} = source;


add_back_all{4} = add_back;
original_IC_numbers_all{4} = original_IC_numbers;
raw_ICs_all{4} = raw_ICs;

close all;
%%
% Mouse 1106
% Notes:
% (all ICs labeled from pre-artifacts removed)
% (done) IC 1 should be re-thresholded
% (done) IC 6 has a corresponding IC on the left (retrosplenial)
% 10 & 11 are same source
% 16 is the corresponding IC to 15
% 24 has a left domain that aligns well with IC 10
% (done) 26 may have corresponding IC on left
% (done) 27 has corresponding IC on right (lateral rostral M2)
% (done) 28 (mid parietal) has strong but small corresponding IC on left
% 30 might be a blood vessel
% (done) Rethreshold 38 (left visual was cut off by blood vessel) 
% (done) Better caudal M2s
add_back = cell(1);
original_IC_numbers = cell(1);
raw_ICs = cell(1);
load('Y:\Sarah\Analysis\Experiments\Random Motorized Treadmill\spatial segmentation\500 SVD components\raw ICs\1106\sources100.mat', 'sources');
load('Y:\Sarah\Analysis\Experiments\Random Motorized Treadmill\preprocessing\masks\masks_m1106.mat');
% Fill masks, convert to absolte value
sources = abs(FillMasks(sources', indices_of_mask, 256, 256));

% IC 1 rethreshold
source =  sources(:, :,1);
figure; imagesc(source);
holder = NaN(256);
holder1 = NaN(256);
holder1(:, 1:125) = source(:, 1:125); % Limit area
inds = find(holder1 > 2.5); % Threshold
holder(inds) = source(inds);
figure; imagesc(holder); colorbar; caxis([0 4]);
add_back{1} = holder;
original_IC_numbers{1} = 1;
raw_ICs{1} = source;

% Left side of IC 6 (left retrosplenial)
source =  sources(:, :,5);
figure; imagesc(source);
holder = NaN(256);
holder1 = NaN(256);
holder1(110:225, 60:120) = source(110:225, 60:120); % Limit area
inds = find(holder1 > 2.5); % Threshold
holder(inds) = source(inds);
figure; imagesc(holder); colorbar; caxis([0 4]);
add_back{2} = holder;
original_IC_numbers{2} = 5;
raw_ICs{2} = source;

% Left side of IC 26 (posterior parietal)
source =  sources(:, :,25);
figure; imagesc(source);
holder = NaN(256);
holder1 = NaN(256);
holder1(168:217,17:80) = source(168:217,17:80); % Limit area
inds = find(holder1 > 2.0); % Threshold
holder(inds) = source(inds);
figure; imagesc(holder); colorbar; caxis([0 4]);
add_back{3} = holder;
original_IC_numbers{3} = 25;
raw_ICs{3} = source;

% Right lateral rostral M2
source =  sources(:, :,26);
figure; imagesc(source);
holder = NaN(256);
holder1 = NaN(256);
holder1(1:80, 172:256) = source(1:80, 172:256); % Limit area
inds = find(holder1 > 2.5); % Threshold
holder(inds) = source(inds);
figure; imagesc(holder); colorbar; caxis([0 4]);
add_back{4} = holder;
original_IC_numbers{4} = 26;
raw_ICs{4} = source;

% left part of 28 (mid parietal)
source =  sources(:, :,27);
figure; imagesc(source);
holder = NaN(256);
holder1 = NaN(256);
holder1(150:180, 1:125) = source(150:180, 1:125); % Limit area
inds = find(holder1 > 2.5); % Threshold
holder(inds) = source(inds);
figure; imagesc(holder); colorbar; caxis([0 4]);
add_back{5} = holder;
original_IC_numbers{5} = 27;
raw_ICs{5} = source;

% rethreshold 38 (left visual) 
source =  sources(:, :,42);
figure; imagesc(source);
holder = NaN(256);
holder1 = NaN(256);
holder1(171:end, 14:125) = source(171:end, 14:125); % Limit area
inds = find(holder1 > 3.5); % Threshold
holder(inds) = source(inds);
figure; imagesc(holder); colorbar; caxis([0 4]);
add_back{6} = holder;
original_IC_numbers{6} = 42;
raw_ICs{6} = source;

% Better left caudal M2
source =  nansum(sources(:, :, [11 12]), 3)./2; % 12
figure; imagesc(source);
holder = NaN(256);
holder1 = NaN(256);
holder1(45:133, 80:125) = source(45:133, 80:125); % Limit area
inds = find(holder1 > 2.0); % Threshold
holder(inds) = source(inds);
figure; imagesc(holder); colorbar; caxis([0 4]);
add_back{7} = holder;
original_IC_numbers{7} = source;
raw_ICs{7} = source;

% Better right caudal M2
source =  nansum(sources(:, :, [11 12]), 3)./2; 
figure; imagesc(source);
holder = NaN(256);
holder1 = NaN(256);
holder1(36:135, 125:160) = source(36:135, 125:160); % Limit area
inds = find(holder1 > 2.5); % Threshold
holder(inds) = source(inds);
figure; imagesc(holder); colorbar; caxis([0 4]);
add_back{8} = holder;
original_IC_numbers{8} = source;
raw_ICs{8} = source;

add_back_all{5} = add_back;
original_IC_numbers_all{5} = original_IC_numbers;
raw_ICs_all{5} = raw_ICs;

close all;
%%
% Mouse 1107
% Notes:
% (all ICs labeled from pre-artifacts removed)
% 1 & 2 are from the same source
% 5 has corresponding left visual, also a right medial motor in it
% 8 is part of 7
% (done) 17 might need to be re-thresholded as a retrosplenial.
% (done) 21 has a left corresponding IC (medial rostral M2) + 26
% (done) 30 has a corresponding rightthat might work as left retrosplenial
% 39 is same source as 37
% 47 might just be gunk
% (done) Find right lateral rostral M2
% (done) re-threshold 48
% (done) re-threshold 52
% (done) Look for caudal M2s (?)
% IC 29 might just be a blood vessel
% (done) Add back posterior-lateral parietal/visual (raw IC 22)
% Right medial M1-type ? 
% Get top of raw IC 1 back 
add_back = cell(1);
original_IC_numbers = cell(1);
raw_ICs = cell(1);

load('Y:\Sarah\Analysis\Experiments\Random Motorized Treadmill\spatial segmentation\500 SVD components\raw ICs\1107\sources100.mat', 'sources');
load('Y:\Sarah\Analysis\Experiments\Random Motorized Treadmill\preprocessing\masks\masks_m1107.mat');
% Fill masks, convert to absolte value
sources = abs(FillMasks(sources', indices_of_mask, 256, 256));

% Left visual
source =  sources(:, :,4);
figure; imagesc(source);
holder = NaN(256);
holder1 = NaN(256);
holder1(170:215, 22:125) = source (170:215, 22:125); % Limit area
inds = find(holder1 > 2.0); % Threshold
holder(inds) = source(inds);
figure; imagesc(holder); colorbar; caxis([0 4]);
add_back{1} = holder;
original_IC_numbers{1} = 4;
raw_ICs{1} = source;

% medial motor -type (corresponding to IC 1)
source =  sources(:, :,4);
figure; imagesc(source);
holder = NaN(256);
holder1 = NaN(256);
holder1(64:142, 125:170) = source (64:142, 125:170); % Limit area
inds = find(holder1 > 2.0); % Threshold
holder(inds) = source(inds);
figure; imagesc(holder); colorbar; caxis([0 4]);
add_back{2} = holder;
original_IC_numbers{2} = 4;
raw_ICs{2} = source;

% Left retrosplenial
source =  nansum(sources(:, :, [15 28]), 3)./2;
figure; imagesc(source);
holder = NaN(256);
holder1 = NaN(256);
holder1(138:243, 74:125) = source(138:243, 74:125); % Limit area
inds = find(holder1 > 2.0); % Threshold
holder(inds) = source(inds);
figure; imagesc(holder); colorbar; caxis([0 4]);
add_back{3} = holder;
original_IC_numbers{3} = source;
raw_ICs{3} = source;

% Left medial rostral M2 
source =  nansum(sources(:, :, [20 26]), 3)./2;
figure; imagesc(source);
holder = NaN(256);
holder1 = NaN(256);
holder1(1:85, 29:92) = source(1:85, 29:92); % Limit area
inds = find(holder1 > 1.5); % Threshold
holder(inds) = source(inds);
figure; imagesc(holder); colorbar; caxis([0 4]);
add_back{4} = holder;
original_IC_numbers{4} = source;
raw_ICs{4} = source;

% Right lateral rostral M2 
source =  sources(:, :,17);
figure; imagesc(source);
holder = NaN(256);
holder1 = NaN(256);
holder1(1:80, 195:end) = source(1:80, 195:end); % Limit area
inds = find(holder1 > 1.5); % Threshold
holder(inds) = source(inds);
figure; imagesc(holder); colorbar; caxis([0 4]);
add_back{5} = holder;
original_IC_numbers{5} = 17;
raw_ICs{5} = source;

% Put back raw source 22
source =  sources(:, :,22);
figure; imagesc(source);
holder = NaN(256);
holder1 = NaN(256);
holder1(187:256, 1:65) = source(187:256, 1:65); % Limit area
inds = find(holder1 > 3.5); % Threshold
holder(inds) = source(inds);
figure; imagesc(holder); colorbar; caxis([0 4]);
add_back{6} = holder;
original_IC_numbers{6} = 22;
raw_ICs{6} = source;

% Re-threshold 48;
source =  sources(:, :,54);
figure; imagesc(source);
holder = NaN(256);
holder1 = NaN(256);
holder1(109:172, 211:end) = source(109:172, 211:end); % Limit area
inds = find(holder1 > 2.5); % Threshold
holder(inds) = source(inds);
figure; imagesc(holder); colorbar; caxis([0 4]);
add_back{7} = holder;
original_IC_numbers{7} = 54;
raw_ICs{7} = source;

% Re-threshold 52
source =  sources(:, :,74);
figure; imagesc(source);
holder = NaN(256);
holder1 = NaN(256);
holder1(134:198 ,1:45) = source(134:198 ,1:45); % Limit area
inds = find(holder1 > 3); % Threshold
holder(inds) = source(inds);
figure; imagesc(holder); colorbar; caxis([0 4]);
add_back{8} = holder;
original_IC_numbers{8} = 74;
raw_ICs{8} = source;

% Left caudal M2
% source =  nansum(sources(:, :, [12 19 26]), 3)./3; 
% figure; imagesc(source);
% holder = NaN(256);
% holder1 = NaN(256);
% holder1(43:114, 66:125) = source(43:114, 66:125); % Limit area
% inds = find(holder1 >2); % Threshold
% holder(inds) = source(inds);
% figure; imagesc(holder); colorbar; caxis([0 4]);
% add_back{9} = holder;
% original_IC_numbers{9} = source;
%raw_ICs{9} = source;

% Right caudal M2
% source =  sources(:,:,39);
% figure; imagesc(source);
% holder = NaN(256);
% holder1 = NaN(256);
% holder1(1:97, 125:190) = source(1:97, 125:190); % Limit area
% inds = find(holder1 > 4); % Threshold
% holder(inds) = source(inds);
% figure; imagesc(holder); colorbar; caxis([0 6]);
% add_back{10} = holder;
% original_IC_numbers{10} = source;
% raw_ICs{10} = source;

% Right medial M1-ish
source =  nansum(sources(:, :, [1 ]), 3)./1; 
figure; imagesc(source);
holder = NaN(256);
holder1 = NaN(256);
holder1(120:200, 125:175) = source(120:200, 125:175); % Limit area
inds = find(holder1 >2.5); % Threshold
holder(inds) = source(inds);
figure; imagesc(holder); colorbar; caxis([0 4]);
add_back{11} = holder;
original_IC_numbers{11} = source;
raw_ICs{11} = source;

% Get top of raw IC 1 back 
source = sources(:, :, 1); 
figure; imagesc(source);
holder = NaN(256);
holder1 = NaN(256);
holder1(50:88, 93:125) = source(50:88, 93:125); % Limit area
inds = find(holder1 >2.5); % Threshold
holder(inds) = source(inds);
figure; imagesc(holder); colorbar; caxis([0 4]);
add_back{12} = holder;
original_IC_numbers{12} = 1;
raw_ICs{12} = source;

add_back_all{6} = add_back;
original_IC_numbers_all{6} = original_IC_numbers;
raw_ICs_all{6} = raw_ICs;

%%
% Mouse 539
% Notes:
% (all ICs labeled from pre-artifacts removed)
% 2 is part of left part of 1 (node 17 & 18)
% 3 is part of left part of 1 (node 17 & 18)
% 4 is part of 1 (node 17 & 18)
% 5 is part of 1 (node 17 & 18)
% 6 is either a blood vessel artifact or retrosplenial
% There's more of 12 in raw
% 15 is part of 14 (node 6 or 14)
% 16 (node 25) might have a right side in it (node 26)
% 17 (node 21) has a strong right side in it (node 22)
% There's more of 25 (maybe node 3)
% There's more of 29
% 33 could maybe be divided into node 5 & 1
% There's a right part of 36 

%%
% Mouse 1099
% Notes:
% (all ICs labeled from pre-artifacts removed)

% (done) IC 14 has a matching raw on the other side
% 16 goes along with 17 on other side
% 19 and 20 are same source
% (done) 31 has a right side version (retrosplenial/node 32)
% (done) There's more of IC 45 in raw
% (done) There's a right side of 51 (also retrosplenial)
% (done) There's more of 55. Put it with node 5
% IC 32 is node 1

add_back = cell(1);
original_IC_numbers = cell(1);
raw_ICs = cell(1);

load('Y:\Sarah\Analysis\Experiments\Random Motorized Treadmill\spatial segmentation\500 SVD components\raw ICs\1099\sources100.mat', 'sources');
load('Y:\Sarah\Analysis\Experiments\Random Motorized Treadmill\preprocessing\masks\masks_m1099.mat');
% Fill masks, convert to absolte value
sources = abs(FillMasks(sources', indices_of_mask, 256, 256));

source =  sources(:, :,13);
figure; imagesc(source);
holder = NaN(256);
holder1 = NaN(256);
holder1(140:167, 178:204) = source(140:167, 178:204); % Limit area
inds = find(holder1 > 2.5); % Threshold
holder(inds) = source(inds);
figure; imagesc(holder); colorbar; caxis([0 4]);
add_back{1} = holder;
original_IC_numbers{1} = 13;
raw_ICs{1} = source;

source =  sources(:, :,26);
figure; imagesc(source);
holder = NaN(256);
holder1 = NaN(256);
holder1(124:204, 135:167) = source(124:204, 135:167); % Limit area
inds = find(holder1 > 2.0); % Threshold
holder(inds) = source(inds);
figure; imagesc(holder); colorbar; caxis([0 4]);
add_back{2} = holder;
original_IC_numbers{2} = 26;
raw_ICs{2} = source;

source =  sources(:, :,39);
figure; imagesc(source);
holder = NaN(256);
holder1 = NaN(256);
holder1(120:170, 7:32) = source(120:170, 7:32); % Limit area
inds = find(holder1 > 3.5); % Threshold
holder(inds) = source(inds);
figure; imagesc(holder); colorbar; caxis([0 4]);
add_back{3} = holder;
original_IC_numbers{3} = 39;
raw_ICs{3} = source;

source =  sources(:, :,48);
figure; imagesc(source);
holder = NaN(256);
holder1 = NaN(256);
holder1(172:223, 138:168) = source(172:223, 138:168); % Limit area
inds = find(holder1 > 1.5); % Threshold
holder(inds) = source(inds);
figure; imagesc(holder); colorbar; caxis([0 4]);
add_back{4} = holder;
original_IC_numbers{4} = 48;
raw_ICs{4} = source;

% More of node 9
source =  sources(:, :,14);
figure; imagesc(source);
holder = NaN(256);
holder1 = NaN(256);
holder1(46:106, 82:120) = source(46:106, 82:120); % Limit area
inds = find(holder1 > 3.5); % Threshold
holder(inds) = source(inds);
figure; imagesc(holder); colorbar; caxis([0 4]);
add_back{5} = holder;
original_IC_numbers{5} = 14;
raw_ICs{5} = source;

add_back_all{8} = add_back;
original_IC_numbers_all{8} = original_IC_numbers;
raw_ICs_all{8} = raw_ICs;

%% Actually add back these ICs to artifacts removed matrices
% Applies CLeanClust function
% load mice_all
load('Y:\Sarah\Analysis\Experiments\Random Motorized Treadmill\mice_all.mat');
% App
% For each mouse, 
for mousei =  1:size(mice_all,2)
    mouse = mice_all(mousei).name;

    % If add_back_all isn't empty for that mouse
    if ~isempty(add_back_all{mousei})
        
        % Load artifacts removed.
        load(['Y:\Sarah\Analysis\Experiments\Random Motorized Treadmill\spatial segmentation\500 SVD components\artifacts removed conditional thresholding\pre addback\' mouse '\sources.mat']);

        % holder for overlay.
        overlay = sources.overlay;

        % Get maximum value already in overlay.
        number_ics = max(max(overlay));

        % For each IC added back in that mouse
        for ici = 1: numel(add_back_all{mousei})
            
            % If not empty
            if ~isempty(add_back_all{mousei}{ici})
    
                % Update number of elements in indices to remove, artifact masks
                sources.indices_to_remove = [sources.indices_to_remove; cell(1)];
                sources.artifact_masks = [sources.artifact_masks; cell(1)];
        
                % Update original IC number. 
                % If empty or has more than 1 number, 
                if isempty(original_IC_numbers_all{mousei}{ici}) || numel(original_IC_numbers_all{mousei}{ici}) > 1
                    sources.originalICNumbers = [sources.originalICNumbers; NaN];
                else 
                    sources.originalICNumbers = [sources.originalICNumbers; original_IC_numbers_all{mousei}{ici}];
                end 

                % Threshold & clean
                new_source = add_back_all{mousei}{ici} > 0;
                [new_source, ~] = CleanClust(new_source);

                % Update overlay. 
            
                    % Add new source to overlay
                    overlay(find(new_source > 0)) = ici + number_ics;

                % Add the source itself in. 
                   
                    % Multiply mask by raw source
                    try
                    color_mask = new_source .* raw_ICs_all{mousei}{ici};
                    catch 
                       disp('here');
                    end 
                    % Concatenate
                    sources.sources = cat(3, sources.sources, color_mask);

            end
        end 

        % Rename overlay
        sources.overlay = overlay;

        % Save
        mkdir(['Y:\Sarah\Analysis\Experiments\Random Motorized Treadmill\spatial segmentation\500 SVD components\artifacts removed conditional thresholding\post addback staging\' mouse '\']);
        save(['Y:\Sarah\Analysis\Experiments\Random Motorized Treadmill\spatial segmentation\500 SVD components\artifacts removed conditional thresholding\post addback staging\' mouse '\sources.mat'], 'sources');
    
        % *** Now do it all for the regularized folder
        load(['Y:\Sarah\Analysis\Experiments\Random Motorized Treadmill\spatial segmentation\500 SVD components\regularized ICs 150 amp 3.5 two conditionals small 2.5 large 1\' mouse '\sources100.mat']);
        
        % For each IC added back in that mouse
        for ici = 1: numel(add_back_all{mousei})
            
            % If not empty
            if ~isempty(add_back_all{mousei}{ici})
    
                % Update original IC number. 
                % If empty or has more than 1 number, 
                if isempty(original_IC_numbers_all{mousei}{ici}) || numel(original_IC_numbers_all{mousei}{ici}) > 1
                    sources.originalICNumber_domainsSplit = [sources.originalICNumber_domainsSplit NaN];
                else 
                    sources.originalICNumber_domainsSplit = [sources.originalICNumber_domainsSplit original_IC_numbers_all{mousei}{ici}];
                end 
                
                % Threshold and clean. 
                new_source = add_back_all{mousei}{ici} > 0;
                [new_source, ~] = CleanClust(new_source);
                
                % Update colormask
                color_mask = new_source .* raw_ICs_all{mousei}{ici};
                
                % Concatenate
                sources.color_mask_domainsSplit = cat(3, sources.color_mask_domainsSplit, color_mask);
            end 
        end

        
        mkdir(['Y:\Sarah\Analysis\Experiments\Random Motorized Treadmill\spatial segmentation\500 SVD components\regularized ICs 150 amp 3.5 two conditionals small 2.5 large 1\post addback\' mouse '\']);
        save(['Y:\Sarah\Analysis\Experiments\Random Motorized Treadmill\spatial segmentation\500 SVD components\regularized ICs 150 amp 3.5 two conditionals small 2.5 large 1\post addback\' mouse '\sources100.mat'], 'sources');
    end
end 

%% Make overlay plots of these new additions to artifact- removed overlays to
% check progress (DON'T SAVE YET)
% close all;
% 
% % load mice_all
% load('Y:\Sarah\Analysis\Experiments\Random Motorized Treadmill\mice_all.mat');
% 
% mice_all = mice_all([1:3, 5, 6]);
% figure; 
% % For each mouse, 
% for mousei = 1:size(mice_all,2)
%     mouse = mice_all(mousei).name;
% 
%     % If add_back_all isn't empty for that mouse
%     if ~isempty(add_back_all{mousei})
%         
%         % Load artifacts removed overlay. 
%         load(['Y:\Sarah\Analysis\Experiments\Random Motorized Treadmill\spatial segmentation\500 SVD components\artifacts removed conditional thresholding\post addback staging\' mouse '\sources.mat']);
% 
%         % holder for overlay.
%         overlay = sources.overlay;
% 
%         % Get maximum value already in overlay.
%         number_ics = max(max(overlay));
% 
%         % For each IC added back in that mouse
%         for ici = 1: numel(add_back_all{mousei})
%             
%             % If not empty
%             if ~isempty(add_back_all{mousei}{ici})
% 
%             
% 
%                 % Threshold & overlay 
%                 new_source = add_back_all{mousei}{ici} > 0;
%                 [new_source, ~] = CleanClust(new_source);
%                
%                 % Add new source to overlay
%                 overlay(find(new_source > 0)) = ici + number_ics;
%             end
% 
%         end 
% 
%         % Plot new overlay
%         subplot(2,3, mousei); imagesc(overlay); axis square; title(mouse);
%         xticks([]); yticks([]); colorbar;
%     end
% end
% 
% % Save figure;
% savefig('Y:\Sarah\Analysis\Experiments\Random Motorized Treadmill\spatial segmentation\added back.fig');
% 

%% Plot overlays from post addback staging

% load mice_all
load('Y:\Sarah\Analysis\Experiments\Random Motorized Treadmill\mice_all.mat');

mice_all = mice_all([1:6, 8]);
figure; 
% For each mouse, 
for mousei = 1:size(mice_all,2)
    mouse = mice_all(mousei).name;
    
    % Load artifacts removed overlay. 
    load(['Y:\Sarah\Analysis\Experiments\Random Motorized Treadmill\spatial segmentation\500 SVD components\artifacts removed conditional thresholding\post addback staging\' mouse '\sources.mat']);
    subplot(3,3,mousei); 
    imagesc(sources.overlay); title(mouse);
end 