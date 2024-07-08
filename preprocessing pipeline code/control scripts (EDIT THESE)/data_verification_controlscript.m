% data_verification_controlscript.m
% Sarah West
% 5/28/22

clear all; 

% Create the experiment name. This is used to name the output folder. 
parameters.experiment_name='Random Motorized Treadmill';

% Output directory name bases
parameters.dir_base='Y:\Sarah\Data\';
parameters.dir_exper=[parameters.dir_base parameters.experiment_name '\'];

% Load in ROIs
load('Y:\Sarah\Analysis\Experiments\Random Motorized Treadmill\data_verification_rois.mat');

% Set up iterators
loop_variables.days = {'051222', '051322', '051822', '051922', '052022', '052522', '052622', '052722'};
loop_variables.mice = {'m539', 'm1099'};

loop_list.iterators = {'day', {'loop_variables.days'}, 'day_iterator';
                       'mouse', {'loop_variables.mice'}, 'mouse_iterator'};

stacks_directory_format = {parameters.dir_exper, 'day', '\', 'mouse', '\stacks\'};

looping_output_list = LoopGenerator(loop_list, loop_variables);

for itemi = 1:size(looping_output_list, 1)
    
    day = looping_output_list(itemi).day;
    mouse = looping_output_list(itemi).mouse;
    disp([day ', ' mouse]);

    directory_name = CreateStrings(stacks_directory_format, {'day', 'mouse'}, {day, mouse});

    % Run
    Data_verification_v3_short_Automated(directory_name, roi_1_mask, roi_2_mask);

end 
