% motor_FindBehaviorPeriods_all.m
% Sarah West
% 12/6/21

% Cycles through mice, days, stacks and calls motorFindBehaviorPeriods
% function. Saves the behavior periods.

function [] = motor_FindBehaviorPeriods_all(parameters)
  
   % Establish base input directory
    dir_in_base=[parameters.dir_exper 'behavior\extracted motor data\'];

    % Establish base output directory
    dir_out_base=[parameters.dir_exper 'behavior\period instances\'];
    
    % Tell user where data is being saved
    disp(['Data saved in ' dir_out_base]); 
    
    % For each mouse 
    for mousei=1:size(parameters.mice_all,2)
        mouse=parameters.mice_all(mousei).name;
        
        % For each day
        for dayi=1:size(parameters.mice_all(mousei).days, 2)
            
            % Get the day name.
            day=parameters.mice_all(mousei).days(dayi).name; 
            
            % Display
            disp(['mouse ' mouse ', day ' day]);
            
            % Find input directory and cleaner output directory. 
            parameters.dir_in = [dir_in_base  mouse '\' day '\'];
            parameters.input_data_name = {'trial', 'stack number', '.mat'};
            dir_out=[dir_out_base mouse '\' day '\']; 
            mkdir(dir_out); 
            
            % Get the stack list
            [stackList]=GetStackList(mousei, dayi, parameters);
            
            % Cycle through the stack files. 
            for stacki=1:size(stackList.filenames,1)
                %disp(stacki);
                
                % Get the stack number and filename for the stack.
                stack_number = stackList.numberList(stacki, :);
                filename = stackList.filenames(stacki, :);
                
                 %disp(['mouse ' mouse ', day ' day ', stack ' stack_number]);
                % Load corresponding data
                load([parameters.dir_in filename]);
                
                % Run motor periods function
                [all_periods] = motorFindBehaviorPeriods(trial, parameters);
                
                % Save
                save([dir_out 'all_periods_' stack_number '.mat'], 'all_periods'); 
            end 
        end
    end
end