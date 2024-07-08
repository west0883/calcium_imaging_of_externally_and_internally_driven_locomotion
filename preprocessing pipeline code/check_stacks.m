% check_stack.m
% Sarah West
% 8/27/21
% Checks if a preprocessed stack is corrupt

function [] = check_stacks(parameters) 
    
    % Assign parameters their original names
    input_data_name = parameters.input_data_name;
    dir_exper = parameters.dir_exper; 
    mice_all = parameters.mice_all; 
    digitNumber = parameters.digitNumber; 
   
    % Set input (and re-preprocessing output) directory base.
    dir_preprocessed_base = parameters.dir_in_base; 
    
    % Create output file name
    output_filename='corrupt_stacks.mat';

    % Tell user what you're doing.
    disp(['Finding corrupt stacks. Saving data in ' dir_preprocessed_base]);
    
    % Initiate an empty cell array
    corrupt_files={};
    
    % Cycle through mice. 
    for mousei=1:size(mice_all,2)
        
        % Get the mouse name.
        mouse=mice_all(mousei).name; 

        % Cycle through days.
        for dayi=1:size(mice_all(mousei).days, 2)
            
            % Get the day name.
            day=mice_all(mousei).days(dayi).name; 
            
            % Tell user where the code is.
            disp(['mouse ' mouse ', ' day ]); 
            
            % Establish input directory. 
            dir_preprocessed=[dir_preprocessed_base '\' mouse '\' day '\'];
            
            % Find the correct stack list entry of mice_all. 
            stackList=mice_all(mousei).days(dayi).stacks; 
            
            % If stackList is a character string (to see if 'all')
            if ischar(stackList)
        
               % If it is a character string, check to see if it's the string
               % 'all'. 
               if strcmp(stackList, 'all')
                   
                   % If it is the character string 'all', list stacks from
                   % the day directory. 
                   list=dir([dir_preprocessed 'data*.mat']);
                   
                   % Assign a flag for marking if this happened. 
                   all_flag=1; 
               end
               
               % If stackList is not a character string, assume it's a
               % vector of integer stacknumbers.  
            else
                list=ListStacks(stackList, digitNumber); 
                
                % Assign a flag for marking if the list was already made
                all_flag=0;
            end 
            
            % For each stack 
            for stacki=1:size(list,1)
                
                % Depending on the value of the all_flag (the format of the
                % stack list), figure out the name of the stacks
                switch all_flag
                    
                    % If the user said use 'all' stacks,
                    case 1
                        
                        % Get the stack number assuming the filename is
                        % 'data' [stack_number] '.mat'.
                        stack_number=list(stacki).name(5:5+digitNumber-1); 
                          
                        % Get the filename of the stack. 
                        filename=[dir_preprocessed list(stacki).name];
                    
                    % If the user said use a specifc list of stacks, use
                    % the list generated from the ListStacks function. 
                    case 0
                        % Get the stack number.
                        stack_number=list(stacki, :); 
                        
                        % Get the filename.
                        stackname=CreateFileStrings(input_data_name, [], [], stack_number, [],  false); 
                        filename=[dir_preprocessed stackname];
                end 
                
               % Run filename in CheckIfCorruptFast 
               [isCorrupt]=CheckIfCorruptFast(filename);
              
               % If the file was corrupt, 
               if isCorrupt==true
                   % Add to list of corrupt stacks.
                   corrupt_files=[corrupt_files; {mouse} {day} {stack_number}]; 

                   % Tell user there was a corrupt stack. 
                   disp('Found a corrupt stack.');
                   disp(['Preprocesing ' mouse ', ' day ' stack ' stack_number ]); 
               
                   % Run the re-do preprocessing here. 
                   
                   % Make a new parameters structure with only the needed
                   % stacks in it. 
                   
                   mice_all_redo(1).name=mouse;
                   mice_all_redo(1).days(1).name=day;
                   mice_all_redo(1).days(1).stacks=str2num(stack_number);
                   
                   parameters_redo=parameters;
                   parameters_redo.mice_all=mice_all_redo;
                   Preprocessing(parameters_redo);
               end
            end
        end 
    end 

% Save list of corrupt stacks, just in case.  
save([dir_preprocessed_base output_filename], 'corrupt_files');

end 