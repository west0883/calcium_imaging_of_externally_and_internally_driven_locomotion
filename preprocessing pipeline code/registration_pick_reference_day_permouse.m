% registration_pick_reference_day_permouse.m
% Sarah West
% 8/23/21

% Cyles through all the mice, calls pick_reference_day, saves the day the
% user picks for registration. Called by main pipeline code

function []=registration_pick_reference_day_permouse(parameters)
   
    % Assign parameters their original names
    dir_exper = parameters.dir_exper; 
    plot_sizes = parameters.plot_sizes; 
    
    % Load the ORIGINAL mice_all 
    load([dir_exper 'mice_all.mat']); 
    mice_all_original = mice_all;
    
    % make the output directory
    dir_out=[dir_exper 'representative images\']; 
    mkdir(dir_out); 
    disp(['data saved in ' dir_out]); 
    
    % See if the reference days have already been started and saved
    if exist([dir_out 'reference_days.mat'])
       
        % if it does, load it
        load([dir_out 'reference_days.mat']); 
    else  
        % if not, create a new variable for holding the reference days 
        reference_days.mouse=cell(size(mice_all_original,2),1);
        reference_days.day=cell(size(mice_all_original,2),1);
        
    end     
    % for each mouse 

    % Use only the mice & days you've been instructed to use, but keep
    % mouse and day indexes accurate to the original mice_all list.


    for mousei=1:size(parameters.mice_all,2)
        mouse=parameters.mice_all(mousei).name;
        
        % Get the "original" mouse index 
        mousei_original = find(strcmp({mice_all.name}, mouse)==1);
        
        % Use a flag for determining if the "pick_reference_day"
        % function should be run.
        pick_flag=0;         
 
        % see if it's been found before
        if size(reference_days.mouse, 1) >= mousei_original && ~isempty(reference_days.mouse{mousei_original}) % if it HAS been found before
            
            % ask user if they want to redo-it
            user_answer = inputdlg(['Would you like to re-find the reference day for mouse' mouse '? (y = yes, n = no)']);
            % get rid of cell formatting
            user_answer=user_answer{1};
            % if they don't want to,
            if user_answer=='n'
                pick_flag=0; % don't run functon
            % if they DO want to,     
            elseif user_answer=='y'
                pick_flag=1;  % run function          
            end
            
        % if it hasn't been found before    
        else
            pick_flag=1; % run function 
        end
        
        if pick_flag==1 
            % find the directories and images for the pick_reference_day function
            % holding variable for the images to show
            file_paths=cell(size(parameters.mice_all(mousei).days,1),1); 

            % put the file paths of each background image into the list
            days_list=vertcat(parameters.mice_all(mousei).days(:).name);
            for dayi=1:size(days_list,1)
                day=parameters.mice_all(mousei).days(dayi).name;
                file_paths{dayi}=[dir_exper 'representative images\' mouse '\' day '\bRep.mat'];

            end 
        
            % Run the pick_reference_day function
            [dayi_output]=PickReferenceDay(mouse, file_paths, plot_sizes);

            % Convert the dayi user input to the day name and put in the
            % variable
            reference_days.mouse{mousei_original}=mouse; 
            reference_days.day{mousei_original}=parameters.mice_all(mousei).days(dayi_output).name;
       
            % Load that bRep so you can save it.
            load(file_paths{dayi_output});

            reference_image = bRep; 

            % Save that image matrix
            save([dir_out mouse '/reference_image.mat'], 'reference_image')
        end 
        
        % Save the user outputs/reference days as a variable. Save each
        % time you look at a mouse, because we don't want something to go
        % wrong in the middle and have to start over.
        save([dir_out 'reference_days.mat'], 'reference_days'); 

    end     
    close all;
end 