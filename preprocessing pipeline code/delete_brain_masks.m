% delete_brain_masks.m
% Sarah West
% 12/4/21

% Called from preprocessing pipeline. Calls on DeleteMasks.m to delete
% brain masks. 

function [] = delete_brain_masks(parameters)

    
    % Establish input and output folders 
    dir_in_rep_base=[parameters.dir_exper 'representative images\'];
    dir_in=[parameters.dir_exper 'masks\']; 

    % Display where data is being saved for user
    disp(['data saved in ' dir_in]); 
    
    % Load reference days
    load([dir_in_rep_base '\reference_days.mat']); 
        
    % Cycle through mice based on the willingness of the user
    mousei=1; 
    while mousei <= size(parameters.mice_all,2) 
        
        % Find the mouse name
        mouse=parameters.mice_all(mousei).name;
        
        % Display which mouse you're working on
        disp(['working on mouse ' mouse]); 
        
        % Find the day you're supposed to register to with this mouse 
        ind = NaN(1,size(reference_days.mouse,1)); 
        for i=1:size(reference_days.mouse,1)
           ind(i)=strcmp(mouse, reference_days.mouse{i}); 
        end
        refdayi=find(ind); 
        reference_day=reference_days.day{refdayi};
        
        % Define input folder based on reference day
        dir_in_rep=[dir_in_rep_base mouse '\' reference_day '\'];
         
        % Load that mouse's Reference bRep
        load([dir_in_rep '\bRep.mat']);
        
        % Check the size of the bRep, cut to size if needed. 
        bRep = FixImageSize(bRep, parameters.pixels); 

        yDim = parameters.pixels(1);
        xDim = parameters.pixels(2);
        
        %  % Determine if a mask file for this mouse already exists.
        existing_mask_flag=isfile([dir_in 'masks_m' mouse '.mat']); 
        
        % If it does exist already, load the mask file
        if existing_mask_flag==1 
           load([dir_in 'masks_m' mouse '.mat']); 
            
        % If it doesn't exist, announce that & go to next mouse. 
        elseif existing_mask_flag==0 
            disp(['No masks for mouse ' mouse ]); 
             mousei = mousei + 1;
            continue 
        end 
        
        % Run deletion protocol.
        [masks, indices_of_mask] = DeleteMasks(masks, indices_of_mask, bRep);
        
        % Save whatever additions you've made to the mask file 
        save([dir_in 'masks_m' mouse '.mat'], 'masks', 'indices_of_mask');
       
        % clear things for next mouse 
        close all; 
        
        % Ask if the user wants to keep working
        user_answer1= inputdlg('Do you want to work on the next mouse? 1=Y, 0=N'); 
        answer1=str2num(user_answer1{1});
        
        % If the user says yes,
        if answer1==1
            % Increase the valuse of mousei and continue 
            mousei=mousei+1; 
        else
            % If the user says anything else, break the while loop so
            % another mouse isn't started 
            break
        end
        
    end 

end 