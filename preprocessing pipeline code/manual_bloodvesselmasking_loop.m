% manual_bloodvesselmasking_loop.m
% Sarah West
% 9/14/21
% Modified from SCA_manual_masking.m, just the manual masking parts. Will
% also put in code from create_bloodvessel_masks_withbackgroundmasks.m . 
% Let's you keep adding masks until you say you're done. You can return to
% it. 

% This has the exact same code as manual_masking_loop, just changes the
% output directories. 

function []=manual_bloodvesselmasking_loop(parameters)
    
   % Assign parameters their original names
   dir_exper = parameters.dir_exper; 
   mice_all = parameters.mice_all; 
   channelNumber = parameters.channelNumber; 
   
   % Make a flag that says if you should continue or not, default to "don't
   % continue." --> This flag is mostly to help make sure I didn't overlook
   % anything in my logic.
   continue_flag=0; 

   % Be sure the user knows this is for blood vessels, not brain masks.
    user_answer1= inputdlg('WARNING: This is for BLOOD VESSEL MASKS. Do you want to continue? 1=Y, 0=N'); 
   
     answer1=str2num(user_answer1{1});
             
     % If user says yes
     if answer1==1
         
         % Do nothing, continue function as normal. 
         
     % If user answers anything else    
     else
         % Leave function
         return 

     end 
    
    switch channelNumber
           case 2
            
            % If there are two channels, warn user that they probably don't
            % need this step. Ask them if they still want to do it anyway.  
             user_answer1= inputdlg('WARNING: You have a violet channel and don''t have to find blood vessel masks...Would you like to draw some anyway? 1=Y, 0=N'); 
             answer1=str2num(user_answer1{1});
             
             % If user says yes
             if answer1==1
                 
                 % Make function run normally by making "continue_flag" = 1
                 continue_flag=1; 
                 
             % If user answers anything else    
             else
                 % Leave function
                 return 
                 
             end 
        
        case 1
            continue_flag=1; 
            
    end
    
    % If there's the right number of channels/the user said make blood
    % vessel masks anyway,.
    if continue_flag==1
        % Establish input and output folders 
        dir_in_base=[dir_exper 'representative images\'];
        dir_out=[dir_exper 'blood vessel masks\']; 
        mkdir(dir_out); 

        % Display where data is being saved for user
        disp(['data saved in ' dir_out]); 

        % Load reference days
        load([dir_in_base '\reference_days.mat']); 

        % Cycle through mice based on the willingness of the user
        mousei=1; 
        while mousei <= size(mice_all,2) 

            % Find the mouse name
            mouse=mice_all(mousei).name;

            % Display which mouse you're working on
            disp(['working on mouse ' mouse]); 

            % Find the reference day
            reference_day=reference_days.day{mousei};

            % Define input folder based on reference day
            dir_in=[dir_in_base mouse '\' reference_day '\'];

            % Load that mouse's Reference bRep
            load([dir_in '\bRep.mat']);
            
            % Get dimensions
            yDim=size(bRep,1);
            xDim=size(bRep,2); 
            
            % Add the brain mask to the bRep image to make it clear not to
            % draw masks outside of that. Then apply the brain mask to the
            % vessel masks to get the correct number of indices.
            
            % Load brain mask 
            load([dir_exper 'masks\masks_m' mouse '.mat'], 'indices_of_mask'); 
            
            % Rename brain mask indices to avoid confusion/overwriting. 
            brain_mask_indices=indices_of_mask; 

            % Get the inverse of the brain mask.
            inverse_brain_mask_indices=true(yDim, xDim); 
            inverse_brain_mask_indices(brain_mask_indices)=false;
            
            % Apply brain masks to the bRep image 
            bRep(inverse_brain_mask_indices)=NaN;             
            
            % Determine if a mask file for this mouse already exists.

            % If it does exist already, load the mask file
            if isfile([dir_out 'bloodvessel_masks_m' mouse '.mat'])
               load([dir_out 'bloodvessel_masks_m' mouse '.mat']); 
            % If it doesn't exist, 
            else
                % Make a starting masks variable that's empty
                vessel_masks=[];
            end 

            % Rename existing masks so they're not confused with the new ones
            % that will be drawn.
            existing_masks=vessel_masks;
            
            % ***Run the function that runs the masking itself***
            [vessel_masks, ~]=ManualMasking(bRep, existing_masks);     

            % Save whatever additions you've made to the mask file 
            save([dir_out 'bloodvessel_masks_m' mouse '.mat'], 'vessel_masks');

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