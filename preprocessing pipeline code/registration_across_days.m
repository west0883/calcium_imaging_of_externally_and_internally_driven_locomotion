% registration_across_days.m
% Sarah West
% 8/23/21
% Modified from SCA_manual_masking.m, just the registration across days
% parts

function []=registration_across_days(parameters)
    
    % Assign parameters their original names
    dir_exper = parameters.dir_exper; 
    mice_all = parameters.mice_all; 
    transformation = parameters.transformation; 
    configuration = parameters.configuration; 
    max_step_length = parameters.max_step_length;
    max_iterations = parameters.max_iterations; 
    
    % Establish input and output directories
    dir_in_base=[dir_exper 'representative images\'];
    dir_out_base=[dir_exper 'tforms across days\'];
    
    % Display to user where everything is saved.
    disp(['data saved in ' dir_out_base]); 
    
    % Load the list of what days should be used to register everything else
    % to. 
    load([dir_in_base 'reference_days.mat']);
    
    % Find parameters for the registration you want to do
    [optimizer, metric] = imregconfig(configuration);
    
    % Change the default max step length and max iterations. 
    optimizer.MaximumStepLength=max_step_length;
    optimizer.MaximumIterations=max_iterations;
  
    
    % **Compute t-forms** 
    % For each mouse 
    for mousei=1:size(mice_all,2)
        mouse=mice_all(mousei).name;
        
        % Get the list of days for that mouse
        days_list=vertcat(mice_all(mousei).days(:).name); 
        
        % Find the day you're supposed to register to with this mouse 
        ind = NaN(1,size(reference_days.mouse,1)); 
        for i=1:size(reference_days.mouse,1)
           ind(i)=strcmp(mouse, reference_days.mouse{i}); 
        end
        refdayi=find(ind); 
        reference_day=reference_days.day{refdayi};
        
        % Load the reference bRep image, rename it
        load([dir_in_base '\' mouse '\' reference_day '\bRep.mat']);
        Reference_bRep=bRep;
        
        % for each day
        for dayi=1:size(days_list,1)
            day=days_list(dayi,:); 
            
            % make a dir_in and dir_out folder name for each day
            dir_in = [dir_in_base mouse '\' day '\']; 
            dir_out= [dir_out_base mouse '\' day '\']; 

            % Create the output folder 
            mkdir(dir_out); 

            % See if this day is the reference day
            if strcmp(day, reference_day)

                % If this is the reference day, make tform empty
                tform=[];

            else
                % If this is NOT the reference day, perform the
                % registration

                % Load day's bRep
                load([dir_in 'bRep.mat']);

                % Perform registration.                 
                 tform = imregtform(bRep, Reference_bRep, transformation, optimizer, metric);

                % Perform a check and save it in the folder
                    % Apply the transform to the bRep
                    result=imwarp(bRep,tform,'OutputView',imref2d(size(Reference_bRep))); 

                    % Plot both images together before and after registration 
                    figure; 
                    subplot(1,2,1); imshowpair(bRep, Reference_bRep); title('before')
                    subplot(1,2,2); imshowpair(result,Reference_bRep); title('after')
                    sgtitle([mouse ', ' day])
                    % Save the check figure 
                    savefig([dir_out 'before_and_after.fig']);  
            end
            % Save the tform for each day (including empy tform variables, 
            % which makes the logic easier in later steps) 
            save([dir_out 'tform.mat'], 'tform');         
            
        end 
    end
end