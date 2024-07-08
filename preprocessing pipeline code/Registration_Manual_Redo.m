% Registration_Manual_Redo.m
% Sarah West
% 8/24/21

% Takes the across-day registrations that weren't done very well by the automatic
% registration and let's you re-do them manually. Adapted from
% "lines_register_multiple_forbackwards.m"


function []=Registration_Manual_Redo(redo, parameters)
    
    % Assign parameters their original names
    dir_exper = parameters.dir_exper;
    
    % Set up bases of input and output directories
    dir_in_base=[dir_exper 'representative images\'];
    dir_out_base=[dir_exper 'tforms across days\'];
   
    % Display to user where everything is saved.
    disp(['data saved in ' dir_out_base]); 
    
    % Load the list of what days should be used to register everything else
    % to. 
    load([dir_in_base 'reference_days.mat']);
    
    % For each day that needs to be re-done
    for redoi=1:size(redo,1)
        
        % Find the mouse and day 
        mouse=redo{redoi,1};
        day=redo{redoi,2}; 
        
        % Create dir_out name  
        dir_out=[dir_out_base mouse '\' day '\']; 
        
        % Find the day you're supposed to register to with this mouse
        mousei=find(contains(reference_days.mouse, mouse));
        reference_day=reference_days.day{mousei};
        
        % Load+rename Reference_bRep
        load([dir_in_base mouse '\' reference_day '\bRep.mat']); 
        Reference_bRep=bRep;
        
        % Load bRep for current day 
        load([dir_in_base mouse '\' day '\bRep.mat']); 
        
        % Plot each image
        fig1=figure; imagesc(Reference_bRep); 
        fig2=figure; imagesc(bRep)
        
        % Use the function "getpts" to slelect 3 or more points from the first image 
        % (I like to use blood vessel branch points), then the same 3 points 
        [x1, y1]=getpts(fig1);
        [x2, y2]=getpts(fig2); 
        
        % Calculate the matrix (tform) that, when multiplied with the image you're 
        % registering, aligns it with the reference image
        tform = fitgeotrans([x2 y2], [x1 y1], 'affine'); 
        
        % Save the new tform.
        save([dir_out 'tform.mat'], 'tform'); 

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
end 