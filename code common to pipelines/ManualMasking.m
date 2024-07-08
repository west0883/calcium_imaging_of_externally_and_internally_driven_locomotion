% ManualMasking.m
% Sarah West
% 8/23/21

% Is a function that lets you draw masks over mice brains manually. Is
% called by manual_masking_loop.m 

% Inputs:
% existing_masks-- 
% image_to_mask-- 
% iterator_name -- just for asking user if they're done.

% Outputs: 
% masks-- a 

function[masks, indices_of_mask]=ManualMasking(image_to_mask, existing_masks, axis_for_drawing, flip)
    
    % If no input for axes, set default axes to current axes.
    if nargin < 3
        axis_for_drawing = gca;
    end

    % If no input for flipping or not, set flip to false (because should
    % now be plotting on the axis given by higher-level function).
    if nargin < 4
        flip = false;
    end
    % Apply any existing masks to image_to_mask
    indices_of_mask=[];
    
    if isempty(existing_masks)
       % If no previous masks, don't need to add anything to image_to_mask 
    else
        % If there are pre-existing masks, apply them to
        % image_to_mask
        for i=1:size(existing_masks,3)
            mask_flat=existing_masks(:,:,i);
            indices_of_mask=[indices_of_mask; find(mask_flat)]; 
        end
        image_to_mask(indices_of_mask)=NaN;
    end 

    % Set the axes to draw on as the passed axes handle.
    axes(axis_for_drawing); 

    % Display image_to_mask. Displays masks as gray (not black, can't see
    % crosshairs).
    mymap=[0.5 0.5 0.5; parula(512)];
    imagesc(image_to_mask); colormap(axis_for_drawing, mymap);
    xticks([]); yticks([]); axis square; 
    
    % Set input dialogue options -- allow for user to still interact with
    % figures.
    opts.WindowStyle = 'normal';

    % Ask user if they want to add a mask
    user_answer1= inputdlg(['Do you want to draw additional masks on this image? y = yes, n = no'], 'User input', 1,{'n'}, opts); 

    %Convert the user's answer into a value
    answer1=user_answer1{1};

    % Include existing_masks in the list of masks you're making 
    masks=existing_masks; 

    % If the user said yes,
    while strcmp(answer1, 'y')      

        % Run function "PolyDraw" on the image with previous masks; will 
        % output the coordinates of the ROI drawn
        ROI1=PolyDraw(axis_for_drawing, 1); 
        
        % Make a mask of the ROI drawn, depends on if call said to flip.
        if flip
            mask1=flipud(poly2mask(ROI1(1,:),ROI1(2,:),size(image_to_mask,1), size(image_to_mask,2)));    
        else 
            mask1 = poly2mask(ROI1(1,:),ROI1(2,:),size(image_to_mask,1), size(image_to_mask,2));
        end 

        % Add the new mask to the matrix of masks that have been drawn so far
        masks=cat(3, masks, mask1); 

        % Find the indices of the new mask
        indices_new=find(mask1);

        % Set the new mask to NaN for display
        image_to_mask(indices_new)=NaN; 
        
        % Add these new indices to the list of mask indices
        indices_of_mask=[indices_of_mask; indices_new];
        
        % Set the axes to draw on as the passed axes handle (again).
        axes(axis_for_drawing);

        % Draw the representative image with the new masks on it
        imagesc(image_to_mask); colormap(axis_for_drawing, mymap); 
        xticks([]); yticks([]); axis square;

        % Repeat
        user_answer1= inputdlg(['Do you want to draw additional masks on this image? y = yes, n = no'], 'User input', 1,{'n'}, opts); 
        answer1=user_answer1{1};
    end

    % leaves while loop if user answers anything other than "1"

end