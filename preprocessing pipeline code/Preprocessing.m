% preprocessing.m
% Sarah Wests
% 8/24/21

% This is happening all in the same code to reduce the amount of data being 
% saved in intermediate processing steps. THIS CODE DOES:

% 1. Reads the data tiffs in as matrics.
% 2. Separates data into blue and violet channels. 
% 3. Registers within-stack/across stacks within a day. 
% 4. Applies the pre-calculated across-day tforms.
% 5. Applies the pre-calculated mask per mouse.
% 6. Corrects hemodynamics via the user's chosen method.
% 7. Apples filtering. 
% 8. Saves preprocessed stacks. 

function [parameters] = Preprocessing(parameters)
    
    % Assign parameters their original names 
    channelNumber = parameters.channelNumber;
    skip = parameters.skip;
    pixel_rows = parameters.pixel_rows;
    pixel_cols = parameters.pixel_cols; 
    filter_flag = parameters.filter_flag;
    b = parameters.b;
    a = parameters.a;
    usfac = parameters.usfac; 
    frames_for_spotchecking = parameters.frames_for_spotchecking;
    minimum_frames = parameters.minimum_frames;
    mask_flag = parameters.mask_flag;
    correction_method = parameters.correction_method;
    im_list = parameters.im_list;
    tform = parameters.tform;

    % Make a cell array of false (default) "don't save" values, one for
    % each output
    parameters.dont_save = repmat({false}, numel(fieldnames(parameters.loop_list.things_to_save)), 1); 
    
    % Establish a holding cell array to keep track of what days need to be
    % checked on.
    if ~isfield(parameters, 'bad_trials')
        bad_trials = {};
    else 
        bad_trials = parameters.bad_trials;
    end

    % Tell user what iteration function is on
    MessageToUser('Preprocessing ', parameters);

    % *** 1. Read in tiffs.***
    disp('Reading tiffs');  
    
    % Find the sizes of the data
    yDim = size(im_list(1).data,1);
    xDim = size(im_list(1).data,2);
    
    % Get total number of images
    nim = size(im_list,2); 

    % Figure out if this frames number is long enough for
    % further processing. If not, quit this stack. 
    if (nim - skip)/channelNumber < minimum_frames
        warning('This stack is too short-- will not be processed.');
        
        bad_trials = [bad_trials; {[dir_in filename], 'too short'}];
        parameters.bad_trials = bad_trials;

        % Make all "don't save" flags true...
        dont_save = repmat({true}, numel(fieldnames(parameters.loop_list.things_to_save)), 1); 

        % .. but still save the bad_trials
        dont_save{find(strcmp(fieldnames(parameters.loop_list.things_to_save), 'bad_trials'))} = false; 
      
        parameters.dont_save = dont_save;

        % Go to next stack
        return 
    end 
    
    % ***2. Separate Channels***
    
    % Only if number of channels is 2. 
    switch channelNumber
        case 2
        disp('Separating channels'); 

        % Pick 2 images after the skip to compare 
        im1 = im_list(skip+1).data; 
        im2 = im_list(skip+2).data;

        % Figure out which is what channel
        [first_image_channel] = DetermineChannel(parameters.blue_brighter, im1, im2, pixel_rows, pixel_cols);

        % Make two lists of which images are what channel.
        switch first_image_channel

            % If the first image is blue
            case 'b'
                % Then assign every other frame starting with "skip"
                % to the blue list, the others to the violet list.
                sel470 = skip+1:2:nim;
                sel405 = skip+2:2:nim;

            % If the first image is violet. 
            case'v' 
                % Then assign every other frame starting with
                % "skip"+1 to the blue list, the others to the violet list.
                sel470 = skip+2:2:nim; 
                sel405 = skip+1:2:nim;
        end

        % Find the minimum stack length of the two channels; make this the "frames" number 
        frames = min(length(sel470),length(sel405));

        % Limit the frame indices for each color stack to the 
        % minimum number of indices (takes care of uneven image 
        % numbers by making them same length).
        sel470 = sel470(1:frames); 
        sel405 = sel405(1:frames);
        
        % Put respective channels into own data matrics
        bData = cell2mat({im_list(sel470).data});
        bData = reshape(bData, yDim, xDim, []);

        vData = cell2mat({im_list(sel405).data});
        vData = reshape(vData, yDim, xDim, []);
       
        % Set aside images for spotcheck 
        spotcheck_data.initial.blue = bData(:,:, frames_for_spotchecking);
        spotcheck_data.initial.violet = vData(:,:, frames_for_spotchecking);
        
    case 1 
       % If only one channel
       
       % Get list of frames after the skip
       frames_list = skip+1:nim; 
       
       % Get the number of frames after the skip
       frames = length(frames_list); 
       
        bData = cell2mat({im_list(frames_list).data});
        bData = reshape(bData, yDim, xDim, []);
        
        % Set aside images for spotcheck 
        spotcheck_data.initial.blue = bData(:,:, frames_for_spotchecking);
    end

    clear im_list;
    
    % ***3. Register within-stack/across stacks within a day.***
    disp('Registering within days'); 

    % Run the within-day registration function; overwrite bData
    % so you don't take up as much memory. 
    [tforms_forblueandviolet] = RegisterStackWithDFT(parameters.bRep, bData, usfac);

    % Apply the calculated tforms to the blue stack. Overwrite bData
    % so you don't take up as much memory.  
    [bData] = RegisterStack_WithPreviousDFTShifts(tforms_forblueandviolet, bData, usfac, yDim, xDim); 

    % Set aside images for spotcheck 
    spotcheck_data.withindayregistered.blue = bData(:,:, frames_for_spotchecking);

    % If more than one channel
    if channelNumber==2
        % Apply the calculated tforms to the violet stack. Overwrite vData
        % so you don't take up as much memory.  
        [vData] = RegisterStack_WithPreviousDFTShifts(tforms_forblueandviolet, vData, usfac, yDim, xDim); 

        % Also set aside image for spotcheck
        spotcheck_data.withindayregistered.violet = vData(:,:, frames_for_spotchecking);
    end 
  

    % *** 4. Apply registration across days ***
    
    % If the tform's empty, then you don't need to register
    if isempty(tform)==1 
        % Do nothing
    else
        % Else (the tform isn't empty) perform the registration/warp. 
        % Use imwarp to tranform the current image to align with the 
        % reference image using the tranform stored in the tform variable. 
        % Should be able to apply to all images in the 3rd dimension at the same time 
        disp('Applying registration across days');  
        bData = imwarp(bData,tform,'nearest', 'OutputView',imref2d([yDim xDim]));
        
        % Set aside images for spotcheck 
        spotcheck_data.registrationacrossdays.blue = bData(:,:, frames_for_spotchecking);
        
        % If more than 1 channel, do for violet channel as well
        if channelNumber == 2
            vData = imwarp(vData,tform,'nearest', 'OutputView',imref2d([yDim xDim]));
            spotcheck_data.registrationacrossdays.violet = vData(:,:, frames_for_spotchecking);
        end 
    end
    
    % ***Check if image of stack are the right size. Use just
    % the first frame. ***
    bData = FixImageSize(bData, parameters.pixels); 

    % If two channels, repeat with violetdata
    if parameters.channelNumber ==2
       vData = FixImageSize(vData, parameters.pixels); 
    end 
    
    % Make sure dimensions are updated to desired dimensions.
    yDim = parameters.pixels(1);
    xDim = parameters.pixels(2);
    
    % *** Reshape data into a 2D matrix (total pixels x frames) for
    % applying the mask, regressions, and the lowpass filter. Overwrite the variable
    % so you don't take up excess memory. ***
    bData = reshape(bData, yDim*xDim, frames);
    
    % If more than 1 channel, do for violet channel as well
    if channelNumber==2
        vData = reshape(vData, yDim*xDim, frames);
    end 
     
    % *** 5. Apply mask *** 
    % Keep only the indices that belong to the mask; Don't rename
    % the variable, because that will take up extra memory/time.
    
    % Only if user said to use mask (if mask_flag =  true)
    if mask_flag
        
        % Tell user what's happening.
        disp('Applying mask')
        
        % Apply mask (keep only pixels included in the mask).
        bData = bData(parameters.indices_of_mask,:);  

        % Set aside images for spotcheck 
        spotcheck_data.masked.blue = bData(:, frames_for_spotchecking);

        % If more than 1 channel, do for violet channel as well
        if channelNumber == 2
            vData = vData(parameters.indices_of_mask,:);
            spotcheck_data.masked.violet = vData(:, frames_for_spotchecking);
        end 
    end
    
     % ** *6. Filter***
    % Filter data.
    
    % Only if the user said they wanted to (if
    % filter_flag = true).
    if filter_flag
        disp('Filtering');

        % filtfilt treats each column as its own channel. Flip 
        % data as you put it into the filter so it's filtered
        % in temporal dimension. (frames x pixesl). 
        bData = filtfilt(b,a, bData'); 

        bData = bData'; 
        
        % Set aside images for spotcheck 
        spotcheck_data.filtered.blue = bData(:, frames_for_spotchecking);
        
        % If 2 channels, repeat for violet channel.
        if parameters.channelNumber == 2
            
            vData = filtfilt(b,a, vData'); 

            vData = vData'; 
        
            % Set aside images for spotcheck 
            spotcheck_data.filtered.violet = vData(:, frames_for_spotchecking);
            
        end 
    end 
    
    % *** 7. Correct hemodynamics. ***
    % Run HemoRegression function; 
    disp('Correcting hemodynamics');
    
    % Depending on the method desired by user
    switch correction_method
    
        case 'regression' 
            % Run regressions. 
            [data, intercepts_data] = HemoRegression(bData, vData);
            
        case 'scaling'
            % Run detrend-rescale version of hemo correction
            % (Laurentiu's version)
            data = HemoCorrection(bData, vData);
            
        case 'vessel regression'
            
            % Convert vessel masks into 2D matrix
            vessel_masks = reshape(vessel_masks, yDim*xDim, size(vessel_masks, 3));
            
            % If user said to use a brain mask,
            if mask_flag
                % Mask each blood vessel mask with the brain mask.
                vessel_masks = vessel_masks(indices_of_mask, :);
            end 
            
            % Run regression against extractions from blood
            % vessel masks. 
            data = VesselRegression(bData, vessel_masks); 
    end

    % Calculate blue channel stack mean. Convert to single
    % precision, put in parameters output structure.
    if isfield(parameters, 'save_stack_mean') && parameters.save_stack_mean

        if channelNumber == 1
            
            parameters.data_mean = single(mean(bData, 2, 'omitnan'));

        elseif channelNumber == 2
            
            parameters.data_mean_blue = single(mean(bData, 2, 'omitnan'));...
            parameters.data_mean_violet = single(mean(vData, 2, 'omitnan'));
            parameters.hemo_corrected_mean = single(mean(data + intercepts_data, 2, 'omitnan'));
        
        end 

        % Also calculate std
        parameters.data_std = single(std(data,[], 2, 'omitnan'));
    end

    % Set aside images for spotcheck 
    spotcheck_data.hemodynamicscorrected = data(:, frames_for_spotchecking);
    
    % *** 8. Save preprocessed stacks***
    disp('Saving');
    
    % Convert data to single precision to take up less space';
    % put in output structure
    parameters.data = single(data);

    % Also save blue channel from before hemo correction
    parameters.data_blue = single(bData);

    % Put spotcheck data into parameters structure
    parameters.spotcheck_data = spotcheck_data;

    % Bad trials info
    parameters.bad_trials = bad_trials;

end
