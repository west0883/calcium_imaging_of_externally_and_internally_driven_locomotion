% registration_SaveRepresentaiveImages.m
% Sarah West
% 8/25/21

% Loads individual images for differentiating between light channels, then
% saves a representaive image for the day. 
           
function [] = registration_SaveRepresentativeImages(parameters)
     
    % Assign parameters their original names
    dir_dataset_name = parameters.dir_dataset_name; 
    input_data_name = parameters.input_data_name;
    dir_exper = parameters.dir_exper;  
    skip = parameters.skip;
 
    rep_stacki = parameters.rep_stacki; 
    rep_framei = parameters.rep_framei; 
    digitNumber = parameters.digitNumber; 
    channelNumber = parameters.channelNumber;
    
    dir_out_base=[dir_exper 'representative images\'];
    disp(['Data saved in ' dir_out_base]); 
    
    % for each mouse
    for mousei=1:size(parameters.mice_all,2)
        % Get the mouse name
        mouse=parameters.mice_all(mousei).name;
        
        % Get the mouse's dataset days
        days_list=vertcat(parameters.mice_all(mousei).days(:).name); 
        
        % For each day of the mouse 
        for dayi=1:size(days_list,1)
            day=days_list(dayi,:);
            
            % get the day name
            day=days_list(dayi,:); 
            
            % make the output directory by day
            dir_out=[dir_out_base '\' mouse '\' day '\']; 
            mkdir(dir_out);
            
            parameters.dir_in = dir_dataset_name;
           
            % Get the stack list
            [stackList]=GetStackList(mousei, dayi, parameters);
            
            % Get the stack number to use for the rep image based on the
            % given stack index. 
            stack_number = stackList.numberList(rep_stacki, :);
            
            %Create the input directory
            dir_in=CreateFileStrings(dir_dataset_name, mouse, day, stack_number, [], false);
            
            % Use the stack number to make an input filename
            stackname=stackList.filenames(rep_stacki, :);
            input_filename=[dir_in stackname];
                
            % If there are 2 channels, 
            if channelNumber == 2 
                
                % Select 2 sequential images to read; Determined by the 
                % skip and the ref_framei given by the user. Need 2 to 
                % confirm which channel is which.
                image_indices=[skip + rep_framei, skip + rep_framei + 1]; 

                % Read those two images
                im_list=tiffreadAltered_SCA(input_filename, image_indices, 'ReadUnknownTags',1);

                im1=im_list(1).data;
                im2=im_list(2).data; 

                % Compare the brightness of the two images; 
                first_image_channel = DetermineChannel(parameters.blue_brighter, im1, im2, parameters.pixel_rows, parameters.pixel_cols);

                if first_image_channel=='b'
                   bRep=im1; 
                elseif first_image_channel=='v'
                   bRep=im2;
                end
            
            else
                % Else, there's only 1 channel.
                % Read in only 1 image
                image_indices=skip + rep_framei; 
                im_list=tiffreadAltered_SCA(input_filename, image_indices, 'ReadUnknownTags',1);
                
                % Assign bRep to be the image you read in. 
                bRep=im_list(1).data;
            end
            
            % Convert bRep to single precision.
            bRep=single(bRep);
            
            % Save the representative image
            save([dir_out 'bRep.mat'], 'bRep');
        
        end 
    end
end 
