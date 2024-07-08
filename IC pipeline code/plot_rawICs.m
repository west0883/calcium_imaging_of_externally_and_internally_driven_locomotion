% plot_rawICs.m
% Sarah West
% 9/16/21
% Plots raw ICs from calculate_ICs.m for looking at/inspection.

function []=plot_rawICs(parameters)
    
    % Return parameters to individual names.
    mice_all = parameters.mice_all;
    num_sources = parameters.num_sources;
    yDim = parameters.yDim;
    xDim = parameters.xDim;
    masked_flag = parameters.masked_flag; 
    
    % Tell user where data is being saved. 
    disp(['data saved in ' parameters.dir_output_base{1}]); 
    
    % For each mouse
    for mousei=1:size(mice_all,2)  
        
        % Get the mouse name and display to user.
        mouse=mice_all(mousei).name;

        % Establish input folder name you're working with. 
        filename = CreateFileStrings([parameters.dir_input_base parameters.input_filename], mouse, [], [], [], false);
        input_variable = CreateFileStrings(parameters.input_variable, mouse, [], [], [], false);
        
        % Load the raw sources, convert variable to something generic.
        load(filename, input_variable); 
        eval(['sources = ' input_variable ';']);
        
        % Reshape and permute for plotting 
        
        % If the data was masked, load mask and reshape specially
        if masked_flag==1
            
            % Find file name of masks
            filename_mask=CreateFileStrings([parameters.dir_input_mask parameters.mask_filename], mouse, [], [], false);
            mask_variable = CreateFileStrings(parameters.mask_variable, mouse, [], [], false);

            % Load mask indices, convert variable to something generic. 
            load(filename_mask, mask_variable); 
            eval(['indices_of_mask = ' mask_variable ';']);
            
            % Flip sources for inputting into FillMasks
            sources=sources';
            
            % Fill in masks with FillMasks.m function
            sources_permute=FillMasks(sources, indices_of_mask, yDim, xDim);
        
        % If no masks, just reshape and permute sources    
        else
            sources_permute=permute(reshape(sources, num_sources, yDim, xDim),[2 3 1]);
        end 
        
        % Get the number of subplots to use, if user put it in.
        if isfield(parameters, 'plot_sizes')
            subplot_rows=parameters.plot_sizes(1);
            subplot_columns=parameters.plot_sizes(2); 
        else
            % Otherwise, get the optimized version.
            [subplot_rows, subplot_columns] = OptimizeSubplotNumbers(num_sources);
        end

        % Plot sources
        figure; 
        for i=1:num_sources 
            subplot(subplot_rows, subplot_columns,i); 
            imagesc(abs(sources_permute(:,:,i))); 
            caxis([0 10]); 
            xticks([]);
            yticks([]);
        end
        sgtitle(['m' mouse ', ' num2str(num_sources) ' absolute value of sources']); 
         
        % Get output names.
        filename_output = CreateFileStrings([parameters.dir_output_base parameters.output_filename], mouse, [], [], [], false);

        % Save
        savefig(filename_output); 
    
    end 
end 