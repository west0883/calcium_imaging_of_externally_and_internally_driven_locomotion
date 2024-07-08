% plot_spotcheck.m
% Sarah West
% 12/1/21

% Plots a subset of the saved "spotcheck" frames from preprocessing and
% saves the plots. Makes a plot for each stack, also a plot for each day.

function [] = plot_spotcheck(parameters)
    
    dir_in_base = parameters.dir_in_base;
 
    % For each mouse
    for mousei = 1:size(parameters.mice_all,2) 
        
        % Get mouse name
        mouse = parameters.mice_all(mousei).name; 
        
        % If masked, load masks
        if parameters.mask_flag
            load([parameters.dir_exper 'masks\masks_m' mouse '.mat']);
        end 
        
        % For each day
        for dayi = 1:size(parameters.mice_all(mousei).days,2)
            
            % Get day name
            day = parameters.mice_all(mousei).days(dayi).name;
            
            % Get input directory for that day
            parameters.dir_in = [parameters.dir_in_base  mouse '\' day '\'];
            
            % Change what the input name you want is
            parameters.input_data_name = {'spotcheck_data', 'stack number', '.mat'};
            
            % Get the stack list
            [stackList]=GetStackList(mousei, dayi, parameters);
            
            % Make holder matrix for all final (hemodynamics corrected)
            % plots
            
            hemo_all = NaN(size(stackList.filenames,1), numel(parameters.frames_to_plot), 256,256);
            
            % For each stack
            for stacki=1:size(stackList.filenames,1)
                
                fig = figure;
                fig.WindowState = 'maximized';
                
                % Get the stack number and filename for the stack.
                stack_number = stackList.numberList(stacki, :);
                filename = stackList.filenames(stacki, :);
        
                % Load the spotcheck data.
                load([parameters.dir_in filename]); 
                
                % Get number of subplots, depending on if this stack needed
                % to be registered across days, and on how many channels
                % there are.
                a = (numel(fieldnames(spotcheck_data)));  
                b =  numel(parameters.frames_to_plot)*parameters.channelNumber;
                number_of_subplots = [a , b];
                    
                % Plot each step, in preprocessing order

                % For each channel (for some plots, need both)
                for channeli = 1:parameters.channelNumber
                    
                    % Keep track of what row you're plotting in
                    row = 1;
                    
                    % Assign the channel name, starting column index
                    if channeli == 1
                       channel = 'blue';
                       column_start = 1;
                    else
                       channel = 'violet'; 
                       column_start = 2; 
                    end
                    
                    % First the initial image
                    for i = 1:numel(parameters.frames_to_plot)
                         
                        subplot(a, b, (row -1) * b + column_start + (i-1) *2); 
                        eval(['imagesc(spotcheck_data.initial.' channel '(:, :, i));']);    
                       if column_start + (i-1) *2 ==1
                            ylabel('initial');
                       end
                       caxis(parameters.color_range);
                    end 
                  
                    % Within day registered
                    row =row + 1;
                    for i = 1:numel(parameters.frames_to_plot)
                         
                        subplot(a, b, (row -1) * b + column_start + (i-1) *2); 
                        eval(['imagesc(spotcheck_data.withindayregistered.' channel '(:, :, i));']);    
                       if column_start + (i-1) *2 ==1
                            ylabel('within-day');
                       end
                       caxis(parameters.color_range);
                    end 
                    
                    % Across day registered
                    if isfield(spotcheck_data,'registrationacrossdays')
                        row =row +1;
                        
                        for i = 1:numel(parameters.frames_to_plot)
                            
                            subplot(a, b, (row -1) * b + column_start + (i-1) *2); 
                            eval(['imagesc(spotcheck_data.registrationacrossdays.' channel '(:, :, i));']); 
                            if column_start + (i-1) *2 ==1
                                ylabel('across-day');
                            end
                            caxis(parameters.color_range);
                        end 
                      
                    end

                    % Masked
                    if parameters.mask_flag
                        row =row + 1;
                        
                        % Fill masks 
                        eval(['data_matrix = spotcheck_data.masked.' channel ';']);
                        [data_matrix_filled]=FillMasks(data_matrix, indices_of_mask, parameters.pixels(1), parameters.pixels(2));
                        
                        for i = 1:numel(parameters.frames_to_plot)
                            
                            subplot(a, b, (row -1) * b + column_start + (i-1) *2); 
                            imagesc(data_matrix_filled(:, :, i)); 
                           if column_start + (i-1) *2 ==1
                                ylabel('masked');
                           end
                           caxis(parameters.color_range);
                        end 
                        
                    end 
                    
                    % Any filtering
                    if parameters.filter_flag
                        row =row + 1;
                        
                        eval(['data = spotcheck_data.filtered.' channel ';']);    
                        if parameters.mask_flag   
                            [data_matrix_filled] = FillMasks(data, indices_of_mask, parameters.pixels(1), parameters.pixels(2));
                        else 
                            data_matrix_filled = data;
                        end        
                        
                        for i = 1:numel(parameters.frames_to_plot) 
                         
                            subplot(a, b, (row -1) * b + column_start + (i-1) *2); 
                            imagesc(data_matrix_filled(:,:, i));
                            if column_start + (i-1) *2 ==1
                                ylabel('filtered');
                            end
                            caxis(parameters.color_range);
                        end
                    end
                    
                end
                
                % Reset column start to 1 for the post-hemo correction
                % parts
                column_start = 1;
                
                % Hemodynamics corrected
                row = row + 1; 
                    
                % Fill masks for hemo 
                if parameters.mask_flag   
                   [data_matrix_filled]=FillMasks(spotcheck_data.hemodynamicscorrected, indices_of_mask, parameters.pixels(1), parameters.pixels(2));
                else 
                    data_matrix_filled = spotcheck_data.hemodynamicscorrected;
                end        
                
                for i = 1:numel(parameters.frames_to_plot)
                    subplot(a, b, (row -1) * b + column_start + (i-1) *2); 
                    imagesc(data_matrix_filled(:, :, i)); 
                    if column_start + (i-1) *2 ==1
                        ylabel('hemo corrected');
                    end
                    caxis(parameters.color_range_hemocorrected);
                    hemo_all(stacki, i, :, :) = data_matrix_filled(:,:,i); 
                end
                
                % Add global title.
                sgtitle(['mouse ' mouse ', day ' day ', stack' stack_number]);
            
                % Save figure
                savefig([parameters.dir_in 'check_bystack_' stack_number '.fig']);
                
            end
            
            % Plot final preprocessed images togther per day.
            
            % Make a new count of suplots
            a = size(stackList.filenames,1);
            b = numel(parameters.frames_to_plot); 
             
            fig = figure;
            fig.WindowState = 'maximized';
            
            % Plot each frame
            % For each stack
            for stacki=1:size(stackList.filenames,1)
                
                for i = 1:numel(parameters.frames_to_plot)
                    subplot( b, a, (i-1)*a + stacki); 
                    imagesc(squeeze(hemo_all(stacki,i,:,:)));
                    caxis(parameters.color_range_hemocorrected);
                    xticks([]);
                    xticklabels([]);
                    yticks([]);
                    yticklabels([]);
                    
                end
            end
            sgtitle(['mouse ' mouse ', day ' day]);
            savefig([parameters.dir_in 'all_hemo.fig']);
            % Close all figures at end of a day
            close all; 
        end 
    end
end