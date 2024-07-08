% PickReferenceDay.m
% Sarah West
% 8/23/21

% Loads the "bRep.mat" file for each day of each mouse (which is the
% reference image for the within- and across-stack registrations) and plots them
% so you can pick what day to use. Called by "registration_pick_reference_day_permouse.m"

function [dayi_output]=PickReferenceDay(mouse, file_paths, plot_sizes)
    
    % Set the dimensions/layout of the subplots you want to use
    plot_rows=plot_sizes(1); 
    plot_columns=plot_sizes(2);
    
    % Initialize the figure to plot in
    figure; sgtitle(['mouse ' mouse]);
    
    % Enlarge figure to full screen.
    set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
    
    % for each day
    for dayi=1:size(file_paths,2)
        
        % load the background image
        load(file_paths{dayi}); 
        
        % plot the background image in the proper subplot
        subplot(plot_rows, plot_columns, dayi); imagesc(bRep);
        
        % give the image a title corresponding to the dayi
        title(num2str(dayi));
    end 
    
    % As the user which image/day they want to use as the reference
     
    dayi_output = inputdlg('Which image would you like to use as the reference image? Type the number below.');
    
    % make the output not a cell array anymore and from str to number
    dayi_output=str2num(dayi_output{1});
end 