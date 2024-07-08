% PlotMiceStrideAveragess.m
% Sarah West 
% 7/18/23

% For a given behavior period plots the average, resampled stride
% (normalized to be 10 data points long) across mice. Uses the standard 
% error of the mean as error bars.
% Runs with RunAnalysis

function [parameters] = PlotMiceStrideAverages(parameters)
    
    % Inputs
    % parameters.average -- a vector, one per mouse; the mean stride
    % parameters.std_dev -- a vector, one per mouse; the standard deviation
    % of that mouse's stride
    % parameters.concatenated_data -- a matrix, one per mouse; all the strides for
    % that period in that mouse; is for getting the number of samples to
    % plot the standard error of the mean
    % parameters.instancesDim -- integer; the dimension of
    % 'concatenated_data' that has different instances 
    % parameters.ylimits -- pair of scalars; is the yaxis limit you want to
    % plot


    % Outputs
    % parameters.fig -- the figure of overlaid strides

    % Tell user what iteration you're on
    MessageToUser('Plotting ', parameters);

    % Pull out inputs 
    average = parameters.average;
    std_dev = parameters.std_dev;
    concatenated_data = parameters.concatenated_data;
    instancesDim = parameters.instancesDim;
    ylimits = parameters.ylimits;
    errorType = parameters.errorType;
    colors = parameters.colors;
    lineTypes = parameters.lineTypes;
        
    fig = figure; 
 
    % Get the current period & iterator
    if any(strcmp(parameters.keywords, 'period'))
    period = parameters.values{strcmp(parameters.keywords, 'period')};
    period_iterator = parameters.values{strcmp(parameters.keywords, 'period_iterator')};
    else 
        period = [];
        period_iterator = 1;
    end
    
    % if no data in average, skip to next 
    if isempty(average)
        parameters.fig = fig;
        parameters.SEM = [];
        return
    end

    % Calculate standard error of the mean (SEM)
    SEM = std_dev./sqrt(size(concatenated_data, instancesDim));

    % Calculate confidence interval (CI)
    ts = tinv([0.025  0.975], size(concatenated_data, instancesDim) - 1);      % T-Score
    CI = ts(2).*SEM;  

    if strcmp(errorType, 'SEM')

        errorRange = SEM;
    elseif strcmp(errorType, 'std_dev')
        errorRange = std_dev;
    else 
        errorRange = CI;
    end 
    % Plot the average and shaded error bars
    s = shadedErrorBar(1:10, average, errorRange, 'patchSaturation', 0.2);
    delete(s.edge); % delete shaded area edges
    set(s.mainLine, 'LineWidth', 1);
   
    % Set the y axis limits
    ylim(ylimits);

    % Make figure title 
    if any(strcmp(parameters.keywords, 'paw'))
        paw = parameters.values{strcmp(parameters.keywords, 'paw')};
        paw_section = [ paw ', '];
    else
        paw_section = [];
    end 
    
    if any(strcmp(parameters.keywords, 'velocity_direction'))
        velocity_direction = parameters.values{strcmp(parameters.keywords, 'velocity_direction')};
        velocity_direction_section = [ velocity_direction ', '];
    else
        velocity_direction_section = [];
    end 
    title(['mean with ' errorType  ', ' paw_section velocity_direction_section period ' ' num2str(period_iterator)], 'Interpreter', 'none');


    % Plot all the averages per paw together 
    if parameters.do_allPeriods 
        if ~isfield(parameters, 'fig_allPeriods')
            fig_allPeriods = figure;
        else 
            fig_allPeriods = parameters.fig_allPeriods;
        end 
        set(0, 'currentfigure', fig_allPeriods);  % set to current figure
        hold on;
    
        % Get colors for this plot
        period_iterator = parameters.values{strcmp(parameters.keywords, 'period_iterator')};
        this_color = colors(period_iterator, :);
        this_lineType = lineTypes{period_iterator};
        s = shadedErrorBar(1:10, average, errorRange, 'patchSaturation', 0.2, 'lineProps',{'Color', this_color, 'LineStyle', this_lineType});
        delete(s.edge); % delete shaded area edges
        set(s.mainLine, 'LineWidth', 2);
    
        % Set title & whatnot 
        title(['mean with ' errorType  ', ' paw_section velocity_direction_section], 'Interpreter', 'none');
        ylim(ylimits);

        parameters.fig_allPeriods = fig_allPeriods;
    end 

    % Put the output figure into outputs
    parameters.fig = fig;
    parameters.SEM = SEM;
    

end 
