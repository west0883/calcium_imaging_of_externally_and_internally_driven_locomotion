% PlotMiceStrideAveragess.m
% Sarah West 
% 7/18/23

% For a given behavior period plots the average, resampled stride
% (normalized to be 10 data points long) across mice. Uses the standard 
% error of the mean as error bars.
% Runs with RunAnalysis

function [parameters] = PlotMiceStrideOneMotorAndSpon(parameters)
    
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
    average_motor = parameters.average_motor;
    std_dev_motor = parameters.std_dev_motor;
    concatenated_data_motor = parameters.concatenated_data_motor;
    average_spon = parameters.average_spon;
    std_dev_spon = parameters.std_dev_spon;
    concatenated_data_spon = parameters.concatenated_data_spon;
    instancesDim = parameters.instancesDim;
    ylimits = parameters.ylimits;
    errorType = parameters.errorType;
    colors = parameters.colors;
    lineTypes = parameters.lineTypes;
        
    fig = figure; 

    % Calculate standard error of the mean (SEM)
    SEM_motor = std_dev_motor./sqrt(size(concatenated_data_motor, instancesDim));
    SEM_spon = std_dev_spon./sqrt(size(concatenated_data_spon, instancesDim));

    % Calculate confidence interval (CI)
    ts_motor = tinv([0.025  0.975], size(concatenated_data_motor, instancesDim) - 1);      % T-Score
    CI_motor = ts_motor(2).*SEM_motor;  
    ts_spon = tinv([0.025  0.975], size(concatenated_data_spon, instancesDim) - 1);      % T-Score
    CI_spon = ts_spon(2).*SEM_spon;

    if strcmp(errorType, 'SEM')
        errorRange_motor = SEM_motor;
        errorRange_spon = SEM_spon;
    elseif strcmp(errorType, 'std_dev')
        errorRange_motor = std_dev_motor;
        errorRange_spon = std_dev_spon;
    else 
        errorRange_motor = CI_motor;
        errorRange_spon = CI_spon;
    end 

    % Plot the average and shaded error bars
    % spontaneous
    this_color = colors(1, :);
    s = shadedErrorBar(1:10, average_spon, errorRange_spon, 'patchSaturation', 0.2, 'lineProps',{'Color', this_color});
    delete(s.edge); % delete shaded area edges
    set(s.mainLine, 'LineWidth', 2);
   
    % motor
    hold on; 
    this_color = colors(2, :);
    s = shadedErrorBar(1:10, average_motor, errorRange_motor, 'patchSaturation', 0.2, 'lineProps',{'Color', this_color});
    delete(s.edge); % delete shaded area edges
    set(s.mainLine, 'LineWidth', 2);
   
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
    title(['mean with ' errorType  ', ' paw_section velocity_direction_section], 'Interpreter', 'none');

    % Put the output figure into outputs
    parameters.fig = fig;

end 
