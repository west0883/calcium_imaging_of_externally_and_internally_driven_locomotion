% PlotDurationHistogram.m
% Sarah West
% 11/4/23

function [parameters] = PlotDurationHistogram(parameters)

    MessageToUser('Plotting ', parameters)

    rest = parameters.rest;
    walk = parameters.walk;
    edges = parameters.edges;

    % get mouse for title 
    mouse = parameters.values{strcmp(parameters.keywords, 'mouse')};
    
    fig = figure; 

    histogram(rest, edges);

    hold on; 
    histogram(walk, edges);

    legend({'rest', 'walk'});
    title(mouse)

    parameters.fig = fig;
end 


