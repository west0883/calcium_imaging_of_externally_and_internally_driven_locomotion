% HistogramOfRatios.m
% Sarah West
% 6/23/22

function [parameters] = HistogramOfRatios(parameters)

 % Turn into a numeric array instead of cell array.
 holder = [parameters.data{:}];

 % Plot histogram.
 parameters.fig = figure; 
 histogram(holder, [0:0.05:1]);

 % Get this mouse number for title.
 index = strcmp(parameters.keywords, 'mouse');
 mouse = parameters.values{index};

 title(['NaN ratios for pupil diamters in m' mouse]);
 
end