% OptimizeSubplotNumbers.m
% Sarah West
% 3/18/22

% A small function to determined the optimized subplot dimensions for a
% given number of plots to fit in one figure. Ratio is row:column ratio,
% depending on user's input (optional, default = 2/3). 

function [subplot_rows, subplot_columns] = OptimizeSubplotNumbers(number_of_plots, ratio)

    % If no ratio input, default to 2/3
    if nargin < 2
        ratio = 2/3; 
    end 
    
    % If there is only one plot, skip all this and make just 1 subplot.
    if number_of_plots == 1
        subplot_rows = 1;
        subplot_columns = 1;

    else

        % Do the math for the columns (unrounded)
        columns = (1/ratio * number_of_plots) ^.5;
        
        % Find number of columns by rounding up to nearest integer (favors more
        % columns than rows). 
        subplot_columns = ceil(columns);
    
        % Find number of rows from unrounded columns number, then round down to
        % nearest integer (favors more columns than rows).
        subplot_rows = floor(ratio * columns);
    
        % Make sure the optimization worked. If not, increase the number of
        % rows (usually looks weird if you try increasing columns here), or if that isn't enough, increase the number of columns
        while subplot_columns * subplot_rows < number_of_plots
            
            % Try increasing number of rows
            if subplot_columns * (subplot_rows + 1) >= number_of_plots
                subplot_rows = subplot_rows +1; 
    
            else
                % If that doesn't work, try also increasing number of columns & try
                % again.
                subplot_rows = subplot_rows + 1; 
                subplot_columns = subplot_columns + 1; 
            end
        end
    end
end 