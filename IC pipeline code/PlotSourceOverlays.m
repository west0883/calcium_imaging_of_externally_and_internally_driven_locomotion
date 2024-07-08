% PlotCleanICs.m
% Sarah West
% 3/21/22

function [parameters] = PlotCleanICs(parameters)

    % If the figure handle doesn't exist yet, create empty figure handle.
    if ~isfield(parameters.sources_artifacts_removed_figure)
        
        fig_overlay = figure; 
        fig_colormasks = figure; 
       
    else 
        % If it does exist, iteratively add to it.
        fig_overlay= parameters.sources_artifacts_overlay; 
        fig_colormasks= parameters.sources_artifacts_colormasks;
    end 

   


    % At the end, rename figure handle for passing.
    parameters.sources_artifacts_removed_overlay = fig_overlay;
    parameters.sources_artifacts_removed_colormasks = fig_colormasks;
    
end 