% PlotRegionAssignments.m
% Sarah West
% 4/4/22

% Plots a subplot of the atlas, subplot of the sources color-coded by their
% assigned regions.

function [parameters] = PlotRegionAssignments(parameters)
    
    % Make a blank image to add the sources to.
    region_assigned_sources = zeros(size(parameters.atlas));

    % For each source, 
    for sourcei = 1:size(parameters.manual_assignments,1)

        % Put the new color in.
        region_assigned_sources(parameters.sources_artifacts_removed.sources(:,:, sourcei) > 0) = parameters.manual_assignments(sourcei, 2);

    end

    % Get the colormap.
    number_of_regions = 72;
    mymap = repelem(jet(number_of_regions/2), 2, 1);

    % Plot atlas.
    parameters.figure_region_assignments = figure;
    subplot(1,2,1);
    imagesc(parameters.atlas); colormap([1 1 1; mymap]); caxis([0 number_of_regions]);
   
    % Plot color-coded sources.
    subplot(1,2,2);
    imagesc(region_assigned_sources); colormap([1 1 1; mymap]); caxis([0 number_of_regions]);

    sgtitle(['m' parameters.values{1}]);
end 