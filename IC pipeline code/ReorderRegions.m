% ReorderRegions.m
% Sarah West
% 1/31/22

% Starting with the LocaNMF map you want, reeorders regions into a logical 
% order based on ant-post location of subdivided atlas. 

function [] = ReorderRegions(parameters)

    % Load subdivided basline atlas. 
    load(parameters.dir_subdivided_atlas);

    % Grab list of area names;
    all_regions = fieldnames(areanames);

    % Get only the left regions. 
    left_regions = all_regions([1:31 63:67]);

    % Create holding matrix of centroids.
    centroids = NaN(numel(left_regions), 2);

    % For each region 
    for regioni = 1:numel(left_regions)
        
        % Get the name of the region and corresponding atlas value.
        region_name = left_regions{regioni};
        eval(['value = areanames.' region_name ';']);

        % Get indices of atlas region.
        indices = find(atlas == value);
        
        % Put into new atlas.
        holder = zeros(size(atlas));
        holder(indices) = 1;
        
        %new_atlas(:, :, regioni) = holder;

        % Make new region names variable
        regions(regioni).name = region_name;
        regions(regioni).value = value;

        % Find centroids/centers of mass
        stats = regionprops(holder,'Centroid');
        centroids(regioni, :) = round(stats.Centroid);

    end
    
    % Sort left centroids by ant-pos.
    [sorted_centroids, sorted_indices] = sort(centroids(:,2));

    % Use sort positions to put in new sort names and values
    sorted_areanames = regions(sorted_indices);

end 