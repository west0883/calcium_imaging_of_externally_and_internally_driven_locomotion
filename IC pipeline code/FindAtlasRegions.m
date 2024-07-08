% FindAtlasRegions.m
% Sarah West
% 4/1/22

% Uses center-of-mass and percent area to find the best atlas region for
% each spatial source.

function [parameters] = FindAtlasRegions(parameters)
    
    % If there's a mask for this mouse & the atlas hasn't been masked yet,
    % do that now. 
    if isfield(parameters, 'indices_of_mask') && ~isfield(parameters, 'atlas_masked')
        
        % Set up atlas_masked as a matrix of zeroes the size of atlas.
        parameters.atlas_masked = zeros(size(parameters.atlas));

        % Set values inside mask to corresponding values of the atlas.
        parameters.atlas_masked(parameters.indices_of_mask) = parameters.atlas(parameters.indices_of_mask);
    
    % If no masking needed, make the atlas_masked the atlas 
    elseif ~isfield(parameters, 'indices_of_mask') && ~isfield(parameters, 'atlas_masked')
        parameters.atlas_masked = parameters.atlas;

    end
    
    % If atlas metrics haven't been calculated yet, do that now.
    if ~isfield(parameters, 'atlas_metrics')
           
        % Calculate metrics for each region name. Get list of all the
        % regions
        if isfield(parameters, 'region_order')
           all_regions = parameters.region_order(:,2);
        else
            all_regions = fieldnames(parameters.region_names);
        end

        % Get number of atlas region
        number_of_regions = numel(all_regions);

        % Make a holder of calulations. Make it a normal numeric array, 
        % because trying to do calculations on cells or structures is a nightmare.
        % This will just be center of mass x and y.
        parameters.atlas_metrics = NaN(number_of_regions, 2);

        % For each region, 
        for regioni = 1:number_of_regions
            
            % Get the value at of region so you can isolate the atlas
            % region
            if isfield(parameters, 'region_order')
                region_value = parameters.region_order{regioni,1};

            else
                region_value = getfield(parameters.region_names, all_regions{regioni}); 
            end 
           
            % Grab the atlas region. 
            region = parameters.atlas_masked == region_value;
    
            % Calculated center of mass (will be empty if mouse doesn't
            % have that region)>
            COM = centerOfMass(double(region));

            % Put it into atlas metrics 
            if ~isempty(COM)
                parameters.atlas_metrics(regioni, :) = COM; 
            end 
        end 

    end 
    
    % Begin comparisons

    % Get number of sources.
    number_of_sources = size(parameters.sources_artifacts_removed.sources, parameters.sourcesDim); 

    % Make a holder for each source's center of mass. 
    sources_center_of_mass = NaN(number_of_sources, 2); 

    % Make a holder for the comparison matrix (n regions x n sources x 2).
    % Third dimension is distance between center-of-mass, weighted area
    % overlap. 
    parameters.metrics.comparison_matrix = NaN(number_of_regions, number_of_sources, 2);

    % Get list of what the sources would've been called in regularized ICs
    %IC_old_list = 1:size(parameters.sources_artifacts_removed.sources, 3);
   % IC_old_list(parameters.sources_artifacts_removed.sources_removed) = [];

    % For each source, 
    for sourcei = 1:number_of_sources

        % Convert to nubmer that it would've been in the regularized ICs
        %source_index = IC_old_list(sourcei);
        
   
        % Pull out source
            % Set up abstractable dimensions
            S = repmat({':'},1, ndims(parameters.sources_artifacts_removed.sources));
            S{parameters.sourcesDim} = sourcei; 

            % Get out source
            source = parameters.sources_artifacts_removed.sources(S{:});

            % Make any NaNs in source (usually outside the mask) into 0 to
            % avoid match issues.
            source(isnan(source)) = 0; 

        % Calculate center of mass. 
        region_COM = centerOfMass(source);

        % Put into storage matrix
        sources_center_of_mass(sourcei, :) = region_COM;

        % Compare to each atlas region. 
        for regioni = 1:number_of_regions

            % *** Calculate distances between centers of mass***
            % Compare distance from center of mass to each atlas region 
            atlas_COM = parameters.atlas_metrics(regioni, :);

            % Calculate Euclidean distance
            parameters.metrics.comparison_matrix(regioni, sourcei, 1) = pdist([region_COM; atlas_COM], 'euclidean');
           
            % *** Calculate weighted overlap***
            % Get the value at of region so you can isolate the atlas
            % region
            region_value = getfield(parameters.region_names, all_regions{regioni}); 

            % Grab the atlas region. 
            region = parameters.atlas_masked == region_value;
           
            % Compare area overlapping with each atlas region (weighted);
            overlap = sum(sum(region .* source ))./(sum(sum(source))); 

            % Put into storage matrix.
            parameters.metrics.comparison_matrix(regioni, sourcei, 2) = overlap; 

        end
    end 

    % Find best region fit.

    % For each source,
    for sourcei = 1:number_of_sources
        
        % Find best COM distance (smallest), note the index 
        holder_COM = parameters.metrics.comparison_matrix(:, sourcei,1);
        index_COM = find(holder_COM == min(holder_COM, [], 'all', 'omitnan'));
        parameters.metrics.best_fit.indices(sourcei).COM = index_COM; 
        parameters.metrics.best_fit.names(sourcei).COM = all_regions{index_COM}; 

        % Find best overlap (largest), note the index & name 
        holder_overlap = parameters.metrics.comparison_matrix(:, sourcei,2);
        index_overlap = find(holder_overlap == max(holder_overlap, [], 'all', 'omitnan'));
        parameters.metrics.best_fit.indices(sourcei).overlap = index_overlap; 
        parameters.metrics.best_fit.names(sourcei).overlap = all_regions{index_overlap} ; 
        
       % If those indices match, assign to a best match field
       if index_COM == index_overlap 
          parameters.metrics.best_fit.indices(sourcei).best = index_COM; 
          parameters.metrics.best_fit.names(sourcei).best = all_regions{index_COM};

       else 
          parameters.metrics.best_fit.indices(sourcei).best = NaN; 
          parameters.metrics.best_fit.names(sourcei).best = [];
       end 
    end

    % Make a figure of 2 subplots: the COM and volume overlap metrics.
    figure_metrics = figure; 
    
    % Take center of mass metric matrix, make NANs -1 for plotting
    holder = parameters.metrics.comparison_matrix(:,:,1);
    holder(isnan(holder)) = -1;

    % Plot center of mass metric 
    subplot(1,2,1); 
    imagesc(holder);
    colormap([1 1 1; flipud(parula(256))]); colorbar;
    caxis([0 50]); % Limit the upper limit of the coloraxis, we don't care about specifics when it's very far away
    xlabel('source number'); ylabel('atlas region, left right');
    
    % Reformat region names for better tick labels. 
    tick_names = all_regions(1:2:end);
    for i = 1:numel(tick_names)
        tick_names{i} = replace(tick_names{i}, {'_L', '_', '1'}, {'', ' ', ''});
    end
    yticks([1:2:size(all_regions,1)]); yticklabels(tick_names);
    title('COM distance');

    % Plot overlap metric
    subplot(1, 2, 2); 
    imagesc(parameters.metrics.comparison_matrix(:,:,2));
    colormap(gca, [1 1 1; parula(256)]); colorbar; 
    xlabel('source number'); ylabel('atlas region, left right');
    yticks([1:2:size(all_regions,1)]); yticklabels(tick_names);
    title('overlap');
    sgtitle(['mouse ' parameters.values{1}]);

    % Make a figure with 3 subplots: the color-coded atlas, the ICs that
    % have a best fit with their color codes, the ICs that didn't fit
    % with any. Use the "new_order" list of regions 
    
    % Make color-coded atlas. 

    % Make a new holding image. 
    atlas_color_coded = zeros(size(parameters.atlas_masked)); 

    % Change value to appropriate color-value.
 
    for regioni = 1:number_of_regions
            
        % Get the value at of region so you can isolate the atlas
        % region
        if isfield(parameters, 'region_order')
            region_value = parameters.region_order{regioni,1};
        else
            region_value = getfield(parameters.region_names, all_regions{regioni}); 
        end 

        % Convert
        atlas_color_coded(parameters.atlas_masked == region_value) = regioni;
    end

    % Put color coded atlas into parameters.
    parameters.atlas_color_coded = atlas_color_coded;

    figure_best_fit = figure;
     
    % Get the colormap.
    mymap = repelem(jet(number_of_regions/2), 2, 1);
 
    % Draw atlas in first subplot; 
    subplot(1,2,1); 
    imagesc(atlas_color_coded); colormap([1 1 1; mymap]); caxis([0 number_of_regions]);
    axis square;

    % Draw agreeing & disagreeing sources with color-coding.
    best_regions_color_coded = zeros(size(parameters.atlas));
    not_fit_regions = atlas_color_coded;

    for sourcei = 1:number_of_sources
        
        % Set up abstractable dimensions
        S = repmat({':'},1, ndims(parameters.sources_artifacts_removed.sources));
        S{parameters.sourcesDim} = sourcei; 

        % Get out source
        source = parameters.sources_artifacts_removed.sources(S{:});

        % If there's a best fit region (not NaN)
        if ~isnan(parameters.metrics.best_fit.indices(sourcei).best)
   
             % Add to color coding, from best maching atlas
             % Get out best matching index
             best_regions_color_coded(source > 0) = parameters.metrics.best_fit.indices(sourcei).best; 

        else
            % If no best fit, plot over the atlas as gray (store as -1 for
            % now); 
           % best_regions_color_coded(source > 0) = -1; 

           % Plot with the overlap metric (so right now is just basically
           % using the overlap metric)
           best_regions_color_coded(source > 0) = parameters.metrics.best_fit.indices(sourcei).overlap; 


        end 

    end 

    % Plot agreeing regions.
    subplot(1, 2, 2); 
    imagesc(best_regions_color_coded); caxis([-1 number_of_regions]);
    colormap([0.5 0.5 0.5; 1 1 1; jet(number_of_regions)]); 
    axis square;
    sgtitle(['mouse ' parameters.values{1}]);
% 
%     % Plot not- agreeing regions.
%     subplot(1, 3, 3); 
%     imagesc(not_fit_regions); caxis([0 number_of_regions]);
%     colormap([0.5 0.5 0.5; 1 1 1; jet(number_of_regions)]);
%     axis square;

    % Put figure handles into parameters structure.
    parameters.figure_metrics = figure_metrics;
    parameters.figure_best_fit = figure_best_fit;

end 

