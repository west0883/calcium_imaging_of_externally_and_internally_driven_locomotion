% ReorderSources.m
% Sarah West
% 4/4/22

function [parameters] = ReorderSources(parameters)

    % Tell user what you're doing.
    MessageToUser('Reordering', parameters);

    % Make a holding matrix
    parameters.sources_reordered = NaN(parameters.yDim, parameters.xDim, parameters.num_new_sources);

    % For each (new) source/node in the manual assignments list
    for sourcei = 1:size(parameters.assigned_region_order,1) 
        
        % Get new place/source number
        new_place = parameters.assigned_region_order(sourcei, 2); 

        % Get original place/source number
        original_place = parameters.assigned_region_order(sourcei, 1);

        % Put source in new place. 
        
        % If user says so, do this by adding so multiple sources that
        % belong to the same (new) source/node so they're both represented.
        % Default is to over-write previous entries. 
        if isfield(parameters, 'add_sources') && parameters.add_sources 
           
            holder = cat(3, parameters.sources_reordered(:,:, new_place), parameters.sources(:,:, original_place));
            parameters.sources_reordered(:,:, new_place) =  sum(holder, 3, 'omitnan');
        
        % Otherwise, overwrite holder entry. 
        else
            parameters.sources_reordered(:,:, new_place) = parameters.sources(:,:, original_place); 
        end 
    end
end