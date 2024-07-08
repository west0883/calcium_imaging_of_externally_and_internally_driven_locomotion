% UseThisMouse.m
% Sarah West
% 11/14/23

function [parameters] = UseThisMouse(parameters)
   
    data = parameters.data; 

    % Find the location & value of the comparison iterator in
    % parameters.values
    iterator_location = strcmp(parameters.keywords, 'comparison_iterator'); 
    comparison_iterator = parameters.values{iterator_location};

    % *** Deal with mice that should be skipped for this comparison ***
    
    % If comparison_type is a field 
    if any(strcmp(parameters.keywords, 'comparison_type'))

        comparison_type = parameters.values{strcmp(parameters.keywords, 'comparison_type')};
        this_comparison_set = parameters.loop_variables.(['comparisons_' comparison_type]); 
        
    elseif isfield(parameters, 'comparison_type')
        comparison_type = parameters.comparison_type; 

        this_comparison_set = parameters.loop_variables.(['comparisons_' comparison_type]); 
        
    % If "mice_not_to_use" is a field in this comparison set,
    else 
        
        this_comparison_set = parameters.this_comparison_set;

    end 

    if isfield(this_comparison_set(comparison_iterator), 'mice_not_to_use')
       
        % Get the list of mice not to use in this comparison (if any)
        mice_not_to_use = this_comparison_set(comparison_iterator).mice_not_to_use;

        % Get the mouse of this iteration
        mouse_location = strcmp(parameters.keywords, 'mouse'); 
        mouse = parameters.values{mouse_location};
    
        % If mice_not_to_use is not empty
        if ~isempty(mice_not_to_use)
    
            % Check if current mouse is included in the list of mice not to
            % use.
            if any(cellfun(@strcmp, mice_not_to_use, repmat({mouse}, size(mice_not_to_use)))) 
                
                % If you don't use this mouse, make data empty
                data = [];
               
            end 
        end
    end

    % Put into output structure
    parameters.data = data; 
end 