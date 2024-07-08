

function [parameters] = ncomponents_std_correlate(parameters)

       MessageToUser('Running ', parameters);

       if ~isfield(parameters, 'ncomponents_concatenated')
          parameters.ncomponents_concatenated = [];
          parameters.sigmas_concatenated = [];

       end 

       % Get the list of mice not to use in this comparison (if any)
        mice_not_to_use = parameters.this_comparison_set(parameters.values{3}).mice_not_to_use;
    
        % Get the mouse of this iteration
        mouse_location = strcmp(parameters.keywords, 'mouse'); 
        mouse = parameters.values{mouse_location};
    
        % If mice_not_to_use is not empty
        if ~isempty(mice_not_to_use)
    
            % Check if current mouse is included in the list of mice not to
            % use. If it is, exit function.
            if any(cellfun(@strcmp, mice_not_to_use, repmat({mouse}, size(mice_not_to_use)))) 
                
                % Tell RunAnalysis not to save anything this iteration.
                parameters.dont_save = true;
    
                % Leave this function-- will move on to next comparison
                return
            end 
        end

       if parameters.values{end} == size(parameters.mice_all,2)
           data = parameters.data.responseVariables;
           sigma = std(data,[], 2);
           parameters.sigmas_concatenated = [parameters.sigmas_concatenated; sigma];
       end
     
       
       ncomponents = parameters.ncomponents_used; 
       parameters.ncomponents_concatenated = [parameters.ncomponents_concatenated; ncomponents]; 

       if parameters.values{3} == size(parameters.this_comparison_set, 2)
           if parameters.values{4} == size(parameters.mice_all,2)
           [R,P] = corrcoef(parameters.sigmas_concatenated, parameters.ncomponents_concatenated);
           parameters.correlation = [R(1,2), P(1,2)];

           parameters = rmfield(parameters, 'sigmas_concatenated');
           parameters = rmfield(parameters, 'ncomponents_concatenated');
           else
               parameters.correlation = [];
           end
       else
           parameters.correlation = [];
       end 
end