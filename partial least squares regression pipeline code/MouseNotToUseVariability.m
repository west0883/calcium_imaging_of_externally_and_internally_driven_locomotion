
function [parameters] = MouseNotToUseVariability(parameters)

    MessageToUser('Checking mouse for ', parameters)
    
    % Get current mice not to use 
    dir_stem_iterator = parameters.values{strcmp(parameters.keywords, 'dir_stem_iterator')};
    comparison_iterator = parameters.values{strcmp(parameters.keywords, 'comparison_iterator')};
    mice_not_to_use = parameters.loop_variables.dir_stems(dir_stem_iterator).comparisons(comparison_iterator).mice_not_to_use;
    
    % If it's not empty,
    if ~isempty(mice_not_to_use)

        % get the current mouse
        mouse = parameters.values{strcmp(parameters.keywords, 'mouse')};

        % If current mouse is in mice not to use, 
        if strcmp(mice_not_to_use, mouse)

            % make data out NaN
            data_out = NaN(size(parameters.data));
        
        % If current mouse isn't in mice not to use
        else
            % make data_out equal to inputted data
            data_out = parameters.data;
         
        end 

    % If empty, make data_out equal to inputted data
    else 
        data_out = parameters.data;
    end 

    parameters.data_out = data_out; 
end 