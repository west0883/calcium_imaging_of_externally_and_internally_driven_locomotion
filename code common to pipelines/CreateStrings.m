% CreateStrings.m
% Sarah West
% 2/25/22

% Takes an input cell array of character strings, and inputs any variable
% to any keyword. Outputs a string.

% Input: 
% string_format_cell--cell array. Establish the format of the file names of compressed data. Each piece
% needs to be a separate entry in a cell array. 
% keywords - a cell array of strings that you look for inside the
% string_format_cell. 
% variables - a cell array of variables (with string values) that you put in to replace the
% corresponding keywords.

% Output: 
% new_string -- string (character array); the new string with keywords
% inserted.

function [new_string]=CreateStrings(string_format_cell, keywords, variables)
    
    % Make a new cell array to manipulate. 
    new_string= string_format_cell;
    
    % For each keyword,
    for keywordi = 1:numel(keywords)
        
        keyword = keywords{keywordi};
    
        % See if there is an entry for mouse number, find where it's located
        keyword_index=find(strcmp(string_format_cell, keyword));
        
        % If there is a position, but there's no corresponding keyword
        if ~isempty(keyword_index) && isempty(variables{keywordi})

            error(['No variable given for ''' keyword '''.']);

        elseif ~isempty(keyword_index)
            % Put the variable in place of the keyword
            for subi = 1:numel(keyword_index)
                new_string(keyword_index) = {num2str(variables{keywordi})};
            end
        end
    end 
    
    % Now concatenate everything into a single string.
    if iscell(new_string)
        new_string=horzcat(new_string{:}); 
    end
end 