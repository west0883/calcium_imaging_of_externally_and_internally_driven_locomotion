% CreateFileStrings.m
% Sarah West
% 9/1/21

% Takes an input cell array of character strings and outputs the file
% string with the correct mouse and day in them. 

% Input: 
% file_format_cell--cell array. Establish the format of the file names of compressed data. Each piece
% needs to be a separate entry in a cell array. Put the string 'mouse', 'day',
% or 'stack number' where the mouse, day, or stack number will be. If you 
% concatenated this as a sigle string, it should create a file name, with the 
% correct mouse/day/stack name inserted accordingly. 
% 
% searching -- true/false Boolean. Inserts * into the place of each entry
% instead of a changing number. If true, all other entries should be empty

function [file_string]=CreateFileStrings(file_format_cell, mouse, day, stack_number, period, searching)
    
    % Make a new cell array to manipulate. 
    file_format_output_cell=file_format_cell;
    
    % See if there is an entry for mouse number, find where it's located
    mouse_index=find(strcmp(file_format_cell,'mouse number'));
    
    % If not empty
    if isempty(mouse_index)==0
         
         % If there is, make sure the mouse entry isn't empty 
         if isempty(mouse)==0

             % Put the mouse number in place of the mouse number tag
             file_format_output_cell(mouse_index)={mouse}; 

         % If the mouse input is empty, and you're searching
         elseif searching
            % We're creating a file string for searching directories, put a * in this place.
            file_format_output_cell(mouse_index)={'*'}; 

         % If mouse input is empty, but you're not searching throw an error
         else 
             error('no mouse number was given'); 
         end
         
    end
     
    % See if there is an entry for day, find where it's located 
     day_index=find(strcmp(file_format_cell,'day'));
     
     if isempty(day_index)==0
          
         % If there is, make sure the day entry isn't empty 
         if isempty(day)==0

             % Put the day number in place of the day number tag
             file_format_output_cell(day_index)={day}; 
         
         % If it's empty, but we're searching 
         elseif searching
             
            % We're creating a file string for searching directories, put a * in this place.
            file_format_output_cell(day_index)={'*'}; 
         
         % If it's empty & we're not creating this for searching, throw an
         % error
         else 
            error('no day was given'); 
             
         end
     end 
     
     % See if there is an entry for stack number, find where it is 
     stack_index=find(strcmp(file_format_cell,'stack number'));
     
     if isempty(stack_index)==0
         
         % If we're creating a file string for searching directories, put a * in this place.
         file_format_output_cell(stack_index)={'*'};  

         % If there is, make sure the stack number entry isn't empty 
         if isempty(stack_number)==0

             % Put the mouse number in place of the mouse number tag
             file_format_output_cell(stack_index)={stack_number}; 

         % If the stack input is empty and we're searching
         elseif searching

            % We're creating a file string for searching directories, put a * in this place.
            file_format_output_cell(stack_index)={'*'}; 

         % If the stack is empty, but we're not searching
         else
             error('no stack number was given'); 
         end 
        
     end 
     
     % See if there is an entry for period, find where it is 
     period_index=find(strcmp(file_format_cell,'period name'));
     
     if isempty(period_index)==0
    
         % If there is, make sure the stack number entry isn't empty 
         if isempty(period)==0

             % Put the mouse number in place of the mouse number tag
             file_format_output_cell(period_index)={period}; 

         % If the period input is empty and we're searching, 
         elseif searching
            
             % We're creating a file string for searching directories, put a * in this place.
             file_format_output_cell(period_index)={'*'}; 
         
         % If period input is empty, but we're not searching, throw an
         % error
         else 
             error('no period name was given'); 
         end 
     end 
    
    % Now concatenate everything into a single string.
    file_string=horzcat(file_format_output_cell{:}); 
end 