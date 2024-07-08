% FixImageSize.m
% Sarah West
% 12/3/21

% Removes the extra pixel in each dimension sometimes captured in
% microManager on 472 rig (no clue how that can happen). 

% Input: 
% data = the matrix of images you want to adjust. First 2 dimensions must
% be pixels by pixels.
% pixels = the y and x dimensions that you want, in a 2-element array. (ie,
% [256, 256]).

function [data] = FixImageSize(data, pixels)
 
   % ***Check if image of stack are the right size. Use just
    % the first frame. ***
    image_dimensions = size(data); 
    
    for i = 1:numel(pixels)
        % Get a list of dimensions, repeat the ':'.
        inds = repmat({':'}, 1, ndims(data));

        % Check fist dimension
        if image_dimensions(i) ~= pixels(i) 

            % Only do this step if the pixels are 1 more than
            % expected, else display a message and continue to next
            % stack.
            if image_dimensions(i) - pixels(i) == 1

                disp(['Fixing image dimension '  num2str(i)]); 

                % Remove very last pixel
                inds{i} = image_dimensions(i); 
                data(inds{:}) = [];

            else 
                error(['Dimension ' num2str(i) ' is wrong size']);  
            end

        end
    end 
end