% DetermineChannel.m
% Sarah West
% 8/25/21

% Takes two images (must be 2 sequential images from a 2-channel
% recording and determines which belongs to which channel depending on the
% brightness in a region of interest over the brain. 

% blue_brighter -- a boolean. Is true if the blue channel is expected to be
% brighter than the 2nd, violet channel. Rarely, the blue channel is dimmer
% because of an LED setting error. 
function [first_image_channel] = DetermineChannel(blue_brighter, im1, im2, pixel_rows, pixel_cols)
       
        % get the mean gray level from a brain ROI 
        lev1=mean2(im1(pixel_rows, pixel_cols));
        lev2=mean2(im2(pixel_rows, pixel_cols));

        if lev1~=lev2 % if lev1 doesn't equal lev2 (i.e. the mean gray levels
            %are not equal between im1 and im2
            
            % If the blue should be brighter
            if blue_brighter
                % If lev1 is greater than lev2, then it's the blue channel
                if lev1 > lev2 
                    first_image_channel='b'; 

                else
                   % If lev1 is less than lev2, then it's the violet channel
                    first_image_channel='v'; 
                end
            
            % If the violet should be brighter.
            else
                % If lev1 is less than lev2, then it's the blue channel
                % (Oppositete from above) 
                if lev1 < lev2 
                    first_image_channel='b'; 

                else
                   % If lev1 greater than lev2, then it's the violet channel
                    first_image_channel='v'; 
                end
                
            end 
        else
           error('Stack illumination does not differentiate hemo signal.'); 
        end
end