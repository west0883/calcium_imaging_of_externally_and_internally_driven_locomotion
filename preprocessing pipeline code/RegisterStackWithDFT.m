% RegisterStackWithDFT
% Sarah West
% 8/26/21

% Runs registration within stacks and across stacks within the same day to
% a representative image of that day. Uses dftregistration function. 

function [all_tforms]=RegisterStackWithDFT(reference_image, stack_to_register, usfac)

    % Inputs: 
    % reference_image -- the representative image (of that recording day) that
    % you're registering all frames of all stacks to. MUST BE BLUE CHANNEL
    % stack_to_register -- the data matrix of a stack to register. 3D (pixels,
    % pixels, frames). MUST BE
    % BLUE CHANNEL
    % usfac -- "upsampling factor" for dftregistration function (for within stack & 
    % day registration); determines the sub-pixel resolution of the registration; 

    % Outputs:
    % registered_stack -- the stack with each frame registered; 3D (pixels,
    % pixels, frames)
    % all_tforms -- all the shifts that were applied; 3D (rows shifted, cols 
    % shifted, frames)

    % Run 2D Fourier transform on reference image  
    fRef=fft2(reference_image); 

    % Find number of frames of stack
    frames=size(stack_to_register,3); 
    
    % Initialize output matrices
    % The "all_tforms" will hold the diffphase, row shift, and column shift for each frame 
    all_tforms=NaN(3,frames);
  
    % Apply the 2D Fourier transform to the frame/image.
    fim=fft2(stack_to_register);
    
    % Use dftregistration function to align the fourier transform of
    % the current BLUE image with the fourier transform of the
    % background/reference image of the reference BLUE image
    
    %for each frame of the stack as indexed by variable t
    parfor t=1:frames
        [output ] = dftregistration(fRef, fim(:,:,t), usfac);

        % Grab variables from "output" that you'll need to find the registered VIOLET image
        % From dftregistration code: output=[error,diffphase,row_shift,col_shift];
        all_tforms(:, t)=output(2:4); 
    end 
end 