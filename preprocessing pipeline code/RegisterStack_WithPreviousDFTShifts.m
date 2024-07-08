% RegisterStack_WithPreviousDFTShifts.m
% Sarah West
% 8/26/21

% Takes the frame-by-frame transforms previously calculated by a
% dftregistration and applies it to a stack. (dftregistration will give
% this to you as Greg variable as you're calculating transforms, but this
% code allows you to apply the transforms to a new stack. This is really
% only important when you're aligning violet channel images based on the
% transforms calculated on the blue images--> for accurate hemodynamics
% correction.) The calculations were copied from dftregistration. 

% Inputs:
% tforms -- Is an output from RegisterStackWithDFT. Is a 2D matrix, with
    % (first dim) diffphase,row_shift,col_shift for each frame (2nd dim).
% stack_to_register -- a (violet channel) stack matrix to register. 3D
    % matrix (pixels, pixels,frames)
% usfac -- "upsampling factor"; Determines the sub-pixel resolution of the registration; 

% Outputs:
% registered_stack -- the registered stack. 3D matrix (pixels, pixels,
    % frames)   
    
function [registered_stack] =RegisterStack_WithPreviousDFTShifts(tforms, stack_to_register, usfac, yDim, xDim) 
    
    % Find number of frames of stack
    frames=size(stack_to_register,3);

    % Create a Fourier-transformed version of the image you want to
    % register, keeping with the name of the variable used in the
    % dftregistration code to try to keep comparisons simple
    buf2ft=fft2(stack_to_register); 

    clear stack_to_register;

    % calculate violet registered image
    if (usfac > 0)

        % find dimensions of images needed for the calulations
        % (from dftregistration.m)
        [nr,nc]=size(buf2ft, [1 2]);
        Nr = ifftshift(-fix(nr/2):ceil(nr/2)-1);
        Nc = ifftshift(-fix(nc/2):ceil(nc/2)-1);
        [Nc,Nr] = meshgrid(Nc,Nr);
        Nc = repmat(Nc, [1 1 frames]);
        Nr = repmat(Nr, [1 1 frames]);

        row_shift = permute(repmat(tforms(2,:), [yDim 1 yDim]), [1 3 2]);
        col_shift = permute(repmat(tforms(3,:), [xDim 1 yDim]), [1 3 2]);

        % First step
        registered_stack = buf2ft .* exp(1i*2*pi*(-row_shift.*Nr/nr-col_shift.*Nc/nc));

        clear row_shift col_shift Nr Nc;
        
        % Second step
        diffphase = permute(repmat(tforms(1,:), [yDim 1 yDim]), [1 3 2]);
        registered_stack = registered_stack.*exp(1i.*diffphase);

    elseif (usfac == 0)
        diffphase = permute(repmat(tforms(1,:), [yDim 1 yDim]), [1 3 2]);
        registered_stack = buf2ft .* exp(1i*diffphase);
    end
    clear diffphase;

    %Get the absolute value of the inverse fourier transform of the
    %registered images.
    registered_stack = abs(ifft2(registered_stack));
        
end 