% HemoRegression.m
% Sarah West
% 9/12/21

% Runs hemodynamic correction on a 2-channel stack by regressing each blue
% pixel against the corresponding violet pixel and keeping only the
% residuals.

function [residual_data, intercepts_data]=HemoRegression(bData, vData)
   % Inputs: 
   % bData--  a 2D matrix of all the blue channel data in the stack; pixels x frames
   % vData--  a 2D matrix of all the violet channel data in the stack; pixels x frames
  
   % Outputs: 
   % residual_data-- a 2D matrix of the residuals left over from the
   %    regression; Is the hemodynamics-corrected data. pixels x frames. 

   % Initialized the holding matrix for residuals. 
   residual_data=NaN(size(bData)); 
   intercepts_data=NaN(size(bData)); 

   % Find lenght of each channel for regression 
   frames=size(vData,2); 
   
   % For each frame,
   parfor i=1:size(bData,1) % parfor
         % The "regress" function needs everything in column format. 
         
         % Create the column of dependent data; Take a pixel from blue
         % data, then flip it into a column. 
         y=bData(i,:)';
         
         % Create the column of predictor data; Take the corresponding
         % pixel from the violet data, then flip it into a column.
         X=vData(i,:)'; 
         
         % Create an intercept vector of ones for the predictor input. Has
         % to be the length of frames. 
         intercepts=ones(frames,1); 
         
         % Combine the X and interceps data into the predictor, run
         % regression.
         [b,~,r] = regress(y, [X intercepts]);

         % Put each residual trace into the corresponding pixel location. 
         residual_data(i,:)=r;
         intercepts_data(i, :) = b(2);
   end
end 
