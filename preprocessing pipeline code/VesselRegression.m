% VesselRegression.m
% Sarah West
% 9/14/21

% Runs hemodynamic correction by regressing pixels against blood vessel
% timeseries.
% RIGHT NOW THE BLOOD VESSEL MASKS AND THE NOW-MASKED bDATA ARE DIFFERENT
% SIZES 
function [corrected_data]= VesselRegression(bData, vessel_masks)
% Inputs:
% bData-- a 2D matrix of all the data in the stacks 

    % Extract blood vessel timecourse averages  
    all_vessel_timecourses=bData'*vessel_masks ./sum(vessel_masks,1);           

    % Establish empty output matrix
    corrected_data=NaN(size(bData));

    % Make a column of ones for intercept (length is equal to number of
    % frames)
    intercept=ones(size(bData,2),1); 

    % For every pixel in the blue data,
    for i=1:size(bData,1)

         % Regress the data matrix against the vessel timecourses.
         [~,~,r] = regress(bData(i,:)', [all_vessel_timecourses intercept]);

         % Keep the residuals of the regression. 
         corrected_data(i,:)=r;
    end
end 

