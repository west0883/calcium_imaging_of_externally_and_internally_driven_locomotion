% CountRolls.m
% Sarah West
% 5/28/22

% Returns the number of data "rolls" (sliding/rolling window windows) that
% can come out of a data period, based on "duration."

% Inputs:
% duration -- the number of time points in a time period.
% windowSize -- the number of time points the window is across.
% stepSize -- the number of time points the window is moved between
    % rolls/counts.

% Outputs:
% roll_number -- the total number of sliding window instances that can fit
    % in this duration.

function [roll_number] = CountRolls(duration, windowSize, stepSize)

    roll_number = (duration - windowSize)/stepSize + 1;

    % If not an integer, round down to nearest integer (preserves any
    % existing integers)
    roll_number = floor(roll_number);

    % Any roll number below 1 (including any weird negative numbers) will be 1. 
    if roll_number < 1
       roll_number = 1;
    end
    
end
