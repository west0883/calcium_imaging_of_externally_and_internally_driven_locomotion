% CheckIfCorruptFast.m
% Sarah West
% 9/2/21
% Checks if a file/stack is corrupt WITHOUT HAVING TO LOAD IT :D
% Basis code is from https://stackoverflow.com/questions/38325121/check-if-mat-file-is-corrupt-without-load
% matfile doesn't need to load in the whole file, and it issues a warning
% if the file is corrupt. 

% Warning ID is : 'MATLAB:whos:UnableToRead'

function [isCorrupt]=CheckIfCorruptFast(filename)

% Inputs:
% filename -- the name of the file you want to check for corruption. Must
% be a .mat file. 

% Outputs:
% isCorrupt-- a true/false flag that says if the file is corrupt (true) or
% not (false).
    
    % Clear the last warning 
    lastwarn(''); 

    % Attempt to use the matfile function on the file.
    matfile(filename);
    
    % Get the last warning.
    [~, warnId] = lastwarn;

    % If the last warning is the 'UnableToRead' warning, assume file is
    % corrupt
    if strcmp(warnId, 'MATLAB:whos:UnableToRead')
        % File is corrupt, mark it 
        isCorrupt=true;
    else
        % If the last warning is anything but the 'UnableToRead' (including
        % the blank cleared warning), assume is not corrupt and mark it.
        isCorrupt=false;
    end
end 