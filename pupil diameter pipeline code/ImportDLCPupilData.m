% ImportDLCPupilData.m
% Sarah West
% 6/14/22

% Imports the pupil edges data from DeepLabCut .csv files. Run with RunAnalysis.
% Is really just a holder function for RunAnalysis to call while it does
% the actual importing.
function [parameters] = ImportDLCPupilData(parameters)

    MessageToUser('Importing ', parameters);

    parameters.import_out = parameters.import_in;


end 