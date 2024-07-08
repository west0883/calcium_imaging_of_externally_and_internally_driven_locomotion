% SearchForData.m
% Sarah West
% 6/21/22

% Function that searches through your iterators for the inputted file
% names. Does NOT use RunAnalysis, but has similar input structure. 

function [] = SearchForData(parameters)

    % Make a loop list from iterators.
    looping_output_list = LoopGenerator(parameters.loop_list, parameters.loop_variables);

    % Initialize missing data cell array. 
    missing_data = cell(1, size(parameters.loop_list.iterators,1) + 1);
    
    % For each item in looping output list,
    for itemi = 1:size(looping_output_list,1)
    
        % Get this list of loading and saving string-creating parameters.keywords and
        % variables
        
        % Keywords should be the names of each iterator, which are in the
        % first column of iterators cell. Also include the iterator names.
        parameters.keywords = [parameters.loop_list.iterators(:,1); parameters.loop_list.iterators(:,3)];

        % Values are the corresponding values in the looping output list
        % for each keyword's field.
        parameters.values = cell(size(parameters.keywords));
        for i = 1: numel(parameters.keywords)
            parameters.values{i} = looping_output_list(itemi).(cell2mat(parameters.keywords(i)));
        end
        
        % Get the file names of the files you're checking for
        dir_cell = parameters.loop_list.things_to_check.dir;
        filename_cell = parameters.loop_list.things_to_check.filename;
   
        input_dir = CreateStrings(dir_cell, parameters.keywords, parameters.values);
        filename = CreateStrings(filename_cell, parameters.keywords, parameters.values);
      
        % Search for data file.
        filename_structure = dir([input_dir filename]);

        % If it could not find the data, or if there was more than one file
        % matching,
        if isempty(filename_structure) 
            missing_data = [missing_data; parameters.values(1:end/2)' {'missing'}];

        elseif size(filename_structure, 1) > 1
            missing_data = [missing_data; parameters.values(1:end/2)' {'duplicates'}];
        end

    end
    % Put missing data into output structure (so you can change the
    % variable name, if you want). 
    parameters.missing_data = missing_data;

    % Get output file name strings.
    dir_cell = parameters.loop_list.missing_data.dir;
    filename_cell = parameters.loop_list.missing_data.filename;
   
    output_dir = CreateStrings(dir_cell, parameters.keywords, parameters.values);
    filename = CreateStrings(filename_cell, parameters.keywords, parameters.values);
  
    % Save
    save([output_dir filename], 'missing_data', '-v7.3');
end 