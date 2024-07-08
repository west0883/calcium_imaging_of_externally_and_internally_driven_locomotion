% CheckSizes.m
% Sarah West
% 6/21/22

% Function that searches through your iterators for the inputted file
% names. Does NOT use RunAnalysis, but has similar input structure. 

function [] = CheckSizes(parameters)

    % Make a loop list from iterators.
    looping_output_list = LoopGenerator(parameters.loop_list, parameters.loop_variables);

    % Initialize missing data cell array. 
    mismatched_data = cell(1, size(parameters.loop_list.iterators,1) + 1);
    
   % try 
    % For each item in looping output list,
    for itemi = 1:size(looping_output_list,1)
    
        % Get this list of loading and saving string-creating parameters.keywords and
        % variables
        
        % Keywords should be the names of each iterator, which are in the
        % first column of iterators cell. Also include the iterator names.
        parameters.keywords = [parameters.loop_list.iterators(:,1); parameters.loop_list.iterators(:,3)];

        % Values are the corresponding parameters.values in the looping output list
        % for each keyword's field.
        parameters.values = cell(size(parameters.keywords));
        for i = 1: numel(parameters.keywords)
            parameters.values{i} = looping_output_list(itemi).(cell2mat(parameters.keywords(i)));
        end

        parameters.RunAnalysis_flag = true;
        MessageToUser('Checking ', parameters);
        
        % Get the file names of the files you're checking for
        dir_cell = parameters.loop_list.things_to_check.dir;
        filename_cell = parameters.loop_list.things_to_check.filename;
        variable_cell = parameters.loop_list.things_to_check.variable;
   
        input_dir = CreateStrings(dir_cell,parameters.keywords, parameters.values);
        filename = CreateStrings(filename_cell,parameters.keywords, parameters.values);
        variable_string_checking = CreateStrings(variable_cell,parameters.keywords, parameters.values);

        % Create object for data file.
        if isfile([input_dir filename])
            if parameters.use_substructure
                file_object_checking = load([input_dir filename]);
            else 

                file_object_checking = matfile([input_dir filename]);
            end 
        else
            if isfield(parameters, 'values')
                holder = strjoin(parameters.values(1:numel(parameters.values)/2), ', ');
            else 
                holder = 'this iteration.';
            end 
            disp(['No file for ' holder]);
            continue
        end

        % Get filenames & object for data to check against
        dir_cell = parameters.loop_list.check_against.dir;
        filename_cell = parameters.loop_list.check_against.filename;
        variable_cell = parameters.loop_list.check_against.variable;
   
        input_dir = CreateStrings(dir_cell,parameters.keywords, parameters.values);
        filename = CreateStrings(filename_cell,parameters.keywords, parameters.values);
        variable_string_check_against = CreateStrings(variable_cell,parameters.keywords, parameters.values);

        if isfile([input_dir filename])
            file_object_check_against = matfile([input_dir filename]);
        else
            continue
        end

        % Get the condition name 
        location = strcmp(parameters.keywords, {'condition'});
        condition = parameters.values{location};

        if strcmp(condition, 'motorized') 
            periods = parameters.periods_motorized;

        else
            periods = parameters.periods_spontaneous;
        end

        % Pull out the variable you want from file_object_checking into
        eval(['variable_checking = file_object_checking.' variable_string_checking ';']);

        % For each period you want to look at, 
        for periodi = 1:size(periods,1)

            % if these files don't exist, skip 
            
            subholder_checking = variable_checking(periodi,1);
            subholder_check_against = file_object_check_against.segmented_timeseries(periodi,1);

            % If they're not empty,
            if ~isempty(subholder_check_against{1})
                % Compare the sizes. If they don't match,
                if size(subholder_checking{1}, 1) ~= size(subholder_check_against{1}, 1) || size(subholder_checking{1}, 2) ~= size(subholder_check_against{1}, 3) 
    
                    disp('Found mismatch');
    
                    % Add to list of mismatched 
                    mismatched_data = [mismatched_data; parameters.values(1:end/2)' {periodi}];
    
                end
            end
        end
    end
   
%     catch
% 
%         % Get output file name strings.
%         dir_cell = parameters.loop_list.mismatched_data.dir;
%         filename_cell = parameters.loop_list.mismatched_data.filename;
%        
%         output_dir = CreateStrings(dir_cell,parameters.keywords, parameters.values);
%         filename = CreateStrings(filename_cell,parameters.keywords, parameters.values);
%       
%         % Save
%         save([output_dir filename], 'mismatched_data', '-v7.3');
%     end


        % Get output file name strings.
        dir_cell = parameters.loop_list.mismatched_data.dir;
        filename_cell = parameters.loop_list.mismatched_data.filename;
       
        output_dir = CreateStrings(dir_cell,parameters.keywords, parameters.values);
        filename = CreateStrings(filename_cell,parameters.keywords, parameters.values);
      
        % Save
        save([output_dir filename], 'mismatched_data', '-v7.3');

end 