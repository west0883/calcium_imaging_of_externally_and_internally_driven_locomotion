% RunAnalysis.m 
% Sarah West
% 3/11/21

% A function that abstracts & generalizes all for loops, loading,
% calculations, and saving for analysis done by the Ebner lab. Information
% and parameters need to be formatted properly into the "parameters"
% structure. 
% Functions are given as a cell array of function handles

% Example data to use
% loop_list_1.iterators = {'mouse', {'mice_all(:).name'}, 'mouse_iterator'; 'day', {'mice_all(', 'mouse_iterator', ').days(:).name'}, 'day_iterator'; 'stack', {'mice_all(', 'mouse_iterator', ').days(', 'day_iterator', ').stacks'}, 'stack_iterator'};
% loop_list_1.things_to_load.mask.load_level = 'mouse';
% loop_variables.mice_all = mice_all;
% loop_list_1.save_level = 'stack';
% loop_list_1.things_to_load.data.load_level = 'stack';
% loop_variables.mice_all = mice_all;

function [] = RunAnalysis(functions, parameters)

    % *** Generate list of information to loop through. ***
    
    % Because this is running analysis, a load and save field is required.
    % Give errors if those are not present.
    % load level field
    if ~isfield(parameters.loop_list, 'things_to_load') || isempty(parameters.loop_list.things_to_load)
       error('A non-empty field in loop list called "things_to_load" is required.');
    end 
    % save level field
    if ~isfield(parameters.loop_list, 'things_to_save') || isempty(parameters.loop_list.things_to_save)
       error('A non-empty field in loop list called "things_to_save" is required.');
    end 

    % Make a flag telling any lower-down functions if it was called with
    % RunAnalysis
    parameters.RunAnalysis_flag = true; 
    
    % Grab the number of digits user wants to use for iterating numbers in
    % filenames in LoopGenerator (like with stacks). Otherwise, default to 3 digist.
    if isfield(parameters, 'digitNumber')
       parameters.loop_list.digitNumber = parameters.digitNumber;
    else 
        parameters.loop_list.digitNumber = 3; 
    end

    % Generate list of things to loop through.
    [looping_output_list, maxIterations] = LoopGenerator(parameters.loop_list, parameters.loop_variables);

    % Initialize iterator for looping through looping_output_list
    itemi = 1; 

    % Make a variable holder.
    load_fields = fieldnames(parameters.loop_list.things_to_load);
    variable_in = cell(numel(load_fields),1);

    save_fields = fieldnames(parameters.loop_list.things_to_save);
    %variable_out = cell(numel(save_fields),1);

    % For each item in the list of information to loop through, (Use a
    % while loop so you can have functions skip iterations if needed).
    while itemi <= size(looping_output_list, 1)

        % Make a default abort flag, which allows continuing to next item
        abort = false; 

        % Make a default don't save flag
        parameters.dont_save = false;

        % Get this list of loading and saving string-creating parameters.keywords and
        % variables
        
        if ~iscell(parameters.loop_list.iterators) && strcmp(parameters.loop_list.iterators, 'none')
            parameters.keywords = {'none'};
            parameters.values = {NaN};
        else
            % Keywords should be the names of each iterator, which are in the
            % first column of iterators cell. Also include the iterator names.
            parameters.keywords = [parameters.loop_list.iterators(:,1); parameters.loop_list.iterators(:,3)];

            % Values are the corresponding values in the looping output list
            % for each keyword's field.
            parameters.values = cell(size(parameters.keywords));
            for i = 1: numel(parameters.keywords)
                parameters.values{i} = looping_output_list(itemi).(cell2mat(parameters.keywords(i)));
            end
        end

        % *** Loading ***
        % Check each potential thing to load
        retrieved_value = cell(numel(load_fields),1); 
        
        for loadi = 1:numel(load_fields)
            
            % Figure out if that item should be loaded
            load_flag = looping_output_list(itemi).([load_fields{loadi} '_load']);
            
            % If it should be loaded, start loading process
            if load_flag
             
                % Clear any value still hanging out in this cell of
                % variable_in
                variable_in{loadi} = [];

                % Get the filename & input variable name formatting cells

                dir_cell = parameters.loop_list.things_to_load.(load_fields{loadi}).dir;
                filename_cell = parameters.loop_list.things_to_load.(load_fields{loadi}).filename;
                variable_cell = parameters.loop_list.things_to_load.(load_fields{loadi}).variable;

                input_dir = CreateStrings(dir_cell, parameters.keywords, parameters.values);
                filename = CreateStrings(filename_cell, parameters.keywords, parameters.values);
                variable_string = CreateStrings(variable_cell, parameters.keywords, parameters.values);
                
                % If there's a * in the filename, use a dir function to
                % search for it. 
                if contains(filename, '*')

                    filename_structure = dir([input_dir filename]);

                    % If it could not find the data,
                    if isempty(filename_structure)
                        warning(['Aborting analysis. Could not find file with name ' input_dir filename ]);
                        abort = true;
                        break
                    else
                        filename = filename_structure(1).name;
                    end
                end
            
                % Load 
                
                % Make sure file exists
                if isfile([input_dir filename])

                    % If there is a load function for this field, use that
                    this_load_item = parameters.loop_list.things_to_load.(load_fields{loadi}); 
                   
                    if isfield(this_load_item, 'load_function')
                        
                        load_function = this_load_item.load_function; 

                        % If there are inputs to put in,
                        if isfield(this_load_item, 'load_function_additional_inputs')
                           
                            %eval(['retrieved_value{loadi} = load_function(''' input_dir filename ''',' this_load_item.load_function_additional_inputs ');']);
                            retrieved_value{loadi} = load_function([input_dir filename], this_load_item.load_function_additional_inputs{:}); 

                        else 
                            retrieved_value{loadi} = load_function([input_dir filename]); 
                        end 
                    else 
                        % Check if "variable" has a period in it (and therefore
                        % would become a structure).
                        if contains(variable_string, '.')
                            
                            % Find location of the period & only use the string
                            % up to that point to load -- because that's how
                            % loading structures works in Matlab. 
                            index = find(variable_string == '.', 1);
                            variable_in{loadi} = load([input_dir filename], variable_string(1:index-1)); 
                        elseif contains(variable_string, '{')
                            
                            % Find location of the period & only use the string
                            % up to that point to load -- because that's how
                            % loading structures works in Matlab. 
                            index = find(variable_string == '{', 1);
                            variable_in{loadi} = load([input_dir filename], variable_string(1:index-1)); 

                         elseif contains(variable_string, '(')
                            
                            % Find location of the period & only use the string
                            % up to that point to load -- because that's how
                            % loading structures works in Matlab. 
                            index = find(variable_string == '(', 1);
                            variable_in{loadi} = load([input_dir filename], variable_string(1:index-1)); 
    
                        else
                            variable_in{loadi}= load([input_dir filename], variable_string); 
                           
                        end 
                    end
                    
                else
                    % If there was no existing file & user didn't specify
                    % what to do, abort this iteration (default is to
                    % abort).
                    if isfield(parameters, 'load_abort_flag') && ~parameters.load_abort_flag

                        % If no file, report (sometimes we want this).
                        warning(['No file for ' load_fields{loadi} ' found at ' input_dir filename]);
 
                    else
                        % Say to abort this item and continue to next item
                        abort = true; 
                        error(['Aborting analysis: No file for ' load_fields{loadi} ' found at ' input_dir filename ]);
                        
                        % Break out of loading loop
                       % break 
                    end
                end 
            end
        end
        
        % If a load abort flag was given, skip to next item
        if abort 
           % Go to the next item i value where there's loading for the
           % failed load field. Loadi is whatever it was when the load look
           % broke.
           holder = [looping_output_list(itemi + 1:end).([load_fields{loadi} '_load'])];
           
           % If there's no next one, finish the loop 
           if isempty (find(holder, 1))
              break
           else
               % If there is a next one, make itemi equal that.
               itemi = itemi + find(holder, 1);
           end
           continue 
        end
        
        % Pull out loaded variable here, so you can still iterate below the
        % load level. 
        for loadi = 1:numel(load_fields)        

            % Get this set of load fields.
            this_load_item = parameters.loop_list.things_to_load.(load_fields{loadi}); 
            
             % If this doesn't need to be loaded in specially 
             if ~isfield(this_load_item, 'load_function')  
                
                % Only if the item was loaded & this load field  of variable_in isn't empty.
                if  ~isempty(variable_in{loadi}) % load_flag &&
                    % Skip if there was a special load function because retrieved
                    % value was already defiened. 
                    this_load_item = parameters.loop_list.things_to_load.(load_fields{loadi});
        
                    variable_cell = parameters.loop_list.things_to_load.(load_fields{loadi}).variable;
                    variable_string = CreateStrings(variable_cell, parameters.keywords, parameters.values);
                    
                    %retrieved_value{loadi} = variable_in{loadi}.(variable_string);
                    % Assign. (Keeping "eval" here because doing dynamic
                    % field names doesn't work for '(', '.', or '{'
                    % sub-indexing.)
                    eval(['retrieved_value{loadi} = variable_in{loadi}.' variable_string ';']);
                end 
             end
           
            % Assign to the specific name in parameters structure 
            parameters.(load_fields{loadi}) = retrieved_value{loadi};
         
        end

        % Create continue flags for each level of iterator, so functions
        % run by user can tell RunAnalysis to skip only certain levels.
        parameters.continue_flag = repmat({true}, size(parameters.loop_list.iterators,1), 1);

        % Pull out the max iteration values that can be reached within
        % each loop level at this loop (i.e, in this mouse, how many days
        % are there). Makes it easier for called functions to throw
        % errors/know where they are in analysis without accessing loop_list.  
        if ~iscell(parameters.loop_list.iterators) && strcmp(parameters.loop_list.iterators, 'none')
        else 
            parameters.maxIterations = [];
            parameters.maxIterations.numbers_only = [];
            for i = 1:size(parameters.loop_list.iterators,1)
                
                parameters.maxIterations.(parameters.loop_list.iterators{i,3}) = maxIterations(itemi, i);
                parameters.maxIterations.numbers_only = [parameters.maxIterations.numbers_only maxIterations(itemi, i)];
            end
        end

        % *** Run functions chosen by user. *** 
        
     
        % For each function, (run recursively on "parameters" structure).
        for functioni = 1:numel(functions)
            
            % Also pass functioni to functions, if needed.
            parameters.functioni = functioni;
            
            % If this is the second or more function & there are things to hold onto, 
            if functioni > 1 
                
                % If there are things to hold onto,
                if isfield(parameters.loop_list, 'things_to_hold')
                    holder_data = struct; 
                    % For each of the items in this level of 'things_to_hold'
                    for holdi = 1:numel(parameters.loop_list.things_to_hold(functioni - 1,:))
                        
                        hold_name = parameters.loop_list.things_to_hold{functioni - 1,holdi};
                        holder_data.(hold_name) = parameters.(hold_name); 
                       
                    end 
                end 
                
                % Rename data that's supposed to be renamed.
                if isfield(parameters.loop_list, 'things_to_rename')
                    for renamei = 1:size(parameters.loop_list.things_to_rename{functioni - 1},1)
                        
                        parameters.(parameters.loop_list.things_to_rename{functioni - 1}{renamei, 2}) = parameters.(parameters.loop_list.things_to_rename{functioni - 1}{renamei, 1});
            
                    end 
                end
            end 

            % Assign the function for this step
            F = functions{functioni};

            % Run the function
            parameters = F(parameters); 

            % If this is the second or more function & there are things to hold onto, 
            % put held data back into position where it was before this
            % function. 
            if functioni > 1 && isfield(parameters.loop_list, 'things_to_hold')
                
                % For each of the items in this level of 'things_to_hold'
                for holdi = 1:numel(parameters.loop_list.things_to_hold(functioni - 1,:))
                    
                    hold_name = parameters.loop_list.things_to_hold{functioni - 1,holdi};
                    parameters.(hold_name) = holder_data.(hold_name);
                   
                end 
            end 
        end

        % *** Save  ***
        
        % For each potential thing to save
        for savei = 1:numel(save_fields)

            % Assign to variable strings now (so you can iterate between
            % saves). 

            % A called function can make the dont_save field into a cell
            % array, one entry for each save field.
            if iscell(parameters.dont_save)
                dont_save = parameters.dont_save{savei};
            else
                % Otherwise, just pull out dont_save
                dont_save = parameters.dont_save;
            end
           
            % Find out if you save at this level.
            save_flag = looping_output_list(itemi).([save_fields{savei} '_save']);

            % If you save at this level, or if the called function sends
            % out the "save_now flag", 
            if ~dont_save && (save_flag || (isfield(parameters, 'save_now') && parameters.save_now))
                
                % Get data out of parameters structure
                variable_out = parameters.(save_fields{savei});
                variable_cell = parameters.loop_list.things_to_save.(save_fields{savei}).variable;
                variable_string = CreateStrings(variable_cell, parameters.keywords, parameters.values);

                % Convert to non-generic variable name
                eval([variable_string ' = variable_out;']);

                clear_variable_out = true;
                
                % Create strings for all saving info
                dir_cell = parameters.loop_list.things_to_save.(save_fields{savei}).dir;
                filename_cell = parameters.loop_list.things_to_save.(save_fields{savei}).filename;
               
                output_dir = CreateStrings(dir_cell, parameters.keywords, parameters.values);
                filename = CreateStrings(filename_cell, parameters.keywords, parameters.values);
                
                % Make output directory, if it doesn't already exist.
                if ~exist(output_dir, 'dir')
                    mkdir(output_dir);
                end
         
                
                % Check if variable is a figure handle. First have to make
                % sure it's a graphics object handle specifically or else 
                % Matlab freaks out. Making sure it's specifically a figure
                % is good to be sure someone didn't make a typo.
                if  strcmp(class(variable_out), 'matlab.ui.Figure')
                    % ~ismatrix(variable) && isgraphics(variable) 
                    % If it is a figure, save as figure. Uses the "compact"
                    % fig format, apparently makes it faster & smaller
                    % file.
                    savefig(variable_out, [output_dir filename],'compact'); 

                    % If user wants to save a fancy version (in addition to
                    % .fig type)

                    % Get out the loop_list info 
                    fig_list_variable = parameters.loop_list.things_to_save.(save_fields{savei});
                    
                    % Check if user put in a saveas type. 
                    if isfield(fig_list_variable, 'saveas_type')

                        % Get the save type
                        fig_type = fig_list_variable.saveas_type;
                        
                        % Save the fancy figure, too. 
                        saveas(variable_out, [output_dir filename], fig_type);

                    end 

                   
                else
                    % If not a figure, save as variable. 
                    
                    % If variable is a structure, see if there's a field in loop list to save
                    % everything as variables (default is as structure)
                    
                    % Get the relevant sub-structure (have to do it like
                    % this because isfield doesn't nest)
                    loop_list_variable = parameters.loop_list.things_to_save.(save_fields{savei});
                  
                    if isstruct(variable_out) && isfield(loop_list_variable, 'save_as_variables') && loop_list_variable.save_as_variable

                            clear_variable_out = false;

                            % Save with structure-to-variables feature
                            save([output_dir filename], '-struct', 'variable_out'); 
                            

                            % To make it even more flexible... but I'll
                            % leave this for later.
                            % Get out the list of variables they want to
                            % save these as. Do it as just a cell array,
                            % I'm tired. 
%                           % subvariables  = getfield(list_variable);
%                           %for i = 1:size(subvariables,1)
%                           %      eval([subvariables{i,1} ' = ' subvariables{i,2} ';']);
%                           % end     
                    % If the variable name shows it's a cell array
                    elseif contains(variable_string, '{')
                            
                            % Find location of the period & only use the string
                            % up to that point to load -- because that's how
                            % loading structures works in Matlab. 
                            index = find(variable_string == '{', 1);
                            variable_string = variable_string(1:index-1); 

                            clear_variable_out = false;
                        
                        % Save
                        save([output_dir filename], variable_string, '-v7.3'); 

                    % If the variable name shows it's a normal array
                    elseif contains(variable_string, '(')
                            
                            % Find location of the period & only use the string
                            % up to that point to load -- because that's how
                            % loading structures works in Matlab. 
                            index = find(variable_string == '(', 1);
                            variable_string = variable_string(1:index-1); 
                            clear_variable_out = false;

                        % Save
                        save([output_dir filename], variable_string, '-v7.3'); 

                    % Or if the variable name shows it's a structure
                    elseif contains(variable_string, '.')
                            
                            % Find location of the period & only use the string
                            % up to that point to load -- because that's how
                            % loading structures works in Matlab. 
                            index = find(variable_string == '.', 1);
                            variable_string = variable_string(1:index-1); 
                            clear_variable_out = false;

                        % Save
                        save([output_dir filename], variable_string, '-v7.3');    

                    else 
                       % Save 
                        save([output_dir filename], variable_string, '-v7.3'); 

                    end 
                end

                % Once you've saved it, remove that variable from parameters 
                % structure so called functions don't over write or add to it.
                parameters = rmfield(parameters,save_fields{savei});

                % Also clear the varible itself (relevant if you dealt with
                % elements of an array or cell array).
                if clear_variable_out 
                    clear(variable_string);
                end 
                
            end
        end
        
        % If user said they want to close figures, 
        if isfield(parameters, 'closeFigures') && parameters.closeFigures
            
            % Close any open figures
            close all

        end 

        % If the run functions return a continue flag & the continue flag
        % on any iterator level is false
        if isfield(parameters, 'continue_flag') && any(~[parameters.continue_flag{:}])
           
            % Get the first (highest-level) iterator-level that's false
            iterator_level_to_skip = find(~[parameters.continue_flag{:}],1);

            % Go to next iteration of level ABOVE that. 
            next_iterator_level = iterator_level_to_skip - 1; 

            % If there are no higher levels (iterator level is 0), then
            % you're done! Break. 
            if next_iterator_level <= 0
                break
            end

            % Get the current value of that iterator 
            current_iterator_to_skip = looping_output_list(itemi).(parameters.loop_list.iterators{next_iterator_level, 3});

            % In looping list, find when that iterator increases by 1
            holder = looping_output_list(:).(parameters.loop_list.iterators{next_iterator_level, 3}); 
            next_iteration = find(holder == current_iterator_to_skip +1,1); 
            
            % If "next_iteration" isn't empty
            if ~isempty(next_iteration)
                
                % Make itemi equal that location in the looping list.
                itemi = next_iteration; 
            
            % If it is empty, break out of item loop & end RunAnalysis
            else
                break
            end
        else
            % If no skipping, just increase the item iterator by 1.
            itemi = itemi + 1;
        end
    end 
end 