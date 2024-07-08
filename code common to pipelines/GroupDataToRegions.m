
function [] = GroupDataToRegions (parameters) 
  
    mice_all = parameters.mice_all; 
    periods_all = parameters.periods_all; 
    dir_groups = parameters.dir_in_groups; 
    groups_file_name = parameters.file_name_groups;
    groups_variable_name = parameters.groups_variable_name;
    dir_input_base = parameters.dir_input_base;
    input_file_name = parameters.input_file_name;
    input_variable_name = parameters.input_variable_name;
    dir_out_base = parameters.dir_out_base;
    groupedDim = parameters.groupedDim;
    output_file_name = parameters.output_file_name;
    output_variable_name = parameters.output_variable_name; 
    
    % Load list of groups
    load([dir_groups groups_file_name]);

    % For each mouse, 
    for mousei=1:size(mice_all,2)
        mouse=mice_all(mousei).name;
    
        % Get the input directory
        dir_in = [dir_input_base '\' mouse '\'];

        % Get the data that corresponds to that period, assign it
        % a generic name.
        eval(['mouse_instances = ' variable_name_input ';']);  

        % Establish output directory.
        dir_out = [dir_out_base '\' mouse '\']; 
        mkdir(dir_out);

        % Get the right groups list and make it generic.
         eval(['group_list=' groups_variable_name '.m' mouse ';']); 
        
        % For each period, 
        for periodi=1:size(periods_all,1)
            period=periods{periodi};
           
            % Get specific name of input filename.
            file_name_input = CreateFileStrings(input_file_name,[], [], [], period, false);

            % Load the data.
            load([dir_in file_name_input]); 

            % Get specific name of input variable
            variable_name_input = input_variable_name; %CreateFileStrings(input_variable_name,[], [], [], period, false);
          
            % Convert to a general name. 
            eval(['data_all=' variable_name_input ';']);
          
            % Set up matrix to hold newly grouped data, depending on size
            % of the data being grouped and the dimensions being grouped
            dimensions_list = NaN(1,ndim(data_all));
            for dimi = 1:ndim(data_all)
                
                % If this is a dimension to be grouped, 
                if any(groupedDim == dimi)
                   
                    % This dimension will be the length of the group list.
                    dimensions_list(dimi) = size(group_list,1);
                
                else
                    % Otherwise, this dimension will be the same length as
                    % the corresponding dimension of the data.
                    dimensions_list(dimi) = size(data_all, dimi); 
                
                end
            end 
                     
            % Now initialize a matrix with the dimensions you just found.
            data_matrix = NaN(dimensions_list); 
            
            % Complicated-looking stuff that makes sure you're
            % using the correct dimension. (From old matrix).
            inds = repmat({':'},1,ndims(data_all));
            
            % Carry out code based on number of entered dimensions to be
            % grouped. 
            switch length(groupedDim) 
                
                case 1  % Only 1 dimension to be grouped.
                    
                    for groupi1=1:size(group_list,1)
                        group1=group_list{groupi1,2}; 
                        
                        % If the mouse has a region in this group,
                        if isnan(group1)== 0 
                            
                            holder =[];
                       
                            for i=1:numel(group1) 
                                
                                % Change the index in the correct dimension
                                % (from old matrix).
                                inds{groupedDim(1)} = group1(i);
                                
                                % Get the data 
                                holder = [holder data_all(inds{:})];
                         
                            end
                             
                            % Change the index in the correct dimension
                            % (to new matrix).
                            inds{groupedDim(1)} = groupi1;
                            
                            % Take the mean. 
                            data_matrix(groupi1)=nanmean(holder);
                        end
                    end
           
                case 2  % 2 dimensions to be grouped. 
                    % For each group in the group list, 
                    for groupi1=1:size(group_list,1)
                        group1=group_list{groupi1,2};   

                        % In a second dimension of grouping, now search if there's a
                        % group in all the 2nd dimensions 
                        for groupi2 = 1:size(group_list,1)
                            group2 = group_list{groupi2,2};

                            % Only if both are not empty. 
                            if any(isnan(group1)) == 0 && any(isnan(group2)) == 0

                                % Initialize an empty matrix for holding data to
                                % average into the same group. 
                                holder=[];
                              
                                % See if there are multiple entries in this group
                                for i = 1:numel(group1)
                                    
                                    % Get the correct index for the
                                    % dimension.
                                    inds{groupedDim(1)} = group1(i);
                                    
                                    % See if there are multiple entries in this group.
                                    for ii=1:numel(group2)
                                        
                                        % Get the correct index for this
                                        % dimension.
                                        inds{groupedDim(2)} = group2(ii); 
                                        
                                        % Put them all together for
                                        % averaging .
                                        holder=[holder; data_all(inds{:})];
                                    end  
                                end 
                                
                                % Complicated-looking stuff that makes sure you're
                                % using the correct dimension. (To new matrix).
                                inds = repmat({':'},1,ndims(data_matrix));
                                inds{groupedDim(1)} = groupi1;
                                inds{groupedDim(2)} = groupi2;

                                % Average the multiple entries.
                                data_matrix(inds{:})=nanmean(holder,1);
                            end
                        end 
                    end
            end 
            
            % Get specific name of variable for saving.
            variable_name_output = CreateFileStrings(output_variable_name,[], [], [], period, false);
             
            % Change matrix name to output-specific name.
            eval([variable_name_output ' = data_matrix;']); 
             
            % Get specific name of file for saving.
            filename_output = CreateFileStrings(output_file_name,[], [], [], period, false);
             
            % Save
            save([dir_out filename_output], variable_name_output);
      
        end
    end
end