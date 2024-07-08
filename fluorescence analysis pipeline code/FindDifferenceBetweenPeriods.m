% FindDifferenceBetweenPeriods.m
% Sarah West
% 2/3/22

function [] = FindDifferenceBetweenPeriods(parameters)

    % Say where data is being saved
    disp(['Data saved in '  parameters.dir_out_base{1}]); 

    % Make colormap for difference plots.
    mymap=flipud(cbrewer('div', 'RdBu', 256, 'linear')); 

    % Load first list of periods
    filename = CreateFileStrings([parameters.dir_in_periods1_base parameters.input_periods1_filename], [], [], [], [], false);
    load(filename);
    variable = CreateFileStrings([parameters.input_periods1_variable], [], [], [], [], false);
    eval(['Periods1 = ' variable ';']);

    % If comparing all to one period, make a Periods2 list
    if parameters.all_from_one_period 
        
        % Make a repeating list of the one period that's the same length as
        % Periods1
        Periods2 = repmat({parameters.one_period}, numel(Periods1), 1); 
    
    else
        % If not comparing all to one period, load the Periods2 list.
        filename = CreateFileStrings([parameters.dir_in_periods1_base parameters.input_periods1_filename], [], [], [], [], false);
        load(filename);
        variable = CreateFileStrings([parameters.input_periods1_variable], [], [], [], [], false);
        eval(['Periods2 = ' variable ';']);
    end 

    % For each mouse 
    for mousei=1:size(parameters.mice_all,2)
        mouse=parameters.mice_all(mousei).name;
      
        % For each comparison
        for compi = 1:numel(Periods1)

            % Get the first and second period names
            period1 = Periods1{compi};
            period2 = Periods2{compi}; 

            % If either period name entries are empty, skip to next loop
            % iteration.
            if isempty(period1) || isempty(period2)
                continue
            end
            
            % Get the file names for each.
            filename1 = CreateFileStrings([parameters.dir_in_base parameters.input_filename], mouse, [], [], period1, false);
            filename2 = CreateFileStrings([parameters.dir_in_base parameters.input_filename], mouse, [], [], period2, false);
            
            % If either file doesn't exist, skip to next loop iteration
            if ~isfile(filename1) || ~isfile(filename2)
                continue
            end
           
            % Load period data, convert variable to something generic.

            % Period 1
            load(filename1);
            variable = CreateFileStrings([parameters.input_variable], mouse, [], [], period1, false);
            eval(['PeriodData1 = ' variable ';']);

            % Period 2
            load(filename2);
            variable = CreateFileStrings([parameters.input_variable], mouse, [], [], period2, false);
            eval(['PeriodData2 = ' variable ';']);

            % If either period data is empty, skip to next loop iteration.
            if isempty(PeriodData1) || isempty(PeriodData2)
                continue
            end

            % Take difference 
            difference = PeriodData1 - PeriodData2;
           
            % Get the output variable name, convert to specific
            % variable name.
            output_variable_name = CreateFileStrings(parameters.output_variable, mouse, [], [], [], false);
            eval([output_variable_name '= difference;']); 

            % Make ouput directory name, make if doesn't already exist.
            dir_out = CreateFileStrings(parameters.dir_out_base, mouse, [], [], period1, false);
            if ~isfolder(dir_out)
                mkdir(dir_out);
            end

            % Get the right names for saving per period.
            % [right now only saves with period1]
            saving_filename = CreateFileStrings(parameters.output_filename, mouse , [], [], period1, false);
            
            % Save per period. 
            save([dir_out saving_filename], output_variable_name);
           
            % Plot and save 
            
            % Close any previous figures.
            close all;
            
            % Plot diff per period.
            figure; imagesc(difference); colormap(mymap); 
            caxis(parameters.color_range); colorbar; 
            
            title(['grouped regions, m' mouse ', difference of ' period1 ' from ' period2]);
            
            % Make a figure filename
            % [right now only saves with period1]
            figure_saving_filename = saving_filename;
            if all(figure_saving_filename([end-2:end]) == 'mat')
               figure_saving_filename([end-2:end]) = 'fig';
            end

            % Save figure;
            savefig([dir_out figure_saving_filename ]);

        end
    end
end