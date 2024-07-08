% create_RandomSubset_mice_all_RunAnalysis.m
% Sarah West
% 2/22/22

% Takes mice_all and gets a list of a subset of stacks distributed across
% all days per mouse, randomly chosen within days. Is so you don't have to
% use as much resources for SVD compression on MSI.

% Parameters for directories

experiment_name='Random Motorized Treadmill';
dir_base='Y:\Sarah\Analysis\Experiments\';
dir_exper=[dir_base experiment_name '\']; 
dir_out = [dir_exper 'preprocessing\stack means\']; 
mkdir(dir_out);

% Will want to double-check the existance of and size of the files (which
% could cause concatenation errors), so give the location of the data. 
file_format_cell = {dir_exper, 'preprocessing\fully preprocessed stacks\', 'mouse number', '\', 'day', '\', 'data', 'stack number', '.mat'};
input_variable = 'data';

% Load mice_all
load([dir_exper 'mice_all.mat']);

% Adjust here if you want to use only some mice.
mice_all = mice_all;

% For each mouse you want, enter the pixel number and frame number you want
% per stack (to check them before putting into the big MSI matrix).
% stack_parameters = [38942, 6000;
%                     43031, 6000;
%                     46335, 6000];

% Paramters for randomizing--Fields you want representation from, and amount of stacks you want (per
% mouse) represented from each.
fields = {'stacks', 15 ;
          'spontaneous', 8};

% Number of digits you want in your stack number name
digitNumber = 2; 

% Make list of stacks available for using.
available_mice_all = mice_all;

% Make empty holder of randomized stacks. 
for mousei = 1:size(mice_all,2)
    random_mice_all(mousei).name = mice_all(mousei).name;
    for dayi = 1:size(mice_all(mousei).days, 2)
        random_mice_all(mousei).days(dayi).name = mice_all(mousei).days(dayi).name;
        for fieldi = 1:size(fields,1)
            field = fields{fieldi, 1};
            random_mice_all(mousei).days(dayi).(field) = {};
        end
    end
end

% For now, prevent overwriting previous randomizations
if isfile([dir_out 'mice_all_random.mat'])
   error('mice_all_RandomSubset.mat already exists!')
end    

% For each field
for fieldi = 1:size(fields,1)
    field = fields{fieldi, 1};

    % For each mouse 
    for mousei = 1:size(mice_all,2)
        mouse = mice_all(mousei).name;
        disp(mouse);
        % Create entry in random_mice_all
        random_mice_all(mousei).name = mouse;
        
        % Will go through each day, initialize counter.
        dayi = 1;

        % Until you hit the desired number of stacks for the field
        stack_total = 0;
       
        while stack_total < fields{fieldi, 2}

            % Get day name
            day = mice_all(mousei).days(dayi).name; 

            % If the field exists in this day and is neither NaN nor empty
            if ~isfield(available_mice_all(mousei).days(dayi), field) || isempty(getfield(available_mice_all(mousei).days(dayi), field))                 dayi = dayi +1;
                dayi = dayi + 1;
                if dayi > size(mice_all(mousei).days,2)
                    dayi = 1; 
                end
                continue
            end
            % If NaN, set the corresponding random output to NaN, continue 
           % if any(isnan(getfield(available_mice_all(mousei).days(dayi), field)))
%                 eval(['random_mice_all(mousei).days(dayi).' field '= NaN;']);
%                 dayi = dayi + 1;
%                 if dayi > size(mice_all(mousei).days,2)
%                     dayi = 1; 
%                 end
%                 continue
            %end    

            % Randomly pull a stack from the listed day & field. Put it
            % into random_mice_all, remove it from available_mice_all.
            stack_list = available_mice_all(mousei).days(dayi).(field);
            
            % If there are no stacks left in this day, skip the day.
            if isempty(stack_list)
                dayi = dayi +1;
                if dayi > size(mice_all(mousei).days,2)
                    dayi = 1; 
                end
                continue
            end    

            % Randomly choose a stack
            index = randsample(numel(stack_list),1);
            
            % Get the stack number as a string.
            stack_number = stack_list{index};

            % Check the existance and size of the stack-- get file name
            filename = CreateFileStrings(file_format_cell, mouse, day, stack_number, [], false);
            
            % Check existance of the stack
            if ~isfile(filename)
                dayi = dayi +1;
                if dayi > size(mice_all(mousei).days,2)
                    dayi = 1; 
                end
                continue
            end

            % Check sizes of the stack.
%             matObj = matfile(filename);
%             pixels = size(matObj,input_variable, 1);
%             frames = size(matObj, input_variable,2);
%             if pixels ~= stack_parameters(mousei,1) || frames ~= stack_parameters(mousei,2)
%                 dayi = dayi +1;
%                 if dayi > size(mice_all(mousei).days,2)
%                     dayi = 1; 
%                 end
%                 continue
%             end

            random_mice_all(mousei).days(dayi).name = mice_all(mousei).days(dayi).name;
            random_mice_all(mousei).days(dayi).(field)(end+1)= stack_list(index);
            available_mice_all(mousei).days(dayi).(field) = stack_list([1:index-1 index+1:end]);
            
            % Increase stack counter
            stack_total = stack_total + 1;

            % Increase day iterator
            dayi = dayi + 1;

            % If this increase to day iterator is greater than the number
            % of days, put back to 1.
            if dayi > size(mice_all(mousei).days,2)
                dayi = 1; 
            end
        end     
    end
end 
mice_all = random_mice_all;
save([dir_out 'mice_all_random.mat'], 'mice_all');
