% NullDistribution
% Sarah West
% 5/19/22

% To be run with RunAnalysis. Makes shufflings between two datasets (parameters.data1 &
% parameters.data2). Also takes the mean of each distribution and subtracts
% the mean of dataset 2 from the mean of dataset 1, for convenience. 

function [parameters] = NullDistribution(parameters)

    % Tell user what's happening.
    MessageToUser('Shuffling ', parameters);

    data1 = parameters.data1; 
    data2 = parameters.data2; 
    shuffleDim = parameters.shuffleDim; 
    shuffle_along_dimension = parameters.shuffle_along_dimension;
    shuffleNumber = parameters.shuffleNumber;
    saveShuffled = parameters.saveShuffled;
    mix_distributions = parameters.mix_distributions;
    reorder_one_distribution = parameters.reorder_one_distribution;
    
    % Set up matrix to hold shuffled data, depending on size
    % of the data being compared and the dimensions being shuffled.
    % Only one dimension should be shuffled. The last dimension is the the
    % number of shuffles.
    dimensions_list_1 = [size(data1) shuffleNumber];
    dimensions_list_2 =[size(data2) shuffleNumber];

    % Complicated-looking stuff that lets you select the correct
    % dimension.
    inds_data = repmat({':'}, 1, ndims(data1) + 1);
    
    % Set up holders, clear any existing values. First entry of 
    % null_distributions will be shuffled dataset 1, second will 
    % be shuffled dataset 2. Keep the distributions together in same 
    % file/variable/context, because they're meaningless on their own. 
    
    if saveShuffled
        distributions1 = NaN(dimensions_list_1);
        distributions2 = NaN(dimensions_list_2);
    end

    % If you're mixing two distributions together,
    if mix_distributions

        % Dimensions of differences will be the same as the datasets, but with
        % the shuffle dimension (instances) removed (will have been averaged
        % across that dimension).
        differences_dimensions_list = [size(data1) shuffleNumber];
        differences_dimensions_list(shuffleDim) = [];
        inds_difference = repmat({':'}, 1, ndims(data1));    
        differences_of_distributions = NaN(differences_dimensions_list);
    
        % Begin shufflings.
        
        % For each repetition, (can't do parfor because of the unknown
        % dimensions)
        for shufflei = 1: shuffleNumber
    
            progress = rem(shufflei, 100);
            if progress == 0
                disp([num2str(shufflei) 'th permutation'])
            end
            
            % Change the index in the correct dimension.                
            inds_data{end} = shufflei;
            inds_difference{end} = shufflei;
            
            % Concatenate all data of both periods.
            all_data = cat(shuffleDim, data1, data2); 
            
            % Make a mixing vector that's made up of a random
            % permutation of the number of total instances of both
            % periods.
            vect_mix = randperm(size(all_data, shuffleDim));
            
            % Apply this mixing vector to the concatenated data and split the
            % resulting shuffled data in two new shuffled datasets that
            % have the same number of instances as the two periods that 
            % were inputted.
            inds_mix_1 = repmat({':'}, 1, ndims(data1));
            inds_mix_2 = repmat({':'}, 1, ndims(data2));
    
            inds_mix_1{shuffleDim} = vect_mix(1:size(data1, shuffleDim));
            inds_mix_2{shuffleDim} = vect_mix((size(data1, shuffleDim) + 1):end);
            
            holder1 = all_data(inds_mix_1{:});
            holder2 = all_data(inds_mix_2{:});
            
            if saveShuffled
                distributions1(inds_data{:}) = holder1;
                distributions2(inds_data{:}) = holder2;
            end 
    
            % Take the difference of the mean of the two shuffled
            % datasets and put it in the right place.
            differences_of_distributions(inds_difference{:}) = mean(holder2, shuffleDim, 'omitnan') - mean(holder1, shuffleDim, 'omitnan');
    
        end

    % If you're permuting one distribution in relation to the other;
    % (permuting data2 in relation to data1)
    elseif reorder_one_distribution

        % for each shuffle, 
        for shufflei = 1: shuffleNumber
    
            % report every 100th shuffle to command line
            progress = rem(shufflei, 100);
            if progress == 0
                disp([num2str(shufflei) 'th permutation'])
            end
    
            % get indices for where you'll put the shuffles 
            indices_out = repmat({':'}, 1, numel(dimensions_list_1));
            indices_out{shuffleDim}= shufflei;
    
            % distributions1 is data1 unaltered
            distributions1(indices_out{:}) = data1; 

            % Get random mixing of indices, for data2
            vect_mix = randperm(size(data2, shuffle_along_dimension));

            % get indices to pull out of data2
            indices_in = repmat({':'}, 1, numel(dimensions_list_2));
            indices_in{shuffle_along_dimension} = vect_mix;
           
            distributions2(indices_out{:}) = data2(indices_in{:});

        end 

    % If no method has been specified
    else
        error('Specify null distribution type.');
    end

    % Put into output structure.
    if saveShuffled
        parameters.null_distributions{1} = distributions1;
        parameters.null_distributions{2} = distributions2; 
    end 
    if mix_distributions
        parameters.differences_of_distributions = differences_of_distributions;
    end
        
end