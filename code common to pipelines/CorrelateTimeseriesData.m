% CorrelateTimeseriesData.m
% Sarah West
% 2/2/21

% A function that correlates timeseries data.
% Is general, so you can use with all kinds of timeseries data.

% need inputs: parameters.sourceDim, parameters.timeDim

function [parameters] = CorrelateTimeseriesData(parameters)
   
    % Display progress message to user.
    MessageToUser('Correlating ', parameters);

    % If no segment ranges, 
    if isempty(parameters.data)
        % Do nothing. Put correlation as empty as well.
        parameters.correlation =[]; 
       
    else
        % Number of sources (for convenience)
        num_sources = size(parameters.data, parameters.sourceDim);

        % If there are more dimensions than just the 2 necessary ones, 
        number_of_dimensions = ndims(parameters.data);
        if number_of_dimensions > 2

            extra_dimensions = setdiff([1:number_of_dimensions], [parameters.sourceDim parameters.timeDim]);
           
            % Get size of each dimension
            extra_dimension_sizes = size(parameters.data, extra_dimensions);

            % Make a list of dimensions to grab. 
            loop_list.iterators = cell(number_of_dimensions-2, 1); 

            % Cycle through each dimension more than just the 2 necessary ones, 
            for dimi = extra_dimensions

                loop_list.iterators(dimi-2, 1:3) = {['dim' num2str(dimi)], {['1:' num2str(extra_dimension_sizes(dimi-2))]}, ['dim' num2str(dimi) '_iterator']};

                % Grab the necessary blocks of timeseries. 
            end 

            looping_output_list = LoopGenerator(loop_list, []);

            % Make a list of grabbign dimensions. 
            dimensions_extraction = repmat({':'},1, number_of_dimensions);

            fields =  fieldnames(looping_output_list); 

            % Make a holder correlation matrix. 
            correlations = NaN([num_sources, num_sources, extra_dimension_sizes]);

            % For each looping output,
            for itemi = 1:size(looping_output_list,1)
                
                for dimi = 1:numel(extra_dimensions)
                    dim =  getfield(looping_output_list, {itemi}, fields{dimi *2});
                    dimensions_extraction{extra_dimensions(dimi)} = dim;
                end

                % Grab needed timeseries.
                data = parameters.data(dimensions_extraction{:});
              
                % Permute to put dimensions in best place, if needed (time
                % in first dimension)
                if parameters.sourceDim < parameters.timeDim

                    data = permute(data, [2 1]);

                end

                % Correlate
                [corrs] = corrcoef(data);

                % Put into a holding matrix of all corrs. 
                correlations(dimensions_extraction{:}) = corrs;

            end 

            % Pass to parameters object.
            parameters.correlation = correlations;

   
        % Or else just 2 dimensions (the minimum)
        else 

             % permute dimensions so time is in correct place (first
             % dimention)
             data = permute(parameters.data, [parameters.timeDim parameters.sourceDim]); 

             % Correlate
             parameters.correlation = corrcoef(data);

        end 
    end
end 
