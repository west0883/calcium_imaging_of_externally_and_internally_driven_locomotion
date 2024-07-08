% SegmentTimeseriesData.m
% Sarah West
% 10/13/21

% A function that segments timeseries data using a list of start and end points for each behvior. 
% Is general, so you can use with all kinds of timeseries data. 


function [parameters] = SegmentTimeseriesData(parameters)
        
        % Display progress message to user.
        MessageToUser('Segmenting ', parameters); 

        % If parameters.time_ranges has more than one entry, in cell form 
        if iscell(parameters.time_ranges)
           
            parameters.segmented_timeseries = cell(size(parameters.time_ranges));

            for celli = 1:numel(parameters.time_ranges)
                if ~isempty(parameters.time_ranges{celli})
                    parameters.segmented_timeseries{celli} = SubSegmenter(parameters.time_ranges{celli}, parameters.timeseries, parameters.segmentDim, parameters.concatDim, celli, parameters); 
    
                end
            end 

        else 
            if ~isempty(parameters.time_ranges)
                celli = [];
                parameters.segmented_timeseries = SubSegmenter(parameters.time_ranges, parameters.timeseries, parameters.segmentDim, parameters.concatDim, celli, parameters); 
            else
                parameters.segmented_timeseries = [];
            end
        end 


        % make this it's own sub-function
end 

function [segmented_timeseries] = SubSegmenter(time_ranges, timeseries, segmentDim, concatDim, celli, parameters)
        
        % Make an empty matrix. 
        if ~isfield(parameters, 'uniformSegments') | parameters.uniformSegments 
            segmented_data = []; 
        else 
            segmented_data = {}; 
        end
       
        % Take the ranges using a flexible number of dimensions
        % C is a holder of as many ':' indices as we need.
        C = repmat({':'},1, ndims(timeseries));

         % For each instance 
        for instancei = 1 : size(time_ranges, 1)
            
            % Convert the ranges to time points. 
            all_ranges = time_ranges(instancei,1):time_ranges(instancei,2);
          
            % Put ranges into list of dimensions.
            C(segmentDim) = {all_ranges}; 
            
            % If user says the final segments should be the same size, try
            % to concatenate
            if ~isfield(parameters, 'uniformSegments') | parameters.uniformSegments 
                % Try to concatenate; if one doesn't fit the
                % dimensions, skip it for now.
                try 
                    segmented_data =cat(concatDim, segmented_data, timeseries(C{:})); 
    
                     % Put into output 
                     segmented_timeseries = segmented_data; 
                catch 
                    if nargin > 3
                    disp(['Dimension error in ' num2str(celli)]);
    
                    else 
                        disp(['Dimension error.']);
                    end
    
                    segmented_timeseries = [];
                    continue
                end

            % If not uniform segments, put into cell array
            else
                segmented_data = [segmented_data; {timeseries(C{:})} ];

                segmented_timeseries = segmented_data;
            
            end 
        end
end 